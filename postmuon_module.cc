////////////////////////////////////////////////////////////////////////
/// \brief   This module writes out a simple ntuple of information about
///          hits that follow stopping muons.
/// \author  M. Strait
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "RecoBase/CellHit.h"
#include "RecoBase/RecoHit.h"
#include "RecoBase/Track.h"
#include "Calibrator/Calibrator.h"

#include "Simulation/ParticleNavigator.h"
#include "RawData/RawTrigger.h"

#include "GeometryObjects/PlaneGeo.h"

#include <string>
#include <algorithm>

#include <signal.h>

#include "art/Framework/Core/EDAnalyzer.h"

static FILE * OUT = NULL;

struct pm{
  int nhit; // number of hits in delayed cluster
            // TODO: seperate nhit counts for the various energy sums.
  int trk; // index of the track
  int first_accepted_time; // time of the first hit in the cluster
  int last_accepted_time; // time of the last hit in the cluster
  float mindist;  // minimum hit distance from track end
  float dist2sum; // summed hit distance from track end
  int tsum;       // summed TDC of a cluster of delayed hits
  int asum;       // summed ADC of a cluster of delayed hits
  float esum;     // summed calibrated energy, in MeV
  float esum_ex;  // Same, but without hits in any module where the track was,
                  // in order to attempt to get an unbiased, sag-free measurement
  float esum_ex2; // Same, but also excluding any module that had a hit coincident
                  // with the track, where anything within 0.4us counts as
                  // coincident.  This is 2-3 sigma in the time resolution and
                  // should catch almost all hits that are prompt, either muon
                  // hits that didn't get reconstructed as part of the track,
                  // or brems, or x-rays from muon atomic capture, or whatever else.
  float esum_xt;  // Same, but without hits that are part of any tracks,
                  // in an attempt to beat down uncorrelated background (doesn't seem
                  // to have much effect -- maybe 5-10%.)
  int nuncal;     // number of uncalibrated hits
  int cluster_i;
};

static pm mkpm()
{
  pm res;
  res.nhit = 0;
  res.trk = 0;
  res.first_accepted_time = -1000000;
  res.last_accepted_time  = -1000000;
  res.mindist = 1000000;
  res.dist2sum = 0;
  res.tsum = 0;
  res.asum = 0;
  res.esum = 0;
  res.esum_ex = 0;
  res.esum_ex2 = 0;
  res.esum_xt = 0;
  res.cluster_i = 0;
  res.nuncal = 0;
  return res;
}

static void reset_pm(pm & res)
{
  res.nhit = 0;
  // do not change trk
  res.first_accepted_time = -1000000;
  res.last_accepted_time  = -1000000;
  res.mindist = 1000000;
  res.dist2sum = 0;
  res.tsum = 0;
  res.asum = 0;
  res.esum = 0;
  res.esum_ex = 0;
  res.esum_ex2 = 0;
  res.esum_xt = 0;
  res.nuncal = 0;
  // do not change cluster_i
}

namespace PostMuon {

  class PostMuon : public art::EDAnalyzer {

    public:

    explicit PostMuon(fhicl::ParameterSet const& pset);
    virtual ~PostMuon();

    void analyze(const art::Event& evt);

    void endJob();

    private:

    int         fRemoveBadChans; ///< whether to remove bad channels
    std::string fRawDataLabel;   ///< label of where to find RawData
    float       fMaxDistInCells; ///< Maximum distance of cluster hits
    bool        fTracksAreDown;  ///< true: down. false: north

  }; // class PostMuon
}


static double USEC_PER_TDC = 1./64.;     // yes, really it is
static double USEC_PER_MICROSLICE = 0.5; // XXX right?

// At least in runs 19107 and 19108, TDCs seem to only come in multiples
// of 4.  Not true later, maybe, because of multipoint readout?
const int TDC_GRANULARITY = 4;

// Geometrically about correct, but perhaps should be scaled by density or
// radiation length or neutron cross section or something.  Or not, since
// which of those is right depends on what you're looking at.
const double planes_per_cell = 76./39.;

/*
 Return the average hit time of the latest 40 hits (or all of them),
 in integer (but not multiple of anything) TDC counts.

 This is meant to provide a robust measure of the track time, even if
 an early Michel decay gets reconstructed as part of the track.

 This also smears out weird timing effects like those shown in
 doc-16889-v2 which may or may not be a good thing.
*/
static double mean_late_track_time(const rb::Track & trk)
{
  std::vector<int> tdcs;
  for(unsigned int i = 0; i < trk.NCell(); i++)
    tdcs.push_back(trk.Cell(i)->TDC());

  std::sort(tdcs.begin(), tdcs.end());

  const int max_for_avg = 40;

  double acc = 0;
  for(unsigned int i = std::max(0, (int)tdcs.size() - max_for_avg);
      i < tdcs.size(); i++)
    acc += tdcs[i];

  return int(acc / std::min(max_for_avg, (int)tdcs.size()) + 0.5);
}

// Returns true if two cells (given that they are in the
// same plane) are in the same module.
static bool same_module(const int cell1, const int cell2)
{
  // Cells are numbered 0-31, 32-63, etc, doc-11570.
  // True in both Near and Far.  So this is really easy.
  return cell1/32 == cell2/32;
}

/*
  Returns true if the given hit is in the same module as any hit coincident
  (within the time given below) with the given track.  Along with excluding the
  track itself, this is meant to do better at excluding cells that are
  inefficient.
*/
static bool hit_in_track_coincident_module(const rb::CellHit & chit,
                                           const rb::Track & trk,
                                           const std::vector<rb::CellHit> & sorted_hits)
{
  const double timecut = 0.4/USEC_PER_TDC;

  const double tracktime = mean_late_track_time(trk);

  for(unsigned int i = 0; i < sorted_hits.size(); i++){
    const rb::CellHit & coin_chit = sorted_hits[i];

    if(tracktime - coin_chit.TDC() > +timecut) continue;
    if(tracktime - coin_chit.TDC() < -timecut) break;

    if(coin_chit.Plane() == chit.Plane() &&
       same_module(coin_chit.Cell(), chit.Cell())) return true;
  }
  return false;
}

/*
  Returns true if the given hit is in the same module as any hit in the
  given track.  This is meant to exclude cells that are inefficient for
  detecting Michel hits because of APD sag (see doc-12802, etc.) and/or
  the lack of second rising edge, so there is no time requirement.
*/
static bool hit_in_track_module(const rb::CellHit & chit,
                                const rb::Track & trk)
{
  for(unsigned int i = 0; i < trk.NCell(); i++){
    const rb::CellHit & trk_chit = *(trk.Cell(i));

    if(trk_chit.Plane() == chit.Plane() &&
       same_module(trk_chit.Cell(), chit.Cell())) return true;
  }
  return false;
}


/* Returns true if this cell hit is in any track in the event (or, really
 * if it is contained in the list of hits passed in).
 *
 * This is a dumb implementation.  Could sort the list of track hits and use a
 * stdlib function to see if it is in the list.  However, this uses only a tiny
 * fraction of the time, so don't bother. Tested for both regular and DDsnews
 * triggers.
 */
static bool hit_in_any_track(const rb::CellHit & chit,
                             const std::vector<rb::CellHit> & trkhits)
{
  for(unsigned int i = 0; i < trkhits.size(); i++)
    if(trkhits[i] == chit) return true;
  return false;
}

/*
  Returns true if we need to reverse a track based on the assumed
  direction (down or north).
*/
static bool shall_we_flip_it(const rb::Track & trk, const bool tracks_are_down)
{
  if(tracks_are_down) return trk.Stop().Y() > trk.Start().Y();
  else                return trk.Stop().Z() < trk.Start().Z();
}

/*
 Returns true if the track goes in the +z direction, assuming that it is
 downward-going (i.e. that it goes in the -y direction).
*/
static bool is_increasing_z(const rb::Track & trk, const bool tracks_are_down)
{
  const bool claims_increasing_z = trk.Stop().Z() > trk.Start().Z();
  const bool needs_flip = shall_we_flip_it(trk, tracks_are_down);
  return claims_increasing_z ^ needs_flip;
}

/*
  Takes a track and a hit and determines the distance between the
  end of the track and the hit.  lasthiti_{even,odd} are the indices
  of the last hits in each view assuming that the track is going in the
  -y direction (down).
*/
static float dist_trackend_to_cell(const rb::Track & __restrict__ trk,
                                   const rb::CellHit & __restrict__ chit,
                                   const int lasthiti_even,
                                   const int lasthiti_odd,
                                   const bool tracks_are_down)
{
  const int lastplane_even = trk.Cell(lasthiti_even)->Plane();
  const int lastplane_odd =  trk.Cell(lasthiti_odd) ->Plane();
  const int lastcell_even =  trk.Cell(lasthiti_even)->Cell();
  const int lastcell_odd =   trk.Cell(lasthiti_odd) ->Cell();

  const bool increasing_z = is_increasing_z(trk, tracks_are_down);

  const int lastplane = increasing_z?std::max(lastplane_even, lastplane_odd)
                                    :std::min(lastplane_even, lastplane_odd);
  return
    chit.Plane()%2 == 0 ?
      sqrt(pow(planes_per_cell*(chit.Plane() - lastplane)   , 2) +
           pow(                 chit.Cell()  - lastcell_even, 2))
   :  sqrt(pow(planes_per_cell*(chit.Plane() - lastplane)   , 2) +
           pow(                 chit.Cell()  - lastcell_odd , 2));
}

/*
  If the hit is near the track and after it in time, return the distance
  in cells.  Otherwise, return -1 if the hit preceeds the track, or -2
  if it is too far away.
*/
static double hit_near_track(const rb::Track & __restrict__ trk,
  const double tracktime,
  const int lasthiti_even, const int lasthiti_odd,
  const rb::CellHit & __restrict__ chit,
  const float maxdist_in_cells,
  const bool tracks_are_down)
{
  const int aftertdc = chit.TDC();
  if(aftertdc < tracktime) return -1;

  // Accept the hit even if it is in the track!  Because if a Michel
  // hit gets swept up into the track, this is the only way to see it.
  //if(hit_is_in_track(chit, trk)) return -1;

  const double dist = dist_trackend_to_cell(trk, chit, lasthiti_even,
                                            lasthiti_odd, tracks_are_down);

  if(dist > maxdist_in_cells) return -2;

  return dist;
}

/*
  Passes back the last hit in each plane, assuming that the track is
  downwards-going (i.e. in the -y direction).
*/
static void last_hits(int & __restrict__ lasthiti_even,
                      int & __restrict__ lasthiti_odd,
                      const rb::Track & __restrict__ trk,
                      const bool tracks_are_down)
{
  float latest_even = 1e30, latest_odd = 1e30;

  const bool increasing_z = is_increasing_z(trk, tracks_are_down);
  const bool decreasing_z = !increasing_z;

  for(int c = 0; c < (int)trk.NCell(); c++){
    const rb::CellHit & chit = *(trk.Cell(c));
    const rb::RecoHit rhit = trk.RecoHit(c);

    if(chit.Plane()%2 == 0){
      if((decreasing_z && rhit.Z() < latest_even) ||
         (increasing_z && rhit.Z() > latest_even)){
        lasthiti_even = c;
        latest_even = rhit.Z();
      }
    }
    else{
      if((decreasing_z && rhit.Z() < latest_odd) ||
         (increasing_z && rhit.Z() > latest_odd)){
        lasthiti_odd = c;
        latest_odd = rhit.Z();
      }
    }
  }
}

static void print_ntuple_line(const art::Event & __restrict__ evt,
                              const rb::Track & __restrict__ trk,
                              const double tracktime,
                              const double eventlength,
                              const pm answer,
                              const bool tracks_are_down)
{
  const double timeleft = eventlength - tracktime*USEC_PER_TDC;

  const bool needs_flip = shall_we_flip_it(trk, tracks_are_down);

  const double tsx = needs_flip? trk.Stop().X(): trk.Start().X();
  const double tsy = needs_flip? trk.Stop().Y(): trk.Start().Y();
  const double tsz = needs_flip? trk.Stop().Z(): trk.Start().Z();

  const double tx = needs_flip? trk.Start().X(): trk.Stop().X();
  const double ty = needs_flip? trk.Start().Y(): trk.Stop().Y();
  const double tz = needs_flip? trk.Start().Z(): trk.Stop().Z();

  fprintf(OUT, "%d %d ", evt.run(), evt.event());
  fprintf(OUT, "%d ", answer.trk);
  fprintf(OUT, "%f ", tracktime*USEC_PER_TDC);
  fprintf(OUT, "%d ", answer.cluster_i);
  if(answer.nhit)
    fprintf(OUT, "%f ", (float(answer.tsum)/answer.nhit-tracktime)*USEC_PER_TDC);
  else
    fprintf(OUT, "-1 ");

  fprintf(OUT, "%.1f %.1f %.1f ", tsx, tsy, tsz);
  fprintf(OUT, "%.1f %.1f %.1f %.3f ", tx, ty, tz, answer.mindist);

  if(answer.dist2sum == 0) fprintf(OUT, "0 ");
  else fprintf(OUT, "%.3f ", answer.dist2sum/answer.nhit);

  fprintf(OUT, "%d %.3f %.3f %.3f %.3f ",
    answer.asum, answer.esum, answer.esum_ex, answer.esum_ex2, answer.esum_xt);

  fprintf(OUT, "%f %d ", timeleft, answer.nhit);
  fprintf(OUT, "%d ", answer.nuncal);
  fprintf(OUT, "%d ", answer.last_accepted_time - answer.first_accepted_time);
  fprintf(OUT, "\n");
}

namespace PostMuon{

PostMuon::PostMuon(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset), fRemoveBadChans(pset.get<bool>("RemoveBadChans")),
  fRawDataLabel(pset.get< std::string >("RawDataLabel")),
  fMaxDistInCells(pset.get<float>("MaxDistInCells")),
  fTracksAreDown(pset.get<bool>("TracksAreDown"))
{

}

void PostMuon::endJob()
{
  fclose(OUT);
}

// Compare tracks by their time
static bool compare_track(const rb::Track & __restrict__ a,
                          const rb::Track & __restrict__ b)
{
  return mean_late_track_time(a) < mean_late_track_time(b);
}

// Compare cellhits by their time.  Note that CellHit has the operator<
// defined, but that compares by plane, then by cell, then by time, which
// is not what we need here.
static bool compare_cellhit(const rb::CellHit & __restrict__ a,
                            const rb::CellHit & __restrict__ b)
{
  return a.TDC() < b.TDC();
}

static std::vector<rb::CellHit> make_trkhits(const std::vector<rb::Track> & trks)
{
  std::vector<rb::CellHit> answer;
  for(unsigned int i = 0; i < trks.size(); i++)
    for(unsigned int j = 0; j < trks[i].NCell(); j++)
      answer.push_back(*(trks[i].Cell(j)));
  return answer;
}

PostMuon::~PostMuon() { }

void PostMuon::analyze(const art::Event& evt)
{
  // I'm so sorry that I have to do this.  And, my goodness, doing
  // it in the constructor isn't sufficient.
  signal(SIGPIPE, SIG_DFL);

  art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
  evt.getByLabel("daq", rawtrigger);

  art::Handle< std::vector<rb::CellHit> > cellcol; // get hits
  evt.getByLabel(fRawDataLabel, cellcol);

  art::Handle< std::vector<rb::Track> > tracks;
  evt.getByLabel("kalmantrackmerge", tracks);

  {
    if(OUT == NULL){
      OUT = fopen(Form("postmuon_%d_%d.20170120.ntuple", evt.run(), evt.subRun()), "w");
      if(OUT == NULL){
         fprintf(stderr, "Could not open output ntuple file\n");
         exit(1);
      }
    }
    static int NOvA = fprintf(OUT,
      "run:event:trk:trktime:"
      "i:"
      "t:"
      "trkstartx:trkstarty:trkstartz:"
      "trkx:trky:trkz:mindist:"
      "dist2:"
      "adc:e:eex:eex2:ext:"
      "timeleft:nhit:"
      "nuncal:"
      "tdclen"
      "\n");
    NOvA = NOvA;
  }

  art::ServiceHandle<calib::Calibrator> calthing;

  if(rawtrigger->empty()) return;

  const double triggerlength = (*rawtrigger)[0].
    fTriggerRange_TriggerLength*USEC_PER_MICROSLICE;

  // CellHits do *not* come in time order
  std::vector<rb::CellHit> sorted_hits;
  if(!tracks->empty()){
    for(int c = 0; c < (int)cellcol->size(); c++)
      sorted_hits.push_back((*cellcol)[c]);
    std::sort(sorted_hits.begin(), sorted_hits.end(), compare_cellhit);
  }

  std::vector<rb::Track> sorted_tracks;
  if(!tracks->empty()){
    for(int c = 0; c < (int)tracks->size(); c++)
      sorted_tracks.push_back((*tracks)[c]);
    std::sort(sorted_tracks.begin(), sorted_tracks.end(), compare_track);
  }

  std::vector<rb::CellHit> trkhits = make_trkhits(sorted_tracks);

  int first_hit_to_consider = 0;

  for(unsigned int t = 0; t < sorted_tracks.size(); t++){
    const rb::Track & trk = sorted_tracks[t];

    // Badly tracked tracks are not useful
    if(trk.Stop().X() == 0 || trk.Stop().Y() == 0) continue;

    int lasthiti_even = 0, lasthiti_odd = 0;
    last_hits(lasthiti_even, lasthiti_odd, trk, fTracksAreDown);
    const double tracktime = mean_late_track_time(trk);

    pm answer = mkpm();
    answer.trk = t;
    const bool needs_flip = shall_we_flip_it(trk, fTracksAreDown);
    for(int c = first_hit_to_consider; c < (int)sorted_hits.size(); c++){
      const rb::CellHit & chit = sorted_hits[c];

      const double dist =
        hit_near_track(trk, tracktime, lasthiti_even, lasthiti_odd,
                       chit, fMaxDistInCells, fTracksAreDown);

      if(dist < 0){
        // This hit is before this track, so it will also be before
        // all further tracks (since we sorted the tracks above). ~5%
        // speed bump from this optimization.
        if(dist == -1) first_hit_to_consider = c+1;

        continue;
      }

      // Hits are time ordered.  Only report on hits that are separated
      // in time from other accepted hits by at least 1 quiet TDC tick.
      const bool newcluster = chit.TDC() > answer.last_accepted_time
                                           + 1*TDC_GRANULARITY;

      if(newcluster && answer.nhit){
        print_ntuple_line(evt, trk, tracktime, triggerlength, answer, fTracksAreDown);
        answer.cluster_i++;
        reset_pm(answer);
      }

      // Add everthing to the cluster *after* the print of the
      // previous cluster (or no-op).

      answer.nhit++;

      const rb::RecoHit rhit = calthing->MakeRecoHit(chit,
         // If the hit is in X, it needs a Y plane to provide W
         chit.View() == geo::kX? (needs_flip?trk.Start().Y():trk.Stop().Y()):
                                 (needs_flip?trk.Start().X():trk.Stop().X()));

      answer.dist2sum += dist*dist;
      if(dist < answer.mindist) answer.mindist = dist;

      answer.last_accepted_time = chit.TDC();
      if(answer.first_accepted_time < 0)
        answer.first_accepted_time = chit.TDC();
      answer.tsum += chit.TDC();
      answer.asum += chit.ADC();

      if(rhit.IsCalibrated()){
        answer.esum += rhit.GeV()*1000;
        if(!hit_in_track_module(chit, trk)){
          answer.esum_ex += rhit.GeV()*1000;
          if(!hit_in_track_coincident_module(chit, trk, sorted_hits))
            answer.esum_ex2 += rhit.GeV()*1000;
        }
        if(!hit_in_any_track(chit, trkhits)) answer.esum_xt += rhit.GeV()*1000;
      }
      else answer.nuncal++;
    }

    // print last cluster, or track info if no cluster
    print_ntuple_line(evt, trk, tracktime, triggerlength, answer, fTracksAreDown);
  }
}

DEFINE_ART_MODULE(PostMuon);

} // end namespace PostMuon
//////////////////////////////////////////////////////////////////////////
