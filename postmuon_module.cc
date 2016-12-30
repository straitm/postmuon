////////////////////////////////////////////////////////////////////////
/// \brief   This module is named PostMuon and does bar.
/// \author  M. Strait
/// \date
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

#include <string>
#include <algorithm>

#include <signal.h>

#include "art/Framework/Core/EDAnalyzer.h"

struct pm{
  int cluster_nhit;
  int last_accepted_time;
  int last_accepted_i;
  float mindist;  // minimum hit distance from track end
  float dist2sum; // summed hit distance from track end
  int asum;       // sumed ADC of a cluster of delayed hits
  float esum;     // summed calibrated energy, in MeV
  int cluster_i;
};


/// Calibrating RawData to Produce CellHits
namespace PostMuon {

  class PostMuon : public art::EDAnalyzer {

  public:

    explicit PostMuon(fhicl::ParameterSet const& pset);
    virtual ~PostMuon();

    void analyze(const art::Event& evt);

  private:

    int         fRemoveBadChans; ///< whether to remove bad channels
    std::string fRawDataLabel;   ///< label of where to find RawData

  }; // class PostMuon
}


static double kUSEC_PER_TDC = 1./64.;
static double kUSEC_PER_MICROSLICE = 0.5; // XXX right?

// Geometrically about correct, but perhaps should be scaled by density or
// radiation length or neutron cross section or something.  Or not, since
// which of those is right depends on what you're looking at.
const double planes_per_cell = 76./39.;

/*
 Returns true if the track goes in the +z direction, assuming that it is
 downward-going (i.e. that it goes in the -y direction).
*/
static bool is_increasing_z(const rb::Track & trk)
{
  const bool claims_increasing_z = trk.Stop().Z() > trk.Start().Z();
  const bool needs_flip = trk.Stop().Y() > trk.Start().Y();
  return claims_increasing_z ^ needs_flip;
}

/*
  Takes a track and a hit and determines the distance between the
  end of the track and the hit.  lasthiti_{even,odd} are the indices
  of the last hits in each view assuming that the track is going in the
  -y direction (down).
*/
static float dist_trackend_to_cell(const rb::Track & trk,
                                   const rb::CellHit & chit,
                                   const int lasthiti_even,
                                   const int lasthiti_odd)
{
  const int lastplane_even = trk.Cell(lasthiti_even)->Plane();
  const int lastplane_odd =  trk.Cell(lasthiti_odd) ->Plane();
  const int lastcell_even =  trk.Cell(lasthiti_even)->Cell();
  const int lastcell_odd =   trk.Cell(lasthiti_odd) ->Cell();

  const bool increasing_z = is_increasing_z(trk);

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
 Return the average hit time of the latest 40 hits (or all of them).

 This is meant to provide a robust measure of the track time, even if
 an early Michel decay gets reconstructed as part of the track.
*/
static int mean_late_track_time(const rb::Track & trk)
{
  std::vector<int> tdcs;
  for(unsigned int i = 0; i < trk.NCell(); i++)
    tdcs.push_back(trk.Cell(i)->TDC());

  std::sort(tdcs.begin(), tdcs.end());

  const int max_for_avg = 40;

  float acc = 0;
  for(unsigned int i = std::max(0, (int)tdcs.size() - max_for_avg);
      i < tdcs.size(); i++)
    acc += tdcs[i];

  return int(acc / std::min(max_for_avg, (int)tdcs.size()) + 0.5);
}

/*
 If the hit is near the track and after it in time, return the distance
 in cells (where a plane is ~2 cells).  Otherwise, return -1.
*/
static double hit_near_track(const rb::Track & trk,
  const int lasthiti_even, const int lasthiti_odd,
  const art::Handle< std::vector<rb::CellHit> > & cellcol,
  const int celli)
{
  const rb::CellHit & chit = (*cellcol)[celli];

  const int tracktime = mean_late_track_time(trk);

  const int aftertdc = chit.TDC();
  if(aftertdc < tracktime) return -1;

  // Accept the hit even if it is in the track!  Because if a Michel
  // hit gets swept up into the track, this is the only way to see it.
  //if(hit_is_in_track(chit, trk)) return -1;

  // Enough to catch ~99% of gammas from neutron capture according to my toy MC
  const double maxdist_in_cells = 20.1;

  const double dist = dist_trackend_to_cell(trk, chit, lasthiti_even,
                                            lasthiti_odd);

  if(dist > maxdist_in_cells) return -1;

  return dist;
}

/*
  Passes back the last hit in each plane, assuming that the track is
  downwards-going (i.e. in the -y direction).
*/
static void last_hits(int & lasthiti_even,
                      int & lasthiti_odd,
                      const rb::Track & trk)
{
  float latest_even = 1e30, latest_odd = 1e30;

  const bool increasing_z = is_increasing_z(trk);
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

static void print_ntuple_line(const art::Event & evt,
                              const rb::Track & trk,
                              const double eventlength,
                              const int hittime_tdc,
                              const pm answer)
{
  const int tracktime = mean_late_track_time(trk);
  const double timeleft = eventlength - tracktime *kUSEC_PER_TDC;

  const bool needs_flip = trk.Stop().Y() > trk.Start().Y();

  const double tx = needs_flip? trk.Start().X(): trk.Stop().X();
  const double ty = needs_flip? trk.Start().Y(): trk.Stop().Y();
  const double tz = needs_flip? trk.Start().Z(): trk.Stop().Z();

  const double tsx = needs_flip? trk.Stop().X(): trk.Start().X();
  const double tsy = needs_flip? trk.Stop().Y(): trk.Start().Y();
  const double tsz = needs_flip? trk.Stop().Z(): trk.Start().Z();

  printf("ntuple: %d %d %d %f "
                 "%f %f %f "
                 "%f %f %f %f "
                 "%f %d %f %d %f\n",
    evt.run(), evt.event(),
    answer.cluster_i,
    (hittime_tdc - tracktime)*kUSEC_PER_TDC,
    tsx, tsy, tsz,
    tx, ty, tz,
    answer.mindist,
    answer.dist2sum/answer.cluster_nhit,
    answer.asum, answer.esum, answer.cluster_nhit, timeleft);
}

namespace PostMuon{

PostMuon::PostMuon(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset), fRemoveBadChans(pset.get<bool>("RemoveBadChans")),
  fRawDataLabel(pset.get< std::string >("RawDataLabel"))
{

}

static bool compare_cellhit(const rb::CellHit & a, const rb::CellHit & b)
{
  return a.TDC() < b.TDC();
}

PostMuon::~PostMuon() { }

static pm mkpm()
{
  pm res;
  res.cluster_nhit = 0;
  res.last_accepted_time = -1;
  res.last_accepted_i = -1;
  res.mindist = 1000000;
  res.dist2sum = 0;
  res.asum = 0;
  res.esum = 0;
  res.cluster_i = 0;
  return res;
}

static void reset_pm(pm & answer)
{
  answer.cluster_nhit = 0;

  answer.mindist = 1000000;
  answer.dist2sum = 0;
  answer.asum = 0;
  answer.esum = 0;
}

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
    static int NOvA = printf(
      "ntuple: run:event:i:t:trkstartx:trkstarty:trkstartz:"
      "trkx:trky:trkz:mindist:dist2:adc:e:nhit:timeleft\n");
    NOvA = NOvA;
  }

  art::ServiceHandle<calib::Calibrator> calthing;

  if(rawtrigger->empty()) return;

  const double triggerlength = (*rawtrigger)[0].
    fTriggerRange_TriggerLength*kUSEC_PER_MICROSLICE;

  for(unsigned int t = 0; t < tracks->size(); t++){

    const rb::Track & trk = (*tracks)[t];
    int lasthiti_even = 0, lasthiti_odd = 0;
    last_hits(lasthiti_even, lasthiti_odd, trk);

    pm answer = mkpm();

    // CellHits do *not* come in time order
    std::vector<rb::CellHit> sorted_hits;
    for(int c = 0; c < (int)cellcol->size(); c++)
      sorted_hits.push_back((*cellcol)[c]);
    std::sort(sorted_hits.begin(), sorted_hits.end(), compare_cellhit);

    for(int c = 0; c < (int)cellcol->size(); c++){
      const double dist = hit_near_track(trk, lasthiti_even,
         lasthiti_odd, cellcol, c);

      if(dist < 0) continue;

      answer.dist2sum += dist*dist;
      if(dist < answer.mindist) answer.mindist = dist;

      const rb::CellHit & chit = (*cellcol)[c];

      const bool needs_flip = trk.Stop().Y() > trk.Start().Y();
      const rb::RecoHit rhit = calthing->MakeRecoHit(chit,
         chit.Plane()%2 == 0? (needs_flip?trk.Start().X():trk.Stop().X()): // doc-11570
                              (needs_flip?trk.Start().Y():trk.Stop().Y()));

      // Hits are time ordered.  Only report on hits that are separated
      // in time from other accepted hits by at least 1 quiet TDC tick.
      const bool newcluster = chit.TDC() > answer.last_accepted_time + 1;

      if(newcluster && answer.cluster_nhit){
        answer.cluster_nhit++;
        print_ntuple_line(evt, trk, triggerlength,
                          chit.TDC(), answer);
        answer.cluster_i++;
        reset_pm(answer);
      }
      else{
        answer.cluster_nhit++;
      }

      answer.asum += chit.ADC();
      // For now, just throw away uncalibrated hits,
      // XXX is there a better approach?
      if(rhit.IsCalibrated()) answer.esum += rhit.GeV()*1000;

      answer.last_accepted_time = chit.TDC();
      answer.last_accepted_i = c;
    }

    // print last cluster
    if(answer.cluster_nhit)
      print_ntuple_line(evt, trk, triggerlength,
                        (*cellcol)[answer.last_accepted_i].TDC(), answer);

  }
}

DEFINE_ART_MODULE(PostMuon);

} // end namespace PostMuon
//////////////////////////////////////////////////////////////////////////
