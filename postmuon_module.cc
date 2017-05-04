////////////////////////////////////////////////////////////////////////
/// \brief   This module writes out a simple ntuple of information about
///          hits that follow stopping muons.
/// \author  M. Strait
////////////////////////////////////////////////////////////////////////

// This has to appear *before* the DAQDataFormats includes, or else
// you get compiler errors that don't seem to have anything to do with
// anything. This solution was found by pure trial and error.
#include "art/Framework/Core/EDAnalyzer.h"

#include "DAQDataFormats/RawEvent.h"
#include "DAQDataFormats/RawTrigger.h"
#include "DAQDataFormats/RawTriggerMask.h"
#include "DAQDataFormats/RawDataBlock.h"
#include "DAQDataFormats/RawMicroBlock.h"
#include "DAQDataFormats/RawMicroSlice.h"
#include "DAQDataFormats/RawNanoSlice.h"
#include "DAQChannelMap/DAQChannelMap.h"
#include "RawData/RawSumDropMB.h"

#include "Simulation/ParticleNavigator.h"

#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Core/FindOneP.h"
#include "Utilities/AssociationUtil.h"


#include "StandardRecord/StandardRecord.h"
#include "MCCheater/BackTracker.h"

#include "RecoBase/CellHit.h"
#include "RecoBase/RecoHit.h"
#include "RecoBase/Track.h"
#include "Calibrator/Calibrator.h"

#include "NumuEnergy/NumuE.h"

#include "RawData/FlatDAQData.h"
#include "RawData/RawTrigger.h"

#include "ReMId/ReMId.h"
#include "ReMId/classes.h"

#include "GeometryObjects/PlaneGeo.h"

#include <string>
#include <algorithm>

#include <signal.h>

static FILE * OUT = NULL;
static int NhitTrackTimeAveraging = 1; // to be overwritten in PostMuon()
static bool TracksAreDown = true; // to be overwritten in PostMuon()
static double MaxDistInCells = 0.123456; // dummy as well


struct evtinfo{
  int run;
  int subrun;
  int event;
  int nslc; // number of slices
  double triggerlength;
  double starttime;
};

struct trkinfo{
  rb::Track trk;
  int i; // index of track in the track array
  int ntrk; // total number of tracks in this event
  int slice; // the index of the slice for this track
  bool contained_slice;
  float slice_energy; // Energy of the slice holding this track
  int true_nupdg; // True PDG code of neutrino making this slice, or 0 for data
  int true_pdg; // True PDG code of particle making this track, or 0 for data
  int true_nucc; // 1 if this slice is MC and true CC, 0 otherwise
  int true_atom_cap; // 1 if this is a mu-, pi- or K- that stops, else 0
  bool primary_in_slice; // Is this the longest track in the slice?
  double time;
  int lasthiti_even, lasthiti_odd;
  float sx, sy, sz, ex, ey, ez; // start and end position, after flip correction
  double remid;
};

struct cluster{
  int i; // which cluster for the given track
  int type; // defines the set of cuts used for this cluster
  float first_accepted_time; // time of the first hit in the cluster
  float last_accepted_time; // time of the last hit in the cluster
  float mindist;  // minimum hit distance from track end
  float maxdist;  // maximum hit distance from track end
  float dist2sum; // summed hit distance from track end
  float tsum;     // summed times, in ns, of a cluster of delayed hits
  float previous_cluster_t; // what it says
  float esum;     // summed calibrated energy, in MeV
  float pesum;    // summed PE
  int   adcsum;   // summed ADC
  int nhit;       // number of hits in full delayed cluster
  int nhitx;      // number of x hits in full delayed cluster
  int nhity;      // number of y hits in full delayed cluster
  int nhitover35pe; // as it says
  int nhitover35pex; //
  int nhitover35pey; //
};

const float INVALID_TIME = -1000000;

static void resetcluster(cluster & res)
{
  res.first_accepted_time = INVALID_TIME;
  res.last_accepted_time  = INVALID_TIME;
  res.mindist = 1000000;
  res.maxdist = 0;
  res.dist2sum = 0;
  res.tsum = 0;
  res.esum = 0;
  res.pesum = 0;
  res.adcsum = 0;
  res.nhit = 0;
  res.nhitx = 0;
  res.nhity = 0;
  res.nhitover35pe = 0;
  res.nhitover35pex = 0;
  res.nhitover35pey = 0;
  res.previous_cluster_t = 0;
}

static cluster mkcluster()
{
  cluster res;
  resetcluster(res);
  res.i = 0;
  res.type = 0;
  return res;
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
    int         fNhitTrackTimeAveraging; ///< N hits for track time average

  }; // class PostMuon
}


static int TDC_PER_US = 64;
static int US_PER_MICROSLICE = 50; // I hope this is always true

// Geometrically about correct, but perhaps should be scaled by density
// or radiation length or neutron cross section or something. Or not,
// since which of those is right depends on what you're looking at.
const double planes_per_cell = 76./39.;

/*
  Return the average hit time of the latest 'NhitTrackTimeAveraging'
  hits, in nanoseconds.

  With a largish 'NhitTrackTimeAveraging', this is meant to provide a
  robust measure of the track time, even if an early Michel decay gets
  reconstructed as part of the track. This also smears out weird timing
  effects like those shown in doc-16889-v2 which may or may not be a
  good thing.
*/
static double mean_late_track_time(const rb::Track & trk)
{
  std::vector<double> times;
  for(unsigned int i = 0; i < trk.NCell(); i++)
    times.push_back(trk.Cell(i)->TNS());

  std::sort(times.begin(), times.end());

  double acc = 0;
  for(unsigned int i = std::max(0, (int)times.size()-NhitTrackTimeAveraging);
      i < times.size(); i++)
    acc += times[i];

  return acc / std::min(NhitTrackTimeAveraging, (int)times.size());
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
  Returns true if the given hit is in the same module as any hit
  coincident (within the time given below) with the given track. Along
  with excluding the track itself, this is meant to do better at
  excluding cells that are inefficient.
*/
static bool hit_in_track_coincident_module(
  const rb::CellHit & __restrict__ chit,
  const trkinfo & __restrict__ tinfo,
  const std::vector<rb::CellHit> & __restrict__ sorted_hits)
{
  const double timecut = 400; // ns

  for(unsigned int i = 0; i < sorted_hits.size(); i++){
    const rb::CellHit & coin_chit = sorted_hits[i];

    if(tinfo.time - coin_chit.TNS() > +timecut) continue;
    if(tinfo.time - coin_chit.TNS() < -timecut) break;

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
static bool shall_we_flip_it(const rb::Track & trk)
{
  if(TracksAreDown) return trk.Stop().Y() > trk.Start().Y();
  else              return trk.Stop().Z() < trk.Start().Z();
}

/*
 Returns true if the track goes in the +z direction, after we
 correct the direction to be downward-going or forward-going.
*/
static bool is_increasing_z(const rb::Track & trk)
{
  if(!TracksAreDown) return true;
  const bool claims_increasing_z = trk.Stop().Z() > trk.Start().Z();
  const bool needs_flip = shall_we_flip_it(trk);
  return claims_increasing_z ^ needs_flip;
}

/*
  If the hit is not in the same view as the last hit of the track,
  instead of using the cell number of the last hit in the wrong
  view as the track's ending cell position, correct it by angle
  assuming the track goes exactly one plane further.
*/
static double cell_number_correction(const bool same_view,
                                     const int hit_in_x,
                                     const rb::Track & trk)
{
  if(same_view) return 0;

  if(trk.StopDir().Z() == 0) return 0; // Yes, this happens

  const double correction =
    planes_per_cell * (hit_in_x?trk.StopDir().X():trk.StopDir().Y())
         /trk.StopDir().Z();

  return correction;
}

// Cells are staggered in the same view from one plane to the next.
// Return the position in cell number including this offset.
static double cell_coord_off(const int plane)
{
  return (plane%2) == ((plane/2)%2)? -0.5: 0.0;
}

// There is extra material at each extrusion boundary;
// count how many we cross (typically zero or one) so that
// it can be corrected for.
static int n_extrusions_boundaries_crossed(const int cell1,
                                           const int cell2)
{
  const int cells_per_extrusion = 16; // 2 extrusions per module
  return abs(cell1/cells_per_extrusion - cell2/cells_per_extrusion);
}

/*
  Takes a track and a hit and determines the distance between the
  end of the track and the hit.  lasthiti_{even,odd} are the indices
  of the last hits in each view.
*/
static float dist_trackend_to_cell(const rb::Track & __restrict__ trk,
                                   const rb::CellHit & __restrict__ chit,
                                   const int lasthiti_even,
                                   const int lasthiti_odd)
{
  const int last_tplane_even = trk.Cell(lasthiti_even)->Plane();
  const int last_tplane_odd  = trk.Cell(lasthiti_odd) ->Plane();
  const int last_tcell_even  = trk.Cell(lasthiti_even)->Cell();
  const int last_tcell_odd   = trk.Cell(lasthiti_odd) ->Cell();

  const bool increasing_z = is_increasing_z(trk);

  const int lastplane = increasing_z?
                        std::max(last_tplane_even, last_tplane_odd)
                        :std::min(last_tplane_even, last_tplane_odd);

  const double hit_cc = chit.Cell() + cell_coord_off(chit.Plane());

  // It is really easy to get confused here.  Checked the results
  // with a number of events of different cases in the event display.

  const double track_cc =
    (chit.Plane()%2 == 0? last_tcell_even: last_tcell_odd)

    + cell_number_correction(chit.Plane()%2 == lastplane%2,
                             chit.View() == geo::kX, trk)

    + cell_coord_off(chit.Plane()%2 == 0?last_tplane_even :last_tplane_odd);

  /*                scint length  & density, PVC length & density
    standardcell = 3.5565 * 1.005 * 0.8530 + 0.368      *  1.49
    endcell      = 3.5565 * 0.964 * 0.8530 + 0.51*2     *  1.49
    endcell/standardcell - 1 = 0.235

    This only does a correction in the view that the hit is in.
    Of course, it is also possible to cross an extrusion boundary
    in the other view.  It's a small correction...
  */
  const double addition_for_extrusion_boundaries =
    0.235 * n_extrusions_boundaries_crossed(chit.Cell(),
              (chit.Plane()%2 == 0? last_tcell_even: last_tcell_odd));
  return sqrt(
    pow(planes_per_cell*(chit.Plane() - lastplane), 2) +
    pow(addition_for_extrusion_boundaries + hit_cc - track_cc, 2));
}

/*
  If the hit is near the track, return the distance in cells. Otherwise,
  return a negative number if it is too far away.
*/
static double hit_near_track(const trkinfo & __restrict__ tinfo,
  const rb::CellHit & __restrict__ chit)
{
  // Accept the hit even if it is in the track!  Because if a Michel
  // hit gets swept up into the track, this is the only way to see it.
  // We might exclude it later, but not here.
  //
  // Accept the hit even if it is before the track!  Because this gets
  // us the background level in an unbiased way.

  const double dist = dist_trackend_to_cell(tinfo.trk, chit,
    tinfo.lasthiti_even, tinfo.lasthiti_odd);

  if(dist > MaxDistInCells) return -1;

  return dist;
}

/*
  Passes back the last hit in each plane, assuming that the track is
  going the direction we have decided it is going (down or forward).
*/
static void last_hits(int & __restrict__ lasthiti_even,
                      int & __restrict__ lasthiti_odd,
                      const rb::Track & __restrict__ trk)
{
  const bool increasing_z = is_increasing_z(trk);
  const bool decreasing_z = !increasing_z;

  // ug.
  float latest_even = 1e30, latest_odd = 1e30;
  if(increasing_z){
    latest_even = -1e30; latest_odd = -1e30;
  }

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

static void print_ntuple_line(const evtinfo & __restrict__ einfo,
                              const trkinfo & __restrict__ tinfo,
                              const cluster & __restrict__ cluster,
                              const int totalclusters)
{
  const double timeleft = einfo.triggerlength - (tinfo.time - einfo.starttime);
  const double timeback =                        tinfo.time - einfo.starttime;

  const double tsx = tinfo.sx;
  const double tsy = tinfo.sy;
  const double tsz = tinfo.sz;

  const double tx = tinfo.ex;
  const double ty = tinfo.ey;
  const double tz = tinfo.ez;

  fprintf(OUT, "%d %d %d ", einfo.run, einfo.subrun, einfo.event);
  fprintf(OUT, "%f ", sqrt(pow(tinfo.sx - tinfo.ex, 2)
                         + pow(tinfo.sy - tinfo.ey, 2)
                         + pow(tinfo.sz - tinfo.ez, 2)));
  fprintf(OUT, "%d ", tinfo.i);
  fprintf(OUT, "%d ", tinfo.ntrk);
  fprintf(OUT, "%.1f ", einfo.triggerlength/1000);
  fprintf(OUT, "%f ", tinfo.time/1000);
  fprintf(OUT, "%d ", cluster.i);
  fprintf(OUT, "%d ", totalclusters);
  fprintf(OUT, "%.1f %.1f %.1f ", tsx, tsy, tsz);
  fprintf(OUT, "%.1f %.1f %.1f ", tx, ty, tz);
  fprintf(OUT, "%f ", timeleft/1000);
  fprintf(OUT, "%f ", timeback/1000);
  fprintf(OUT, "%f ", tinfo.remid);
  fprintf(OUT, "%d ", tinfo.primary_in_slice);
  fprintf(OUT, "%f ", tinfo.slice_energy);
  fprintf(OUT, "%d ", einfo.nslc);
  fprintf(OUT, "%d ", tinfo.contained_slice);
  fprintf(OUT, "%d ", tinfo.true_nupdg);
  fprintf(OUT, "%d ", tinfo.true_pdg);
  fprintf(OUT, "%d ", tinfo.true_nucc);
  fprintf(OUT, "%d ", tinfo.true_atom_cap);


  fprintf(OUT, "%d ", cluster.type);

  if(cluster.nhit)
    fprintf(OUT, "%f ",
            (float(cluster.tsum)/cluster.nhit-tinfo.time)/1000);
  else
    fprintf(OUT, "-1 ");

  if(cluster.nhit)
    fprintf(OUT, "%f ",
      (float(cluster.tsum)/cluster.nhit-cluster.previous_cluster_t)/1000);
  else
    fprintf(OUT, "0 ");

  fprintf(OUT, "%.3f ", cluster.mindist);
  fprintf(OUT, "%.3f ", cluster.maxdist);

  if(cluster.dist2sum == 0) fprintf(OUT, "0 ");
  else fprintf(OUT, "%.3f ", cluster.dist2sum/cluster.nhit);

  fprintf(OUT, "%.3f ", cluster.esum);
  fprintf(OUT, "%.3f ", cluster.pesum);
  fprintf(OUT, "%d ", cluster.adcsum);

  fprintf(OUT, "%d %d %d ", cluster.nhit, cluster.nhitx, cluster.nhity);
  fprintf(OUT, "%d %d %d ", cluster.nhitover35pe,
          cluster.nhitover35pex, cluster.nhitover35pey);
  fprintf(OUT, "%.3f ",
          (cluster.last_accepted_time - cluster.first_accepted_time)/1000);
  fprintf(OUT, "\n");
}

namespace PostMuon{

PostMuon::PostMuon(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset), fRemoveBadChans(pset.get<bool>("RemoveBadChans")),
  fRawDataLabel(pset.get< std::string >("RawDataLabel")),
  fMaxDistInCells(pset.get<float>("MaxDistInCells")),
  fTracksAreDown(pset.get<bool>("TracksAreDown")),
  fNhitTrackTimeAveraging(pset.get<int>("NhitTrackTimeAveraging"))
{
  NhitTrackTimeAveraging = fNhitTrackTimeAveraging;
  TracksAreDown = fTracksAreDown;
  MaxDistInCells = fMaxDistInCells;
}

void PostMuon::endJob()
{
  if(OUT != NULL) fclose(OUT);
  OUT = NULL;
}

// Compare tracks by their time
static bool compare_track(const trkinfo & __restrict__ a,
                          const trkinfo & __restrict__ b)
{
  return mean_late_track_time(a.trk) < mean_late_track_time(b.trk);
}

// Compare cellhits by their time.  Note that CellHit has the operator<
// defined, but that compares by plane, then by cell, then by time, which
// is not what we need here.
static bool compare_cellhit(const rb::CellHit & __restrict__ a,
                            const rb::CellHit & __restrict__ b)
{
  return a.TNS() < b.TNS();
}

static std::vector<rb::CellHit> make_trkhits(const std::vector<trkinfo> & trks)
{
  std::vector<rb::CellHit> answer;
  for(unsigned int i = 0; i < trks.size(); i++)
    for(unsigned int j = 0; j < trks[i].trk.NCell(); j++)
      answer.push_back(*(trks[i].trk.Cell(j)));
  return answer;
}

/*
all: all hits passing basic cuts

ex: Without hits in any module where the track was, in order to attempt
to get an unbiased, sag-free measurement

ex2: Same, but also excluding any module that had a hit coincident with
the track, where anything within 0.4us counts as coincident. This is 2-3
sigma in the time resolution and should catch almost all hits that are
prompt, either muon hits that didn't get reconstructed as part of the
track, or brems, or x-rays from muon atomic capture, or whatever else.

xt: Same as ex, but without hits that are part of any tracks, in an attempt to
beat down uncorrelated background (doesn't seem to have much effect --
maybe 5-10%.)
*/
enum clustertype { all, ex, ex2, xt };

struct printinfo{
  evtinfo einfo;
  trkinfo tinfo;
  cluster clu;
};

static printinfo make_printinfo(const evtinfo & __restrict__ einfo_,
                                const trkinfo & __restrict__ tinfo_,
                                const cluster & __restrict__ clu_)
{
  printinfo ans;
  ans.einfo = einfo_;
  ans.tinfo = tinfo_;
  ans.clu= clu_;
  return ans;
}

static void cluster_search(const int type,
  const evtinfo & __restrict__ einfo,
  const std::vector<rb::CellHit> & __restrict__ sorted_hits,
  const std::vector<rb::CellHit> & __restrict__ trkhits,
  const trkinfo & __restrict__ tinfo)
{
  art::ServiceHandle<calib::Calibrator> calthing;

  cluster clu = mkcluster();
  clu.i = 0;
  clu.type = type;
  clu.previous_cluster_t = tinfo.time;

  std::vector<printinfo> toprint;

  for(unsigned int c = 0; c < sorted_hits.size(); c++){
    const rb::CellHit & chit = sorted_hits[c];

    const double dist = hit_near_track(tinfo, chit);

    if(dist < 0) continue;

    const rb::RecoHit rhit = calthing->MakeRecoHit(chit,
       // If the hit is in X, it needs a Y plane to provide W
       chit.View() == geo::kX? tinfo.ey: tinfo.ex);

    if(!rhit.IsCalibrated()
       ||
       (type == ex &&hit_in_track_module(chit, tinfo.trk))
       ||
       (type == ex2&&hit_in_track_coincident_module(chit, tinfo, sorted_hits))
       ||
       (type == xt &&hit_in_any_track(chit, trkhits))
      )
      continue;

    // Hits are time ordered.  Only report on hits that are separated in time
    // from other accepted hits by at least 500ns, which is a few TDC ticks
    // (1/16us) and about twice the gaussian resolution of fine timing for very
    // small hits. It seems like the right amount of time by inspection.
    const bool newcluster = chit.TNS() > clu.last_accepted_time + 500.0;

    if(newcluster && clu.nhit){
      toprint.push_back(make_printinfo(einfo, tinfo, clu));
      clu.i++;
      const float thistime = float(clu.tsum)/clu.nhit;
      resetcluster(clu);
      clu.previous_cluster_t = thistime;
    }

    // Add everything to the cluster *after* the print of the
    // previous cluster (or no-op).

    clu.nhit++;
    if(chit.View() == geo::kX) clu.nhitx++;
    else                       clu.nhity++;
    if(chit.PE() > 35.0){
      clu.nhitover35pe++;
      if(chit.View() == geo::kX) clu.nhitover35pex++;
      else                       clu.nhitover35pey++;
    }

    clu.dist2sum += dist*dist;
    if(dist < clu.mindist) clu.mindist = dist;
    if(dist > clu.maxdist) clu.maxdist = dist;

    clu.last_accepted_time = chit.TNS(); // because they're sorted by TNS
    if(clu.first_accepted_time == INVALID_TIME)
      clu.first_accepted_time = chit.TNS();
    clu.tsum += chit.TNS();

    clu.esum += rhit.GeV()*1000;
    clu.pesum += chit.PE();
    clu.adcsum += chit.ADC();
  }

  // print last cluster, or track info if no cluster
  toprint.push_back(make_printinfo(einfo, tinfo, clu));
  for(unsigned int i = 0; i < toprint.size(); i++)
    print_ntuple_line(toprint[i].einfo,
                      toprint[i].tinfo,
                      toprint[i].clu,
                      toprint.size());
}

PostMuon::~PostMuon() { }

/*
  Get the length of the event in TDC ticks, typically 550*64, and
  "delta_tdc", the time between the trigger time and the time the event
  starts. You can subtract this off of the time that the offline gives
  each hit to get the time since the beginning of the readout, and with
  the event length, the time until the end of the readout.

  delta_tdc is a signed 64 bit integer, even though it should always be
  a small positive number, just in case. Ditto for the event length.

  Returns whether this information was successfully extracted.
*/
static bool delta_and_length(int64_t & event_length_tdc,
  int64_t & delta_tdc,
  const art::Handle< std::vector<rawdata::FlatDAQData> > & flatdaq,
  const art::Handle< std::vector<rawdata::RawTrigger> > & rawtrigger)
{
  daqdataformats::RawEvent raw;
  if(flatdaq->empty()) return false;

  raw.readData((*flatdaq)[0].getRawBufferPointer());
  if(raw.getDataBlockNumber() == 0) return false;

  raw.setFloatingDataBlock(0);
  daqdataformats::RawDataBlock& datablock = *raw.getFloatingDataBlock();

  uint64_t event_start_time = 0xffffffffffffffff;
  uint64_t event_end_time   = 0x0000000000000000;

  for(unsigned int di = 0; di < raw.getDataBlockNumber(); di++){
    raw.setFloatingDataBlock(di);
    datablock = (*raw.getFloatingDataBlock());

    if(datablock.getHeader()->getMarker() ==
         daqdataformats::datablockheader::SummaryBlock_Marker ||
       !datablock.getHeader()->checkMarker()) continue;

    for(unsigned int mi = 0; mi < datablock.getNumMicroBlocks(); mi++){
      datablock.setFloatingMicroBlock(mi);
      daqdataformats::RawMicroBlock * ub = datablock.getFloatingMicroBlock();

      // The time is always in the second and third words of the
      // microslice, which follows two words of microblock header, so
      // just get it. Justin says you can also get it from getTime(),
      // but this already works and I'm not touching it.
      const uint32_t t_marker_low  = ((uint32_t *)(ub->getBuffer()))[3];
      const uint32_t t_marker_high = ((uint32_t *)(ub->getBuffer()))[4];

      uint64_t time_marker = t_marker_low;
      time_marker |= (uint64_t)t_marker_high << 32;
      if(time_marker < event_start_time) event_start_time = time_marker;
      if(time_marker > event_end_time  ) event_end_time   = time_marker;
    }
  }

  delta_tdc = (int64_t)((*rawtrigger)[0].fTriggerTimingMarker_TimeStart
                        - event_start_time);

  // Assume that microblocks are always 50us. I hope that's true for all
  // relevant data.
  event_length_tdc = ((int64_t)(event_end_time - event_start_time))
                     + US_PER_MICROSLICE*TDC_PER_US;
  return true; // ok
}

static void ntuple_header(const art::Event & evt)
{
  if(OUT == NULL){
    OUT = fopen(Form("postmuon_%d_%d.20170308.ntuple",
                     evt.run(), evt.subRun()), "w");
    if(OUT == NULL){
       fprintf(stderr, "Could not open output ntuple file\n");
       exit(1);
    }
    fprintf(OUT,
      "run/I:"
      "subrun/I:"
      "event/I:"
      "trklen/F:"
      "trk/I:"
      "ntrk/I:"
      "triggerlength/F:"
      "trktime/F:"
      "i/I:"
      "nclu/I:"
      "trkstartx/F:"
      "trkstarty/F:"
      "trkstartz/F:"
      "trkx/F:"
      "trky/F:"
      "trkz/F:"
      "timeleft/F:"
      "timeback/F:"
      "remid/F:"
      "primary/I:"
      "slce/F:"
      "nslc/I:"
      "contained/I:"
      "true_nupdg/I:"
      "true_pdg/I:"
      "true_nucc/I:"
      "true_atom_cap/I:"

      "type/I:"
      "t/F:"
      "dt/F:"
      "mindist/F:"
      "maxdist/F:"
      "dist2/F:"
      "e/F:"
      "pe/F:"
      "adc/I:"
      "nhit/I:"
      "nhitx/I:"
      "nhity/I:"
      "nhitover35pe/I:"
      "nhitover35pex/I:"
      "nhitover35pey/I:"
      "ctlen/F"
      "\n");
  }
}

// return the index in the slice array that the given track is in
static int which_slice_is_this_track_in(
  const trkinfo & t,
  const art::Handle< std::vector<rb::Cluster> > & slice)
{
  // I'm sure this is not the best way, but I have had it with trying to
  // figure out what the best way is, and I'm just going to do it *some*
  // way.
  const art::Ptr<rb::CellHit> ahit =
    t.trk.Cell(0); // some random hit on the track

  // Could probably skip slice 0 since it is the noise slice, but let's not
  // in case that convention changes.
  for(unsigned int i = 0; i < slice->size(); i++){
    const rb::Cluster & slc = (*slice)[i];
    for(unsigned int j = 0; j < slc.NCell(); j++){
      const art::Ptr<rb::CellHit> shit = slc.Cell(j);
      if(*ahit == *shit) return i;
    }
  }
  return -1;
}

// Fill in for each track whether it is the one with the highest remid
// in the slice
static void fill_primary_track_info(std::vector<trkinfo> & ts, const int nslc)
{
  // For each slice...
  for(int i = 0; i < nslc; i++){
    double maxremid = -1e40;
    unsigned int best = 0;
    // ... find the best track
    for(unsigned int j = 0; j < ts.size(); j++){
      if(ts[j].slice == i && ts[j].remid > maxremid){
        maxremid = ts[j].remid;
        best = j;
      }
    }

    // ... and label it as such and the rest as not such
    for(unsigned int j = 0; j < ts.size(); j++)
      if(ts[j].slice == i)
        ts[j].primary_in_slice = (best == j);
  }
}

// Return true iff this slice is numu-contained
static bool containedND(const art::Ptr<caf::StandardRecord> sr)
{
  // Lifted from CAFAna/Cuts/NumuCuts.h, development 2017-03-28 I've
  // removed the trk.kalman.tracks lines because they make it seg
  // fault, even if I check that the vector is not empty, why?! But
  // since I output the start and stop position anyway, and don't allow
  // tracks that include the muon catcher, track position checks are
  // not really needed. Ditto for sr->energy.numu.ndhadcalcatE and
  // ndhadcaltranE, which I found to be always nan. Alex R says "I think
  // the calibration in development is still the one that's broken for
  // the nd muon catcher." I don't need them since I disallow the muon
  // catcher anyway.
  const bool a = sr->trk.kalman.ntracks > sr->trk.kalman.idxremid,
             b = sr->slc.ncellsfromedge > 1,
             c = sr->slc.firstplane > 1,   // skip 0 and 1
             d = sr->slc.lastplane  < 212, // skip 212 and 213
             e = sr->sel.contain.kalfwdcellnd > 4,
             f = sr->sel.contain.kalbakcellnd > 8;

  return a && b && c && d && e && f;
}

void PostMuon::analyze(const art::Event& evt)
{
  // I'm so sorry that I have to do this.  And, my goodness, doing
  // it in the constructor isn't sufficient.  If this isn't done,
  // it responds to PIPE by going into an endless loop.
  signal(SIGPIPE, SIG_DFL);

  art::Handle< std::vector<rawdata::FlatDAQData> > flatdaq;
  evt.getByLabel("daq", flatdaq);

  art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
  evt.getByLabel("daq", rawtrigger);

  art::Handle< std::vector<rb::CellHit> > cellcol;
  evt.getByLabel(fRawDataLabel, cellcol);

  art::Handle< std::vector<rb::Cluster> > slice;
  evt.getByLabel("slicer", slice);

  art::Handle< std::vector<rb::Track> > tracks;
  evt.getByLabel("kalmantrackmerge", tracks);

  // works better for cosmics (bizzarely), but has no remid
  //evt.getByLabel("windowtrack", tracks);

  art::FindOneP<remid::ReMId> track2remid(tracks, evt, "remid");
  if(!track2remid.isValid()) { puts("No track2remid"); return; }

  art::FindOneP<numue::NumuE> slice2numue(slice, evt, "numue");
  if(!slice2numue.isValid()) { puts("No slice2numue"); return; }

  art::FindOneP<caf::StandardRecord> slice2caf(slice, evt, "cafmaker");
  if(!slice2caf.isValid()) { puts("No slice2caf"); return; }

  ntuple_header(evt);

  if(rawtrigger->empty()) return;
  if(tracks    ->empty()) return;

  art::ServiceHandle<cheat::BackTracker> backtracker_thing;
  const sim::ParticleNavigator& pnav = backtracker_thing->ParticleNavigator();

  // Is this an ok way to test for data vs. MC?
  const bool is_data = !backtracker_thing->HaveTruthInfo();

  int64_t event_length_tdc = 0, delta_tdc = 0;
  if(is_data){
    art::Handle< std::vector<rawdata::FlatDAQData> > flatdaq;
    evt.getByLabel("daq", flatdaq);

    art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
    evt.getByLabel("daq", rawtrigger);

    if(rawtrigger->empty()) return;

    if(!delta_and_length(event_length_tdc, delta_tdc, flatdaq, rawtrigger))
      return;
  }
  else{ // XXX
    event_length_tdc = 500 * 64;
    delta_tdc = 224 * 64;
  }

  evtinfo einfo;
  einfo.triggerlength = event_length_tdc * 1000. / TDC_PER_US;
  einfo.starttime = -(delta_tdc * 1000. / TDC_PER_US);
  einfo.run = evt.run();
  einfo.subrun = evt.subRun();
  einfo.event = evt.event();
  einfo.nslc = slice->size();

  // CellHits do *not* come in time order
  std::vector<rb::CellHit> sorted_hits;
  for(int c = 0; c < (int)cellcol->size(); c++)
    sorted_hits.push_back((*cellcol)[c]);
  std::sort(sorted_hits.begin(), sorted_hits.end(), compare_cellhit);

  std::vector<trkinfo> sorted_tracks;

  for(unsigned int c = 0; c < tracks->size(); c++){
    // Badly tracked tracks are not useful
    if((*tracks)[c].Stop().X() == 0 || (*tracks)[c].Stop().Y() == 0) continue;

    trkinfo t;
    t.trk = (*tracks)[c];
    t.remid = track2remid.at(c)->Value();

    if(0 > (t.slice = which_slice_is_this_track_in(t, slice))) return;

    t.slice_energy = slice2numue.at(t.slice)->E();

    t.contained_slice = containedND(slice2caf.at(t.slice));

    t.true_pdg = t.true_nupdg = t.true_nucc = t.true_atom_cap = 0;
    if(!is_data){
      // Horrible. I cannot figure out how to mash what I have into
      // any of the acceptable types for SliceToNeutrinoInteractions
      // in a cheap way, so I am just constructing my own goddam
      // vector of CellHits. I think this is very slow...
      std::vector<rb::CellHit> this_slc;
      for(unsigned int h = 0; h < (*slice)[t.slice].NCell(); h++)
        this_slc.push_back(*((*slice)[t.slice].Cell(h)));

      std::vector<cheat::NeutrinoEffPur> truths =
        backtracker_thing->SliceToNeutrinoInteractions(this_slc, sorted_hits);
      if(!truths.empty()){
        t.true_nupdg = truths[0].neutrinoInt->GetNeutrino().Nu().PdgCode();
        t.true_nucc  = !truths[0].neutrinoInt->GetNeutrino().CCNC();
      }

      // Ditto above complaint. If there's a way to directly pass a
      // track in, I don't see it, so I have to copy the whole track
      // into a vector. I'm supposed to at least be able to pass in
      // a vector of pointers, except when I try to construct such a
      // thing, I get horribly tangled up in art::Ptrs when real pointers
      // are wanted and conflicting sets of where "const" goes, and I
      // give up. I have other work to do.
      std::vector<rb::CellHit> trackhits;
      for(unsigned int i = 0; i < (*tracks)[c].NCell(); i++)
        trackhits.push_back(*(*tracks)[c].Cell(i));

      const std::vector<const sim::Particle *> particles = backtracker_thing
        ->HitsToParticle(trackhits);

      // Can be empty if all the hits are noise
      if(!particles.empty()){
        t.true_pdg = particles[0]->PdgCode();

        // As of 2017-05-03, sim::Particle::EndE() seems to always
        // return the particle's mass, no matter what happens to it.
        // sim::Particle::EndProcess() is also not useful for finding
        // out whether a particle has stopped, since it always holds an
        // empty string. So to find out what happens to a particle, you
        // have to laboriously look at its daughters...
        for(int d = 0; d < particles[0]->NumberDaughters(); d++){
          sim::ParticleNavigator::const_iterator it =
            pnav.find(particles[0]->Daughter(d));

          if(it == pnav.end()) continue;

          const std::string dproc = it->second->Process();
          if(dproc == "muMinusCaptureAtRest" ||
             dproc == "hBertiniCaptureAtRest"){
            t.true_atom_cap = 1;
            break;
          }
        }
      }
    }

    sorted_tracks.push_back(t);
  }

  fill_primary_track_info(sorted_tracks, slice->size());

  std::sort(sorted_tracks.begin(), sorted_tracks.end(), compare_track);

  std::vector<rb::CellHit> trkhits = make_trkhits(sorted_tracks);

  for(unsigned int t = 0; t < sorted_tracks.size(); t++){
    trkinfo & tinfo = sorted_tracks[t];
    rb::Track & trk = tinfo.trk;

    const bool needs_flip = shall_we_flip_it(trk);
    tinfo.sx = needs_flip? trk.Stop ().X(): trk.Start().X();
    tinfo.sy = needs_flip? trk.Stop ().Y(): trk.Start().Y();
    tinfo.sz = needs_flip? trk.Stop ().Z(): trk.Start().Z();
    tinfo.ex = needs_flip? trk.Start().X(): trk.Stop ().X();
    tinfo.ey = needs_flip? trk.Start().Y(): trk.Stop ().Y();
    tinfo.ez = needs_flip? trk.Start().Z(): trk.Stop ().Z();

    // It's a waste of time to search for clusters around tracks that almost
    // certainly represent exiters.  Note that this does almost nothing to near
    // detector events
    if(fabs(tinfo.ex) >  750) continue;
    if(     tinfo.ey  < -730) continue;
    if(     tinfo.ez  <   25 || tinfo.ez  > 5935) continue;

    // Don't exclude tracks based on their starting positions, since I want this
    // to work for both cosmics and beam
    ;

    tinfo.trk = trk;
    tinfo.i = t;
    tinfo.ntrk = sorted_tracks.size();
    tinfo.remid = sorted_tracks[t].remid;
    tinfo.lasthiti_even = 0, tinfo.lasthiti_odd = 0;
    last_hits(tinfo.lasthiti_even, tinfo.lasthiti_odd, trk);
    tinfo.time = mean_late_track_time(trk);

    cluster_search(all, einfo, sorted_hits, trkhits, tinfo);
    cluster_search(ex,  einfo, sorted_hits, trkhits, tinfo);
    cluster_search(ex2, einfo, sorted_hits, trkhits, tinfo);
    cluster_search(xt,  einfo, sorted_hits, trkhits, tinfo);
  }
}

DEFINE_ART_MODULE(PostMuon);

} // end namespace PostMuon
//////////////////////////////////////////////////////////////////////////
