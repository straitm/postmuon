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
#include "Utilities/AssociationUtil.h"

#include "Metadata/MetadataManager.h"

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

// For MC cross-section weights
#include "StandardRecord/Proxy/SRProxy.h"
#include "StandardRecord/Proxy/CopyRecord.h"
#include "CAFAna/Vars/GenieWeights.h"
#include "CAFAna/Vars/XsecTunes.h"

#include "GeometryObjects/PlaneGeo.h"

#include <IFDH_service.h>
#include "SummaryData/SpillData.h"

#include <string>
#include <set>
#include <algorithm>

#include <signal.h>

#include "TRandom3.h"
static TRandom3 rand_stamp(0);

static FILE * OUT = NULL;
static int NhitTrackTimeAveraging = 1; // to be overwritten in PostMuon()
static bool TracksAreDown = true; // to be overwritten in PostMuon()
static double MaxDistInCells = 0.123456; // dummy as well


struct evtinfo{
  int run;
  int subrun;
  int cycle; // for ND MC - each cycle in a subrun reuses the same rock events
  int event;
  int nslc; // number of slices
  float pot; // Number of 1E12 pot
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
  int true_nuint; // Interaction type: 0 = QE, 1 = Res... , -1 otherwise
  int true_atom_cap; // 1 if this is a mu-, pi- or K- that stops, else 0
  int true_neutrons; // Number of true neutron daughters
  bool primary_in_slice; // Is this the highest remid track in the slice?
  double time;
  int last_plane_even; // plane number -- need to cache these for offspace
  int last_plane_odd; // plane number  -- sample because the cells are no
  int last_cell_even; // cell number   -- longer accessible once we move
  int last_cell_odd; // cell number    -- to the next spill.
  float sx, sy, sz, ex, ey, ez; // start and end position, after flip correction
  float end_dx, end_dy, end_dz; // ending direction
  double remid;
  double cvn;
  float mcweight;
  float dother; // Distance in cells to nearest other slice
};

// A light track that is just used to determine if we've seen this exact
// same track before (because of DAQ problems).  It's enough to be sure
// of uniqueness.
struct unique_track{
  unique_track(const trkinfo & in)
  {
    sx = in.sx;
    sy = in.sy;
    sz = in.sz;
    ex = in.ex;
    ey = in.ey;
    ez = in.ez;
    end_dx = in.end_dx;
    end_dy = in.end_dy;
    end_dz = in.end_dz;
  }

  bool operator <(const unique_track & b) const
  {
    if(sx != b.sx) return sx < b.sx;
    if(sy != b.sy) return sy < b.sy;
    if(sz != b.sz) return sz < b.sz;
    if(ex != b.ex) return ex < b.ex;
    if(ey != b.ey) return ey < b.ey;
    if(ez != b.ez) return ez < b.ez;
    if(end_dx != b.end_dx) return end_dx < b.end_dx;
    if(end_dy != b.end_dy) return end_dy < b.end_dy;
    return end_dz < b.end_dz;
  }

  float sx, sy, sz, ex, ey, ez; // start and end position, after flip correction
  float end_dx, end_dy, end_dz; // ending direction
};

struct cluster{
  int i; // which cluster for the given track
  int type; // defines the set of cuts used for this cluster
  int offspace; // 1 if this is an off-space background sample
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

  // In units of cell widths:
  float x_dz; // In x view, mean difference in z from track end to hit
  float dx;   // mean difference in x in x view
  float y_dz; // Ditto y view
  float dy;   // Ditto y view

  float cosx; // cosine of the mean direction from track end
  float cosy; // wrt track dir in the x and y views, respectively.
              // Definition is a little janky because it uses my
              // plane/cell distance defintions for the cluster, which
              // has some density weighting, but rb::Track's direction
              // for the track, which is in plain old space.
};

// Circular buffers of tracks from previous spills which can be used to
// provide unbiased off-space samples from which to measure pileup.
const int MAXSAVEDTRACKTYPE = 2*2*2;
static std::deque<trkinfo> offspacetrk[MAXSAVEDTRACKTYPE];
const unsigned int savedtrack_size = 1024;

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
  res.x_dz = 0;
  res.y_dz = 0;
  res.dx = 0;
  res.dy = 0;
  res.cosx = -2;
  res.cosy = -2;
}

static cluster mkcluster()
{
  cluster res;
  resetcluster(res);
  res.i = 0;
  res.type = 0;
  res.offspace = 0;
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

/* Returns true if this cell hit is in this track. */
static bool hit_in_track(const rb::CellHit & chit,
                         const rb::Track & trk)
{
  for(unsigned int i = 0; i < trk.NCell(); i++){
    const rb::CellHit & trk_chit = *(trk.Cell(i));

    if(trk_chit == chit) return true;
  }
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
  // Same test for ND and FD.
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
static float dist_trackend_to_cell(const trkinfo & __restrict__ tinfo,
                                   const rb::CellHit & __restrict__ chit,
                                   float & __restrict__ dplane,
                                   float & __restrict__ dcell)
{
  const bool increasing_z = is_increasing_z(tinfo.trk);

  const int lastplane = increasing_z?
                        std::max(tinfo.last_plane_even, tinfo.last_plane_odd)
                        :std::min(tinfo.last_plane_even, tinfo.last_plane_odd);

  const double hit_cc = chit.Cell() + cell_coord_off(chit.Plane());

  // It is really easy to get confused here.  Checked the results
  // with a number of events of different cases in the event display.

  const double track_cc =
    (chit.Plane()%2 == 0? tinfo.last_cell_even: tinfo.last_cell_odd)

    + cell_number_correction(chit.Plane()%2 == lastplane%2,
                             chit.View() == geo::kX, tinfo.trk)

    + cell_coord_off(chit.Plane()%2 == 0?tinfo.last_plane_even:
                                         tinfo.last_plane_odd);

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
              (chit.Plane()%2 == 0? tinfo.last_cell_even: tinfo.last_cell_odd));

  dplane = planes_per_cell*(chit.Plane() - lastplane);
  dcell = addition_for_extrusion_boundaries + hit_cc - track_cc;

  return sqrt(pow(dplane, 2) + pow(dcell, 2));
}

static float make_dother(const trkinfo & __restrict__ tinfo,
  const art::Handle< std::vector<rb::Cluster> > & __restrict__ slice)
{
  float mindist = 10000;
  // start at 1 so as not to get the noise slice
  for(unsigned int i = 1; i < slice->size(); i++){
    // Don't search in the track's own slice
    if(tinfo.slice == (int)i) continue;

    for(unsigned int j = 0; j < (*slice)[i].NCell(); j++){
      float dplane, dcell; // unused
      const float dist = dist_trackend_to_cell(tinfo, *(*slice)[i].Cell(j),
        dplane, dcell);
      if(dist < mindist) mindist = dist;
      if(dist == 0) return 0; // optimization
    }
  }
  return mindist;
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

  // Optimization by quickly rejecting hits by plane alone. Measured to
  // substantially reduce time for FD SNEWS triggers (the worst case)
  if( fabs(chit.Plane()-tinfo.last_plane_even)*planes_per_cell > MaxDistInCells
   && fabs(chit.Plane()-tinfo.last_plane_odd )*planes_per_cell > MaxDistInCells)
    return -1;

  float dplane, dcell; /* unused here */

  const double dist = dist_trackend_to_cell(tinfo, chit,
    /**/ dplane, dcell /**/ );

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

  fprintf(OUT, "%d %d %d %d ", einfo.run,einfo.subrun,einfo.cycle,einfo.event);
  fprintf(OUT, "%f ", einfo.pot);
  fprintf(OUT, "%f ", tinfo.trk.TotalLength());
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
  fprintf(OUT, "%f ", tinfo.cvn);
  fprintf(OUT, "%d ", tinfo.primary_in_slice);
  fprintf(OUT, "%f ", tinfo.slice_energy);
  fprintf(OUT, "%d ", einfo.nslc);
  fprintf(OUT, "%d ", tinfo.contained_slice);
  fprintf(OUT, "%d ", tinfo.true_nupdg);
  fprintf(OUT, "%d ", tinfo.true_pdg);
  fprintf(OUT, "%d ", tinfo.true_nucc);
  fprintf(OUT, "%d ", tinfo.true_nuint);
  fprintf(OUT, "%d ", tinfo.true_atom_cap);
  fprintf(OUT, "%d ", tinfo.true_neutrons);
  fprintf(OUT, "%f %f ", cluster.cosx, cluster.cosy);
  fprintf(OUT, "%f ", tinfo.mcweight);
  fprintf(OUT, "%f ", tinfo.dother);

  fprintf(OUT, "%d ", cluster.type);

  const float tim = (float(cluster.tsum)/cluster.nhit-tinfo.time)/1000;

  if(cluster.nhit)
    fprintf(OUT, "%f ", tim);
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

xt: Without hits that are part of any tracks, in an attempt to
beat down uncorrelated background (doesn't seem to have much effect --
maybe 5-10%.)

x: Without any hits from the track.
*/
enum clustertype { all, ex, ex2, xt, x };

struct printinfo{
  evtinfo einfo;
  trkinfo tinfo;
  cluster clu;
};

static float make_cosine(const int nhit, const float tend_dz,
                         const float tend_dw, float dz, float dw)
{
  if(nhit == 0) return -2;

  dz /= nhit;
  dw /= nhit;

  // clu_mag can be zero if the hits are all exactly at the track end.
  const double clu_mag = sqrt(pow(dz, 2) + pow(dw, 2));
  if(clu_mag == 0) return -2;

  // The track direction is a unit vector in 3D, but we need the
  // magnitude in each view. track_mag should probably not be zero, but
  // probably is for extremely short tracks or something.
  const double track_mag = sqrt(pow(tend_dw, 2) + pow(tend_dz, 2));
  if(track_mag == 0) return -2;

  return (dz * tend_dz + dw * tend_dw)/track_mag/clu_mag;
}

// Given a cluster that has undivided x_dz, etc, produce the cosx and
// cosy variables.
static void make_cosines(cluster & __restrict__ clu,
                         const trkinfo & __restrict__ tinfo)
{
  clu.cosx = make_cosine(clu.nhitx, tinfo.end_dz,tinfo.end_dx, clu.x_dz,clu.dx);
  clu.cosy = make_cosine(clu.nhity, tinfo.end_dz,tinfo.end_dy, clu.y_dz,clu.dy);
}

static printinfo make_printinfo(const evtinfo & __restrict__ einfo_,
                                const trkinfo & __restrict__ tinfo_,
                                      cluster & __restrict__ clu_)
{
  printinfo ans;

  make_cosines(clu_, tinfo_);

  ans.einfo = einfo_;
  ans.tinfo = tinfo_;
  ans.clu= clu_;
  return ans;
}

static int which_saved_track_array(const trkinfo & real_tinfo)
{
  return real_tinfo.primary_in_slice * 0x1
       + real_tinfo.contained_slice  * 0x2
       + (real_tinfo.ez > 1300)      * 0x4;
}

// Return a saved track of the same sort that this track is.
// if no tracks of that sort have been saved yet, return NULL.
static trkinfo * select_offspace_track(const trkinfo & real_tinfo)
{
  const int which = which_saved_track_array(real_tinfo);

  if(offspacetrk[which].empty()) return NULL;

  // Is using my own TRandom3 the right thing to do?
  const int i = rand_stamp.Integer(offspacetrk[which].size());

  return & offspacetrk[which][i];
}

static void cluster_search(const int type,
  const bool offspace,
  const evtinfo & __restrict__ einfo,
  const std::vector<rb::CellHit> & __restrict__ sorted_hits,
  const std::vector<rb::CellHit> & __restrict__ trkhits,
  const trkinfo & __restrict__ real_tinfo)
{
  art::ServiceHandle<calib::Calibrator> calthing;

  cluster clu = mkcluster();
  clu.i = 0;

  // So, e.g. type 3 is xt and type 13 is the xt off-space sample
  // Avoids accidentially selecting both together.
  clu.type = type + 10*offspace;

  const trkinfo * tinfo = &real_tinfo;
  if(offspace) tinfo = select_offspace_track(real_tinfo);
  if(tinfo == NULL) return;

  clu.previous_cluster_t = tinfo->time;

  std::vector<printinfo> toprint;

  for(unsigned int c = 0; c < sorted_hits.size(); c++){
    const rb::CellHit & chit = sorted_hits[c];

    const double dist = hit_near_track(*tinfo, chit);

    if(dist < 0) continue;

    // Don't try to optimize this out of the loop.  We have to recalibrate
    // for each track end point.
    const rb::RecoHit rhit = calthing->MakeRecoHit(chit,
       // If the hit is in X, it needs a Y plane to provide W
       chit.View() == geo::kX? tinfo->ey: tinfo->ex);

    if(!rhit.IsCalibrated()
       ||
       (type == ex &&hit_in_track_module(chit, tinfo->trk))
       ||
       (type == ex2&&hit_in_track_coincident_module(chit, *tinfo, sorted_hits))
       ||
       (type == xt &&hit_in_any_track(chit, trkhits))
       ||
       (type == x  &&hit_in_track(chit, tinfo->trk))
      )
      continue;

    // Hits are time ordered.  Only report on hits that are separated in time
    // from other accepted hits by at least 500ns, which is a few TDC ticks
    // (1/16us) and about twice the gaussian resolution of fine timing for very
    // small hits. It seems like the right amount of time by inspection.
    const bool newcluster = chit.TNS() > clu.last_accepted_time + 500.0;

    if(newcluster && clu.nhit){
      toprint.push_back(make_printinfo(einfo, *tinfo, clu));
      clu.i++;
      const float thistime = float(clu.tsum)/clu.nhit;
      resetcluster(clu);
      clu.previous_cluster_t = thistime;
    }

    // Add everything to the cluster *after* the print of the
    // previous cluster (or no-op).

    float dplane, dcell;

    dist_trackend_to_cell(*tinfo, chit, dplane, dcell);

    clu.nhit++;
    if(chit.View() == geo::kX){
      clu.nhitx++;
      clu.x_dz += dplane;
      clu.dx   += dcell;
    }
    else{
      clu.nhity++;
      clu.y_dz += dplane;
      clu.dy   += dcell;
    }

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

  // Print last cluster, or track info if no cluster. But only print clusters
  // for off-space, since we already have track info from the on-space sample.
  // Affects tag in stagehotelcurlknotgoat in analysis.
  if(clu.nhit > 0 || !offspace)
    toprint.push_back(make_printinfo(einfo, *tinfo, clu));
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

static void ntuple_header(const evtinfo & einfo)
{
  if(OUT == NULL){
    OUT = fopen(Form("postmuon_%d_%d_%d.ntuple",
                     einfo.run, einfo.subrun, einfo.cycle), "w");
    if(OUT == NULL){
       fprintf(stderr, "Could not open output ntuple file\n");
       exit(1);
    }
    fprintf(OUT,
      "run/I:"
      "subrun/I:"
      "cycle/I:"
      "event/I:"
      "pot/F:"
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
      "cvn/F:"
      "primary/I:"
      "slce/F:"
      "nslc/I:"
      "contained/I:"
      "true_nupdg/I:"
      "true_pdg/I:"
      "true_nucc/I:"
      "true_nuint/I:"
      "true_atom_cap/I:"
      "true_neutrons/I:"
      "cosx/F:"
      "cosy/F:"
      "mcweight/F:"
      "dother/F:"

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

  // Skip slice 0 since it is the noise slice
  for(unsigned int i = 1; i < slice->size(); i++){
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
  for(int i = 1 /* sic */; i < nslc; i++){
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


// Return true if the spill was good
static bool goodspill(const art::Ptr<caf::StandardRecord> sr)
{
  // Only apply to NuMI triggers
  if(sr->spill.trigger != 0) return true;

  // From CAFAna/Cuts/SpillCuts.h, except that I do *not* want to
  // drop low POT spills.
  return fabs(sr->spill.deltaspilltimensec) < 0.5e9 &&
         sr->spill.hornI >= -202 && sr->spill.hornI <= -198 &&
         sr->spill.posx >= -2.00 && sr->spill.posx <= +2.00 &&
         sr->spill.posy >= -2.00 && sr->spill.posy <= +2.00 &&
         sr->spill.widthx >= 0.57 && sr->spill.widthx <= 1.58 &&
         sr->spill.widthy >= 0.57 && sr->spill.widthy <= 1.58;
}


// Return true iff this slice is numu-contained
static bool containedND(const art::Ptr<caf::StandardRecord> sr)
{
  // Lifted from CAFAna/Cuts/NumuCuts2017.h:kNumuContainND2017 on 2017-12-14.
  if(sr->vtx.nelastic < 1) return false;
  // reconstructed showers all contained
  for(unsigned int i = 0; i < sr->vtx.elastic[0].fuzzyk.nshwlid; i++){
    TVector3 start = sr->vtx.elastic[0].fuzzyk.png[i].shwlid.start;
    TVector3 stop  = sr->vtx.elastic[0].fuzzyk.png[i].shwlid.stop;
    if(std::min(start.X(), stop.X()) < -180.0) return false;
    if(std::max(start.X(), stop.X()) >  180.0) return false;
    if(std::min(start.Y(), stop.Y()) < -180.0) return false;
    if(std::max(start.Y(), stop.Y()) >  180.0) return false;
    if(std::min(start.Z(), stop.Z()) <   20.0) return false;
    if(std::max(start.Z(), stop.Z()) > 1525.0) return false;
  }

  // only primary muon track present in muon catcher
  if(sr->trk.kalman.ntracks < 1) return false;
  for(unsigned int i = 0; i < sr->trk.kalman.ntracks; i++){
    if(i == sr->trk.kalman.idxremid ) continue;
    else if(sr->trk.kalman.tracks[i].start.Z() > 1275 ||
            sr->trk.kalman.tracks[i].stop.Z()  > 1275)
      return false;
  }

  return sr->trk.kalman.ntracks > sr->trk.kalman.idxremid
         && sr->slc.firstplane > 1   // skip 0 and 1
         && sr->slc.lastplane  < 212 // skip 212 and 213
         && sr->trk.kalman.tracks[0].start.Z() < 1100
         // vertex definitely outside mC
         && ( sr->trk.kalman.tracks[0].stop.Z() < 1275
              || sr->sel.contain.kalyposattrans < 55 ) // air gap
         && sr->sel.contain.kalfwdcellnd > 5
         && sr->sel.contain.kalbakcellnd > 10;
}


static float getpot(const art::Event & evt)
{
  // Try to get POT spill info if this is real data
  art::Handle<sumdata::SpillData> spillPot;

  if(!evt.isRealData()) return 0;

  evt.getByLabel("ifdbspillinfo", spillPot);

  if(spillPot.failedToGet()) return -1;

  return spillPot->spillpot;
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
  if(!evt.getByLabel("kalmantrackmerge", tracks)){
    fprintf(stderr, "Could not get kalmantrackmerge\n");
    return;
  }

  // works better for cosmics (bizzarely), and also at connecting the last
  // bend of a pion track, but has no remid
  //
  // Could get both kalman and window, use the window end point, and then
  // find the kalman track that shared the most hits and use that remid
  // score.  Errigiifhfhhhghghgggg.
  //evt.getByLabel("windowtrack", tracks);

  art::FindOneP<remid::ReMId> track2remid(tracks, evt, "remid");
  if(!track2remid.isValid()) { fputs("No track2remid\n", stderr); }

  art::FindOneP<numue::NumuE> slice2numue(slice, evt, "numue");
  if(!slice2numue.isValid()) { fputs("No slice2numue\n", stderr); }

  art::FindOneP<caf::StandardRecord> slice2caf(slice, evt, "cafmaker");
  if(!slice2caf.isValid()) { fputs("No slice2caf\n", stderr); }

  evtinfo einfo;
  einfo.run = evt.run();
  einfo.subrun = evt.subRun();
  einfo.pot = getpot(evt);

  std::map<std::string, std::string> metadata =
    meta::MetadataManager::getInstance().GetMetadata();
  einfo.cycle =
    metadata.count("simulated.cycle")?std::stoi(metadata["simulated.cycle"]):0;

  ntuple_header(einfo);

  if(rawtrigger->empty()){ fprintf(stderr, "No raw trigger!\n"); return; }
  if(tracks    ->size() == 0){ fprintf(stderr, "No tracks!\n"); return; }

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

    if(!delta_and_length(event_length_tdc, delta_tdc, flatdaq, rawtrigger)){
      fprintf(stderr, "Could not get delta and length\n");
      return;
    }
  }
  else{ // XXX
    event_length_tdc = 500 * 64;
    delta_tdc = 224 * 64;
  }

  einfo.triggerlength = event_length_tdc * 1000. / TDC_PER_US;
  einfo.starttime = -(delta_tdc * 1000. / TDC_PER_US);

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
    t.remid = track2remid.isValid()?track2remid.at(c)->Value():0;

    if(0 > (t.slice = which_slice_is_this_track_in(t, slice))) return;

    t.slice_energy = slice2numue.isValid()?slice2numue.at(t.slice)->E():0;

    if(slice2caf.isValid()){
      const art::Ptr<caf::StandardRecord> sr = slice2caf.at(t.slice);
      if(!goodspill(sr)){
        printf("Skipping bad spill\n");
        return;
      }
      t.contained_slice = containedND(sr);
      t.cvn = sr->sel.cvn.numuid;

      static caf::SRProxy srProxy(0, "");
      CopyRecord(*sr, srProxy);

      if(sr->mc.nu.empty()){
        t.mcweight = 1;
      }
      else{
        t.mcweight = 0;
        // Take the average for slices with multiple neutrino interactions.
        // But typically there is just one.
        for(unsigned int truthi = 0; truthi < sr->mc.nu.size(); truthi++)
          t.mcweight +=
            // This is the correct weight for the flux as of 2018-02-01.
            // https://neutrino.slack.com/archives/C02FS4L15/p1517526480000156
            sr->mc.nu[truthi].rwgt.ppfx.cv *
            // Chris Backhouse told me on Slack to do this
            // https://neutrino.slack.com/archives/C02FS4L15/p1517528737000298
            ana::kXSecCVWgt2018(&srProxy)/sr->mc.nu.size();
      }
    }

    t.true_pdg = t.true_nupdg = t.true_nucc = t.true_atom_cap
      = t.true_neutrons = 0;
    t.true_nuint = -1;
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
        t.true_nuint = truths[0].neutrinoInt->GetNeutrino().Mode();
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
      if(particles.empty()) continue;

      t.true_pdg = particles[0]->PdgCode();

      // As of 2017-05-03, sim::Particle::EndE() seems to always
      // return the particle's mass, no matter what happens to it.
      // sim::Particle::EndProcess() is also not useful for finding out
      // whether a particle has stopped, since it always holds an empty
      // string. So to find out what happens to a particle, you have to
      // laboriously look at its daughters...
      for(int d = 0; d < particles[0]->NumberDaughters(); d++){
        sim::ParticleNavigator::const_iterator it =
          pnav.find(particles[0]->Daughter(d));

        if(it == pnav.end()) continue;

        // Count neutrons produced sensibly close to the end of the
        // track, i.e. attempt to avoid counting neutrons produced in
        // some inelastic scatter far away, because we won't pick those
        // up in data. Otherwise, we find many seemingly-mysterious
        // cases of mu+ making neutrons. I use a mix of reconstructed
        // and true information here because I'm not sure I can reliably
        // get the position that a particle ends. See above comments on
        // EndE(). Also, in data I will be searching around the ends of
        // reconstructed tracks, not true ones, so probably this is the
        // better definition anyway.
        t.true_neutrons += (
          it->second->PdgCode() == 2112 &&
          sqrt(pow((*tracks)[c].Stop().X() - it->second->Position().X(), 2) +
            pow((*tracks)[c].Stop().Y() - it->second->Position().Y(), 2) +
            pow((*tracks)[c].Stop().Z() - it->second->Position().Z(), 2)) < 20.0
        );

        const std::string dproc = it->second->Process();
        if(dproc == "muMinusCaptureAtRest" ||
           dproc == "hBertiniCaptureAtRest"){
          t.true_atom_cap = 1;
        }
        else if(dproc == "Decay"){
          // I'm stuffing "it decayed" into the same variable because
          // decay is like "really really didn't capture". Just will
          // have to be careful to test for specific values and not
          // treat as boolean.
          t.true_atom_cap = -1;
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
    tinfo.end_dx = trk.StopDir().X();
    tinfo.end_dy = trk.StopDir().Y();
    tinfo.end_dz = trk.StopDir().Z();

    {
      static std::set<unique_track> unique_tracks;
      const std::pair<std::set<unique_track>::iterator, bool> does_it_blend =
        unique_tracks.insert(sorted_tracks[t]);
      if(!does_it_blend.second){ // true if the new track was inserted
        printf("Seen exactly this track before.  Known DAQ problem.  Skipping\n");
        continue;
      }
    }

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
    tinfo.cvn   = sorted_tracks[t].cvn;
    int lasthiti_even = 0, lasthiti_odd = 0;
    last_hits(lasthiti_even, lasthiti_odd, trk);
    tinfo.last_plane_even = trk.Cell(lasthiti_even)->Plane();
    tinfo.last_plane_odd  = trk.Cell(lasthiti_odd) ->Plane();
    tinfo.last_cell_even  = trk.Cell(lasthiti_even)->Cell();
    tinfo.last_cell_odd   = trk.Cell(lasthiti_odd) ->Cell();
    tinfo.time = mean_late_track_time(trk);
    tinfo.dother = make_dother(tinfo, slice);

    cluster_search(xt,  false, einfo, sorted_hits, trkhits, tinfo);
#if 0
    cluster_search(all, false, einfo, sorted_hits, trkhits, tinfo);
    cluster_search(ex,  false, einfo, sorted_hits, trkhits, tinfo);
    cluster_search(ex2, false, einfo, sorted_hits, trkhits, tinfo);
    cluster_search(x,   false, einfo, sorted_hits, trkhits, tinfo);
    for(int O = 0; O < 4; O++)
      cluster_search(xt,  true,  einfo, sorted_hits, trkhits, tinfo);
#endif
  }

  // if(neardet) // how do you test this?
  for(unsigned int t = 0; t < sorted_tracks.size(); t++){
    const int which = which_saved_track_array(sorted_tracks[t]);
    offspacetrk[which].push_back(sorted_tracks[t]);
  }
  for(unsigned int j = 0; j < MAXSAVEDTRACKTYPE; j++)
    while(offspacetrk[j].size() > savedtrack_size)
      offspacetrk[j].pop_front();
}

DEFINE_ART_MODULE(PostMuon)

} // end namespace PostMuon
//////////////////////////////////////////////////////////////////////////
