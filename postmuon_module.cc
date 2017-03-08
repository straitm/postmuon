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

#include "RawData/FlatDAQData.h"
#include "RawData/DAQHeader.h"
#include "RawData/RawTrigger.h"

#include "GeometryObjects/PlaneGeo.h"

#include <string>
#include <algorithm>

#include <signal.h>

static FILE * OUT = NULL;
static int NhitTrackTimeAveraging = 1; // dummy, to be overwritten in PostMuon()
static bool TracksAreDown = true; // to be overwritten in PostMuon()
static double MaxDistInCells = 0.123456; // dummy as well


struct evtinfo{
  int run;
  int event;
  double triggerlength;
  double starttime;
};

struct trkinfo{
  rb::Track trk;
  int i; // index of track in the track array
  double time;
  int lasthiti_even, lasthiti_odd;
  float sx, sy, sz, ex, ey, ez; // start and end position, after flip correction
};

struct cluster{
  int i; // which cluster for the given track
  int type; // defines the set of cuts used for this cluster
  float first_accepted_time; // time of the first hit in the cluster
  float last_accepted_time; // time of the last hit in the cluster
  float mindist;  // minimum hit distance from track end
  float dist2sum; // summed hit distance from track end
  float tsum;     // summed times, in ns, of a cluster of delayed hits
  float previous_cluster_t; // what it says
  float esum;     // summed calibrated energy, in MeV
  float pesum;    // summed PE
  int   adcsum;   // summed ADC
  int nhit;       // number of hits in full delayed cluster
};

static void resetcluster(cluster & res)
{
  res.first_accepted_time = -1000000;
  res.last_accepted_time  = -1000000;
  res.mindist = 1000000;
  res.dist2sum = 0;
  res.tsum = 0;
  res.esum = 0;
  res.pesum = 0;
  res.adcsum = 0;
  res.nhit = 0;
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


static double NS_PER_TDC = 1000/64.;
static double NS_PER_MICROSLICE = 50e3;
static int DCMS_PER_FD = 14 * 12; // Hack.  Should look this up in case they aren't all on for a given run or event
                                  // Also want to use ND data sometimes XXX

// Geometrically about correct, but perhaps should be scaled by density or
// radiation length or neutron cross section or something.  Or not, since
// which of those is right depends on what you're looking at.
const double planes_per_cell = 76./39.;

/*
  Return the average hit time of the latest 'NhitTrackTimeAveraging' hits, in
  nanoseconds.

  With a largish 'NhitTrackTimeAveraging', this is meant to provide a robust
  measure of the track time, even if an early Michel decay gets
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
  for(unsigned int i = std::max(0, (int)times.size() - NhitTrackTimeAveraging);
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
  Returns true if the given hit is in the same module as any hit coincident
  (within the time given below) with the given track.  Along with excluding the
  track itself, this is meant to do better at excluding cells that are
  inefficient.
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
  Takes a track and a hit and determines the distance between the
  end of the track and the hit.  lasthiti_{even,odd} are the indices
  of the last hits in each view.
*/
static float dist_trackend_to_cell(const rb::Track & __restrict__ trk,
                                   const rb::CellHit & __restrict__ chit,
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
  If the hit is near the track and after it in time, return the distance
  in cells.  Otherwise, return -1 if the hit preceeds the track, or -2
  if it is too far away.
*/
static double hit_near_track(const trkinfo & __restrict__ tinfo,
  const rb::CellHit & __restrict__ chit)
{
  const float aftertime = chit.TNS();
  if(aftertime < tinfo.time) return -1;

  // Accept the hit even if it is in the track!  Because if a Michel
  // hit gets swept up into the track, this is the only way to see it.
  //if(hit_is_in_track(chit, trk)) return -1;

  const double dist = dist_trackend_to_cell(tinfo.trk, chit, tinfo.lasthiti_even,
                                            tinfo.lasthiti_odd);

  if(dist > MaxDistInCells) return -2;

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
                              const cluster & __restrict__ cluster)
{
  const double timeleft = einfo.triggerlength - tinfo.time + (NS_PER_MICROSLICE - einfo.starttime);

//  if(tinfo.i == 0)
//    printf("%f %f %f\n", einfo.triggerlength, einfo.starttime, tinfo.time);


  const double tsx = tinfo.sx;
  const double tsy = tinfo.sy;
  const double tsz = tinfo.sz;

  const double tx = tinfo.ex;
  const double ty = tinfo.ey;
  const double tz = tinfo.ez;

  fprintf(OUT, "%d %d ", einfo.run, einfo.event);
  fprintf(OUT, "%d ", tinfo.i);
  fprintf(OUT, "%.1f ", einfo.triggerlength);
  fprintf(OUT, "%f ", tinfo.time);
  fprintf(OUT, "%d ", cluster.i);
  fprintf(OUT, "%.1f %.1f %.1f ", tsx, tsy, tsz);
  fprintf(OUT, "%.1f %.1f %.1f ", tx, ty, tz);
  fprintf(OUT, "%f ", timeleft);

  fprintf(OUT, "%d ", cluster.type);

  if(cluster.nhit)
    fprintf(OUT, "%f ",
            float(cluster.tsum)/cluster.nhit-tinfo.time);
  else
    fprintf(OUT, "-1 ");

  if(cluster.nhit)
    fprintf(OUT, "%f ", float(cluster.tsum)/cluster.nhit-cluster.previous_cluster_t);
  else
    fprintf(OUT, "0 ");

  fprintf(OUT, "%.3f ", cluster.mindist);

  if(cluster.dist2sum == 0) fprintf(OUT, "0 ");
  else fprintf(OUT, "%.3f ", cluster.dist2sum/cluster.nhit);

  fprintf(OUT, "%.3f ", cluster.esum);
  fprintf(OUT, "%.3f ", cluster.pesum);
  fprintf(OUT, "%d ", cluster.adcsum);

  fprintf(OUT, "%d ", cluster.nhit);
  fprintf(OUT, "%.3f ",
          cluster.last_accepted_time - cluster.first_accepted_time);
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

  // If I don't do this, it responds to SEGV by backgrounding itself and
  // stopping, and upon being foregrounded, by claiming to have been
  // aborted and going back to sleep, and upon a second foregrounding,
  // finally actually dying. What.
  signal(SIGSEGV, SIG_DFL);
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
  return a.TNS() < b.TNS();
}

static bool compare_cellhit_TDC(const rb::CellHit & __restrict__ a,
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

static void cluster_search(const int type,
  const evtinfo & __restrict__ einfo,
  const std::vector<rb::CellHit> & __restrict__ sorted_hits,
  const std::vector<rb::CellHit> & __restrict__ trkhits,
  int & first_hit_to_consider, const trkinfo & __restrict__ tinfo)
{
  art::ServiceHandle<calib::Calibrator> calthing;

  cluster clu = mkcluster();
  clu.i = 0;
  clu.type = type;
  clu.previous_cluster_t = tinfo.time;

  for(int c = first_hit_to_consider; c < (int)sorted_hits.size(); c++){
    const rb::CellHit & chit = sorted_hits[c];

    const double dist = hit_near_track(tinfo, chit);

    if(dist < 0){
      // This hit is before this track, so it will also be next time we look
      // and for all further tracks. ~5% speed bump from this optimization.
      if(dist == -1) first_hit_to_consider = c+1;

      continue;
    }

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
      print_ntuple_line(einfo, tinfo, clu);
      clu.i++;
      const float thistime = float(clu.tsum)/clu.nhit;
      resetcluster(clu);
      clu.previous_cluster_t = thistime;
    }

    // Add everything to the cluster *after* the print of the
    // previous cluster (or no-op).

    clu.nhit++;

    clu.dist2sum += dist*dist;
    if(dist < clu.mindist) clu.mindist = dist;

    clu.last_accepted_time = chit.TNS();
    if(clu.first_accepted_time < 0)
      clu.first_accepted_time = chit.TNS();
    clu.tsum += chit.TNS();

    clu.esum += rhit.GeV()*1000;
    clu.pesum += chit.PE();
    clu.adcsum += chit.ADC();
  }

  // print last cluster, or track info if no cluster
  print_ntuple_line(einfo, tinfo, clu);
}

PostMuon::~PostMuon() { }

void PostMuon::analyze(const art::Event& evt)
{
  // I'm so sorry that I have to do this.  And, my goodness, doing
  // it in the constructor isn't sufficient.  If this isn't done,
  // it responds to PIPE by going into an endless loop.
  signal(SIGPIPE, SIG_DFL);

  art::Handle< std::vector<rawdata::FlatDAQData> > flatdaq;
  evt.getByLabel("daq", flatdaq);

  art::Handle<rawdata::DAQHeader> daqheader;
  evt.getByLabel("daq", daqheader);

  art::Handle< std::vector<rawdata::RawTrigger> > rawtrigger;
  evt.getByLabel("daq", rawtrigger);

  art::Handle< std::vector<rb::CellHit> > cellcol; // get hits
  evt.getByLabel(fRawDataLabel, cellcol);

  art::Handle< std::vector<rb::Track> > tracks;
  evt.getByLabel("windowtrack", tracks);

  if(OUT == NULL){
    OUT = fopen(Form("postmuon_%d_%d.20170120.ntuple",
                     evt.run(), evt.subRun()), "w");
    if(OUT == NULL){
       fprintf(stderr, "Could not open output ntuple file\n");
       exit(1);
    }
    fprintf(OUT,
      "run/I:"
      "event/I:"
      "trk/I:"
      "triggerlength/F:"
      "trktime/F:"
      "i/I:"
      "trkstartx/F:"
      "trkstarty/F:"
      "trkstartz/F:"
      "trkx/F:"
      "trky/F:"
      "trkz/F:"
      "timeleft/F:"

      "type/I:"
      "t/F:"
      "dt/F:"
      "mindist/F:"
      "dist2/F:"
      "e/F:"
      "pe/F:"
      "adc/I:"
      "nhit/I:"
      "ctlen/F"
      "\n");
    printf("DEBUG TimeStart/l:TimeStart5/I:TimeStart8/I:TDCT0/l:ExtractionStart/l:GenTime/l:FirstTNS/f:LastTNS/f:nhit/d:FirstTDC/d:LastTDC/d:tdclen/d\n");
  }

  if(rawtrigger->empty()) return;

  daqdataformats::RawEvent raw;
  if(flatdaq->empty()) return;

  raw.readData((*flatdaq)[0].getRawBufferPointer());

  if(raw.getDataBlockNumber() == 0) return;

  raw.setFloatingDataBlock(0);
  daqdataformats::RawDataBlock& datablock = (*raw.getFloatingDataBlock());

  if(datablock.getHeader()->getMarker() == daqdataformats::datablockheader::SummaryBlock_Marker &&
     datablock.getHeader()->checkMarker()) return;

    for(unsigned int mi = 0; mi < datablock.getNumMicroBlocks(); mi++){
      datablock.setFloatingMicroBlock(mi);
      daqdataformats::RawMicroBlock * ub = datablock.getFloatingMicroBlock();

      // For gosh sakes, I can't figure out if there's a nice named function
      // that would provide this, but it is always the second word of the microslice,
      // which follows two words of microblock header, so just get it
      const uint32_t time_marker_low  = ((uint32_t *)(ub->getBuffer()))[3];
      const uint32_t time_marker_high = ((uint32_t *)(ub->getBuffer()))[4];
      printf("microblock %5d %08x %08x\n", mi, time_marker_low, time_marker_high);
    }
  }
  
  for(unsigned int i = 0; i < flatdaq->size(); i++){
    printf("FLATDAQ %d\n", i);
    for(unsigned int j = 0; j < (*flatdaq)[i].fRawBuffer.size(); j++){
      const int jflip = 4*(j/4) + (3 - j%4); // fix endianness

      if(j%4 == 0 ) printf("%08x  ", j);
      printf("%02x ", (unsigned char)(*flatdaq)[i].fRawBuffer[jflip]);
      if(j%4 == 3) printf("\n");
    }
    printf("\n");
  }

  evtinfo einfo;

  einfo.triggerlength =
    daqheader->TotalMicroSlices()*NS_PER_MICROSLICE/DCMS_PER_FD;
  einfo.starttime =
    (*rawtrigger)[0].TDCT0()*NS_PER_TDC;


  einfo.run = evt.run();
  einfo.event = evt.event();

  int firsttdc = 0, lasttdc = 0;

  // CellHits do *not* come in time order
  std::vector<rb::CellHit> sorted_hits;

  for(int c = 0; c < (int)cellcol->size(); c++)
    sorted_hits.push_back((*cellcol)[c]);

  std::sort(sorted_hits.begin(), sorted_hits.end(), compare_cellhit_TDC);

  firsttdc = sorted_hits[0].TDC();
  lasttdc = sorted_hits[sorted_hits.size()-1].TDC();

  std::sort(sorted_hits.begin(), sorted_hits.end(), compare_cellhit);

  printf("DEBUG %12llu %5llu %5llu %10llu %12llu %12llu %lf %lf %d %d %d %d\n",
    (*rawtrigger)[0].fTriggerTimingMarker_TimeStart, 
    (*rawtrigger)[0].fTriggerTimingMarker_TimeStart & 0x1fULL, 
    (*rawtrigger)[0].fTriggerTimingMarker_TimeStart & 0xffULL, 
    (*rawtrigger)[0].TDCT0(),
    (*rawtrigger)[0].fTriggerTimingMarker_ExtractionStart, // always zero
    (*rawtrigger)[0].fTriggerTime_GenTime,   // same as TimeStart
    sorted_hits[0].TNS(),
    sorted_hits[sorted_hits.size()-1].TNS(),
    (int)sorted_hits.size(),
    firsttdc,
    lasttdc, 
    lasttdc-firsttdc);
  fflush(stdout);

  return;

  std::vector<rb::Track> sorted_tracks;
  if(!tracks->empty()){
    for(int c = 0; c < (int)tracks->size(); c++){
      // Badly tracked tracks are not useful
      if((*tracks)[c].Stop().X() == 0 ||
         (*tracks)[c].Stop().Y() == 0) continue;

      sorted_tracks.push_back((*tracks)[c]);
    }
    std::sort(sorted_tracks.begin(), sorted_tracks.end(), compare_track);
  }

  std::vector<rb::CellHit> trkhits = make_trkhits(sorted_tracks);

  int first_hit_to_consider = 0;

  for(unsigned int t = 0; t < sorted_tracks.size(); t++){
    const rb::Track & trk = sorted_tracks[t];

    trkinfo tinfo;
    const bool needs_flip = shall_we_flip_it(trk);
    tinfo.sx = needs_flip? trk.Stop ().X(): trk.Start().X();
    tinfo.sy = needs_flip? trk.Stop ().Y(): trk.Start().Y();
    tinfo.sz = needs_flip? trk.Stop ().Z(): trk.Start().Z();
    tinfo.ex = needs_flip? trk.Start().X(): trk.Stop ().X();
    tinfo.ey = needs_flip? trk.Start().Y(): trk.Stop ().Y();
    tinfo.ez = needs_flip? trk.Start().Z(): trk.Stop ().Z();

    // It's a waste of time to search for clusters around tracks that almost
    // certainly represent through-goers.  Note that this does nothing to near
    // detector events
    if(fabs(tinfo.ex) >  750) continue;
    if(     tinfo.ey  < -730) continue;
    if(     tinfo.ez  <   25) continue;
    if(     tinfo.ez  > 5935) continue;

    // Don't exclude tracks based on their starting positions, since I want this
    // to work for both cosmics and beam
    ;

    tinfo.trk = trk;
    tinfo.i = t;
    tinfo.lasthiti_even = 0, tinfo.lasthiti_odd = 0;
    last_hits(tinfo.lasthiti_even, tinfo.lasthiti_odd, trk);
    tinfo.time = mean_late_track_time(trk);

    cluster_search(all, einfo, sorted_hits, trkhits, first_hit_to_consider, tinfo);
    cluster_search(ex,  einfo, sorted_hits, trkhits, first_hit_to_consider, tinfo);
    cluster_search(ex2, einfo, sorted_hits, trkhits, first_hit_to_consider, tinfo);
    cluster_search(xt,  einfo, sorted_hits, trkhits, first_hit_to_consider, tinfo);
  }
}

DEFINE_ART_MODULE(PostMuon);

} // end namespace PostMuon
//////////////////////////////////////////////////////////////////////////
