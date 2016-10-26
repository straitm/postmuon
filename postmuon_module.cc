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

#include <string>
#include <algorithm>

#include <signal.h>

#include "art/Framework/Core/EDAnalyzer.h"


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


__attribute__((unused)) static bool hit_is_in_track(const rb::CellHit & chit,
                                                    const rb::Track & trk)
{
  for(unsigned int i = 0; i < trk.NCell(); i++){
    const rb::CellHit & trk_chit = *(trk.Cell(i));

    // Maybe not strictly true as there could be more than one hit in the same
    // cell at different times, but let's not worry about that now.
    // Also exclude hits that are directly adjacent to track hits.
    if(abs(trk_chit.Plane() - chit.Plane()) <= 2 &&
       abs(trk_chit. Cell() - chit.Cell ()) <= 1) return true;
  }
  return false;
}


__attribute__((unused)) static double kUSEC_PER_TDC = 1./64.;

// Geometrically about correct, but perhaps should be scaled by density or
// radiation length or neutron cross section or something.  Or not, since
// which of those is right depends on what you're looking at.
const double planes_per_cell = 76./39.;

static float dist_trackend_to_cell(const rb::Track & trk,
                                   const rb::CellHit & chit,
                                   const int lasthiti_even,
                                   const int lasthiti_odd)
{
  const int lastplane_even = trk.Cell(lasthiti_even)->Plane();
  const int lastplane_odd =  trk.Cell(lasthiti_odd) ->Plane();
  const int lastcell_even =  trk.Cell(lasthiti_even)->Cell();
  const int lastcell_odd =   trk.Cell(lasthiti_odd) ->Cell();

  return
    chit.Plane()%2 == 0 ?
      sqrt(pow(planes_per_cell*(chit.Plane() - lastplane_even), 2) + 
           pow(                 chit.Cell()  - lastcell_even  , 2))
   :  sqrt(pow(planes_per_cell*(chit.Plane() - lastplane_odd ), 2) + 
           pow(                 chit.Cell()  - lastcell_odd   , 2));

}

static bool hit_near_track(const rb::Track & trk,
  const int lasthiti_even, const int lasthiti_odd, 
  const art::Handle< std::vector<rb::CellHit> > & cellcol,
  const int celli)
{
  const rb::CellHit & chit = (*cellcol)[celli];

  const int tracktime = (trk.Cell(lasthiti_even)->TDC() +
                         trk.Cell(lasthiti_odd )->TDC())/2;

  const int aftertdc = chit.TDC();
  if(aftertdc < tracktime) return false;

  // Accept the hit even if it is in the track!  Because if a Michel
  // hit gets swept up into the track, this is the only way to see it.
  //if(hit_is_in_track(chit, trk)) return false;

  // Enough to catch ~99% of gammas from neutron capture according to my toy MC
  const double maxdist_in_cells = 20.1;

  const double dist = dist_trackend_to_cell(trk, chit, lasthiti_even,
                                            lasthiti_odd);


  // Hits too close, but otherwise correct, probably signal that we
  // are picking up a Michel decay.  Bail out.
  if(dist > maxdist_in_cells) return false;

  return true;
}

static void last_hits(int & lasthiti_even,
                      int & lasthiti_odd,
                      const rb::Track & trk)
{
  float lowest_even = 1e30, lowest_odd = 1e30;

  for(int c = 0; c < (int)trk.NCell(); c++){
    const rb::CellHit & chit = *(trk.Cell(c));
    const rb::RecoHit rhit = trk.RecoHit(c);

    if(chit.Plane()%2 == 0 && rhit.Y() < lowest_even){
      lasthiti_even = c;
      lowest_even = rhit.Y();
    }

    if(chit.Plane()%2 == 1 && rhit.Y() < lowest_odd){
      lasthiti_odd = c;
      lowest_even = rhit.Y();
    }
  }
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
  for(unsigned int i = std::max(0, (int)tdcs.size() - max_for_avg); i < tdcs.size(); i++)
    acc += tdcs[i];

  return int(acc / std::min(max_for_avg, (int)tdcs.size()) + 0.5); 
}

static void print_ntuple_line(const rb::Track & trk,
                              const rb::CellHit & chit,
                              const int asum,
                              const float esum,
                              const int cluster_i,
                              const int cluster_nhit)
{

  int lasthiti_even = 0, lasthiti_odd = 0;
  last_hits(lasthiti_even, lasthiti_odd, trk);
  const double dist = dist_trackend_to_cell(
    trk, chit, lasthiti_even, lasthiti_odd);

  const int tracktime = mean_late_track_time(trk);

  printf("ntuple: %d %f %f %f %f %f %f %d %f %d\n",
    cluster_i,    tracktime *kUSEC_PER_TDC,
    (chit.TDC() - tracktime)*kUSEC_PER_TDC,
    trk.Stop().X(), trk.Stop().Y(), trk.Stop().Z(),
    dist, asum, esum, cluster_nhit);
}

namespace PostMuon{

PostMuon::PostMuon(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset), fRemoveBadChans(pset.get<bool>("RemoveBadChans")),
  fRawDataLabel(pset.get< std::string >("RawDataLabel"))
{
  
}

PostMuon::~PostMuon() { }

void PostMuon::analyze(const art::Event& evt)
{
  // I'm so sorry that I have to do this.  And, my goodness, doing
  // it in the constructor isn't sufficient.
  signal(SIGPIPE, SIG_DFL);

  art::Handle< std::vector<rb::CellHit> > cellcol; // get hits
  evt.getByLabel(fRawDataLabel, cellcol);

  art::Handle< std::vector<rb::Track> > tracks;
  evt.getByLabel("kalmantrackmerge", tracks);

  {
    static int NOvA = printf(
      "ntuple: i:trktime:t:trkx:trky:trkz:dist:adc:e:nhit\n");
    NOvA = NOvA;
  }

  art::ServiceHandle<calib::Calibrator> calthing;

  for(unsigned int t = 0; t < tracks->size(); t++){

    const rb::Track & trk = (*tracks)[t];
    int lasthiti_even = 0, lasthiti_odd = 0;
    last_hits(lasthiti_even, lasthiti_odd, trk);

    int cluster_nhit = 0;
    int last_accepted_time = -1;
    int last_accepted_i = -1;
    int asum = 0; // sumed ADC of a cluster of delayed hits
    float esum = 0; // summed calibrated energy, in MeV
    int cluster_i = 0;

    for(int c = 0; c < (int)(*cellcol).size(); c++){
      if(!hit_near_track(trk, lasthiti_even,
         lasthiti_odd, cellcol, c)) continue;

      const rb::CellHit & chit = (*cellcol)[c];

      const rb::RecoHit rhit = calthing->MakeRecoHit(chit,
         chit.Plane()%2 == 0? trk.Stop().X(): // doc-11570
                              trk.Stop().Y());

      // Hits are time ordered.  Only report on hits that are
      // separated in time from other accepted hits by at least
      // one quiet TDC tick.
      const bool newcluster = chit.TDC() > last_accepted_time + 1;

      if(newcluster && cluster_nhit){
        print_ntuple_line(trk, (*cellcol)[c-1],
                          asum, esum, cluster_i++, cluster_nhit);
        esum = cluster_nhit = 0;
      }
     

      cluster_nhit++;
      asum += chit.ADC();
      if(rhit.IsCalibrated()){
        esum += rhit.GeV()*1000;
      }
      else{
        // For now, just throw away uncalibrated hits,
        // XXX is there a better approach?
        //printf("Uncalibrated ADC: %d\n", chit.ADC());
      }

      last_accepted_time = chit.TDC();
      last_accepted_i = c;
    }

    // print last cluster
    if(cluster_nhit)
      print_ntuple_line(trk, (*cellcol)[last_accepted_i],
                        asum, esum, cluster_i++, cluster_nhit);

  }
}

DEFINE_ART_MODULE(PostMuon);

} // end namespace PostMuon
//////////////////////////////////////////////////////////////////////////
