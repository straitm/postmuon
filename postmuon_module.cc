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

#include "Simulation/ParticleNavigator.h"

#include <string>

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


static bool hit_is_in_track(const rb::CellHit & chit,
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

// Returns time between a stray hit near the end of a track that might
// be from an x-ray and the time of the end of the track.  If there are
// no such hits, returns -999.
static double hits_near_track_end(const rb::Track & trk,
  const int lasthiti_even, const int lasthiti_odd, 
  const art::Handle< std::vector<rb::CellHit> > & cellcol)
{
  // for each view, check for small stray hits near
  // track end and in time with the track.

  const int lastplane_even = trk.Cell(lasthiti_even)->Plane();
  const int lastcell_even = trk.Cell(lasthiti_even)->Cell();
  const int lastplane_odd = trk.Cell(lasthiti_odd)->Plane();
  const int lastcell_odd = trk.Cell(lasthiti_odd)->Cell();

  const double mev_per_adc_even =
    1000*trk.RecoHit(lasthiti_even).GeV()/trk.Cell(lasthiti_even)->ADC();
  const double mev_per_adc_odd  =
    1000*trk.RecoHit(lasthiti_odd ).GeV()/trk.Cell(lasthiti_odd )->ADC();

  __attribute__((unused)) static double kUSEC_PER_TDC = 1./64.;

  const double lasttime = kUSEC_PER_TDC * (trk.Cell(lasthiti_even)->TDC() +
                                           trk.Cell(lasthiti_odd )->TDC())/2.;

  for(unsigned int i = 0; i < cellcol->size(); i++){
    const rb::CellHit & chit = (*cellcol)[i];
    if(hit_is_in_track(chit, trk)) continue;

    const double maxdist = 6.1; // maybe reasonable?  Mean free path is ~25cm
    // must not be directly adjacent to track end (i.e. a reco failure)
    const double mindist = sqrt(2*2*1.5*1.5 + 1) + 0.01;
    const double max_dt = 0.5; // ok? Dunno.
    const int max_mev = 6.0;

    const double dist = 
      chit.Plane()%2 == 0 ?
        sqrt(pow(1.5*(chit.Plane() - lastplane_even), 2) + 
             pow(     chit.Cell () - lastcell_even  , 2))
     :  sqrt(pow(1.5*(chit.Plane() - lastplane_odd ), 2) + 
             pow(     chit.Cell () - lastcell_odd   , 2));

    const bool right_dist = dist < maxdist && dist > mindist;

    const bool close_in_time = fabs(lasttime - kUSEC_PER_TDC * chit.TDC()) < max_dt;

    const double mev = chit.ADC() * (chit.Plane()%2 == 0? mev_per_adc_even
                                                        : mev_per_adc_odd);

    const bool right_energy = mev < max_mev;

    if(right_dist && close_in_time && right_energy)
      return lasttime - kUSEC_PER_TDC * chit.TDC();

    // Hits too close, but otherwise correct, probably signal that we
    // are picking up a Michel decay.  Bail out.
    if(dist < mindist && close_in_time && right_energy) return -999;
  }

  return -999;
}

namespace PostMuon{

PostMuon::PostMuon(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset), fRemoveBadChans(pset.get<bool>("RemoveBadChans")),
  fRawDataLabel(pset.get< std::string >("RawDataLabel"))
{
  
}

//......................................................................
PostMuon::~PostMuon() { }

//......................................................................

void PostMuon::analyze(const art::Event& evt)
{
  // I'm so sorry that I have to do this.  And, my goodness, doing
  // it in the constructor isn't sufficient.
  signal(SIGPIPE, SIG_DFL);

  art::Handle< std::vector<rb::CellHit> > cellcol; // get hits
  evt.getByLabel(fRawDataLabel, cellcol);

  art::Handle< std::vector<rb::Track> > tracks;
  evt.getByLabel("kalmantrackmerge", tracks);

  for(unsigned int t = 0; t < 1 && t < tracks->size(); t++){
    printf("I have a track\n");
    {
      static int HQNLM = printf(
        "event:truemicheltime:xrayhits:ntrueneutrons:neutrone:pvc:"
        "endx:endy:endz:costheta:heavy"
        ":e1:dedx1:plane1:cell1:truedx1:trueebirks1:truee1"
        ":e0:dedx0:plane0:cell0:truedx0:trueebirks0:truee0\n");
      HQNLM++;
    }

    const rb::Track & trk = (*tracks)[t];

    for(int c = 0; c < (int)trk.NCell(); c++){
      //const rb::CellHit & chit = *(trk.Cell(c));

    }
  }
}

DEFINE_ART_MODULE(PostMuon);

} // end namespace PostMuon
//////////////////////////////////////////////////////////////////////////
