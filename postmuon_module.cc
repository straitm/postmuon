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


__attribute__((unused)) static double kUSEC_PER_TDC = 1./64.;

static bool hit_near_track_end_and_after_it(const rb::Track & trk,
  const int lasthiti_even, const int lasthiti_odd, 
  const art::Handle< std::vector<rb::CellHit> > & cellcol,
  const int celli)
{
  const rb::CellHit & chit = (*cellcol)[celli];

  if(hit_is_in_track(chit, trk)) return false;

  const int afterplane = chit.Plane();
  const int aftercell  = chit.Cell();
  const int aftertdc   = chit.TDC();

  const int lastplane_even = trk.Cell(lasthiti_even)->Plane();
  const int lastplane_odd =  trk.Cell(lasthiti_odd) ->Plane();
  const int lastcell_even =  trk.Cell(lasthiti_even)->Cell();
  const int lastcell_odd =   trk.Cell(lasthiti_odd) ->Cell();

  const int tracktime = (trk.Cell(lasthiti_even)->TDC() +
                         trk.Cell(lasthiti_odd )->TDC())/2.;

  //const bool after = tracktime < aftertdc;

  //if(!after) return false;

  const double maxdist_in_cells = 8.1; // maybe reasonable? 

  // Geometrically about correct, but perhaps should be scaled by density or
  // radiation length or neutron cross section or something.  Or not, since
  // which of those is right depends on what you're looking at.
  const double planes_per_cell = 76./39.;

  const double dist = 
    chit.Plane()%2 == 0 ?
      sqrt(pow(planes_per_cell*(afterplane - lastplane_even), 2) + 
           pow(                 aftercell  - lastcell_even  , 2))
   :  sqrt(pow(planes_per_cell*(afterplane - lastplane_odd ), 2) + 
           pow(                 aftercell  - lastcell_odd   , 2));


  // Hits too close, but otherwise correct, probably signal that we
  // are picking up a Michel decay.  Bail out.
  if(dist > maxdist_in_cells) return false;

  printf("ntuple: %f %f %f %f %f %f %d\n",
         tracktime*kUSEC_PER_TDC,
         (aftertdc - tracktime)*kUSEC_PER_TDC,
         trk.Stop().X(), trk.Stop().Y(), trk.Stop().Z(),
         dist,

         // Would be nice to find the energy as in a RecoHit, but
         // takes some work since no existing code gives it a W position.
         chit.ADC());
  return true;
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
    {
      static int HQNLM = printf(
        "ntuple: trktime:t:trkx:trky:trkz:dist:adc"
        "\n"
        );
      HQNLM++;
    }

    const rb::Track & trk = (*tracks)[t];

    const bool south = trk.Dir().Z() < 0;

    int lasthiti_even = 0;
    int lastplane_even = 0;
    int lasthiti_odd = 0;
    int lastplane_odd = 0;

    for(int c = 0; c < (int)trk.NCell(); c++){
      const rb::CellHit & chit = *(trk.Cell(c));

      if(chit.Plane()%2 == 0 && (south ^ (chit.Plane() > lastplane_even))){
        lasthiti_even = c;
        lastplane_even = chit.Plane();
      }

      if(chit.Plane()%2 == 1 && (south ^ (chit.Plane() > lastplane_odd))){
        lasthiti_odd = c;
        lastplane_odd = chit.Plane();
      }
    }

    for(int c = 0; c < (int)(*cellcol).size(); c++){
      if(hit_near_track_end_and_after_it(
        trk, lasthiti_even, lasthiti_odd, cellcol, c)){
      }
    }
  }
}

DEFINE_ART_MODULE(PostMuon);

} // end namespace PostMuon
//////////////////////////////////////////////////////////////////////////
