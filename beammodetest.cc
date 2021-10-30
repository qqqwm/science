#include <evt/Event.h>
#include <evt/SimEvent.h>
#include <evt/RecEvent.h>
#include <evt/sim/VertexTrack.h>
#include <evt/rec/RecEventConst.h>
#include <evt/rec/Trigger.h>

#include <io/EventFileChain.h>
#include <io/IoCodes.h>
#include <io/EventFile.h>

#include <fwk/CentralConfig.h>

#include <det/Detector.h>
#include <det/DetectorConst.h>
#include <det/TPCConst.h>
#include <det/Target.h>
#include <det/Beam.h>
#include <det/TriggerConst.h>

#include <utl/PDGParticleIds.h>
#include <utl/ShineUnits.h>
#include <utl/Vector.h>
#include <utl/MathConst.h>
#include <utl/PhysicalConst.h>
#include <utl/PDGParticleIds.h>
#include <utl/DatabasePDG.h>
#include <utl/GeometryUtilities.h>
#include <utl/UTCDateTime.h>

#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
//#include <TRandom.h>
//#include <TH2I.h>

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;
//using namespace io;
using namespace utl;
using namespace evt;
//using namespace det::TriggerConst;

int main (int argc, char* argv[])
{
  if (argc < 1)
  {
    cerr << "usage: " << argv[0] << "<input filenames>" << endl;
    return 1;
  }
  vector<string> fileNames(argv + 2, argv + argc);
  EventFileChain eventFileChain(fileNames);
  eventFileChain.SetParsingOptions(false);
  eventFileChain.StopOnOpenProblem(false);
  Event event;
  fwk::CentralConfig::GetInstance("bootstrap.xml");
  det::Detector& detector = det::Detector::GetInstance();
  
  // <runNumber> 26451 </runNumber>   <eventTime> 2016-11-25T06:26:03.0 </eventTime>
  const TimeStamp timeStamp = UTCDateTime(2016, 11, 25, 06, 26, 03).GetTimeStamp();
  detector.Update(timeStamp, 26451);
  const det::Target& target = detector.GetTarget();
  int nEvent = 0, nPassedCuts = 0;
  TH1D* eventCuts = new TH1D("eventCuts", "", 1, 0, 0);
  TH1D* hist_bsvpz  = new TH1D("hist_bsvpz", "BeamStopVertexPosition;z [cm];Entries", 320, -600., 1000.);
  TH1D* hist_bsvpzd  = new TH1D("hist_bsvpzd", "BeamStopVertexPosition;z [cm];Entries", 320, -600., -580.);
  TH2D* hist_bsvpxy = new TH2D("hist_bsvpxy", "BeamStopVertexPosition; x[cm]; y[cm]", 200, -7., 7., 200, -7., 7. );
  //printf("target.GetLength() = %f\ntarget.GetMaterialZ() = %f\ntarget.GetMaterialA() = %f\ntarget.GetMaterialDensity() = %f\ntarget.GetMaterialInteractionLength() = %f\ntarget.GetMaterialMolarMass() = %f\ntarget.GetStatus() = %i\n", target.GetLength(), target.GetMaterialZ(), target.GetMaterialA(), target.GetMaterialDensity(), target.GetMaterialInteractionLength(), target.GetMaterialMolarMass(), static_cast<int> target.GetStatus());
  while (eventFileChain.Read(event) == eSuccess) {
    ++nEvent;

    eventCuts->Fill("Events", 1.);
    //cout << "Event# " << nEvent << endl;
    const RawEvent& rawEvent = event.GetRawEvent();
    const raw::Beam& rawBeam = rawEvent.GetBeam();
    const raw::Trigger& rawTrigger = rawBeam.GetTrigger();
    const RecEvent& recEvent = event.GetRecEvent();
    const rec::Beam& recBeam = recEvent.GetBeam();
    const rec::Trigger& recTrigger = recBeam.GetTrigger();
    const SimEvent& simEvent = event.GetSimEvent();

    const evt::EventHeader& header = event.GetEventHeader();

    if (header.IsInitialized())
      detector.Update(header.GetTime(), header.GetRunNumber());

    bool hasT2prsc = rawTrigger.HasTrigger(det::TriggerConst::eT2, det::TriggerConst::ePrescaled);
    //cout << "HasTrigger(T2 prescaled): " << hasT2prsc;
    if (hasT2prsc) {
      //cout <<  " = " << rawTrigger.IsTrigger(det::TriggerConst::eT2, det::TriggerConst::ePrescaled);
    }
    //cout << endl;

    if (!simEvent.HasMainVertex()) {
      eventCuts->Fill("NOMainVertex", 1.);
      //cout << "Has No MainVertex\n";
    }
    else
    {
      //const Point& p2 = simEvent.GetMainVertex().GetPosition();
      //cout << "MainVertex: " << p2.GetX() << ", " << p2.GetY() << ", " << p2.GetZ() << endl;
    }
    int nVertexTracks = 0;
    for (list<evt::sim::VertexTrack>::const_iterator simTrackIter = simEvent.Begin<evt::sim::VertexTrack>(),
         simTrackEnd = simEvent.End<evt::sim::VertexTrack>(); simTrackIter != simTrackEnd; ++simTrackIter) {
      ++nVertexTracks;
      const evt::sim::VertexTrack& beamVertexTrack = *simTrackIter;
      if (beamVertexTrack.GetType() == sim::VertexTrackConst::eBeam) {
        //cout << "beamVertexTrack.GetType() == sim::VertexTrackConst::eBeam\n";
        if (beamVertexTrack.HasStopVertex()) {
          eventCuts->Fill("HasStopVertex", 1.);
          //cout << "beamVertexTrack.HasStopVertex()\n";
          evt::Index< evt::sim::Vertex > indexBeamStopVertex = beamVertexTrack.GetStopVertexIndex();
          evt::sim::Vertex BeamStopVertex = simEvent.Get(indexBeamStopVertex);
          const utl::Point& BeamStopVertexPosition = BeamStopVertex.GetPosition();
          hist_bsvpz->Fill(BeamStopVertexPosition.GetZ());
          hist_bsvpzd->Fill(BeamStopVertexPosition.GetZ());
          hist_bsvpxy->Fill(BeamStopVertexPosition.GetX(), BeamStopVertexPosition.GetY());
          //cout << "BeamStopVertex: " << BeamStopVertexPosition.GetX() << ", " << BeamStopVertexPosition.GetY() << ", " << BeamStopVertexPosition.GetZ() << endl;
          if (target.IsIn(BeamStopVertexPosition))
            //eventCuts->Fill("BeamStopVertexInTarget", 1.);
          if (BeamStopVertex.GetType() == sim::VertexConst::eHadronicInelastic)
            eventCuts->Fill("BeamStopVertexHadrInelas", 1.);
          //cout << "is BeamStopVertex in target: " << target.IsIn(BeamStopVertexPosition) << endl;
          //cout << "is BeamStopVertex HadronicInelastic: " << (BeamStopVertex.GetType() == sim::VertexConst::eHadronicInelastic) << endl;
          if (!(target.IsIn(BeamStopVertexPosition) && BeamStopVertex.GetType() == sim::VertexConst::eHadronicInelastic)) {
            //cout << "eContinueLoop\n";
          }
          else {
            //++nPassedCuts;
            eventCuts->Fill("passed", 1.);
            //cout << "Cut Passed\n";
          }
        } else {
          //cout << "!beamVertexTrack.HasStopVertex(): eContinueLoop\n";
        }
      }
    }
    //cout << "nVertexTracks: " << nVertexTracks << endl;

  }
  //cout << "nPassedCuts: " << nPassedCuts << endl;
  TFile outFile(Form("t.%s.root", argv[1]), "RECREATE");
  outFile.cd();
  //eventCuts->Write();
  hist_bsvpz->Write();
  hist_bsvpzd->Write();
  hist_bsvpxy->Write();
  
}

// /afs/cern.ch/work/y/yebondar/PbPb30_BeamMode_26.sim.root