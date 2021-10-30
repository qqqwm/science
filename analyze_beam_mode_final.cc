/**************************************************************************************************************
 * Author (April 2021): Daria Prokhorova daria.prokhorova@cern.ch
 **************************************************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>

#ifndef SHINELIB
#define SHINELIB
#include <evt/Event.h>
#include <evt/RecEvent.h>
#include <evt/rec/Trigger.h>
#include <fwk/CentralConfig.h>
#include <utl/ShineUnits.h>
#include <io/EventFile.h>
#include <io/EventFileChain.h>
#include <utl/ShineUnits.h>
#include <evt/rec/RecEventConst.h>
#include <det/TriggerConst.h>
#include <utl/MathConst.h>
#include <utl/PhysicalConst.h>
#include <det/TPCConst.h>
#include <det/Target.h>
#include <det/TargetConst.h>
#include <det/PSD.h>
#include <det/BPD.h>
#include <det/Beam.h>
#include <det/Target.h>
#include <det/Detector.h>
#include <det/BeamCounters.h>
#include <evt/EventHeader.h>
#include <evt/sim/Hit.h>
#endif

#ifndef ROOTLIB
#define ROOTLIB
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH3.h>
#include <TCutG.h>
#include <TProfile.h>
#include <TList.h>
#include <TString.h>
#include <TMath.h>
#include <utl/Vector.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#endif
using namespace std;
using namespace io;
using namespace utl;
using namespace evt;

//********************************************************************************************************
// TO RUN DO:

// 1. make executable with usual SHINE Makefile
// 2. create yourFile.txt with the list of the root files to analyze:
//
// yourProducedFile0.root
// yourProducedFile1.root
// yourProducedFile2.root
// ...
//
// 3. run as ./on_twiki yourFile.txt
//********************************************************************************************************

//***********************************************************************************************************************************************************
// Some explanation:
// In this production we kept events that are either T2 or events in which beam had inelastic interaction in the target. This helps us to estimate during the
// analysis T2 trigger inefficiencies. They are of two types:
//
// - some events don't have label T2, but they are good and beam had inelastic interaction in the target. Therefore, they are MISSED by T2 trigger
// (likely because some produced particle hit S4 counter)
//
// - some events have label T2, but there was NO beam inelastic interaction in the target. Therefore, they are FAKE T2 events (likely because
// beam somehow missed S4 counter: stopped before S4 or was deflected)
//
// Hence, during the analysis one has to search for events containing beam inelastic interaction in the target. This example code does this.
//
// Another thing is that there is no main vertex in SIM event in the production in beam-mode (note that for REC events everything is just normal). Consequently,
// one has to look for the final particles produced in the beam inelastic interaction in the target. However, some of those decay quickly (in the target),
// so it is needed to loop over daughter products of decayed particles from "first generation" in order to find the analogs to "eGeneratorFinal" particles.
// This example code does this. Also one has to be careful when analyzing ion reactions since there might be some remnants of nuclei, hence, a special selection on the PDG number is needed.
//
// Have fun!
//***********************************************************************************************************************************************************

vector<string> readFilelist(string filelist)
{
    ifstream in(filelist.c_str());
    vector<string> files;
    string line;
    while (getline(in, line)) {
        if (line == "//-1" || line == "#-1")
            break;
        else if (line.substr(0, 2) == "//" || line.substr(0, 1) == "#")
            continue;
        files.push_back(line);
    }
    return files;
}

int main(int argc, char* argv[])
{
    cout << "start" << endl;
    if (argc < 2) {
        cerr << "usage: "<< argv[0]<< " <filenames>" << endl;
        return 1;
    }
    
    // -- set chain of input files in a file
    vector<string> fileNames = readFilelist(argv[1]);// get list of file names from file
    
    EventFileChain eventFileChain(fileNames); // set event chain
    cout << "event file chain was set" << endl;
    Event event;
    
    // -- create histograms of interest
    TH1D* MultHist_sim = new TH1D("MultHist_sim", "entries;N", 30, -0.5, 29.5);
    TH1D* MultHist_rec = new TH1D("MultHist_rec", "entries;N", 30, -0.5, 29.5);
    
    TH1D* ptHist_sim = new TH1D("ptHist_sim", "entries;pt", 280, 0, 7);
    TH1D* ptHist_rec = new TH1D("ptHist_rec", "entries;pt", 280, 0, 7);
    TH1D* pxHist_sim = new TH1D("pxHist_sim", "entries;px", 280, 0, 7);
    TH1D* pxHist_rec = new TH1D("pxHist_rec", "entries;px", 280, 0, 7);
    TH1D* pyHist_sim = new TH1D("pyHist_sim", "entries;py", 280, 0, 7);
    TH1D* pyHist_rec = new TH1D("pyHist_rec", "entries;py", 280, 0, 7);
    TH1D* pzHist_sim = new TH1D("pzHist_sim", "entries;pz", 900, 0, 160);
    TH1D* pzHist_rec = new TH1D("pzHist_rec", "entries;pz", 900, 0, 160);
    
    TH1D* etaHist_sim = new TH1D("etaHist_sim", "entries;eta", 100, -5, 15);
    TH1D* etaHist_rec = new TH1D("etaHist_rec", "entries;eta", 100, -5, 15);
    
    TH1D* phiHist_sim = new TH1D("phiHist_sim", "entries;phi", 80, -4, 4);
    TH1D* phiHist_rec = new TH1D("phiHist_rec", "entries;phi", 80, -4, 4);
    
    TH2D* etaptHist_sim = new TH2D("etaptHist_sim", "entries;eta;Pt", 100, -5, 15, 50, 0, 5);
    TH2D* etaptHist_rec = new TH2D("etaptHist_rec", "entries;eta;Pt", 100, -5, 15, 50, 0, 5);
    
    TH1D* vertexTypeHist_sim_in_good_events = new TH1D("vertexTypeHist_sim_in_good_events", "entries;sim vertex type", 52, -1.5, 50.5);
    TH1D* vertexTypeHist_sim_in_target = new TH1D("vertexTypeHist_sim_in_target", "entries;sim vertex type", 52, -1.5, 50.5);
    
    TH1D* MultHist_sim_inelastic_noT2 = new TH1D("MultHist_sim_inelastic_noT2", "entries;N", 30, -0.5, 29.5);
    TH2D* etaptHist_sim_inelastic_noT2 = new TH2D("etaptHist_sim_inelastic_noT2", "entries;eta;Pt", 100, -5, 15, 50, 0, 5);
    TH1D* ptHist_sim_inelastic_noT2 = new TH1D("ptHist_sim_inelastic_noT2", "entries;pt", 280, 0, 7);
    TH1D* pHist_sim_inelastic_noT2 = new TH1D("pHist_sim_inelastic_noT2", "entries;p", 1000, 0, 200);
    TH1D* etaHist_sim_inelastic_noT2 = new TH1D("etaHist_sim_inelastic_noT2", "entries;eta", 100, -5, 15);
    TH1D* phiHist_sim_inelastic_noT2 = new TH1D("phiHist_sim_inelastic_noT2", "entries;phi", 80, -4, 4);
    TH1D* idHist_sim_inelastic_noT2 = new TH1D("idHist_sim_inelastic_noT2", "entries;pt", 20001, -10000.5, 10000.5);
    TH2D* beamXYpositionAtS4Z_T2_didntstopHist = new TH2D("beamXYpositionAtS4Z_T2_didntstopHist", "beam position at S4 Z;beam X AtS4Z;beam Y AtS4Z", 200, -5, 5, 200, -5, 5);
    TH2D* beamXYpositionAtS4Z_T2_stop_downstreamS4Hist = new TH2D("beamXYpositionAtS4Z_T2_stop_downstreamS4Hist", "beam position at S4 Z;beam X AtS4Z;beam Y AtS4Z", 200, -5, 5, 200, -5, 5);
    TH1D* beamStopedOutsideTarget = new TH1D("beamStopedOutsideTarget", "entries;beam stop vertex position Z", 22000, -1000, 10000);
    
    TH1D* hNumberOfBeamElasticVertices_fakeT2_stop_downstreamS4Hist = new TH1D("hNumberOfBeamElasticVertices_fakeT2_stop_downstreamS4Hist", "number of beam elastic vertices if it stopped downstream S4;number_of_beam_elastic_vertices;entries", 51, -0.5, 50.5);
    TH1D* hNumberOfBeamElasticVertices_fakeT2_didntstopHist = new TH1D("hNumberOfBeamElasticVertices_fakeT2_didntstopHist", "number of beam elastic vertices if it didn't stop;number_of_beam_elastic_vertices;entries", 51, -0.5, 50.5);
    
    TH1D* s4_energy_noT2 = new TH1D("s4_energy_noT2","S4 energy deposited, noT2 Trigger (S4 was Hit);S4energy;Entries",100000,0,10);
    
    TH1D* s4ADCinelasticnoT2 = new TH1D("s4ADCinelasticnoT2","S4 ADC, T2 Trigger (S4 Not Hit);ADC;Entries",10000,0,10000);
    TH1D* s4ADCnoT2 = new TH1D("s4ADCnoT2","S4 ADC, T2 Trigger (S4 Not Hit);ADC;Entries",10000,0,10000);
    
    // -- create multiplicity counters
    int mult_sim;
    int mult_rec;
    int mult_sim_inelastic_noT2;
    
    // -- create event counters
    int nOfT2events=0;
    int nOfT2events_ingood=0;
    int nOfBeamParticles=0;
    int nOfBeamParticlesInteractedInsideTarget=0;
    int nOfBeamParticlesInteractedInsideTargetStrongly=0;
    int nInelasticButNotT2=0;
    
    int onlyT2=0;
    int onlyInelastic=0;
    int T2andInelastic=0;
    
    int beamdidntstop=0;
    int nbeamstopedoutsidetarget=0;
    int nbeamstopedoutsidetargetDOWNSTREAM_s4=0;
    
    int allevents=0;
    
    // -- loop over events
    while (eventFileChain.Read(event) == eSuccess) {
        
        allevents++;
        if (!(allevents % 10000))
            cout << "\n ---> processing event # " << allevents << endl;
        
        const evt::EventHeader& eventHeader = event.GetEventHeader();
        int runNum = eventHeader.GetRunNumber();
        int eventNum = eventHeader.GetId();
        const utl::TimeStamp& timeStamp = eventHeader.GetTime();
        
        // -- set up to use det/target
        fwk::CentralConfig::GetInstance("bootstrap.xml"); // -- put bootstrap.xml with the proper GlobalKey in the same directory or specify the path to it here
        det::Detector::GetInstance().Update(timeStamp, runNum);
        det::Detector& detector = det::Detector::GetInstance();
        const det::Target& target = detector.GetTarget();
        const det::MagneticField& magneticField = detector.GetMagneticField();
        const det::MagneticFieldTracker tracker(magneticField,utl::Tracker::eTrackStepper);
        const double targetCenterZ = target.GetCenterPosition().GetZ();
        double s4ZPosition = detector.GetBeamCounters().GetCounter(det::BeamCounterConst::eS4).GetCenterPosition().GetZ();
        
        const SimEvent& simEvent = event.GetSimEvent();
        const RecEvent& recEvent = event.GetRecEvent();
        
        const raw::Trigger& trigger = event.GetRawEvent().GetBeam().GetTrigger();
        
        bool t2event=false;
        
        if (trigger.IsTrigger(det::TriggerConst::eT2, det::TriggerConst::ePrescaled))
        {
            nOfT2events++;
            t2event=true;
        }
        else
        {
            s4ADCnoT2->Fill(trigger.GetADC(det::TriggerConst::eS4));
            
            double energyDepositS4 = 0.0;
            for (SimEvent::HitIterator iter = event.GetSimEvent().HitsBegin(det::Const::eS4), end = event.GetSimEvent().HitsEnd(det::Const::eS4); iter != end; ++iter)
            {
                energyDepositS4 += iter->GetEnergyDeposit();
            }
            
            if(energyDepositS4>0)
            {
                s4_energy_noT2->Fill(energyDepositS4);
            }
        }
        
        mult_sim=0;
        mult_rec=0;
        mult_sim_inelastic_noT2=0;
        
        // -- loop over event sim vertex tracks
        
        // -- now we will select events that had beam inelastic interaction in the target
        
        vector<evt::sim::VertexTrack> FinalTracks;
        // -- it will be a vector with the analogs of "eGeneratorFinal" particles. Note, that we do not decay pi0 (the same way it was done in old EPOS production)
        
        for (list<evt::sim::VertexTrack>::const_iterator simTrackIter = simEvent.Begin<evt::sim::VertexTrack>(), simTrackEnd = simEvent.End<evt::sim::VertexTrack>();simTrackIter != simTrackEnd; ++simTrackIter)
        {
            const evt::sim::VertexTrack& beamVertexTrack = *simTrackIter;
            
            if(beamVertexTrack.GetType()==sim::VertexTrackConst::eBeam) // -- this is beam particle
            {
                nOfBeamParticles++;
                
                if(beamVertexTrack.HasStopVertex()) // -- beam has stop interaction
                {
                    evt::Index< evt::sim::Vertex > indexBeamStopVertex = beamVertexTrack.GetStopVertexIndex();
                    evt::sim::Vertex BeamStopVertex = simEvent.Get(indexBeamStopVertex);
                    const utl::Point& BeamStopVertexPosition = BeamStopVertex.GetPosition();
                    //cout<<"BeamStopVertex.GetType() = "<<BeamStopVertex.GetType()<<endl;
                    
                    if (target.IsIn(BeamStopVertexPosition)) // -- beam has stop interaction is inside the target
                    {
                        nOfBeamParticlesInteractedInsideTarget++;
                        //if(beamVertexTrack.GetNumberOfElasticVertices()==0) // -- beam didn't have elastic vertices before the inelastic interaction - we didn't discuss yet whether we will need this check (probably effect is negligible)
                        //{
                        
                        //cout<<"BeamStopVertex.GetType() = "<<BeamStopVertex.GetType()<<endl;
                        if(BeamStopVertex.GetType()==sim::VertexConst::eHadronicInelastic) // -- beam has stop interaction of eHadronicInelastic type
                        {
                            nOfBeamParticlesInteractedInsideTargetStrongly++;
                            
                            if (!(nOfBeamParticlesInteractedInsideTargetStrongly % 5000))
                                cout << "\n ---> processing GOOD event # " << nOfBeamParticlesInteractedInsideTargetStrongly << endl;
                            
                            if (trigger.IsTrigger(det::TriggerConst::eT2, det::TriggerConst::ePrescaled))
                            {
                                nOfT2events_ingood++;
                                T2andInelastic++;
                            }
                            else
                            {
                                onlyInelastic++;
                                s4ADCinelasticnoT2->Fill(trigger.GetADC(det::TriggerConst::eS4));
                            }
                            
                            for (std::list<evt::sim::Vertex>::const_iterator vIter =
                                 simEvent.Begin<evt::sim::Vertex>();
                                 vIter != simEvent.End<evt::sim::Vertex>(); ++vIter)
                            {
                                evt::sim::Vertex EventVertex_in_good = simEvent.Get(vIter->GetIndex());
                                vertexTypeHist_sim_in_good_events->Fill(EventVertex_in_good.GetType());
                                
                                if(target.IsIn(EventVertex_in_good.GetPosition()))
                                {
                                    vertexTypeHist_sim_in_target->Fill(EventVertex_in_good.GetType());
                                }
                            }
                            
                            
                            for (evt::sim::VertexTrackIndexIterator simvtxTrackIter = BeamStopVertex.DaughterTracksBegin(); simvtxTrackIter != BeamStopVertex.DaughterTracksEnd(); ++simvtxTrackIter)
                            {
                                const evt::sim::VertexTrack& ProbablyPrimaryTrack = simEvent.Get(*simvtxTrackIter);
                                
                                    if((ProbablyPrimaryTrack.HasStopVertex()) && (simEvent.Get(ProbablyPrimaryTrack.GetStopVertexIndex()).GetType()==sim::VertexConst::eDecay)&&(target.IsIn(simEvent.Get(ProbablyPrimaryTrack.GetStopVertexIndex()).GetPosition()))&&(ProbablyPrimaryTrack.GetParticleId()!=111))
                                    {
                                        for (evt::sim::VertexTrackIndexIterator simvtxTrackIter4 = simEvent.Get(ProbablyPrimaryTrack.GetStopVertexIndex()).DaughterTracksBegin(); simvtxTrackIter4 != simEvent.Get(ProbablyPrimaryTrack.GetStopVertexIndex()).DaughterTracksEnd(); ++simvtxTrackIter4)
                                        {
                                            const evt::sim::VertexTrack& ProbablyPrimaryTrackN = simEvent.Get(*simvtxTrackIter4);
                                            if((ProbablyPrimaryTrackN.HasStopVertex()) && (simEvent.Get(ProbablyPrimaryTrackN.GetStopVertexIndex()).GetType()==sim::VertexConst::eDecay)&&(target.IsIn(simEvent.Get(ProbablyPrimaryTrackN.GetStopVertexIndex()).GetPosition()))&&(ProbablyPrimaryTrackN.GetParticleId()!=111))
                                            {
                                                for (evt::sim::VertexTrackIndexIterator simvtxTrackIter5 = simEvent.Get(ProbablyPrimaryTrackN.GetStopVertexIndex()).DaughterTracksBegin(); simvtxTrackIter5 != simEvent.Get(ProbablyPrimaryTrackN.GetStopVertexIndex()).DaughterTracksEnd(); ++simvtxTrackIter5)
                                                {
                                                    const evt::sim::VertexTrack& ProbablyPrimaryTrackNN = simEvent.Get(*simvtxTrackIter5);
                                                    if((ProbablyPrimaryTrackNN.HasStopVertex()) && (simEvent.Get(ProbablyPrimaryTrackNN.GetStopVertexIndex()).GetType()==sim::VertexConst::eDecay)&&(target.IsIn(simEvent.Get(ProbablyPrimaryTrackNN.GetStopVertexIndex()).GetPosition()))&&(ProbablyPrimaryTrackNN.GetParticleId()!=111))
                                                    {
                                                        for (evt::sim::VertexTrackIndexIterator simvtxTrackIter6 = simEvent.Get(ProbablyPrimaryTrackNN.GetStopVertexIndex()).DaughterTracksBegin(); simvtxTrackIter6 != simEvent.Get(ProbablyPrimaryTrackNN.GetStopVertexIndex()).DaughterTracksEnd(); ++simvtxTrackIter6)
                                                        {
                                                            const evt::sim::VertexTrack& ProbablyPrimaryTrackNNN = simEvent.Get(*simvtxTrackIter6);
                                                            if((ProbablyPrimaryTrackNNN.HasStopVertex()) && (simEvent.Get(ProbablyPrimaryTrackNNN.GetStopVertexIndex()).GetType()==sim::VertexConst::eDecay)&&(target.IsIn(simEvent.Get(ProbablyPrimaryTrackNNN.GetStopVertexIndex()).GetPosition()))&&(ProbablyPrimaryTrackNNN.GetParticleId()!=111))
                                                            {
                                                                for (evt::sim::VertexTrackIndexIterator simvtxTrackIter7 = simEvent.Get(ProbablyPrimaryTrackNNN.GetStopVertexIndex()).DaughterTracksBegin(); simvtxTrackIter7 != simEvent.Get(ProbablyPrimaryTrackNNN.GetStopVertexIndex()).DaughterTracksEnd(); ++simvtxTrackIter7)
                                                                {
                                                                    const evt::sim::VertexTrack& ProbablyPrimaryTrackNNNN = simEvent.Get(*simvtxTrackIter7);
                                                                    
                                                                    if((ProbablyPrimaryTrackNNNN.HasStopVertex()) && (simEvent.Get(ProbablyPrimaryTrackNNNN.GetStopVertexIndex()).GetType()==sim::VertexConst::eDecay)&&(target.IsIn(simEvent.Get(ProbablyPrimaryTrackNNNN.GetStopVertexIndex()).GetPosition()))&&(ProbablyPrimaryTrackNNNN.GetParticleId()!=111))
                                                                    {
                                                                        for (evt::sim::VertexTrackIndexIterator simvtxTrackIter8 = simEvent.Get(ProbablyPrimaryTrackNNNN.GetStopVertexIndex()).DaughterTracksBegin(); simvtxTrackIter8 != simEvent.Get(ProbablyPrimaryTrackNNNN.GetStopVertexIndex()).DaughterTracksEnd(); ++simvtxTrackIter8)
                                                                        {
                                                                            const evt::sim::VertexTrack& ProbablyPrimaryTrackNNNNN = simEvent.Get(*simvtxTrackIter8);
                                                                            if((ProbablyPrimaryTrackNNNNN.HasStopVertex()) && (simEvent.Get(ProbablyPrimaryTrackNNNNN.GetStopVertexIndex()).GetType()==sim::VertexConst::eDecay)&&(target.IsIn(simEvent.Get(ProbablyPrimaryTrackNNNNN.GetStopVertexIndex()).GetPosition()))&&(ProbablyPrimaryTrackNNNNN.GetParticleId()!=111))
                                                                            {
                                                                                for (evt::sim::VertexTrackIndexIterator simvtxTrackIter9 = simEvent.Get(ProbablyPrimaryTrackNNNNN.GetStopVertexIndex()).DaughterTracksBegin(); simvtxTrackIter9 != simEvent.Get(ProbablyPrimaryTrackNNNNN.GetStopVertexIndex()).DaughterTracksEnd(); ++simvtxTrackIter9)
                                                                                {
                                                                                    const evt::sim::VertexTrack& ProbablyPrimaryTrackNNNNNN = simEvent.Get(*simvtxTrackIter9);
                                                                                    if((ProbablyPrimaryTrackNNNNNN.HasStopVertex()) && (simEvent.Get(ProbablyPrimaryTrackNNNNNN.GetStopVertexIndex()).GetType()==sim::VertexConst::eDecay)&&(target.IsIn(simEvent.Get(ProbablyPrimaryTrackNNNNNN.GetStopVertexIndex()).GetPosition()))&&(ProbablyPrimaryTrackNNNNNN.GetParticleId()!=111))
                                                                                    {
                                                                                        for (evt::sim::VertexTrackIndexIterator simvtxTrackIter10 = simEvent.Get(ProbablyPrimaryTrackNNNNNN.GetStopVertexIndex()).DaughterTracksBegin(); simvtxTrackIter10 != simEvent.Get(ProbablyPrimaryTrackNNNNNN.GetStopVertexIndex()).DaughterTracksEnd(); ++simvtxTrackIter10)
                                                                                        {
                                                                                            const evt::sim::VertexTrack& ProbablyPrimaryTrack7N = simEvent.Get(*simvtxTrackIter10);
                                                                                            if((ProbablyPrimaryTrack7N.HasStopVertex()) && (simEvent.Get(ProbablyPrimaryTrack7N.GetStopVertexIndex()).GetType()==sim::VertexConst::eDecay)&&(target.IsIn(simEvent.Get(ProbablyPrimaryTrack7N.GetStopVertexIndex()).GetPosition()))&&(ProbablyPrimaryTrack7N.GetParticleId()!=111))
                                                                                            {
                                                                                                for (evt::sim::VertexTrackIndexIterator simvtxTrackIter11 = simEvent.Get(ProbablyPrimaryTrack7N.GetStopVertexIndex()).DaughterTracksBegin(); simvtxTrackIter11 != simEvent.Get(ProbablyPrimaryTrack7N.GetStopVertexIndex()).DaughterTracksEnd(); ++simvtxTrackIter11)
                                                                                                {
                                                                                                    const evt::sim::VertexTrack& ProbablyPrimaryTrack8N = simEvent.Get(*simvtxTrackIter10);
                                                                                                    FinalTracks.push_back(ProbablyPrimaryTrack8N);
                                                                                                }
                                                                                            }
                                                                                            else
                                                                                            {
                                                                                                FinalTracks.push_back(ProbablyPrimaryTrack7N);
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                    else
                                                                                    {
                                                                                        FinalTracks.push_back(ProbablyPrimaryTrackNNNNNN);
                                                                                    }
                                                                                }
                                                                            }
                                                                            else
                                                                            {
                                                                                FinalTracks.push_back(ProbablyPrimaryTrackNNNNN);
                                                                            }
                                                                        }
                                                                    }
                                                                    else
                                                                    {
                                                                        FinalTracks.push_back(ProbablyPrimaryTrackNNNN);
                                                                    }
                                                                }
                                                            }
                                                            else
                                                            {
                                                                FinalTracks.push_back(ProbablyPrimaryTrackNNN);
                                                            }
                                                        }
                                                    }
                                                    else
                                                    {
                                                        FinalTracks.push_back(ProbablyPrimaryTrackNN);
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                FinalTracks.push_back(ProbablyPrimaryTrackN);
                                            }
                                        }
                                    }
                                    else {FinalTracks.push_back(ProbablyPrimaryTrack);}
                                } // -- loop over beam produced particles is over, we selected particles for analysis in FinalTracks
                            // -- this is assumed to be equivalent to the check (primarySimTrack == eGeneratorFinal)
                            
                            int size = FinalTracks.size();
                            
                            for(int i=0;i<size;i++)
                            {
                                const evt::sim::VertexTrack& vtxTrack = FinalTracks.at(i);
                                // -- this is assumed to be equivalent to the check (primarySimTrack == eGeneratorFinal)
                                
                                // -- add more sim track quality selection cuts here ....
                                mult_sim++;
                                
                                double_t p, pX, pY, pZ, pT, phi, eta;
                                Vector vtxMomentumGood = vtxTrack.GetMomentum();
                                pX = vtxMomentumGood.GetX();
                                pY = vtxMomentumGood.GetY();
                                pZ = vtxMomentumGood.GetZ();
                                p = vtxMomentumGood.GetMag();
                                phi = TMath::ATan2(pY, pX);
                                pT = sqrt(pX*pX + pY*pY);
                                eta = -0.5*TMath::Log((p - pZ) / (p + pZ));
                                
                                etaptHist_sim->Fill(eta, pT);
                                ptHist_sim->Fill(pT);
                                pxHist_sim->Fill(pX);
                                pyHist_sim->Fill(pY);
                                pzHist_sim->Fill(pZ);
                                etaHist_sim->Fill(eta);
                                phiHist_sim->Fill(phi);
                                
                                // -- see what happend in this event if it is not labeled as T2
                                if (!trigger.IsTrigger(det::TriggerConst::eT2, det::TriggerConst::ePrescaled))
                                {
                                    for (SimEvent::HitIterator iter = event.GetSimEvent().HitsBegin(det::Const::eS4), end = event.GetSimEvent().HitsEnd(det::Const::eS4); iter != end; ++iter)
                                    {
                                        double energyDepositS4 = 0.0;
                                        energyDepositS4 = iter->GetEnergyDeposit();
                                        const evt::Index<evt::sim::Hit>& hitIndex = iter->GetIndex();
                                        
                                        if (energyDepositS4 > 0.0) {
                                            
                                            const evt::sim::Hit& hitTouched = simEvent.Get(hitIndex);
                                            
                                            const evt::sim::VertexTrack& vertexTrackInS4 = simEvent.Get<sim::VertexTrack>(hitTouched.GetVertexTrackIndex());
                                            idHist_sim_inelastic_noT2->Fill(vertexTrackInS4.GetParticleId());
                                            
                                            double_t pS4, pXS4, pYS4, pZS4, pTS4, phiS4, etaS4;
                                            Vector vtxMomentumS4 = vertexTrackInS4.GetMomentum();
                                            pXS4 = vtxMomentumS4.GetX();
                                            pYS4 = vtxMomentumS4.GetY();
                                            pZS4 = vtxMomentumS4.GetZ();
                                            pS4 = vtxMomentumS4.GetMag();
                                            phiS4 = TMath::ATan2(pYS4, pXS4);
                                            pTS4 = sqrt(pXS4*pXS4 + pYS4*pYS4);
                                            etaS4 = -0.5*TMath::Log((pS4 - pZS4) / (pS4 + pZS4));
                                            
                                            ptHist_sim_inelastic_noT2->Fill(pTS4);
                                            pHist_sim_inelastic_noT2->Fill(pS4);
                                            etaHist_sim_inelastic_noT2->Fill(etaS4);
                                            phiHist_sim_inelastic_noT2->Fill(phiS4);
                                            etaptHist_sim_inelastic_noT2->Fill(etaS4, pTS4);
                                            
                                            mult_sim_inelastic_noT2++;
                                            
                                        }
                                    }
                                }
                                
                                // -- fill other histograms...
                                
                            } // -- end of the loop over event eGeneratorFinal particles
                            
                            MultHist_sim->Fill(mult_sim);
                            
                            if (!trigger.IsTrigger(det::TriggerConst::eT2, det::TriggerConst::ePrescaled))
                            {
                                MultHist_sim_inelastic_noT2->Fill(mult_sim_inelastic_noT2);
                                nInelasticButNotT2++;
                            }
                        } // -- eHadronicElastic stop vertex
                        else
                        {
                            if(t2event)
                            {
                                onlyT2++;
                            }
                        }
                        //} // -- beam doesn't have elastic interactions before
                    } // -- stop vertex is inside target
                    else
                    {
                        if(t2event)
                        {
                            onlyT2++;
                            nbeamstopedoutsidetarget++;
                            beamStopedOutsideTarget->Fill(BeamStopVertexPosition.GetZ());
                            if(BeamStopVertexPosition.GetZ()>-210) // -211 is a Z coordinate of S4 counter that defines the T2 trigger
                            {
                                nbeamstopedoutsidetargetDOWNSTREAM_s4++;
                                const sim::Beam& beamSim = simEvent.GetBeam();
                                const utl::Point beamPositionAtTarget(beamSim.GetOrigin().GetX(),
                                                                      beamSim.GetOrigin().GetY(),
                                                                      beamSim.GetOrigin().GetZ());
                                
                                const utl::Vector& beamMomentum = beamSim.GetMomentum();
                                const int beamCharge = 1; // For protons.
                                
                                utl::Point extrapolatedPosition;
                                utl::Vector extrapolatedMomentum;
                                
                                tracker.TrackToZ(s4ZPosition,
                                                 beamCharge,
                                                 beamPositionAtTarget,
                                                 beamMomentum,
                                                 extrapolatedPosition,
                                                 extrapolatedMomentum);
                                
                                beamXYpositionAtS4Z_T2_stop_downstreamS4Hist->Fill(extrapolatedPosition.GetX(),extrapolatedPosition.GetY());
                                hNumberOfBeamElasticVertices_fakeT2_stop_downstreamS4Hist->Fill(beamVertexTrack.GetNumberOfElasticVertices());
                                
                            }
                        }
                    }
                } // -- beam has stop vertex
                else
                {
                    if(t2event)
                    {
                        onlyT2++;
                        beamdidntstop++;
                        
                        const sim::Beam& beamSim = simEvent.GetBeam();
                        const utl::Point beamPositionAtTarget(beamSim.GetOrigin().GetX(),
                                                              beamSim.GetOrigin().GetY(),
                                                              beamSim.GetOrigin().GetZ());
                        
                        const utl::Vector& beamMomentum = beamSim.GetMomentum();
                        const int beamCharge = 1; // For protons.
                        
                        utl::Point extrapolatedPosition;
                        utl::Vector extrapolatedMomentum;
                        
                        tracker.TrackToZ(s4ZPosition,
                                         beamCharge,
                                         beamPositionAtTarget,
                                         beamMomentum,
                                         extrapolatedPosition,
                                         extrapolatedMomentum);
                        
                        beamXYpositionAtS4Z_T2_didntstopHist->Fill(extrapolatedPosition.GetX(),extrapolatedPosition.GetY());
                        hNumberOfBeamElasticVertices_fakeT2_didntstopHist->Fill(beamVertexTrack.GetNumberOfElasticVertices());
                        
                    }
                }
            } // -- this is beam
        } // -- end of the loop over all sim event vertex tracks
        
        // -- loop over rec tracks, this is a standart procedure
        
                if(recEvent.HasMainVertex())
                {
                    const evt::rec::Vertex &mainVertex = recEvent.GetMainVertex();
        
                    for (evt::rec::VertexTrackIndexIterator vtxTrackIter = mainVertex.DaughterTracksBegin(); vtxTrackIter != mainVertex.DaughterTracksEnd(); ++vtxTrackIter)
                    {
                        const evt::rec::VertexTrack& recVertexTrack = recEvent.Get(*vtxTrackIter);
        
                        // -- add more rec quality selection cuts here ....
                        mult_rec++;
        
                        const Vector& momentum = recVertexTrack.GetMomentum();
                        double_t p, pX, pY, pZ, pT, phi, eta;
                        pX = momentum.GetX();
                        pY = momentum.GetY();
                        pZ = momentum.GetZ();
                        p = momentum.GetMag();
                        phi = TMath::ATan2(pY, pX);
                        pT = sqrt(pX*pX + pY*pY);
                        eta = -0.5*TMath::Log((p - pZ) / (p + pZ));
        
                        etaptHist_rec->Fill(eta, pT);
                        ptHist_rec->Fill(pT);
                        pxHist_rec->Fill(pX);
                        pyHist_rec->Fill(pY);
                        pzHist_rec->Fill(pZ);
                        etaHist_rec->Fill(eta);
                        phiHist_rec->Fill(phi);
        
                        // -- fill other histograms....
                    }
                }
        
        MultHist_rec->Fill(mult_rec);
        
    } // -- end of event loop
    
    // -- create output file
    TFile outFile("histo.root", "RECREATE");
    
    // -- fill output file
    MultHist_sim->Write();
    MultHist_rec->Write();
    pxHist_sim->Write();
    pxHist_rec->Write();
    pyHist_sim->Write();
    pyHist_rec->Write();
    pzHist_sim->Write();
    pzHist_rec->Write();
    ptHist_sim->Write();
    ptHist_rec->Write();
    etaHist_sim->Write();
    etaHist_rec->Write();
    etaptHist_sim->Write();
    etaptHist_rec->Write();
    phiHist_sim->Write();
    phiHist_rec->Write();
    vertexTypeHist_sim_in_good_events->Write();
    vertexTypeHist_sim_in_target->Write();
    
    MultHist_sim_inelastic_noT2->Write();
    etaptHist_sim_inelastic_noT2->Write();
    ptHist_sim_inelastic_noT2->Write();
    pHist_sim_inelastic_noT2->Write();
    etaHist_sim_inelastic_noT2->Write();
    phiHist_sim_inelastic_noT2->Write();
    idHist_sim_inelastic_noT2->Write();
    
    beamStopedOutsideTarget->Write();
    beamXYpositionAtS4Z_T2_didntstopHist->Write();
    beamXYpositionAtS4Z_T2_stop_downstreamS4Hist->Write();
    hNumberOfBeamElasticVertices_fakeT2_stop_downstreamS4Hist->Write();
    hNumberOfBeamElasticVertices_fakeT2_didntstopHist->Write();
    
    s4ADCinelasticnoT2->Write();
    s4ADCnoT2->Write();
    
    s4_energy_noT2->Write();

    TAxis *xAxis(idHist_sim_inelastic_noT2->GetXaxis());
    int nbinsX(xAxis->GetNbins());
    double sum(0);
    
    for (int binx = 1; binx < nbinsX + 1; ++binx)
    {
        double content = idHist_sim_inelastic_noT2->GetBinContent(binx);
        sum+=content;
    }
    cout<<"sum = "<<sum<<endl;
    for (int binx = 1; binx < nbinsX + 1; ++binx)
    {
        double content = idHist_sim_inelastic_noT2->GetBinContent(binx);
        
        if(content!=0)
        {
            cout<<"id = "<<xAxis->GetBinCenter(binx)<<", "<<(content/sum)*100<<" %"<<endl;
        }
        
    }
    
    outFile.Write();
    outFile.Close();
    
    double dnOfBeamParticles = nOfBeamParticles;
    double dnOfBeamParticlesInteractedInsideTargetStrongly = nOfBeamParticlesInteractedInsideTargetStrongly;
    double dnOfT2events = nOfT2events;
    double dnOfT2events_ingood = nOfT2events_ingood;
    
    
    // -- see the event counters and their ratios
    cout<<"nOfT2events = "<<nOfT2events<<endl;
    cout<<"nOfBeamParticles = "<<nOfBeamParticles<<endl; // -- NOTE that in the beam-mode we decided to keep either T2 labeled events or events with beam mode inelastic interaction
    cout<<"nOfBeamParticlesInteractedInsideTarget = "<<nOfBeamParticlesInteractedInsideTarget<<endl;
    cout<<"nOfBeamParticlesInteractedInsideTargetStrongly = "<<nOfBeamParticlesInteractedInsideTargetStrongly<<endl;
    
    cout<<"interaction probability = nOfBeamParticlesInteractedInsideTargetStrongly/nOfBeamParticles = "<< dnOfBeamParticlesInteractedInsideTargetStrongly/dnOfBeamParticles<<endl;
    
    cout<<"T2 trigger efficiency = nOfT2events/nOfBeamParticlesInteractedInsideTargetStrongly = "<< dnOfT2events/dnOfBeamParticlesInteractedInsideTargetStrongly<<endl;
    
    cout<<"T2 trigger efficiency* = nOfT2events_ingood/nOfBeamParticlesInteractedInsideTargetStrongly = "<< dnOfT2events_ingood/dnOfBeamParticlesInteractedInsideTargetStrongly<<endl;

    cout<<" "<<endl;
    cout<<"onlyT2 (fakes) = "<<onlyT2<<endl;
    cout<<"beam didn't stop (but it's T2) = "<<beamdidntstop<<endl;
    cout<<"beam stoped outside the target (but it's T2) = "<<nbeamstopedoutsidetarget<<endl;
    cout<<"nbeamstopedoutsidetargetDOWNSTREAM_s4 (but it's T2) = "<<nbeamstopedoutsidetargetDOWNSTREAM_s4<<endl;
    cout<<" "<<endl;
    cout<<"onlyInelastic (missed by T2) = "<<onlyInelastic<<endl;
    cout<<"T2andInelastic = "<<T2andInelastic<<endl;
    cout<<"T2 or Inelastic = "<<T2andInelastic+onlyT2+onlyInelastic<<endl;
    cout<<" "<<endl;
    cout<<"all T2 = "<< onlyT2+T2andInelastic<<endl;
    cout<<"all inelastic = "<< onlyInelastic+T2andInelastic<<endl;
    cout<<" "<<endl;
    return 0;
}
