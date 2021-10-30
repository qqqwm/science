#include <evt/Event.h>
#include <evt/SimEvent.h>
#include <evt/RecEvent.h>
#include <io/EventFileChain.h>
#include <io/IoCodes.h>
#include <utl/PDGParticleIds.h>
#include <evt/Event.h>
#include <evt/RecEvent.h>
#include <evt/rec/RecEventConst.h>
#include <fwk/CentralConfig.h>
#include <det/Detector.h>
#include <utl/ShineUnits.h>
#include <evt/rec/Trigger.h>
#include <utl/Vector.h>
#include <utl/MathConst.h>
#include <utl/PhysicalConst.h>
#include <utl/PDGParticleIds.h>
#include <utl/DatabasePDG.h>
#include <io/EventFileChain.h>

#include <TTree.h>
#include <TFile.h>
#include <TH1I.h>
#include <TMath.h>
#include <TRandom.h>
#include <TH2I.h>

#include <iostream>
#include <cmath>

using namespace std;
using namespace io;
using namespace utl;
using namespace evt;

Vector TransformToNewBasis(const Vector& nx, const Vector& ny, const Vector& nz, const Vector& p)
{
    return Vector(nx.Dot(p), ny.Dot(p), nz.Dot(p));
}

int main(int argc, char* argv[])
{
    if (argc < 3) {
        cerr << "usage: " << argv[0] << "<output filename> <filenames...>" << endl;
        return 1;
    }
    const DatabasePDG& particleDB = DatabasePDG::GetInstance();
    const double protonMass = particleDB.GetParticle(ParticleConst::eProton).Mass() * GeV;
    const double lambdaMass = particleDB.GetParticle(ParticleConst::eLambda).Mass() * GeV;
    //const double pionMass = particleDB.GetParticle(ParticleConst::ePiMinus).Mass()*GeV;

    // helper counters
    int nEvents = 0;
    int nHasMainVertex = 0;
    int nHasStopVertex = 0;
    int nVertexTrackGeneratorFinal = 0;
    int nHasDecayType = 0;
    int nHasTwoDaughterTracks = 0;
    int nLambdaDecay = 0;
    int nFinal = 0;
    Double_t cos_x, cos_y, cos_z, DeltaZ, y, p_full_lambda, cosPhi;
    Double_t alpha, p_t, p_t_lambda, phi_lambda;
    UInt_t nHits11, nHits12, nHits21, nHits22;
    UInt_t XYTargetPass;
    TTree * t = new TTree("Lambdas", "Lambdas pp_158 MC");
    t->Branch("cos_x", &cos_x, "cos_x/D");
    t->Branch("cos_y", &cos_y, "cos_y/D");
    t->Branch("cos_z", &cos_z, "cos_z/D");
    t->Branch("nHits11", &nHits11, "nHits11/i");
    t->Branch("nHits12", &nHits12, "nHits12/i");
    t->Branch("nHits21", &nHits21, "nHits21/i");
    t->Branch("nHits22", &nHits22, "nHits22/i");
    t->Branch("XYTargetPass", &XYTargetPass, "XYTargetPass/i");
    t->Branch("DeltaZ", &DeltaZ, "DeltaZ/D");
    t->Branch("y", &y, "y/D");
    t->Branch("alpha", &alpha, "alpha/D");
    t->Branch("p_t", &p_t, "p_t/D");
    t->Branch("p_t_lambda", &p_t_lambda, "p_t_lambda/D");
    t->Branch("p_full_lambda", &p_full_lambda, "p_full_lambda/D");
    t->Branch("phi_lambda", &phi_lambda, "phi_lambda/D");
    t->Branch("cosPhi", &cosPhi, "cosPhi/D");
    // ---  set chain of input files
    vector<string> fileNames;
    fileNames.reserve(argc - 2);
    for (int i = 2; i < argc; ++i)
    {
        string s("root://eospublic.cern.ch/");
        s += argv[i];
        fileNames.push_back(s);
    }
    EventFileChain eventFileChain(fileNames);
    eventFileChain.SetParsingOptions(false);
    eventFileChain.StopOnOpenProblem(false);
    // --- loop over events
    Event event;
    while (eventFileChain.Read(event) == eSuccess) {
        ++nEvents;
        if (!(nEvents % 1000))
            cout << "\n ---> processing event # " << nEvents << endl;
        using namespace evt::sim;
        const SimEvent& simEvent = event.GetSimEvent();
        if (!simEvent.HasMainVertex()) {
            cout << "event without simulated main vertex?? skip!" << endl;
            continue;
        }
        ++nHasMainVertex;
        const Beam& beam = simEvent.GetBeam();
        auto p_in = beam.GetMomentum();
        for (list<VertexTrack>::const_iterator
                simTrackIter = simEvent.Begin<VertexTrack>(),
                simTrackEnd = simEvent.End<VertexTrack>();
                simTrackIter != simTrackEnd; ++simTrackIter) {
            const VertexTrack& simTrack = *simTrackIter;
            if (simTrack.GetType() != VertexTrackConst::eGeneratorFinal)
                continue;
            ++nVertexTrackGeneratorFinal;
            if (!simTrack.HasStopVertex())
                continue;
            ++nHasStopVertex;
            const Vertex& stopVertex = simEvent.Get(simTrack.GetStopVertexIndex());
            //is the stop vertex a V0 vertex?
            if (stopVertex.GetType() != VertexConst::eDecay)
                continue;
            ++nHasDecayType;
            const auto nofdt = stopVertex.GetNumberOfDaughterTracks();
            if (nofdt != 2)
                continue;
            ++nHasTwoDaughterTracks;
            if (simTrack.GetParticleId() != ParticleConst::eLambda)
                continue;
            ++nLambdaDecay;

            VertexTrackIndexIterator decayTrackIter = stopVertex.DaughterTracksBegin();
            const VertexTrack& vtxTrack1 = simEvent.Get(*(decayTrackIter));
            const VertexTrack& vtxTrack2 = simEvent.Get(*(++decayTrackIter));

            // different charges ?
            const int charge1 = vtxTrack1.GetCharge();
            const int charge2 = vtxTrack2.GetCharge();
            //cout << "ch1 " << charge1 << " ch2 " << charge2 << " pid1 " << pid1 << " pid2 " << pid2 << endl;
            if (charge1 * charge2 != -1)
                continue;
            const int pid1 = vtxTrack1.GetParticleId();
            const int pid2 = vtxTrack2.GetParticleId();
            Vector protonMomentum, pionMomentum;
            if ((pid1 == ParticleConst::eProton) && (pid2 == ParticleConst::ePiMinus)) {
                protonMomentum = vtxTrack1.GetMomentum();
                pionMomentum = vtxTrack2.GetMomentum();
            }
            else if ((pid2 == ParticleConst::eProton) && (pid1 == ParticleConst::ePiMinus)) {
                protonMomentum = vtxTrack2.GetMomentum();
                pionMomentum = vtxTrack1.GetMomentum();
            }
            else {
                cout << "Not L-> p + pi^- decay\n";
                continue;
            }
            const Vector& p_lambda = simTrack.GetMomentum();
            p_full_lambda = p_lambda.GetMag();
            p_t_lambda = p_lambda.GetRho();
            phi_lambda = p_lambda.GetPhi();
            const double e_lambda = sqrt(p_lambda.GetMag2() + lambdaMass * lambdaMass);
            const Vector nz = Normalized(p_lambda);
            const Vector nx = Normalized(Cross(p_lambda, p_in));
            const Vector ny = Normalized(Cross(p_lambda, Cross(p_lambda, p_in)));

            const double p1l = protonMomentum.Dot(p_lambda) / p_lambda.GetMag();
            const double p2l = pionMomentum.Dot(p_lambda) / p_lambda.GetMag();
            p_t = protonMomentum.Cross(p_lambda).GetMag() / p_lambda.GetMag() / GeV;
            alpha = (p1l - p2l) / (p1l + p2l);

            const Vector newpproton = TransformToNewBasis(nx, ny, nz, protonMomentum);

            const double proton_e = sqrt(protonMomentum.GetMag2() + protonMass * protonMass);
            //next - lorentz boost in z direction
            const double boostedproton_x = newpproton.GetX();
            const double boostedproton_y = newpproton.GetY();
            const double boostedproton_z = (e_lambda * newpproton.GetZ() - p_lambda.GetMag() * proton_e ) / lambdaMass;
            const double boostedproton_p = sqrt(boostedproton_x * boostedproton_x + boostedproton_y * boostedproton_y + boostedproton_z * boostedproton_z);
            //const double boostedproton_e = (e_lambda * proton_e - p_lambda.GetMag() * newpproton.GetZ() ) / lambdaMass;
            cos_x = boostedproton_x / boostedproton_p;
            cos_y = boostedproton_y / boostedproton_p;
            cos_z = boostedproton_z / boostedproton_p;

            double vertexZ = simEvent.GetMainVertex().GetPosition().GetZ() / cm;
            double zv0  = stopVertex.GetPosition().GetZ() / cm;
            DeltaZ = zv0 - vertexZ;
            //double y = 0.5 * TMath::Log((finalStateEnergy + p_lambda.GetZ()) / (finalStateEnergy - p_lambda.GetZ())) - 2.909878177287;
            y = atanh(p_lambda.GetZ() / e_lambda) - 2.909878177287;
            ++nFinal;

            const Plane* vertexPlane = new Plane(simEvent.GetMainVertex().GetPosition(), Vector(0, 0, 1));
			const Line* v0Trajectory = new Line(stopVertex.GetPosition(), protonMomentum + pionMomentum);
			const Point v0AtVertexPlane = GeometryUtilities::Intersection(*vertexPlane, *v0Trajectory);
			const Vector mainToVertexPoint = v0AtVertexPlane - simEvent.GetMainVertex().GetPosition();
			auto targX =  mainToVertexPoint.GetX(),
				targY =  mainToVertexPoint.GetY();
			
			if (pow(0.5 * targX, 2) + pow(targY, 2) >= 1.)
				XYTargetPass = 0;
			else
				XYTargetPass = 1;

            nHits11 = vtxTrack1.GetNumberOfHits(evt::rec::TrackConst::eVTPC1);
            nHits12 = vtxTrack1.GetNumberOfHits(evt::rec::TrackConst::eVTPC2);
            nHits21 = vtxTrack2.GetNumberOfHits(evt::rec::TrackConst::eVTPC1);
            nHits22 = vtxTrack2.GetNumberOfHits(evt::rec::TrackConst::eVTPC2);
            /*
            nHitsPass = 1;
			if ((nHits11 < 15) && (nHits12 < 15)) {
				nHitsPass = 0;
			}
			if ((nHits21 < 15) && (nHits22 < 15)) {
				nHitsPass = 0;
			}
			*/
            Vector y(0, 1, 0); //y axis
			Vector norm = Normalized(protonMomentum.Cross(pionMomentum)); //normal to decay plane
			Vector y_prim = Normalized(p_lambda.Cross(y.Cross(p_lambda))); //perpendicular to V0-momentum, on the plane spanned by y and p
			cosPhi = fabs(norm.Dot(y_prim));

            t->Fill();
        }
    }
    // --- write histogram of correction factors to ROOT file
    TFile outFile(Form("y-sel-sim.%s.root", argv[1]), "RECREATE", "", 0);
    outFile.cd();
    t->Write();
    cout << "Event quality cuts: " << endl;
    cout << "nEvents:\t" << nEvents << endl;
    cout << "nHasMainVertex:\t" << nHasMainVertex << endl;
    cout << "nVertexTrackGeneratorFinal:\t" << nVertexTrackGeneratorFinal << endl;
    cout << "nHasStopVertex:\t" << nHasStopVertex << endl;
    cout << "nHasDecayType:\t" << nHasDecayType << endl;
    cout << "nHasTwoDaughterTracks:\t" << nHasTwoDaughterTracks << endl;
    cout << "nLambdaDecay:\t" << nLambdaDecay << endl;
    cout << "nFinal:\t" << nFinal << endl;

    return 0;
}
