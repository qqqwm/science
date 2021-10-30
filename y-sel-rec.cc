#include <evt/Event.h>
#include <evt/SimEvent .h>
#include <evt/RecEvent.h>
#include <evt/sim/VertexTrack.h>
#include <evt/rec/RecEventConst.h>
#include <evt/rec/Trigger.h>
#include <io/EventFileChain.h>
#include <io/IoCodes.h>
#include <utl/PDGParticleIds.h>
#include <fwk/CentralConfig.h>
#include <det/Detector.h>
#include <det/DetectorConst.h>
#include <det/TPCConst.h>
#include <utl/ShineUnits.h>
#include <utl/Vector.h>
#include <utl/MathConst.h>
#include <utl/PhysicalConst.h>
#include <utl/PDGParticleIds.h>
#include <utl/DatabasePDG.h>

#include <utl/GeometryUtilities.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1I.h>
#include <TMath.h>
#include <TRandom.h>
#include <TH2I.h>

#include <iostream>
#include <cmath>
#include <vector>

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
    //const double lambdaMass = particleDB.GetParticle(ParticleConst::eLambda).Mass()*GeV;
    const double pionMass = particleDB.GetParticle(ParticleConst::ePiMinus).Mass() * GeV;

    // helper counters
    int nEvents = 0;
    int nHasMainVertex = 0;
    int nPerfectFit = 0;
    int nGoodZPosition = 0;
    int nS4cut = 0;
    int nTrackStatus = 0;
    int nHasTrack = 0;
    int nHasStopVertex = 0;
    int nHasDecayType = 0;
    int nHasTwoDaughterTracks = 0;
    int nHasCorrectCharges = 0;
    int nIdentifiedLambdaVertex = 0;
    int nXYTarget = 0;
    int nEnoughClusters = 0;
    int nDeltaz = 0;
    int nFinal = 0;

    Double_t targX = 0, targY = 0;

    Double_t cos_x, cos_y, cos_z, DeltaZ, y, cosPhi, v0Mass;
    Double_t alpha, p_t;
    Double_t p_full_lambda, p_t_lambda, phi_lambda;
    ////UInt_t nHits1, nHits2;
    UInt_t nClusters11, nClusters12, nClusters21, nClusters22;
    UInt_t nClustersPass, XYTargetPass;
    UInt_t nS4;

    TTree * t = new TTree("Lambdas", "Lambdas pp_158 MC REC");
    t->Branch("cos_x", &cos_x, "cos_x/D");
    t->Branch("cos_y", &cos_y, "cos_y/D");
    t->Branch("cos_z", &cos_z, "cos_z/D");
    //t->Branch("nHits1", &nHits1, "nHits1/i");
    //t->Branch("nHits2", &nHits2, "nHits2/i");
    //t->Branch("nClusters1", &nClusters1, "nClusters1/i");
    //t->Branch("nClusters2", &nClusters2, "nClusters2/i");
    t->Branch("DeltaZ", &DeltaZ, "DeltaZ/D");
    t->Branch("y", &y, "y/D");
    t->Branch("nS4", &nS4, "nS4/i");
    t->Branch("cosPhi", &cosPhi, "cosPhi/D");
    t->Branch("nClustersPass", &nClustersPass, "nClustersPass/i");
    t->Branch("XYTargetPass", &XYTargetPass, "XYTargetPass/i");
    //t->Branch("v0Mass", &v0Mass, "v0Mass/D");
    t->Branch("p_pi_mass", &v0Mass, "p_pi_mass/D");
    //t->Branch("flightpath", &flightpath, "flightpath/D");
    t->Branch("alpha", &alpha, "alpha/D");
    t->Branch("p_t", &p_t, "p_t/D");
    t->Branch("p_t_lambda", &p_t_lambda, "p_t_lambda/D");
    t->Branch("p_full_lambda", &p_full_lambda, "p_full_lambda/D");
    t->Branch("phi_lambda", &phi_lambda, "phi_lambda/D");
    

    // ---  set chain of input files
    //const vector<string> fileNames(argv + 2, argv + argc);
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
        if (!(nEvents % 10000))
            cout << "\n ---> processing event # " << nEvents << endl;
        using namespace evt::rec;
        const RecEvent& recEvent = event.GetRecEvent();
        const SimEvent& simEvent = event.GetSimEvent();
        if (!recEvent.HasMainVertex())
            continue;
        ++nHasMainVertex;
        //const Vertex& mainVertex = recEvent.GetMainVertex();

        // 7. // Fit quality
        if (!recEvent.HasPrimaryVertex(evt::rec::VertexConst::ePrimaryFitZ))
            continue;

        const Vertex& fitVertex = recEvent.GetPrimaryVertex(VertexConst::ePrimaryFitZ);
        if (fitVertex.GetFitQuality() != FitQuality::ePerfect)
            continue;
        ++nPerfectFit;

        // 8. // Fitted Vtx Z position
        if ((double)(fitVertex.GetPosition().GetZ() / cm) > -540.3 \
                || (double)(fitVertex.GetPosition().GetZ() / cm) < -620.3)
            continue;
        ++nGoodZPosition;

        auto p_in = recEvent.GetBeam().GetMomentum();
        nS4 = 0;
        for (evt::SimEvent::HitIterator iter = simEvent.HitsBegin(det::Const::eS4),
                end = simEvent.HitsEnd(det::Const::eS4); iter != end; ++iter, ++nS4);
        if (nS4 == 0)
            ++nS4cut;

        for (list<VertexTrack>::const_iterator
                recTrackIter = recEvent.Begin<VertexTrack>(),
                recTrackEnd = recEvent.End<VertexTrack>();
                recTrackIter != recTrackEnd; ++recTrackIter) {
            const VertexTrack& recTrack = *recTrackIter;
            // good track status?
            if (recTrack.GetStatus() != 0)
                continue;
            ++nTrackStatus;
            //if (!recTrack.HasTrack())
            //  continue;
            ++nHasTrack;
            if (!recTrack.HasStopVertex())
                continue;
            ++nHasStopVertex;
            const Vertex& stopVertex = recEvent.Get(recTrack.GetStopVertexIndex());
            //is the stop vertex a V0 vertex?
            if (stopVertex.GetType() != VertexConst::eV0Decay)
                continue;
            ++nHasDecayType;
            const auto nofdt = stopVertex.GetNumberOfDaughterTracks();
            if (nofdt != 2)
                continue;
            ++nHasTwoDaughterTracks;

            VertexTrackIndexIterator decayTrackIter = stopVertex.DaughterTracksBegin();
            const VertexTrack& vtxTrack1 = recEvent.Get(*(decayTrackIter));
            const VertexTrack& vtxTrack2 = recEvent.Get(*(++decayTrackIter));
            // different charges ?
            const int charge1 = vtxTrack1.GetCharge();
            const int charge2 = vtxTrack2.GetCharge();
            //cout << "ch1 " << charge1 << " ch2 " << charge2 << " pid1 " << pid1 << " pid2 " << pid2 << endl;
            if (charge1 * charge2 != -1)
                continue;
            ++nHasCorrectCharges;

            const Track& track1 = recEvent.Get(vtxTrack1.GetTrackIndex());
            const Track& track2 = recEvent.Get(vtxTrack2.GetTrackIndex());

            int pid1 = 0, pid2 = 0;
            int highest_number_of_shared_points = 0;
            sim::VertexTrack simVtxTrackMatched_1, simVtxTrackMatched_2;
            for (sim::VertexTrackIndexIterator simVtxTrackIter  = track1.SimVertexTracksBegin();
                    simVtxTrackIter != track1.SimVertexTracksEnd(); ++simVtxTrackIter)
            {
                const sim::VertexTrack simVtxTrack = simEvent.Get(*simVtxTrackIter);
                if (simVtxTrack.GetRecTrackWithMaxCommonPoints() == track1.GetIndex()) {
                    auto number_of_shared_points = simVtxTrack.GetNumberOfCommonPoints(track1.GetIndex());
                    if ( number_of_shared_points > highest_number_of_shared_points ) {
                        highest_number_of_shared_points = number_of_shared_points;
                        simVtxTrackMatched_1 = simVtxTrack;
                        pid1 = simVtxTrack.GetParticleId();
                    }
                }
            }
            if ((pid1 != ParticleConst::eProton) && (pid1 != ParticleConst::ePiMinus)) {
                continue;
            }

            highest_number_of_shared_points = 0;
            for (sim::VertexTrackIndexIterator simVtxTrackIter  = track2.SimVertexTracksBegin();
                    simVtxTrackIter != track2.SimVertexTracksEnd(); ++simVtxTrackIter)
            {
                const sim::VertexTrack simVtxTrack = simEvent.Get(*simVtxTrackIter);
                if (simVtxTrack.GetRecTrackWithMaxCommonPoints() == track2.GetIndex()) {
                    auto number_of_shared_points = simVtxTrack.GetNumberOfCommonPoints(track2.GetIndex());
                    if ( number_of_shared_points > highest_number_of_shared_points ) {
                        highest_number_of_shared_points = number_of_shared_points;
                        simVtxTrackMatched_2 = simVtxTrack;
                        pid2 = simVtxTrack.GetParticleId();
                    }
                }
            }
            if ((pid2 != ParticleConst::eProton) && (pid2 != ParticleConst::ePiMinus)) {
                continue;
            }
            if (!(simVtxTrackMatched_2.HasStartVertex() && simVtxTrackMatched_1.HasStartVertex()))
                continue;
            if (simVtxTrackMatched_2.GetStartVertexIndex() != simVtxTrackMatched_1.GetStartVertexIndex())
                continue;
            sim::Vertex matchVtx1 = simEvent.Get<sim::Vertex>(simVtxTrackMatched_1.GetStartVertexIndex());
            if (matchVtx1.GetNumberOfParentTracks() != 1)
                continue;
            sim::VertexTrack lambdaTrack = simEvent.Get<sim::VertexTrack>(matchVtx1.GetFirstParentTrackIndex());
            if (lambdaTrack.GetParticleId() != 3122)
                //not a sim lambda!
                continue;
            ++nIdentifiedLambdaVertex;

            const Plane* vertexPlane = new Plane(fitVertex.GetPosition(), Vector(0, 0, 1));
            const Line* v0Trajectory = new Line(stopVertex.GetPosition(), vtxTrack1.GetMomentum() + vtxTrack2.GetMomentum());
            const Point v0AtVertexPlane = GeometryUtilities::Intersection(*vertexPlane, *v0Trajectory);
            const Vector mainToVertexPoint = v0AtVertexPlane - fitVertex.GetPosition();
            targX =  mainToVertexPoint.GetX();
            targY =  mainToVertexPoint.GetY();

            if (pow(0.5 * targX, 2) + pow(targY, 2) >= 1.)
                XYTargetPass = 0;
            else
                XYTargetPass = 1;
            ++nXYTarget;

            //impact distance wr to tracks
            //double Imp1_X, Imp1_Y;
            //double Imp2_X, Imp2_Y;
            //const Point &impact1 = vtxTrack1.GetImpactPoint();
            //const Point &impact2 = vtxTrack2.GetImpactPoint();
            //const Vector Imp1 impact1 - mainVertex.GetPosition();
            //const Vector Imp2 impact2 - mainVertex.GetPosition();
            //Imp1_X = Imp1.GetX();
            //Imp1_Y = Imp1.GetY();
            //Imp2_X = Imp2.GetX();
            //Imp2_Y = Imp2.GetY();

            //nHits1 = simVtxTrackMatched_1.GetNumberOfHits(TrackConst::eVTPC1) + simVtxTrackMatched_1.GetNumberOfHits(TrackConst::eVTPC2);
            //nHits2 = simVtxTrackMatched_2.GetNumberOfHits(TrackConst::eVTPC1) + simVtxTrackMatched_2.GetNumberOfHits(TrackConst::eVTPC2);

            nClusters11 = track1.GetNumberOfClusters(TrackConst::eVTPC1);
            nClusters12 = track1.GetNumberOfClusters(TrackConst::eVTPC2);
            nClusters21 = track2.GetNumberOfClusters(TrackConst::eVTPC1);
            nClusters22 = track2.GetNumberOfClusters(TrackConst::eVTPC2);
            /*
            if ((nClusters11 < 15) && (nClusters12 < 15)) {
                continue;
            }
            if ((nClusters21 < 15) && (nClusters22 < 15)) {
                continue;
            }
            */
            nClustersPass = 1;
            if ((nClusters11 < 15) && (nClusters12 < 15)) {
                nClustersPass = 0;
            }
            if ((nClusters21 < 15) && (nClusters22 < 15)) {
                nClustersPass = 0;
            }
            ++nEnoughClusters;

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
                //cout << "Not V0-> p + pi^- decay\n";
                //cout << "pid1:\t" << pid1 << "pid2:\t" << pid2 << endl;
                continue;
            }

            // Lambda --> p + pi^- hypothesis
            const double massSquared1 = pow(charge1 > 0 ? protonMass : pionMass, 2);
            const double massSquared2 = pow(charge2 > 0 ? protonMass : pionMass, 2);
            const Vector& v0p1 = vtxTrack1.GetMomentum();
            const Vector& v0p2 = vtxTrack2.GetMomentum();
            const double pSquared1 = v0p1.GetMag2();
            const double pSquared2 = v0p2.GetMag2();
            const double energy1 = sqrt(pSquared1 + massSquared1);
            const double energy2 = sqrt(pSquared2 + massSquared2);
            const double finalStateEnergy = energy1 + energy2;

            const Vector& v0Momentum = recTrack.GetMomentum();
            v0Mass = sqrt(pow(finalStateEnergy, 2) - v0Momentum.GetMag2());
            //p_pi_mass = sqrt(pionMass * pionMass + protonMass * protonMass + 2 * energy1 * energy2 - 2 * v0p1.Dot(v0p2));

            const double p1l = v0p1.Dot(v0Momentum) / v0Momentum.GetMag();
            const double p2l = v0p2.Dot(v0Momentum) / v0Momentum.GetMag();
            p_t = v0p1.Cross(v0Momentum).GetMag() / v0Momentum.GetMag() / GeV;
            alpha = charge1 * (p1l - p2l) / (p1l + p2l);

            const Vector& p_lambda = recTrack.GetMomentum();
            p_full_lambda = p_lambda.GetMag();
            p_t_lambda = p_lambda.GetRho();
            phi_lambda = p_lambda.GetPhi();
            //const double e_lambda = sqrt(p_lambda.GetMag2() + lambdaMass*lambdaMass);
            const double e_lambda = finalStateEnergy;
            const Vector nz = Normalized(p_lambda);
            const Vector nx = Normalized(Cross(p_lambda, p_in));
            const Vector ny = Normalized(Cross(p_lambda, Cross(p_lambda, p_in)));

            const Vector newpproton = TransformToNewBasis(nx, ny, nz, protonMomentum);

            const double proton_e = sqrt(protonMomentum.GetMag2() + protonMass * protonMass);
            //next - lorentz boost in z direction
            const double boostedproton_x = newpproton.GetX();
            const double boostedproton_y = newpproton.GetY();
            const double boostedproton_z = (e_lambda * newpproton.GetZ() - p_lambda.GetMag() * proton_e ) / v0Mass;
            const double boostedproton_p = sqrt(boostedproton_x * boostedproton_x + boostedproton_y * boostedproton_y + boostedproton_z * boostedproton_z);
            //const double boostedproton_e = (e_lambda * proton_e - p_lambda.GetMag() * newpproton.GetZ() ) / lambdaMass;
            cos_x = boostedproton_x / boostedproton_p;
            cos_y = boostedproton_y / boostedproton_p;
            cos_z = boostedproton_z / boostedproton_p;

            double vertexZ = recEvent.GetMainVertex().GetPosition().GetZ() / cm;
            double zv0  = stopVertex.GetPosition().GetZ() / cm;
            DeltaZ = zv0 - vertexZ;
            //double y = 0.5 * TMath::Log((finalStateEnergy + v0Momentum.GetZ()) / (finalStateEnergy - v0Momentum.GetZ())) - 2.909878177287;
            y = atanh(p_lambda.GetZ() / e_lambda) - 2.909878177287;
            
            bool deltazy_cut = false;
            if (y < 0.25)
            {
                if (DeltaZ > 10)
                    deltazy_cut = true;
            }
            else if (y >= 0.25 && y < 0.75)
            {
                if (DeltaZ > 15)
                    deltazy_cut = true;
            }
            else if (y >= 0.75 && y < 1.25)
            {
                if (DeltaZ > 40)
                    deltazy_cut = true;
            }
            else if (y >= 1.25)
            {
                if (DeltaZ > 60)
                    deltazy_cut = true;
            }
            if (!deltazy_cut) {
                continue;
            }
            ++nDeltaz;
            
            Vector y(0, 1, 0); //y axis
            Vector norm = Normalized(protonMomentum.Cross(pionMomentum)); //normal to decay plane
            Vector y_prim = Normalized(p_lambda.Cross(y.Cross(p_lambda))); //perpendicular to V0-momentum, on the plane spanned by y and p
            cosPhi = fabs(norm.Dot(y_prim));

            ++nFinal;
            //flightpath = recTrack.GetPathLength();
            t->Fill();
        }
    }

    // --- write histogram of correction factors to ROOT file
    TFile outFile(Form("y-selection-rec.%s.root", argv[1]), "RECREATE", "", 0);
    outFile.cd();
    t->Write();

    cout << "Event quality cuts: " << endl;
    cout << "nEvents:\t" << nEvents << endl;
    cout << "nHasMainVertex:\t" << nHasMainVertex << endl;
    cout << "nPerfectFit:\t" << nPerfectFit << endl;
    cout << "nGoodZPosition:\t" << nGoodZPosition << endl;
    cout << "nS4cut:\t" << nS4cut << endl;
    cout << "nTrackStatus:\t" << nTrackStatus << endl;
    cout << "nHasTrack:\t" << nHasTrack << endl;
    cout << "nHasStopVertex:\t" << nHasStopVertex << endl;
    cout << "nHasDecayType:\t" << nHasDecayType << endl;
    cout << "nHasTwoDaughterTracks:\t" << nHasTwoDaughterTracks << endl;
    cout << "nHasCorrectCharges:\t" << nHasCorrectCharges << endl;
    cout << "nIdentifiedLambdaVertex:\t" << nIdentifiedLambdaVertex << endl;
    cout << "nXYTarget:\t" << nXYTarget << endl;
    cout << "nEnoughClusters:\t" << nEnoughClusters << endl;
    cout << "nDeltaz:\t" << nDeltaz << endl;
    cout << "nFinal:\t" << nFinal << endl;

    return 0;
}
