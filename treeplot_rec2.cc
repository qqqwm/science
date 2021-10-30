typedef ROOT::Math::XYZVector Vector;
Vector TransformToNewBasis(const Vector& nx, const Vector& ny, const Vector& nz, const Vector& p)
{
  return Vector(nx.Dot(p), ny.Dot(p), nz.Dot(p));
}

Vector Normalized(Vector v)
{
  return v.Unit();
}
const Double_t alpha = 0.732;
// Define a function with 1 parameter
Double_t fitf(const Double_t *x, const Double_t *par) {
  return (1. + alpha * par[0] * x[0]) / 2.;
}
void BinLogX(TH1* h) {
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];

  for (int i = 0; i <= bins; i++) {
    new_bins[i] = TMath::Power(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  //delete new_bins;
}
void treeplot_rec2() {
  const int Nbins = 200;
  const Float_t font_size = 0.05, title_offset = 1.2, small_margin = 0.05, big_margin = 0.12, medium_margin = 0.09;
  const Int_t n_divisions = 805;

  const double min_mass = 1.08, max_mass = 1.25;
  const double max_pt = 2, max_y = 2.8, max_deltaz = 200;
  const double lambdaMass = 1.115683, protonMass = 0.93827208816, pionMass = 0.13957039, lambdaWidth = 0.00276;
  const double gamma = 9.204522397440854, gammabeta = -9.150040030786212, p_lambda_max_cms = 8.53795463662479;
  double maxRelErr = 0, maxAbsErr = 0, maxRec = 0;

  Double_t cos_x, cos_y, cos_z, DeltaZ, y;
  //double max_hist, min_hist;
  //vector<double> maxes;
  UInt_t nHits1, nHits2, nClusters1, nClusters2;
  Double_t v0Mass, cosPhi;
  UInt_t nClustersPass, XYTargetPass;
  UInt_t nS4;
  Double_t targX, targY;
  Double_t alpha, p_t;
  Double_t p_full_lambda, p_t_lambda, phi_lambda;
  Double_t pRec[9], pSim[9];
  Double_t recVertex[9], simVertex[6];

  TFile f("pp158simrec.1708.root");

  TTree* simTree = (TTree*) f.Get("LambdasMCSimTracks");

  simTree->SetBranchAddress("cos_x", &cos_x);
  simTree->SetBranchAddress("cos_y", &cos_y);
  simTree->SetBranchAddress("cos_z", &cos_z);
  simTree->SetBranchAddress("p_t_lambda", &p_t_lambda);
  simTree->SetBranchAddress("p_full_lambda", &p_full_lambda);

  TTree* t = (TTree*) f.Get("LambdasMCRecTracks");

  t->SetBranchAddress("nClusters1", &nClusters1);
  t->SetBranchAddress("nClusters2", &nClusters2);
  t->SetBranchAddress("cosPhi", &cosPhi);
  t->SetBranchAddress("targX", &targX);
  t->SetBranchAddress("targY", &targY);
  t->SetBranchAddress("pRec", pRec);
  t->SetBranchAddress("pSim", pSim);
  t->SetBranchAddress("recVertex", recVertex);
  t->SetBranchAddress("simVertex", simVertex);

  TH1D* hist_relerrlambda = new TH1D("RelativeErrorLambda", "Relative Lambda Momentum Reconstruction Error;RelativeError, % ;Entries", Nbins, -100., 100.);
  TH1D* hist_relerrproton = new TH1D("RelativeErrorProton", "Relative Proton Momentum Reconstruction Error;RelativeError, % ;Entries", Nbins, -100., 100.);
  TH1D* hist_relerrpion = new TH1D("RelativeErrorPion", "Relative Pion Momentum Reconstruction Error;RelativeError, % ;Entries", Nbins, -100., 100.);

  TH1D* hist_relerrlambda_cut = new TH1D("RelativeErrorLambda_cut", "Relative Lambda Momentum Reconstruction Error(for cos_z < -0.99);RelativeError, % ;Entries", Nbins, -100., 100.);
  TH1D* hist_relerrproton_cut = new TH1D("RelativeErrorProton_cut", "Relative Proton Momentum Reconstruction Error(for cos_z < -0.99);RelativeError, % ;Entries", Nbins, -100., 100.);
  TH1D* hist_relerrpion_cut = new TH1D("RelativeErrorPion_cut", "Relative Pion Momentum Reconstruction Error(for cos_z < -0.99);RelativeError, % ;Entries", Nbins, -100., 100.);

  TH1D* hist_cos_x_rec = new TH1D("COS_X_rec", ";cos #theta_{x};f(cos #theta_{x})", Nbins, -1., 1.);
  TH1D* hist_cos_y_rec = new TH1D("COS_Y_rec", ";cos #theta_{y};f(cos #theta_{y})", Nbins, -1., 1.);
  TH1D* hist_cos_z_rec = new TH1D("COS_Z_rec", ";cos #theta_{z};f(cos #theta_{z})", Nbins, -1., 1.);

  TH1D* hist_cos_x_sim = new TH1D("COS_X_sim", ";cos #theta_{x};f(cos #theta_{x})", Nbins, -1., 1.);
  TH1D* hist_cos_y_sim = new TH1D("COS_Y_sim", ";cos #theta_{y};f(cos #theta_{y})", Nbins, -1., 1.);
  TH1D* hist_cos_z_sim = new TH1D("COS_Z_sim", ";cos #theta_{z};f(cos #theta_{z})", Nbins, -1., 1.);

  TH1D* hist_cos_x_simrec = new TH1D("COS_X_simrec", "Sim values;cos #theta_{x};f(cos #theta_{x})", Nbins, -1., 1.);
  TH1D* hist_cos_y_simrec = new TH1D("COS_Y_simrec", "Sim values;cos #theta_{y};f(cos #theta_{y})", Nbins, -1., 1.);
  TH1D* hist_cos_z_simrec = new TH1D("COS_Z_simrec", "Sim values;cos #theta_{z};f(cos #theta_{z})", Nbins, -1., 1.);

  TH2D* hist_cos_x_2d = new TH2D("COS_X_2d", "Sim and rec values - comparison;rec cos #theta_{x};sim cos #theta_{x}", Nbins, -1., 1., Nbins, -1., 1.);
  TH2D* hist_cos_y_2d = new TH2D("COS_Y_2d", "Sim and rec values - comparison;rec cos #theta_{y};sim cos #theta_{y}", Nbins, -1., 1., Nbins, -1., 1.);
  TH2D* hist_cos_z_2d = new TH2D("COS_Z_2d", "Sim and rec values - comparison;rec cos #theta_{z};sim cos #theta_{z}", Nbins, -1., 1., Nbins, -1., 1.);

  TH1D* hist_xf_sim = new TH1D("xf_sim", "x_{F} distribution;", 400, -1, 1.);

  Double_t *new_bins = new Double_t[Nbins + 1];
  for (int i = 0; i <= Nbins; i++) {
    new_bins[i] = TMath::Power(10, (i - 200) / 50.);
  }
  TH1D* hist_logxf_sim = new TH1D("logxf_sim", "x_{F} distribution;", Nbins, new_bins);
  //BinLogX(hist_logxf_sim);

  TH1D* hist_logpxf_rec = new TH1D("pxf_rec", "x_{F} distribution;", 200, -4., 0.);
  TH1D* hist_lognxf_rec = new TH1D("nxf_rec", "-x_{F} distribution;", 200, 0., 4.);

  TH1D* hist_nClusters1 = new TH1D("nClusters1", "nClusters1 distribution;# of clusters", 73, -0.5, 72.5);
  TH1D* hist_nClusters2 = new TH1D("nClusters2", "nClusters2 distribution;# of clusters", 73, -0.5, 72.5);

  TH2D* hist_logpxf_recpt = new TH2D("xfpt", "x_{F}-p_{T} distribution;x_{F};p_{T}, [GeV/c]", 200, 1e-4, 1., 150, 0.0, max_pt);
  TH2D* hist_ypt = new TH2D("ypt", "y-p_{T} distribution;y;p_{T}, [GeV/c]", 200, -max_y, max_y, 150, 0.0, max_pt);
  TH1D* hist_z = new TH1D("Deltaz", "#Deltaz distribution;#Deltaz,[cm]", 200, 0.0, max_deltaz);
  TH2D* hist_pxy = new TH2D("pxy", "p_{x}-p_{y} distribution;p_{x}, [GeV/c];p_{y}, [GeV/c]", 100, -max_pt, max_pt, 100, -max_pt, max_pt);

  TH1D* hist_pvz_sim = new TH1D("PrimVertex", "Sim PrimVertex distribution;z,[cm]", 200, -600., -540.);
  TH2D* hist_pvxy_sim = new TH2D("pvxy_sim", "Sim PrimVertex distribution;x [cm];y [cm]", 100, -5., 5, 100, -5., 5.);
  TH1D* hist_pvz_rec = new TH1D("pvz_rec", "Rec PrimVertex-MainVertex distribution;z [cm]", 200, -10., 10.);
  TH2D* hist_pvxy_rec = new TH2D("pvxy_rec", "Rec PrimVertex-MainVertex distribution;x [cm];y [cm]", 101, -0.5, 0.5, 101, -0.5, 0.5);

  TH1D* hist_svz_sim = new TH1D("svz_sim", "Sim DecayVertex distribution;z [cm]", 200, -600., -600. + max_deltaz);
  TH2D* hist_svxy_sim = new TH2D("svxy_sim", "Sim DecayVertex distribution;x [cm];y [cm]", 100, -5., 5, 100, -5., 5.);
  TH1D* hist_svz_rec = new TH1D("svz_rec", "Rec DecayVertex distribution;z,[cm]", 200, -1., 1);
  TH2D* hist_svxy_rec = new TH2D("svxy_rec", "Rec DecayVertex distribution;x [cm];y [cm]", 101, -0.5, 0.5, 101, -0.5, 0.5);

  TH2D* hist_ypt_cut = new TH2D("ypt_0.99cut", "y-p_{T} distribution(cos #theta_{z}<-0.99);y;p_{T}, [GeV/c]", 200, -max_y, max_y, 150, 0.0, max_pt);
  TH1D* hist_z_cut = new TH1D("z_0.99cut", "#Deltaz distribution(cos #theta_{z}<-0.99);#Deltaz,[cm]", 200, 0.0, max_deltaz);
  TH2D* hist_pxy_cut = new TH2D("ypxy_0.99cut", "p_{x}-p_{y} distribution(cos #theta_{z}<-0.99);p_{x}, [GeV/c];p_{y}, [GeV/c]", 100, -max_pt, max_pt, 100, -max_pt, max_pt);
  TH1D* hist_sv_cut = new TH1D("StopVertex_0.99cut", "#Lambda decay point distribution(cos #theta_{z}<-0.99);z,[cm]", 200, -600., -600. + max_deltaz);

  TH1D* hist_nClusters1_cut = new TH1D("nClusters1_cut", "nClusters1 distribution(cos #theta_{z}<-0.99);# of clusters", 73, -0.5, 72.5);
  TH1D* hist_nClusters2_cut = new TH1D("nClusters2_cut", "nClusters2 distribution(cos #theta_{z}<-0.99);# of clusters", 73, -0.5, 72.5);

  TH2D* hist_arm = new TH2D("arm", ";#alpha;p_{T}, [GeV/c]", 200, -1., 1., 150, 0.0, 0.5);
  TH2D* hist_arm_cut = new TH2D("arm_cut", "(cos #theta_{z}<-0.99);#alpha;p_{T}, [GeV/c]", 200, -1., 1., 150, 0.0, 0.5);

  TH2D* hist_armcos = new TH2D("armcos", ";#alpha;cos #theta_{z}", 200, -1., 1., 200, -1., 1.);
  TH2D* hist_armcos_cut = new TH2D("armcos_cut", "(cos #theta_{z}<-0.99);#alpha;cos #theta_{z}", 200, -1., 1., 200, -1., 1.);

  TH2D* hist_targXY = new TH2D("targXY", "Impart Parameters distribution;x [cm];y [cm]", 100, -3., 3., 100, -3., 3.);


  TH1D* hist_invmass = new TH1D("inv", "Inv mass;m_{p#pi};", 200, 1.08, 1.14);
  TH1D* hist_invmass_cut = new TH1D("inv_cut", "Inv mass(cos #theta_{z}<-0.99);m_{p#pi};", 200, 1.08, 1.14);

  TH1D* hist_ppmomangle = new TH1D("angle", ";angle;", 200, 0., 0.2);

  /*
    const Long64_t simTreenentries = simTree->GetEntries();
    for (Long64_t i = 0; i < simTreenentries; ++i) {
      simTree->GetEntry(i);
      hist_cos_x_sim->Fill(cos_x);
      hist_cos_y_sim->Fill(cos_y);
      hist_cos_z_sim->Fill(cos_z);
      Double_t p_z_lambda_sim = sqrt(p_full_lambda * p_full_lambda - p_t_lambda * p_t_lambda);
      Double_t e_lambda_sim = sqrt(p_full_lambda * p_full_lambda + lambdaMass * lambdaMass);
      double x_f = (gammabeta * e_lambda_sim + gamma * p_z_lambda_sim) / p_lambda_max_cms;

      hist_xf_sim->Fill(x_f);
      hist_logxf_sim->Fill(fabs(x_f));
    }
  */
  const Long64_t nentries = t->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    t->GetEntry(i);
    Double_t p_lambda_sim = sqrt(pSim[0] * pSim[0] + pSim[1] * pSim[1] + pSim[2] * pSim[2]);
    Double_t p_lambda_rec = sqrt(pRec[0] * pRec[0] + pRec[1] * pRec[1] + pRec[2] * pRec[2]);

    Double_t e_lambda_rec = sqrt(pRec[0] * pRec[0] + pRec[1] * pRec[1] + pRec[2] * pRec[2] + lambdaMass * lambdaMass);
    Double_t e_lambda_sim = sqrt(pSim[0] * pSim[0] + pSim[1] * pSim[1] + pSim[2] * pSim[2] + lambdaMass * lambdaMass);

    y = atanh(pRec[2] / e_lambda_rec) - 2.909878177287;
    DeltaZ = recVertex[8] - recVertex[2];
    //hist_targXY->Fill(targX, targY);

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

    if ((nClusters1 < 15) || (nClusters2 < 15))
      continue;

    if (pow(0.5 * targX, 2) + pow(targY, 2) > 1.)
      continue;

    bool cosPhi_cut = false;
    if (y < -0.25)
    {
      if (cosPhi < 0.95)
        cosPhi_cut = true;
    }
    else if (y >= -0.25 && y < 0.75)
    {
      if (cosPhi < 0.9)
        cosPhi_cut = true;
    }
    else if (y >= 0.75)
    {
      if (cosPhi < 0.8)
        cosPhi_cut = true;
    }
    if (!cosPhi_cut) {
      continue;
    }
    if (recVertex[2] > -569 || recVertex[2] < -591)
      continue;

    Double_t p_proton_sim = sqrt(pSim[3] * pSim[3] + pSim[4] * pSim[4] + pSim[5] * pSim[5]);
    Double_t p_proton_rec = sqrt(pRec[3] * pRec[3] + pRec[4] * pRec[4] + pRec[5] * pRec[5]);

    Double_t p_pion_sim = sqrt(pSim[6] * pSim[6] + pSim[7] * pSim[7] + pSim[8] * pSim[8]);
    Double_t p_pion_rec = sqrt(pRec[6] * pRec[6] + pRec[7] * pRec[7] + pRec[8] * pRec[8]);

    Double_t e_pion_rec = sqrt(pRec[6] * pRec[6] + pRec[7] * pRec[7] + pRec[8] * pRec[8] + pionMass * pionMass);
    Double_t e_proton_rec = sqrt(pRec[3] * pRec[3] + pRec[4] * pRec[4] + pRec[5] * pRec[5] + protonMass * protonMass);

    Double_t pt_lambda_rec = sqrt(pRec[0] * pRec[0] + pRec[1] * pRec[1]);
    Double_t pt_lambda_sim = sqrt(pSim[0] * pSim[0] + pSim[1] * pSim[1]);

    //REC momenta, point-to-point direction of Lambda
    const Vector p_in(0., 0., 1.);
    Vector lambdadir(recVertex[6] - recVertex[0], recVertex[7] - recVertex[1], recVertex[8] - recVertex[2]);
    Vector p_lambda = lambdadir * (p_lambda_rec / lambdadir.R());
    //Vector p_lambda(pRec[0], pRec[1], pRec[2]);
    Vector protonMomentum(pRec[3], pRec[4], pRec[5]);
    Vector pionMomentum(pRec[6], pRec[7], pRec[8]);
    Vector nz = Normalized(p_lambda);
    Vector nx = Normalized(p_lambda.Cross(p_in));
    Vector ny = Normalized(p_lambda.Cross(p_lambda.Cross(p_in)));

    Vector newpproton = TransformToNewBasis(nx, ny, nz, protonMomentum);

    double proton_e = sqrt(protonMomentum.Mag2() + protonMass * protonMass);
    //next - lorentz boost in z direction
    double boostedproton_x = newpproton.X();
    double boostedproton_y = newpproton.Y();
    double boostedproton_z = (e_lambda_rec * newpproton.Z() - p_lambda_rec * proton_e ) / lambdaMass;
    double boostedproton_p = sqrt(boostedproton_x * boostedproton_x + boostedproton_y * boostedproton_y + boostedproton_z * boostedproton_z);
    //double boostedproton_e = (e_lambda_rec * proton_e - p_lambda_rec * newpproton.Z() ) / lambdaMass;
    cos_x = boostedproton_x / boostedproton_p;
    cos_y = boostedproton_y / boostedproton_p;
    cos_z = boostedproton_z / boostedproton_p;

    p_t = sqrt(protonMomentum.Cross(p_lambda).Mag2() / p_lambda.Mag2());
    //alpha = (protonMomentum - pionMomentum).Dot(p_lambda) / p_lambda.Mag2();
    v0Mass  = sqrt(pionMass * pionMass + protonMass * protonMass + 2 * e_pion_rec * e_proton_rec - 2 * protonMomentum.Dot(pionMomentum));

    if (v0Mass < lambdaMass - 3 * lambdaWidth)
      continue;
    if (v0Mass > lambdaMass + 3 * lambdaWidth)
      continue;
    if (p_t < 0.018)
      continue;
    //REC momenta

    //lambdadir = Vector(simVertex[3] - simVertex[0], simVertex[4] - simVertex[1], simVertex[5] - simVertex[2]);
    //p_lambda = lambdadir * (p_lambda_sim / lambdadir.R());
    p_lambda = Vector(pRec[0], pRec[1], pRec[2]);
    nz = Normalized(p_lambda);
    nx = Normalized(p_lambda.Cross(p_in));
    ny = Normalized(p_lambda.Cross(p_lambda.Cross(p_in)));

    newpproton = TransformToNewBasis(nx, ny, nz, protonMomentum);

    proton_e = sqrt(protonMomentum.Mag2() + protonMass * protonMass);
    //next - lorentz boost in z direction
    boostedproton_x = newpproton.X();
    boostedproton_y = newpproton.Y();
    boostedproton_z = (e_lambda_sim * newpproton.Z() - p_lambda_sim * proton_e ) / lambdaMass;
    boostedproton_p = sqrt(boostedproton_x * boostedproton_x + boostedproton_y * boostedproton_y + boostedproton_z * boostedproton_z);
    //boostedproton_e = (e_lambda_sim * proton_e - p_lambda_sim * newpproton.Z() ) / lambdaMass;
    double cos_x_mom = boostedproton_x / boostedproton_p;
    double cos_y_mom = boostedproton_y / boostedproton_p;
    double cos_z_mom = boostedproton_z / boostedproton_p;

    //lets rotate point-to-point direction of Lambda lambdadir and boost
    newpproton = TransformToNewBasis(nx, ny, nz, lambdadir) * (p_lambda_rec / lambdadir.R());
    boostedproton_x = newpproton.X();
    boostedproton_y = newpproton.Y();
    boostedproton_z = e_lambda_rec * (newpproton.Z() - p_lambda_rec) / lambdaMass;
    boostedproton_p = sqrt(boostedproton_x * boostedproton_x + boostedproton_y * boostedproton_y + boostedproton_z * boostedproton_z);


    //double x_f = (gammabeta * e_lambda_rec + gamma * pRec[2]) / p_lambda_max_cms;

    //hist_arm->Fill(alpha, p_t);
    //hist_armcos->Fill(alpha, cos_z);
    //hist_invmass->Fill(v0Mass);
    hist_ppmomangle->Fill(acos(p_lambda.Dot(lambdadir) / sqrt(p_lambda.Mag2() * lambdadir.Mag2())));
    //hist_ppmomangle->Fill(boostedproton_z / boostedproton_p);

    hist_cos_x_rec->Fill(cos_x);
    hist_cos_y_rec->Fill(cos_y);
    hist_cos_z_rec->Fill(cos_z);

    hist_cos_x_sim->Fill(cos_x_mom);
    hist_cos_y_sim->Fill(cos_y_mom);
    hist_cos_z_sim->Fill(cos_z_mom);

    //hist_cos_x_2d->Fill(cos_x, cos_x_sim);
    //hist_cos_y_2d->Fill(cos_y, cos_y_sim);
    //hist_cos_z_2d->Fill(cos_z, cos_z_sim);
    /*
    if (x_f < 0) hist_lognxf_rec->Fill(-TMath::Log10(-x_f));
    else hist_logpxf_rec->Fill(TMath::Log10(x_f));

    hist_logpxf_recpt->Fill(fabs(x_f), pt_lambda_sim);

    hist_ypt->Fill(y, pt_lambda_rec);
    hist_z->Fill(DeltaZ);
    hist_pxy->Fill(pRec[0], pRec[1]);

    hist_nClusters1->Fill(nClusters1);
    hist_nClusters2->Fill(nClusters2);

    hist_pvz_sim->Fill(simVertex[2]);
    hist_pvxy_sim->Fill(simVertex[0], simVertex[1]);
    hist_pvz_rec->Fill(recVertex[2]);
    hist_pvxy_rec->Fill(recVertex[0], recVertex[1]);

    hist_svz_sim->Fill(simVertex[5]);
    hist_svxy_sim->Fill(simVertex[3], simVertex[4]);
    hist_svz_rec->Fill(recVertex[8]);
    hist_svxy_rec->Fill(recVertex[6], recVertex[7]);

    hist_relerrlambda->Fill(abs(p_lambda_rec - p_lambda_sim) * 100. / p_lambda_sim);
    hist_relerrproton->Fill(abs(p_proton_rec - p_proton_sim) * 100. / p_proton_sim);
    hist_relerrpion->Fill(abs(p_pion_rec - p_pion_sim) * 100. / p_pion_sim);

    if (cos_z < -0.99)
    {
      hist_arm_cut->Fill(alpha, p_t);
      hist_armcos_cut->Fill(alpha, cos_z);
      hist_invmass_cut->Fill(v0Mass);

      hist_relerrlambda_cut->Fill(abs(p_lambda_rec - p_lambda_sim) * 100. / p_lambda_sim);
      hist_relerrproton_cut->Fill(abs(p_proton_rec - p_proton_sim) * 100. / p_proton_sim);
      hist_relerrpion_cut->Fill(abs(p_pion_rec - p_pion_sim) * 100. / p_pion_sim);
      hist_ypt_cut->Fill(y, pt_lambda_rec);
      hist_z_cut->Fill(DeltaZ);
      hist_sv_cut->Fill(recVertex[8]);
      hist_pxy_cut->Fill(pRec[0], pRec[1]);

      hist_nClusters1_cut->Fill(nClusters1);
      hist_nClusters2_cut->Fill(nClusters2);
    }
    */
  }

  cout << "Finished processing trees\n";
  TCanvas *c0 = new TCanvas("c0", "canvas", 0, 0, 1000, 800);

  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);
  gPad->SetLogy(1);
  hist_ppmomangle->Draw("hist");
  gPad->SetLeftMargin(big_margin);
  gPad->SetBottomMargin(big_margin);
  gPad->SetRightMargin(small_margin);
  gPad->SetTopMargin(medium_margin);
  hist_ppmomangle->GetXaxis()->SetLabelSize(font_size);
  hist_ppmomangle->GetYaxis()->SetLabelSize(font_size);
  hist_ppmomangle->GetXaxis()->SetTitleSize(font_size);
  hist_ppmomangle->GetYaxis()->SetTitleSize(font_size);
  hist_ppmomangle->GetYaxis()->SetTitleOffset(title_offset);
  hist_ppmomangle->GetXaxis()->SetNdivisions(n_divisions);
  hist_ppmomangle->GetYaxis()->SetNdivisions(n_divisions);
  c0->SaveAs("hist_ppmomangle_boosted_allcuts.pdf");
  /*
    TF1 *func = new TF1("fit", fitf, -1, 1, 1);
    // total cosines distributions
    //cout << hist_cos_x->GetEntries() << endl;
    func->SetParNames("P_{x}:");
    hist_cos_x_sim->Scale(hist_cos_x_sim->GetNbinsX() / (2 * hist_cos_x_sim->GetEntries()));
    auto fitresultx = hist_cos_x_sim->Fit("fit");
    TF1 *fitx = hist_cos_x_sim->GetFunction("fit");
    //Double_t x_p1 = fitx->GetParameter(0);
    //Double_t x_e1 = fitx->GetParError(0);
    func->SetParNames("P_{y}:");
    hist_cos_y_sim->Scale(hist_cos_y_sim->GetNbinsX() / (2 * hist_cos_y_sim->GetEntries()));
    auto fitresulty = hist_cos_y_sim->Fit("fit");
    TF1 *fity = hist_cos_y_sim->GetFunction("fit");
    //Double_t y_p1 = fity->GetParameter(0);
    //Double_t y_e1 = fity->GetParError(0);
    func->SetParNames("P_{z}:");
    hist_cos_z_sim->Scale(hist_cos_z_sim->GetNbinsX() / (2 * hist_cos_z_sim->GetEntries()));
    auto fitresultz = hist_cos_z_sim->Fit("fit");
    TF1 *fitz = hist_cos_z_sim->GetFunction("fit");
    //Double_t z_p1 = fitz->GetParameter(0);
    //Double_t z_e1 = fitz->GetParError(0);

    double max_hist = 0.8, min_hist = 0.;

    hist_cos_x_sim->SetMaximum(max_hist);
    hist_cos_x_sim->SetMinimum(min_hist);
    hist_cos_x_sim->Draw("hist");
    hist_cos_x_sim->GetListOfFunctions()->FindObject("fit")->Draw("same");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_cos_x_sim->GetXaxis()->SetLabelSize(font_size);
    hist_cos_x_sim->GetYaxis()->SetLabelSize(font_size);
    hist_cos_x_sim->GetXaxis()->SetTitleSize(font_size);
    hist_cos_x_sim->GetYaxis()->SetTitleSize(font_size);
    hist_cos_x_sim->GetYaxis()->SetTitleOffset(title_offset);
    hist_cos_x_sim->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_x_sim->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_cos_x_mom_cuts+567.pdf");

    hist_cos_y_sim->SetMaximum(max_hist);
    hist_cos_y_sim->SetMinimum(min_hist);
    hist_cos_y_sim->Draw("hist");
    hist_cos_y_sim->GetListOfFunctions()->FindObject("fit")->Draw("same");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_cos_y_sim->GetXaxis()->SetLabelSize(font_size);
    hist_cos_y_sim->GetYaxis()->SetLabelSize(font_size);
    hist_cos_y_sim->GetXaxis()->SetTitleSize(font_size);
    hist_cos_y_sim->GetYaxis()->SetTitleSize(font_size);
    hist_cos_y_sim->GetYaxis()->SetTitleOffset(title_offset);
    hist_cos_y_sim->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_y_sim->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_cos_y_mom_cuts+567.pdf");

    hist_cos_z_sim->SetMaximum(max_hist);
    hist_cos_z_sim->SetMinimum(min_hist);
    hist_cos_z_sim->Draw("hist");
    hist_cos_z_sim->GetListOfFunctions()->FindObject("fit")->Draw("same");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_cos_z_sim->GetXaxis()->SetLabelSize(font_size);
    hist_cos_z_sim->GetYaxis()->SetLabelSize(font_size);
    hist_cos_z_sim->GetXaxis()->SetTitleSize(font_size);
    hist_cos_z_sim->GetYaxis()->SetTitleSize(font_size);
    hist_cos_z_sim->GetYaxis()->SetTitleOffset(title_offset);
    hist_cos_z_sim->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_z_sim->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_cos_z_mom_cuts+567.pdf");

    //max_hist = 1, min_hist = 0;
    hist_cos_x_rec->Scale(hist_cos_x_rec->GetNbinsX() / (2 * hist_cos_x_rec->GetEntries()));
    hist_cos_y_rec->Scale(hist_cos_y_rec->GetNbinsX() / (2 * hist_cos_y_rec->GetEntries()));
    hist_cos_z_rec->Scale(hist_cos_z_rec->GetNbinsX() / (2 * hist_cos_z_rec->GetEntries()));

    hist_cos_x_rec->SetMaximum(max_hist);
    hist_cos_x_rec->SetMinimum(min_hist);
    hist_cos_x_rec->Draw("hist");
    //hist_cos_x->GetListOfFunctions()->FindObject("fit")->Draw("same");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_cos_x_rec->GetXaxis()->SetLabelSize(font_size);
    hist_cos_x_rec->GetYaxis()->SetLabelSize(font_size);
    hist_cos_x_rec->GetXaxis()->SetTitleSize(font_size);
    hist_cos_x_rec->GetYaxis()->SetTitleSize(font_size);
    hist_cos_x_rec->GetYaxis()->SetTitleOffset(title_offset);
    hist_cos_x_rec->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_x_rec->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_cos_x_pp_cuts+567.pdf");

    hist_cos_y_rec->SetMaximum(max_hist);
    hist_cos_y_rec->SetMinimum(min_hist);
    hist_cos_y_rec->Draw("hist");
    //hist_cos_y->GetListOfFunctions()->FindObject("fit")->Draw("same");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_cos_y_rec->GetXaxis()->SetLabelSize(font_size);
    hist_cos_y_rec->GetYaxis()->SetLabelSize(font_size);
    hist_cos_y_rec->GetXaxis()->SetTitleSize(font_size);
    hist_cos_y_rec->GetYaxis()->SetTitleSize(font_size);
    hist_cos_y_rec->GetYaxis()->SetTitleOffset(title_offset);
    hist_cos_y_rec->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_y_rec->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_cos_y_pp_cuts+567.pdf");

    hist_cos_z_rec->SetMaximum(max_hist);
    hist_cos_z_rec->SetMinimum(min_hist);
    hist_cos_z_rec->Draw("hist");
    //hist_cos_z->GetListOfFunctions()->FindObject("fit")->Draw("same");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_cos_z_rec->GetXaxis()->SetLabelSize(font_size);
    hist_cos_z_rec->GetYaxis()->SetLabelSize(font_size);
    hist_cos_z_rec->GetXaxis()->SetTitleSize(font_size);
    hist_cos_z_rec->GetYaxis()->SetTitleSize(font_size);
    hist_cos_z_rec->GetYaxis()->SetTitleOffset(title_offset);
    hist_cos_z_rec->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_z_rec->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_cos_z_pp_cuts+567.pdf");*/
  /*
    hist_cos_x_simrec->Scale(hist_cos_x_simrec->GetNbinsX() / (2 * hist_cos_x_simrec->GetEntries()));
    hist_cos_y_simrec->Scale(hist_cos_y_simrec->GetNbinsX() / (2 * hist_cos_y_simrec->GetEntries()));
    hist_cos_z_simrec->Scale(hist_cos_z_simrec->GetNbinsX() / (2 * hist_cos_z_simrec->GetEntries()));

    hist_cos_x_simrec->SetMaximum(max_hist);
    hist_cos_x_simrec->SetMinimum(min_hist);
    hist_cos_x_simrec->Draw("hist");
    //hist_cos_x->GetListOfFunctions()->FindObject("fit")->Draw("same");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_cos_x_simrec->GetXaxis()->SetLabelSize(font_size);
    hist_cos_x_simrec->GetYaxis()->SetLabelSize(font_size);
    hist_cos_x_simrec->GetXaxis()->SetTitleSize(font_size);
    hist_cos_x_simrec->GetYaxis()->SetTitleSize(font_size);
    hist_cos_x_simrec->GetYaxis()->SetTitleOffset(title_offset);
    hist_cos_x_simrec->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_x_simrec->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_cos_x_simrec.pdf");

    hist_cos_y_simrec->SetMaximum(max_hist);
    hist_cos_y_simrec->SetMinimum(min_hist);
    hist_cos_y_simrec->Draw("hist");
    //hist_cos_y->GetListOfFunctions()->FindObject("fit")->Draw("same");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_cos_y_simrec->GetXaxis()->SetLabelSize(font_size);
    hist_cos_y_simrec->GetYaxis()->SetLabelSize(font_size);
    hist_cos_y_simrec->GetXaxis()->SetTitleSize(font_size);
    hist_cos_y_simrec->GetYaxis()->SetTitleSize(font_size);
    hist_cos_y_simrec->GetYaxis()->SetTitleOffset(title_offset);
    hist_cos_y_simrec->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_y_simrec->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_cos_y_simrec.pdf");

    hist_cos_z_simrec->SetMaximum(max_hist);
    hist_cos_z_simrec->SetMinimum(min_hist);
    hist_cos_z_simrec->Draw("hist");
    //hist_cos_z->GetListOfFunctions()->FindObject("fit")->Draw("same");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_cos_z_simrec->GetXaxis()->SetLabelSize(font_size);
    hist_cos_z_simrec->GetYaxis()->SetLabelSize(font_size);
    hist_cos_z_simrec->GetXaxis()->SetTitleSize(font_size);
    hist_cos_z_simrec->GetYaxis()->SetTitleSize(font_size);
    hist_cos_z_simrec->GetYaxis()->SetTitleOffset(title_offset);
    hist_cos_z_simrec->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_z_simrec->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_cos_z_simrec.pdf");

    hist_arm->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_arm->GetXaxis()->SetLabelSize(font_size);
    hist_arm->GetYaxis()->SetLabelSize(font_size);
    hist_arm->GetXaxis()->SetTitleSize(font_size);
    hist_arm->GetYaxis()->SetTitleSize(font_size);
    hist_arm->GetYaxis()->SetTitleOffset(title_offset);
    hist_arm->GetXaxis()->SetNdivisions(n_divisions);
    hist_arm->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_arm.pdf");

    hist_arm_cut->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_arm_cut->GetXaxis()->SetLabelSize(font_size);
    hist_arm_cut->GetYaxis()->SetLabelSize(font_size);
    hist_arm_cut->GetXaxis()->SetTitleSize(font_size);
    hist_arm_cut->GetYaxis()->SetTitleSize(font_size);
    hist_arm_cut->GetYaxis()->SetTitleOffset(title_offset);
    hist_arm_cut->GetXaxis()->SetNdivisions(n_divisions);
    hist_arm_cut->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_arm_cut.pdf");

    hist_armcos->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_armcos->GetXaxis()->SetLabelSize(font_size);
    hist_armcos->GetYaxis()->SetLabelSize(font_size);
    hist_armcos->GetXaxis()->SetTitleSize(font_size);
    hist_armcos->GetYaxis()->SetTitleSize(font_size);
    hist_armcos->GetYaxis()->SetTitleOffset(title_offset);
    hist_armcos->GetXaxis()->SetNdivisions(n_divisions);
    hist_armcos->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_armcos.pdf");

    hist_armcos_cut->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_armcos_cut->GetXaxis()->SetLabelSize(font_size);
    hist_armcos_cut->GetYaxis()->SetLabelSize(font_size);
    hist_armcos_cut->GetXaxis()->SetTitleSize(font_size);
    hist_armcos_cut->GetYaxis()->SetTitleSize(font_size);
    hist_armcos_cut->GetYaxis()->SetTitleOffset(title_offset);
    hist_armcos_cut->GetXaxis()->SetNdivisions(n_divisions);
    hist_armcos_cut->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_armcos_cut.pdf");

    hist_invmass->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    gPad->SetLogy(1);
    hist_invmass->GetXaxis()->SetLabelSize(font_size);
    hist_invmass->GetYaxis()->SetLabelSize(font_size);
    hist_invmass->GetXaxis()->SetTitleSize(font_size);
    hist_invmass->GetYaxis()->SetTitleSize(font_size);
    hist_invmass->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_invmass.pdf");

    hist_invmass_cut->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_invmass_cut->GetXaxis()->SetLabelSize(font_size);
    hist_invmass_cut->GetYaxis()->SetLabelSize(font_size);
    hist_invmass_cut->GetXaxis()->SetTitleSize(font_size);
    hist_invmass_cut->GetYaxis()->SetTitleSize(font_size);
    hist_invmass_cut->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_invmass_cut.pdf");
    gPad->SetLogy(0);

    hist_targXY->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_targXY->GetXaxis()->SetLabelSize(font_size);
    hist_targXY->GetYaxis()->SetLabelSize(font_size);
    hist_targXY->GetXaxis()->SetTitleSize(font_size);
    hist_targXY->GetYaxis()->SetTitleSize(font_size);
    hist_targXY->GetYaxis()->SetTitleOffset(title_offset);
    hist_targXY->GetXaxis()->SetNdivisions(n_divisions);
    hist_targXY->GetYaxis()->SetNdivisions(n_divisions);
    TEllipse *el3 = new TEllipse(0.,0.,2.,1.);
    el3->SetLineColor(2);
    el3->SetFillColor(0);
    el3->SetFillStyle(0);
    el3->Draw();
    c0->SaveAs("hist_targXY.pdf");

    hist_xf_sim->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_xf_sim->GetXaxis()->SetLabelSize(font_size);
    hist_xf_sim->GetYaxis()->SetLabelSize(font_size);
    hist_xf_sim->GetXaxis()->SetTitleSize(font_size);
    hist_xf_sim->GetYaxis()->SetTitleSize(font_size);
    hist_xf_sim->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_xf_sim.pdf");

    hist_logxf_sim->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    gPad->SetLogx(1);
    hist_logxf_sim->GetXaxis()->SetLabelSize(font_size);
    hist_logxf_sim->GetYaxis()->SetLabelSize(font_size);
    hist_logxf_sim->GetXaxis()->SetTitleSize(font_size);
    hist_logxf_sim->GetYaxis()->SetTitleSize(font_size);
    hist_logxf_sim->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_logxf_sim.pdf");
    gPad->SetLogx(0);*/
  /*
    hist_cos_x_2d->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_cos_x_2d->GetXaxis()->SetLabelSize(font_size);
    hist_cos_x_2d->GetYaxis()->SetLabelSize(font_size);
    hist_cos_x_2d->GetXaxis()->SetTitleSize(font_size);
    hist_cos_x_2d->GetYaxis()->SetTitleSize(font_size);
    hist_cos_x_2d->GetYaxis()->SetTitleOffset(title_offset);
    hist_cos_x_2d->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_x_2d->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_cos_x_recsim_vertexposition.pdf");

    hist_cos_y_2d->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_cos_y_2d->GetXaxis()->SetLabelSize(font_size);
    hist_cos_y_2d->GetYaxis()->SetLabelSize(font_size);
    hist_cos_y_2d->GetXaxis()->SetTitleSize(font_size);
    hist_cos_y_2d->GetYaxis()->SetTitleSize(font_size);
    hist_cos_y_2d->GetYaxis()->SetTitleOffset(title_offset);
    hist_cos_y_2d->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_y_2d->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_cos_y_recsim_vertexposition.pdf");

    hist_cos_z_2d->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_cos_z_2d->GetXaxis()->SetLabelSize(font_size);
    hist_cos_z_2d->GetYaxis()->SetLabelSize(font_size);
    hist_cos_z_2d->GetXaxis()->SetTitleSize(font_size);
    hist_cos_z_2d->GetYaxis()->SetTitleSize(font_size);
    hist_cos_z_2d->GetYaxis()->SetTitleOffset(title_offset);
    hist_cos_z_2d->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_z_2d->GetYaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_cos_z_recsim_vertexposition.pdf");
    */
  /*
  TGaxis *axis = new TGaxis(-4,0,0,0, 1e-4,1,504,"G");

  hist_logpxf_rec->Draw("hist");
  gPad->SetLeftMargin(big_margin);
  gPad->SetBottomMargin(big_margin);
  gPad->SetRightMargin(big_margin);
  gPad->SetTopMargin(medium_margin);
  hist_logpxf_rec->GetXaxis()->SetLabelSize(font_size);
  hist_logpxf_rec->GetYaxis()->SetLabelSize(font_size);
  hist_logpxf_rec->GetXaxis()->SetTitleSize(font_size);
  hist_logpxf_rec->GetYaxis()->SetTitleSize(font_size);
  hist_logpxf_rec->GetYaxis()->SetTitleOffset(title_offset);
  hist_logpxf_rec->GetXaxis()->SetTickLength(0.);
  axis->Draw();
  c0->SaveAs("hist_logpxf_rec_recsim.pdf");

  hist_lognxf_rec->Draw("hist");
  gPad->SetLeftMargin(big_margin);
  gPad->SetBottomMargin(big_margin);
  gPad->SetRightMargin(big_margin);
  gPad->SetTopMargin(medium_margin);
  hist_lognxf_rec->GetXaxis()->SetLabelSize(font_size);
  hist_lognxf_rec->GetYaxis()->SetLabelSize(font_size);
  hist_lognxf_rec->GetXaxis()->SetTitleSize(font_size);
  hist_lognxf_rec->GetYaxis()->SetTitleSize(font_size);
  hist_lognxf_rec->GetYaxis()->SetTitleOffset(title_offset);
  hist_lognxf_rec->GetXaxis()->SetTickLength(0.);
  axis->Draw();
  c0->SaveAs("hist_lognxf_rec_recsim.pdf");

  hist_logpxf_recpt->Draw("COLZ");
  gPad->SetLeftMargin(big_margin);
  gPad->SetBottomMargin(big_margin);
  gPad->SetRightMargin(big_margin);
  gPad->SetTopMargin(medium_margin);
  hist_logpxf_recpt->GetXaxis()->SetLabelSize(font_size);
  hist_logpxf_recpt->GetYaxis()->SetLabelSize(font_size);
  hist_logpxf_recpt->GetXaxis()->SetTitleSize(font_size);
  hist_logpxf_recpt->GetYaxis()->SetTitleSize(font_size);
  hist_logpxf_recpt->GetYaxis()->SetTitleOffset(title_offset);
  c0->SaveAs("hist_logpxf_recpt_recsim.pdf");
  //gPad->SetLogx(0);
  */
  /*
  hist_ypt->Draw("COLZ");
  gPad->SetLeftMargin(big_margin);
  gPad->SetBottomMargin(big_margin);
  gPad->SetRightMargin(big_margin);
  gPad->SetTopMargin(medium_margin);
  hist_ypt->GetXaxis()->SetLabelSize(font_size);
  hist_ypt->GetYaxis()->SetLabelSize(font_size);
  hist_ypt->GetXaxis()->SetTitleSize(font_size);
  hist_ypt->GetYaxis()->SetTitleSize(font_size);
  hist_ypt->GetYaxis()->SetTitleOffset(title_offset);
  c0->SaveAs("hist_ypt_rec.pdf");

  hist_z->Draw("hist");
  gPad->SetLeftMargin(big_margin);
  gPad->SetBottomMargin(big_margin);
  gPad->SetRightMargin(small_margin);
  gPad->SetTopMargin(medium_margin);
  hist_z->GetXaxis()->SetLabelSize(font_size);
  hist_z->GetYaxis()->SetLabelSize(font_size);
  hist_z->GetXaxis()->SetTitleSize(font_size);
  hist_z->GetYaxis()->SetTitleSize(font_size);
  hist_z->GetYaxis()->SetTitleOffset(title_offset);
  c0->SaveAs("hist_delta_z_rec.pdf");

  hist_pxy->Draw("COLZ");
  gPad->SetLeftMargin(big_margin);
  gPad->SetBottomMargin(big_margin);
  gPad->SetRightMargin(big_margin);
  gPad->SetTopMargin(medium_margin);
  hist_pxy->GetXaxis()->SetLabelSize(font_size);
  hist_pxy->GetYaxis()->SetLabelSize(font_size);
  hist_pxy->GetXaxis()->SetTitleSize(font_size);
  hist_pxy->GetYaxis()->SetTitleSize(font_size);
  hist_pxy->GetYaxis()->SetTitleOffset(title_offset);
  c0->SaveAs("hist_pxy_rec.pdf");
  *//*
hist_pvz_sim->Draw("hist");
gPad->SetLeftMargin(big_margin);
gPad->SetBottomMargin(big_margin);
gPad->SetRightMargin(small_margin);
gPad->SetTopMargin(medium_margin);
hist_pvz_sim->GetXaxis()->SetLabelSize(font_size);
hist_pvz_sim->GetYaxis()->SetLabelSize(font_size);
hist_pvz_sim->GetXaxis()->SetTitleSize(font_size);
hist_pvz_sim->GetYaxis()->SetTitleSize(font_size);
hist_pvz_sim->GetYaxis()->SetTitleOffset(title_offset);
hist_pvz_sim->GetXaxis()->SetNdivisions(n_divisions);
c0->SaveAs("hist_pvz_sim.pdf");

hist_pvxy_sim->Draw("COLZ");
gPad->SetLeftMargin(big_margin);
gPad->SetBottomMargin(big_margin);
gPad->SetRightMargin(big_margin);
gPad->SetTopMargin(medium_margin);
hist_pvxy_sim->GetXaxis()->SetLabelSize(font_size);
hist_pvxy_sim->GetYaxis()->SetLabelSize(font_size);
hist_pvxy_sim->GetXaxis()->SetTitleSize(font_size);
hist_pvxy_sim->GetYaxis()->SetTitleSize(font_size);
hist_pvxy_sim->GetYaxis()->SetTitleOffset(title_offset);
c0->SaveAs("hist_pvxy_sim.pdf");

hist_pvz_rec->Draw("hist");
gPad->SetLeftMargin(big_margin);
gPad->SetBottomMargin(big_margin);
gPad->SetRightMargin(small_margin);
gPad->SetTopMargin(medium_margin);
hist_pvz_rec->GetXaxis()->SetLabelSize(font_size);
hist_pvz_rec->GetYaxis()->SetLabelSize(font_size);
hist_pvz_rec->GetXaxis()->SetTitleSize(font_size);
hist_pvz_rec->GetYaxis()->SetTitleSize(font_size);
hist_pvz_rec->GetYaxis()->SetTitleOffset(title_offset);
hist_pvz_rec->GetXaxis()->SetNdivisions(n_divisions);
c0->SaveAs("hist_p-mvz_rec.pdf");

hist_pvxy_rec->Draw("COLZ");
gPad->SetLeftMargin(big_margin);
gPad->SetBottomMargin(big_margin);
gPad->SetRightMargin(big_margin);
gPad->SetTopMargin(medium_margin);
hist_pvxy_rec->GetXaxis()->SetLabelSize(font_size);
hist_pvxy_rec->GetYaxis()->SetLabelSize(font_size);
hist_pvxy_rec->GetXaxis()->SetTitleSize(font_size);
hist_pvxy_rec->GetYaxis()->SetTitleSize(font_size);
hist_pvxy_rec->GetYaxis()->SetTitleOffset(title_offset);
c0->SaveAs("hist_p-mvxy_rec.pdf");*/
  /*
    hist_svz_sim->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_svz_sim->GetXaxis()->SetLabelSize(font_size);
    hist_svz_sim->GetYaxis()->SetLabelSize(font_size);
    hist_svz_sim->GetXaxis()->SetTitleSize(font_size);
    hist_svz_sim->GetYaxis()->SetTitleSize(font_size);
    hist_svz_sim->GetYaxis()->SetTitleOffset(title_offset);
    hist_svz_sim->GetXaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_svz_sim.pdf");

    hist_svxy_sim->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_svxy_sim->GetXaxis()->SetLabelSize(font_size);
    hist_svxy_sim->GetYaxis()->SetLabelSize(font_size);
    hist_svxy_sim->GetXaxis()->SetTitleSize(font_size);
    hist_svxy_sim->GetYaxis()->SetTitleSize(font_size);
    hist_svxy_sim->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_svxy_sim.pdf");

    hist_svz_rec->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_svz_rec->GetXaxis()->SetLabelSize(font_size);
    hist_svz_rec->GetYaxis()->SetLabelSize(font_size);
    hist_svz_rec->GetXaxis()->SetTitleSize(font_size);
    hist_svz_rec->GetYaxis()->SetTitleSize(font_size);
    hist_svz_rec->GetYaxis()->SetTitleOffset(title_offset);
    hist_svz_rec->GetXaxis()->SetNdivisions(n_divisions);
    c0->SaveAs("hist_svz_rec.pdf");

    hist_svxy_rec->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_svxy_rec->GetXaxis()->SetLabelSize(font_size);
    hist_svxy_rec->GetYaxis()->SetLabelSize(font_size);
    hist_svxy_rec->GetXaxis()->SetTitleSize(font_size);
    hist_svxy_rec->GetYaxis()->SetTitleSize(font_size);
    hist_svxy_rec->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_svxy_rec.pdf");

    hist_ypt_cut->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_ypt_cut->GetXaxis()->SetLabelSize(font_size);
    hist_ypt_cut->GetYaxis()->SetLabelSize(font_size);
    hist_ypt_cut->GetXaxis()->SetTitleSize(font_size);
    hist_ypt_cut->GetYaxis()->SetTitleSize(font_size);
    hist_ypt_cut->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_ypt_cut_rec.pdf");

    hist_z_cut->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_z_cut->GetXaxis()->SetLabelSize(font_size);
    hist_z_cut->GetYaxis()->SetLabelSize(font_size);
    hist_z_cut->GetXaxis()->SetTitleSize(font_size);
    hist_z_cut->GetYaxis()->SetTitleSize(font_size);
    hist_z_cut->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_delta_z_cut_rec.pdf");

    hist_pxy_cut->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_pxy_cut->GetXaxis()->SetLabelSize(font_size);
    hist_pxy_cut->GetYaxis()->SetLabelSize(font_size);
    hist_pxy_cut->GetXaxis()->SetTitleSize(font_size);
    hist_pxy_cut->GetYaxis()->SetTitleSize(font_size);
    hist_pxy_cut->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_pxy_cut_rec.pdf");

    hist_sv_cut->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_sv_cut->GetXaxis()->SetLabelSize(font_size);
    hist_sv_cut->GetYaxis()->SetLabelSize(font_size);
    hist_sv_cut->GetXaxis()->SetTitleSize(font_size);
    hist_sv_cut->GetYaxis()->SetTitleSize(font_size);
    hist_sv_cut->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_sv_cut_rec.pdf");

    hist_nClusters1->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_nClusters1->GetXaxis()->SetLabelSize(font_size);
    hist_nClusters1->GetYaxis()->SetLabelSize(font_size);
    hist_nClusters1->GetXaxis()->SetTitleSize(font_size);
    hist_nClusters1->GetYaxis()->SetTitleSize(font_size);
    hist_nClusters1->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_nClusters1.pdf");

    hist_nClusters2->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_nClusters2->GetXaxis()->SetLabelSize(font_size);
    hist_nClusters2->GetYaxis()->SetLabelSize(font_size);
    hist_nClusters2->GetXaxis()->SetTitleSize(font_size);
    hist_nClusters2->GetYaxis()->SetTitleSize(font_size);
    hist_nClusters2->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_nClusters2.pdf");

    hist_nClusters1_cut->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_nClusters1_cut->GetXaxis()->SetLabelSize(font_size);
    hist_nClusters1_cut->GetYaxis()->SetLabelSize(font_size);
    hist_nClusters1_cut->GetXaxis()->SetTitleSize(font_size);
    hist_nClusters1_cut->GetYaxis()->SetTitleSize(font_size);
    hist_nClusters1_cut->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_nClusters1_cut.pdf");

    hist_nClusters2_cut->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_nClusters2_cut->GetXaxis()->SetLabelSize(font_size);
    hist_nClusters2_cut->GetYaxis()->SetLabelSize(font_size);
    hist_nClusters2_cut->GetXaxis()->SetTitleSize(font_size);
    hist_nClusters2_cut->GetYaxis()->SetTitleSize(font_size);
    hist_nClusters2_cut->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_nClusters2_cut.pdf");

    hist_relerrlambda->Draw("hist");
    gPad->SetLogy();
    c0->SaveAs("hist_relerrlambda.png");

    hist_relerrproton->Draw("hist");
    //gPad->SetLogy();
    c0->SaveAs("hist_relerrproton.png");

    hist_relerrpion->Draw("hist");
    //gPad->SetLogy();
    c0->SaveAs("hist_relerrpion.png");

    hist_relerrlambda_cut->Draw("hist");
    //gPad->SetLogy();
    c0->SaveAs("hist_relerrlambda_0.99cut.png");

    hist_relerrproton_cut->Draw("hist");
    //gPad->SetLogy();
    c0->SaveAs("hist_relerrproton_0.99cut.png");

    hist_relerrpion_cut->Draw("hist");
    //gPad->SetLogy();
    c0->SaveAs("hist_relerrpion_0.99cut.png");
    gPad->SetLogy(0);

    delete c0;

    c0 = new TCanvas("canvas_pm", "x_{F} distribution", 0, 0, 1366, 768);
    gStyle->SetPalette(103);
    TColor::InvertPalette();
    //gPad->SetRightMargin(0.03);
    //gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetFrameBorderMode(0);
    c0->Divide(2, 1, -0.01, 0.0);
    c0->cd(1);
    gStyle->SetOptStat(0);
    gPad->SetRightMargin(0);
    gPad->SetLeftMargin(0.19);
    gPad->SetFrameLineColor(0);
    hist_lognxf_rec->Draw("hist");
    hist_lognxf_rec->GetYaxis()->SetTitle("Entries");
    //hist_lognxf_rec->GetYaxis()->SetLabelFont(43);
    //hist_lognxf_rec->GetYaxis()->SetTitleOffset(1.3);
    hist_lognxf_rec->GetXaxis()->SetLabelOffset(30.);
    hist_lognxf_rec->GetXaxis()->SetTickLength(0.);
    hist_lognxf_rec->SetMaximum(13000);
    hist_lognxf_rec->SetMinimum(0);
    hist_lognxf_rec->GetYaxis()->SetTitleSize(0.043);
    hist_lognxf_rec->GetYaxis()->SetLabelSize(0.043);
    //hist_lognxf_rec->GetYaxis()->SetBinLabel(1,"1");
    //hist_lognxf_rec->GetYaxis()->SetNdivisions(510);
    double xmin = 1e-4;
    double xmax = 1;
    double ymin = -6;
    double ymax = 12000;
    TGaxis *axis = new TGaxis(-TMath::Log10(xmin), ymin, -TMath::Log10(xmax), ymin, xmin, xmax, 504, "G-L");
    axis->SetLabelFont(43);
    axis->SetLabelSize(30);
    axis->SetLabelOffset(0.055);
    axis->SetTitle("-x_{F}");
    axis->SetTitleFont(43);
    axis->SetTitleSize(33);
    axis->SetTitleOffset(1.15);
    axis->Draw();
    hist_lognxf_rec->SetTitle("x_{F}<0");
    gStyle->SetNumberContours(99);

    c0->cd(2);
    gStyle->SetOptStat(0);
    hist_logpxf_rec->Draw("hist");
    BinLogX(hist_logpxf_rec);
    hist_logpxf_rec->GetXaxis()->SetTitle("x_{F}");
    gStyle->SetPalette(103);
    //TColor::InvertPalette();
    gPad->SetLogx();
    gStyle->SetNumberContours(99);
    hist_logpxf_rec->SetMaximum(13000);
    hist_logpxf_rec->SetMinimum(0);
    hist_logpxf_rec->GetXaxis()->SetLabelFont(43);
    hist_logpxf_rec->GetXaxis()->SetLabelSize(30);
    hist_logpxf_rec->GetXaxis()->SetLabelOffset(0.00);
    hist_logpxf_rec->GetXaxis()->SetTitle("x_{F}");
    hist_logpxf_rec->GetXaxis()->SetTitleFont(43);
    hist_logpxf_rec->GetXaxis()->SetTitleSize(33);
    hist_logpxf_rec->GetXaxis()->SetTitleOffset(1.15);
    hist_logpxf_rec->GetYaxis()->SetLabelOffset(30.);
    hist_logpxf_rec->GetYaxis()->SetTickLength(0.);
    hist_logpxf_rec->GetYaxis()->SetAxisColor(kWhite);
    //hist_logpxf_rec->GetZaxis()->SetLabelFont(43);
    hist_logpxf_rec->SetMinimum(0);
    hist_logpxf_rec->SetTitle("x_{F}>0");
    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0);
    gPad->SetFrameLineColor(0);
    c0->cd();
    c0->SetTitle("x_{F} distribution");
    c0->SaveAs("hist_logpxf_rec_recsim.png");
    TColor::InvertPalette();
    */
}
