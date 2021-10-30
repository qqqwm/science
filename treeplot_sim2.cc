const Double_t alpha = 0.732;
// Define a function with 1 parameter
Double_t fitf(const Double_t *x, const Double_t *par) {
    return (1. + alpha * par[0] * x[0]) / 2.;
}
// Lorentzian Peak function
Double_t lorentzianPeak(Double_t *x, Double_t *par) {
    return (0.5 * par[0] * par[2] / TMath::Pi()) / TMath::Max(1.e-10,
            (x[0] - par[1]) * (x[0] - par[1]) + .25 * par[2] * par[2]);
}
Double_t pol0(const Double_t *x, const Double_t *par) {
    return par[0];
}
Double_t max(Double_t a, Double_t b) {
    return (a > b) ? a : b;
}
void mccorr() {
    TCanvas *c0 = new TCanvas("c0", "canvas", 0, 0, 1000, 800);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(11);

    TFile f("y-selection-rec.first.0308.root");
    double max_hist = 0.8, min_hist = 0.3;
    const Float_t font_size = 0.05, title_offset = 1.2, small_margin = 0.05, big_margin = 0.12, medium_margin = 0.09;
    const Int_t n_divisions = 805;
    const int Nbins = 100;
    const double min_mass = 1.08, max_mass = 1.25;
    const double max_pt = 1.5, max_y = 2.8, max_deltaz = 100;


    TH1D* hist_cos_x = new TH1D("COS_X", ";cos #theta_{x};f(cos #theta_{x})", Nbins, -1., 1.);
    TH1D* hist_cos_y = new TH1D("COS_Y", ";cos #theta_{y};f(cos #theta_{y})", Nbins, -1., 1.);
    TH1D* hist_cos_z = new TH1D("COS_Z", ";cos #theta_{z};f(cos #theta_{z})", Nbins, -1., 1.);

    TH2D* hist_ypt_cut3 = new TH2D("ypt_cut3", "y-p_{T} distribution;y;p_{T}, [GeV/c]", 200, -max_y, max_y, 150, 0.0, max_pt);
    TH1D* hist_y_cut3 = new TH1D("y_cut3", "y distribution;y;", 200, -max_y, max_y);
    TH1D* hist_pt_cut3 = new TH1D("pt_cut3", "p_{T} distribution;p_{T}, [GeV/c];", 150, 0.0, max_pt);
    TH1D* hist_phi_cut3 = new TH1D("phi_cut3", "#phi distribution;#phi", 150, -3.14, 3.14);
    TH1D* hist_z_cut3 = new TH1D("z_cut3", "#Deltaz distribution;#Deltaz,[cm]", 200, 0.0, max_deltaz);
    TH2D* hist_y_z_cut3 = new TH2D("y_z_cut3", "y-#Deltaz distribution;y;#Deltaz,[cm]", 200, -max_y, max_y, 200, 0.0, max_deltaz);

    TH1D* hist_p_pi_mass = new TH1D("hist_p_pi_mass", ";m_{p#pi^{-}} (GeV/c^{2});Entries", 3 * Nbins, min_mass, max_mass);

    Double_t cos_x, cos_y, cos_z, DeltaZ, y;
    UInt_t nClusters1, nClusters2;
    Double_t p_pi_mass, cosPhi, targX, targY;
    Double_t p_t_lambda;

    TTree* st = (TTree*) f.Get("LambdasMCSimTracks");
    TTree* rt = (TTree*) f.Get("LambdasMCRecTracks");

    st->Print();
    rt->Print();

    st->SetBranchAddress("cos_x", &cos_x);
    st->SetBranchAddress("cos_y", &cos_y);
    st->SetBranchAddress("cos_z", &cos_z);
    st->SetBranchAddress("y", &y);
    st->SetBranchAddress("p_t_lambda", &p_t_lambda);

    rt->SetBranchAddress("cos_x", &cos_x);
    rt->SetBranchAddress("cos_y", &cos_y);
    rt->SetBranchAddress("cos_z", &cos_z);
    rt->SetBranchAddress("nClusters1", &nClusters1);
    rt->SetBranchAddress("nClusters2", &nClusters2);
    rt->SetBranchAddress("DeltaZ", &DeltaZ);
    rt->SetBranchAddress("y", &y);
    rt->SetBranchAddress("cosPhi", &cosPhi);
    rt->SetBranchAddress("p_pi_mass", &p_pi_mass);
    rt->SetBranchAddress("p_t_lambda", &p_t_lambda);
    rt->SetBranchAddress("targX", &targX);
    rt->SetBranchAddress("targY", &targY);

    const int MCgen[10][10][10][10] = {}, MCacc[10][10][10][10] = {};

    const Long64_t stentries = st->GetEntries();
    for (Long64_t i = 0; i < stentries; ++i) {
        if (i % 100000 == 0)
            cout << "SIM Entry: " << i << endl;
        st->GetEntry(i);
        int bin1 = p_t_lambda * 5;
        int bin1 = (y + 2.5) * 2;
        int bin3 = (cos_x + 1) * 5;
        int bin4 = atan2(cos_y, cos_z) * 4.9999999 / pi + 5;
        if (bin1 < 0 || bin1 > 9) continue;
        if (bin2 < 0 || bin2 > 9) continue;
        if (bin3 < 0 || bin3 > 9) continue;
        if (bin4 < 0 || bin4 > 9) continue;
        ++MCgen[bin1][bin2][bin3][bin4];
    }

    const Long64_t rtentries = rt->GetEntries();
    for (Long64_t i = 0; i < rtentries; ++i) {
        if (i % 100000 == 0)
            cout << "REC Entry: " << i << endl;
        rt->GetEntry(i);

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
        int bin1 = p_t_lambda * 5;
        int bin1 = (y + 2.5) * 2;
        int bin3 = (cos_x + 1) * 5;
        int bin4 = atan2(cos_y, cos_z) * 4.9999999 / pi + 5;
        if (bin1 < 0 || bin1 > 9) continue;
        if (bin2 < 0 || bin2 > 9) continue;
        if (bin3 < 0 || bin3 > 9) continue;
        if (bin4 < 0 || bin4 > 9) continue;
        ++MCacc[bin1][bin2][bin3][bin4];
        hist_p_pi_mass->FIll(p_pi_mass);
    }
    
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            for (int k = 0; k < 10; ++k) {
                for (int l = 0; l < 10; ++l) {
                    if (MCgen[i][j][k][l] > 3)
                        coeff = MCgen[i][j][k][l] / MCacc[i][j][k][l];
                }
            }
        }
    }
    
    gStyle->SetOptStat(0);
    TF1 *func2 = new TF1("lorentzianPeak", lorentzianPeak, min_mass, max_mass, 3);
    func2->SetParNames ("Norm.Const", "m_{0}", "#Gamma");
    func2->SetParameters(hist_p_pi_mass->GetEntries() / hist_p_pi_mass->GetNbinsX(), 1.11566, 2.59380e-03);
    func2->SetNpx(2000);
    func2->SetLineWidth(1);
    func2->SetLineColor(kBlue);

    auto fitresult_v0Mass = hist_p_pi_mass->Fit("lorentzianPeak");
    TF1 *fit_v0Mass = hist_p_pi_mass->GetFunction("lorentzianPeak");
    //hist_p_pi_mass->SetLineStyle(2);
    //hist_p_pi_mass->SetLineWidth(2);
    //hist_p_pi_mass->Draw("P");
    //hist_p_pi_mass->SetMarkerSize(3);
    hist_p_pi_mass->GetListOfFunctions()->FindObject("lorentzianPeak")->Draw("same");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_p_pi_mass->GetXaxis()->SetLabelSize(font_size);
    hist_p_pi_mass->GetYaxis()->SetLabelSize(font_size);
    hist_p_pi_mass->GetXaxis()->SetTitleSize(font_size);
    hist_p_pi_mass->GetYaxis()->SetTitleSize(font_size);
    hist_p_pi_mass->GetYaxis()->SetTitleOffset(title_offset);
    hist_p_pi_mass->GetXaxis()->SetNdivisions(508);
    c0->SaveAs("hist_p_pi_mass_rec.pdf");

}
