const Double_t alpha = 0.732;
// Define a function with 1 parameter
Double_t fitf(const Double_t *x, const Double_t *par) {
    return (1. + alpha * par[0] * x[0]) / 2.;
}
// Define a function with 1 parameter
Double_t pol0(const Double_t *x, const Double_t *par) {
    return par[0];
}
/*
hist_cos_x - bare data
hist_cos_x_cut - DeltaZ and nHits cuts
hist_cos_x_cut1 - +XYTargetPass cut
hist_cos_x_cut2 - +cosPhi cut
hist_cos_x_cut3 - +XYTargetPass+cosPhi cut
*/
void treeplot_sim() {
    TCanvas *c0 = new TCanvas("c0", "canvas", 0, 0, 1000, 800);
    gStyle->SetOptStat(0);
    TFile f("y-sel-sim.1606.pp158.root");
    //TH1D* hist_invmass = new TH1D("InvMass", ";m_{inv},GeV", 200, 0, 1.5);
    //TH2D* ap = new TH2D("APplot", "APplot;alpha;p_{T}, GeV", Nbins, -1., 1., 200, 0.0, 0.5);
    const int Nbins = 100;
    const double max_pt = 1.5, max_y = 2.8, max_deltaz = 100;

    TH1D* hist_cos_x = new TH1D("COS_X", ";cos #theta_{x};f(cos #theta_{x})", Nbins, -1., 1.);
    TH1D* hist_cos_y = new TH1D("COS_Y", ";cos #theta_{y};f(cos #theta_{y})", Nbins, -1., 1.);
    TH1D* hist_cos_z = new TH1D("COS_Z", ";cos #theta_{z};f(cos #theta_{z})", Nbins, -1., 1.);

    TH1D* hist_cos_x_cut = new TH1D("COS_X_cut", ";cos #theta_{x};f(cos #theta_{x})", Nbins, -1., 1.);
    TH1D* hist_cos_y_cut = new TH1D("COS_Y_cut", ";cos #theta_{y};f(cos #theta_{y})", Nbins, -1., 1.);
    TH1D* hist_cos_z_cut = new TH1D("COS_Z_cut", ";cos #theta_{z};f(cos #theta_{z})", Nbins, -1., 1.);

    TH1D* hist_cos_x_cut1 = new TH1D("COS_X_cut1", ";cos #theta_{x};f(cos #theta_{x})", Nbins, -1., 1.);
    TH1D* hist_cos_y_cut1 = new TH1D("COS_Y_cut1", ";cos #theta_{y};f(cos #theta_{y})", Nbins, -1., 1.);
    TH1D* hist_cos_z_cut1 = new TH1D("COS_Z_cut1", ";cos #theta_{z};f(cos #theta_{z})", Nbins, -1., 1.);

    TH1D* hist_cos_x_cut2 = new TH1D("COS_X_cut2", ";cos #theta_{x};f(cos #theta_{x})", Nbins, -1., 1.);
    TH1D* hist_cos_y_cut2 = new TH1D("COS_Y_cut2", ";cos #theta_{y};f(cos #theta_{y})", Nbins, -1., 1.);
    TH1D* hist_cos_z_cut2 = new TH1D("COS_Z_cut2", ";cos #theta_{z};f(cos #theta_{z})", Nbins, -1., 1.);

    TH1D* hist_cos_x_cut3 = new TH1D("COS_X_cut3", ";cos #theta_{x};f(cos #theta_{x})", Nbins, -1., 1.);
    TH1D* hist_cos_y_cut3 = new TH1D("COS_Y_cut3", ";cos #theta_{y};f(cos #theta_{y})", Nbins, -1., 1.);
    TH1D* hist_cos_z_cut3 = new TH1D("COS_Z_cut3", ";cos #theta_{z};f(cos #theta_{z})", Nbins, -1., 1.);

    TH2D* hist_ypt = new TH2D("ypt", "y-p_{T} distribution;y;p_{T}, [GeV/c]", 200, -max_y, max_y, 150, 0.0, max_pt);
    TH1D* hist_y = new TH1D("y", "y distribution;y;", 200, -max_y, max_y);
    TH1D* hist_pt = new TH1D("pt", "p_{T} distribution;p_{T}, [GeV/c];", 150, 0.0, max_pt);
    TH1D* hist_phi = new TH1D("phi", "#phi distribution;#phi", 150, -3.14, 3.14);
    TH1D* hist_z = new TH1D("z", "#Deltaz distribution;#Deltaz,[cm]", 200, 0.0, max_deltaz);
    TH2D* hist_y_z = new TH2D("y_z", "y-#Deltaz distribution;y;#Deltaz,[cm]", 200, -max_y, max_y, 200, 0.0, max_deltaz);

    TH2D* hist_ypt_cut = new TH2D("ypt_cut", "y-p_{T} distribution;y;p_{T}, [GeV/c]", 200, -max_y, max_y, 150, 0.0, max_pt);
    TH1D* hist_y_cut = new TH1D("y_cut", "y distribution;y;", 200, -max_y, max_y);
    TH1D* hist_pt_cut = new TH1D("pt_cut", "p_{T} distribution;p_{T}, [GeV/c];", 150, 0.0, max_pt);
    TH1D* hist_phi_cut = new TH1D("phi_cut", "#phi distribution;#phi", 150, -3.14, 3.14);
    TH1D* hist_z_cut = new TH1D("z_cut", "#Deltaz distribution;#Deltaz,[cm]", 200, 0.0, max_deltaz);
    TH2D* hist_y_z_cut = new TH2D("y_z_cut", "y-#Deltaz distribution;y;#Deltaz,[cm]", 200, -max_y, max_y, 200, 0.0, max_deltaz);

    TH2D* hist_ypt_cut1 = new TH2D("ypt_cut1", "y-p_{T} distribution;y;p_{T}, [GeV/c]", 200, -max_y, max_y, 150, 0.0, max_pt);
    TH1D* hist_y_cut1 = new TH1D("y_cut1", "y distribution;y;", 200, -max_y, max_y);
    TH1D* hist_pt_cut1 = new TH1D("pt_cut1", "p_{T} distribution;p_{T}, [GeV/c];", 150, 0.0, max_pt);
    TH1D* hist_phi_cut1 = new TH1D("phi_cut1", "#phi distribution;#phi", 150, -3.14, 3.14);
    TH1D* hist_z_cut1 = new TH1D("z_cut1", "#Deltaz distribution;#Deltaz,[cm]", 200, 0.0, max_deltaz);
    TH2D* hist_y_z_cut1 = new TH2D("y_z_cut1", "y-#Deltaz distribution;y;#Deltaz,[cm]", 200, -max_y, max_y, 200, 0.0, max_deltaz);

    TH2D* hist_ypt_cut2 = new TH2D("ypt_cut2", "y-p_{T} distribution;y;p_{T}, [GeV/c]", 200, -max_y, max_y, 150, 0.0, max_pt);
    TH1D* hist_y_cut2 = new TH1D("y_cut2", "y distribution;y;", 200, -max_y, max_y);
    TH1D* hist_pt_cut2 = new TH1D("pt_cut2", "p_{T} distribution;p_{T}, [GeV/c];", 150, 0.0, max_pt);
    TH1D* hist_phi_cut2 = new TH1D("phi_cut2", "#phi distribution;#phi", 150, -3.14, 3.14);
    TH1D* hist_z_cut2 = new TH1D("z_cut2", "#Deltaz distribution;#Deltaz,[cm]", 200, 0.0, max_deltaz);
    TH2D* hist_y_z_cut2 = new TH2D("y_z_cut2", "y-#Deltaz distribution;y;#Deltaz,[cm]", 200, -max_y, max_y, 200, 0.0, max_deltaz);

    TH2D* hist_ypt_cut3 = new TH2D("ypt_cut3", "y-p_{T} distribution;y;p_{T}, [GeV/c]", 200, -max_y, max_y, 150, 0.0, max_pt);
    TH1D* hist_y_cut3 = new TH1D("y_cut3", "y distribution;y;", 200, -max_y, max_y);
    TH1D* hist_pt_cut3 = new TH1D("pt_cut3", "p_{T} distribution;p_{T}, [GeV/c];", 150, 0.0, max_pt);
    TH1D* hist_phi_cut3 = new TH1D("phi_cut3", "#phi distribution;#phi", 150, -3.14, 3.14);
    TH1D* hist_z_cut3 = new TH1D("z_cut3", "#Deltaz distribution;#Deltaz,[cm]", 200, 0.0, max_deltaz);
    TH2D* hist_y_z_cut3 = new TH2D("y_z_cut3", "y-#Deltaz distribution;y;#Deltaz,[cm]", 200, -max_y, max_y, 200, 0.0, max_deltaz);

    Double_t cos_x, cos_y, cos_z, DeltaZ, y, p_full_lambda, cosPhi;
    Double_t alpha, p_t, p_t_lambda, phi_lambda;
    UInt_t nHits11, nHits12, nHits21, nHits22;
    UInt_t XYTargetPass;

    double max_hist;
    vector<double> maxes;

    TTree* t = (TTree*) f.Get("Lambdas");
    t->Print();
    t->SetBranchAddress("cos_x", &cos_x);
    t->SetBranchAddress("cos_y", &cos_y);
    t->SetBranchAddress("cos_z", &cos_z);
    t->SetBranchAddress("nHits11", &nHits11);
    t->SetBranchAddress("nHits12", &nHits12);
    t->SetBranchAddress("nHits21", &nHits21);
    t->SetBranchAddress("nHits22", &nHits22);
    t->SetBranchAddress("XYTargetPass", &XYTargetPass);
    t->SetBranchAddress("DeltaZ", &DeltaZ);
    t->SetBranchAddress("y", &y);
    t->SetBranchAddress("alpha", &alpha);
    t->SetBranchAddress("p_t", &p_t);
    t->SetBranchAddress("p_t_lambda", &p_t_lambda);
    t->SetBranchAddress("p_full_lambda", &p_full_lambda);
    t->SetBranchAddress("phi_lambda", &phi_lambda);
    t->SetBranchAddress("cosPhi", &cosPhi);

    bool deltazy_cut, nhits_cut;
    const Long64_t nentries = t->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        if (i % 100000 == 0)
            cout << "Entry: " << i << endl;
        t->GetEntry(i);
        hist_y->Fill(y);
        hist_ypt->Fill(y, p_t_lambda);
        hist_pt->Fill(p_t_lambda);
        hist_phi->Fill(phi_lambda);
        hist_z->Fill(DeltaZ);
        hist_y_z->Fill(y, DeltaZ);

        hist_cos_x->Fill(cos_x);
        hist_cos_y->Fill(cos_y);
        hist_cos_z->Fill(cos_z);

        nhits_cut = true;
        if ((nHits11 < 15) && (nHits12 < 15)) {
            nhits_cut = false;
        }
        if ((nHits21 < 15) && (nHits22 < 15)) {
            nhits_cut = false;
        }
        deltazy_cut = false;
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
        if (deltazy_cut & nhits_cut)
        {
            hist_cos_x_cut->Fill(cos_x);
            hist_cos_y_cut->Fill(cos_y);
            hist_cos_z_cut->Fill(cos_z);

            hist_y_cut->Fill(y);
            hist_ypt_cut->Fill(y, p_t_lambda);
            hist_pt_cut->Fill(p_t_lambda);
            hist_phi_cut->Fill(phi_lambda);
            hist_z_cut->Fill(DeltaZ);
            hist_y_z_cut->Fill(y, DeltaZ);

            if (XYTargetPass == 1)
            {
                hist_cos_x_cut1->Fill(cos_x);
                hist_cos_y_cut1->Fill(cos_y);
                hist_cos_z_cut1->Fill(cos_z);

                hist_y_cut1->Fill(y);
                hist_ypt_cut1->Fill(y, p_t_lambda);
                hist_pt_cut1->Fill(p_t_lambda);
                hist_phi_cut1->Fill(phi_lambda);
                hist_z_cut1->Fill(DeltaZ);
                hist_y_z_cut1->Fill(y, DeltaZ);
            }
            if (cosPhi_cut)
            {
                hist_cos_x_cut2->Fill(cos_x);
                hist_cos_y_cut2->Fill(cos_y);
                hist_cos_z_cut2->Fill(cos_z);

                hist_y_cut2->Fill(y);
                hist_ypt_cut2->Fill(y, p_t_lambda);
                hist_pt_cut2->Fill(p_t_lambda);
                hist_phi_cut2->Fill(phi_lambda);
                hist_z_cut2->Fill(DeltaZ);
                hist_y_z_cut2->Fill(y, DeltaZ);
            }
            if ((XYTargetPass == 1) && cosPhi_cut)
            {
                hist_cos_x_cut3->Fill(cos_x);
                hist_cos_y_cut3->Fill(cos_y);
                hist_cos_z_cut3->Fill(cos_z);

                hist_y_cut3->Fill(y);
                hist_ypt_cut3->Fill(y, p_t_lambda);
                hist_pt_cut3->Fill(p_t_lambda);
                hist_phi_cut3->Fill(phi_lambda);
                hist_z_cut3->Fill(DeltaZ);
                hist_y_z_cut3->Fill(y, DeltaZ);
            }
        }
    }

    const Float_t font_size = 0.05, title_offset = 1.2, small_margin = 0.05, big_margin = 0.12, medium_margin = 0.09;
    const Int_t n_divisions = 805;

    hist_y->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_y->GetXaxis()->SetLabelSize(font_size);
    hist_y->GetYaxis()->SetLabelSize(font_size);
    hist_y->GetXaxis()->SetTitleSize(font_size);
    hist_y->GetYaxis()->SetTitleSize(font_size);
    hist_y->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_y_sim.pdf");

    hist_pt->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_pt->GetXaxis()->SetLabelSize(font_size);
    hist_pt->GetYaxis()->SetLabelSize(font_size);
    hist_pt->GetXaxis()->SetTitleSize(font_size);
    hist_pt->GetYaxis()->SetTitleSize(font_size);
    hist_pt->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_pt_sim.pdf");

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
    c0->SaveAs("hist_ypt_sim.pdf");

    // TODO: FIT TO STRAIGHT LINE
    TF1 *pol0f = new TF1("pol0f", pol0, -1, 1, 1);
    hist_phi->Fit("pol0f");
    hist_phi->Draw("hist");
    hist_phi->GetListOfFunctions()->FindObject("pol0f")->Draw("same");
    hist_phi->SetMinimum(0);
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_phi->GetXaxis()->SetLabelSize(font_size);
    hist_phi->GetYaxis()->SetLabelSize(font_size);
    hist_phi->GetXaxis()->SetTitleSize(font_size);
    hist_phi->GetYaxis()->SetTitleSize(font_size);
    hist_phi->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_phi_sim.pdf");

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
    c0->SaveAs("hist_delta_z_sim.pdf");

    hist_y_z->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_y_z->GetXaxis()->SetLabelSize(font_size);
    hist_y_z->GetYaxis()->SetLabelSize(font_size);
    hist_y_z->GetXaxis()->SetTitleSize(font_size);
    hist_y_z->GetYaxis()->SetTitleSize(font_size);
    hist_y_z->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_y_z_sim.pdf");


    hist_y_cut->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_y_cut->GetXaxis()->SetLabelSize(font_size);
    hist_y_cut->GetYaxis()->SetLabelSize(font_size);
    hist_y_cut->GetXaxis()->SetTitleSize(font_size);
    hist_y_cut->GetYaxis()->SetTitleSize(font_size);
    hist_y_cut->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_y_cut_sim.pdf");

    hist_pt_cut->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_pt_cut->GetXaxis()->SetLabelSize(font_size);
    hist_pt_cut->GetYaxis()->SetLabelSize(font_size);
    hist_pt_cut->GetXaxis()->SetTitleSize(font_size);
    hist_pt_cut->GetYaxis()->SetTitleSize(font_size);
    hist_pt_cut->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_pt_cut_sim.pdf");

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
    c0->SaveAs("hist_ypt_cut_sim.pdf");

    // TODO: FIT TO STRAIGHT LINE
    hist_phi_cut->Fit("pol0f");
    hist_phi_cut->Draw("hist");
    hist_phi_cut->GetListOfFunctions()->FindObject("pol0f")->Draw("same");
    hist_phi_cut->SetMinimum(0);
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_phi_cut->GetXaxis()->SetLabelSize(font_size);
    hist_phi_cut->GetYaxis()->SetLabelSize(font_size);
    hist_phi_cut->GetXaxis()->SetTitleSize(font_size);
    hist_phi_cut->GetYaxis()->SetTitleSize(font_size);
    hist_phi_cut->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_phi_cut_sim.pdf");

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
    c0->SaveAs("hist_delta_z_cut_sim.pdf");

    hist_y_z_cut->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_y_z_cut->GetXaxis()->SetLabelSize(font_size);
    hist_y_z_cut->GetYaxis()->SetLabelSize(font_size);
    hist_y_z_cut->GetXaxis()->SetTitleSize(font_size);
    hist_y_z_cut->GetYaxis()->SetTitleSize(font_size);
    hist_y_z_cut->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_y_z_cut_sim.pdf");

    //cut1
    
    hist_y_cut1->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_y_cut1->GetXaxis()->SetLabelSize(font_size);
    hist_y_cut1->GetYaxis()->SetLabelSize(font_size);
    hist_y_cut1->GetXaxis()->SetTitleSize(font_size);
    hist_y_cut1->GetYaxis()->SetTitleSize(font_size);
    hist_y_cut1->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_y_cut1_sim.pdf");

    hist_pt_cut1->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_pt_cut1->GetXaxis()->SetLabelSize(font_size);
    hist_pt_cut1->GetYaxis()->SetLabelSize(font_size);
    hist_pt_cut1->GetXaxis()->SetTitleSize(font_size);
    hist_pt_cut1->GetYaxis()->SetTitleSize(font_size);
    hist_pt_cut1->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_pt_cut1_sim.pdf");

    hist_ypt_cut1->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_ypt_cut1->GetXaxis()->SetLabelSize(font_size);
    hist_ypt_cut1->GetYaxis()->SetLabelSize(font_size);
    hist_ypt_cut1->GetXaxis()->SetTitleSize(font_size);
    hist_ypt_cut1->GetYaxis()->SetTitleSize(font_size);
    hist_ypt_cut1->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_ypt_cut1_sim.pdf");

    // TODO: FIT TO STRAIGHT LINE
    hist_phi_cut1->Fit("pol0f");
    hist_phi_cut1->Draw("hist");
    hist_phi_cut1->GetListOfFunctions()->FindObject("pol0f")->Draw("same");
    hist_phi_cut1->SetMinimum(0);
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_phi_cut1->GetXaxis()->SetLabelSize(font_size);
    hist_phi_cut1->GetYaxis()->SetLabelSize(font_size);
    hist_phi_cut1->GetXaxis()->SetTitleSize(font_size);
    hist_phi_cut1->GetYaxis()->SetTitleSize(font_size);
    hist_phi_cut1->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_phi_cut1_sim.pdf");

    hist_z_cut1->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_z_cut1->GetXaxis()->SetLabelSize(font_size);
    hist_z_cut1->GetYaxis()->SetLabelSize(font_size);
    hist_z_cut1->GetXaxis()->SetTitleSize(font_size);
    hist_z_cut1->GetYaxis()->SetTitleSize(font_size);
    hist_z_cut1->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_delta_z_cut1_sim.pdf");

    hist_y_z_cut1->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_y_z_cut1->GetXaxis()->SetLabelSize(font_size);
    hist_y_z_cut1->GetYaxis()->SetLabelSize(font_size);
    hist_y_z_cut1->GetXaxis()->SetTitleSize(font_size);
    hist_y_z_cut1->GetYaxis()->SetTitleSize(font_size);
    hist_y_z_cut1->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_y_z_cut1_sim.pdf");

    //cut2

    hist_y_cut2->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_y_cut2->GetXaxis()->SetLabelSize(font_size);
    hist_y_cut2->GetYaxis()->SetLabelSize(font_size);
    hist_y_cut2->GetXaxis()->SetTitleSize(font_size);
    hist_y_cut2->GetYaxis()->SetTitleSize(font_size);
    hist_y_cut2->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_y_cut2_sim.pdf");

    hist_pt_cut2->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_pt_cut2->GetXaxis()->SetLabelSize(font_size);
    hist_pt_cut2->GetYaxis()->SetLabelSize(font_size);
    hist_pt_cut2->GetXaxis()->SetTitleSize(font_size);
    hist_pt_cut2->GetYaxis()->SetTitleSize(font_size);
    hist_pt_cut2->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_pt_cut2_sim.pdf");

    hist_ypt_cut2->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_ypt_cut2->GetXaxis()->SetLabelSize(font_size);
    hist_ypt_cut2->GetYaxis()->SetLabelSize(font_size);
    hist_ypt_cut2->GetXaxis()->SetTitleSize(font_size);
    hist_ypt_cut2->GetYaxis()->SetTitleSize(font_size);
    hist_ypt_cut2->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_ypt_cut2_sim.pdf");

    // TODO: FIT TO STRAIGHT LINE
    hist_phi_cut2->Fit("pol0f");
    hist_phi_cut2->Draw("hist");
    hist_phi_cut2->GetListOfFunctions()->FindObject("pol0f")->Draw("same");
    hist_phi_cut2->SetMinimum(0);
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_phi_cut2->GetXaxis()->SetLabelSize(font_size);
    hist_phi_cut2->GetYaxis()->SetLabelSize(font_size);
    hist_phi_cut2->GetXaxis()->SetTitleSize(font_size);
    hist_phi_cut2->GetYaxis()->SetTitleSize(font_size);
    hist_phi_cut2->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_phi_cut2_sim.pdf");

    hist_z_cut2->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_z_cut2->GetXaxis()->SetLabelSize(font_size);
    hist_z_cut2->GetYaxis()->SetLabelSize(font_size);
    hist_z_cut2->GetXaxis()->SetTitleSize(font_size);
    hist_z_cut2->GetYaxis()->SetTitleSize(font_size);
    hist_z_cut2->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_delta_z_cut2_sim.pdf");

    hist_y_z_cut2->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_y_z_cut2->GetXaxis()->SetLabelSize(font_size);
    hist_y_z_cut2->GetYaxis()->SetLabelSize(font_size);
    hist_y_z_cut2->GetXaxis()->SetTitleSize(font_size);
    hist_y_z_cut2->GetYaxis()->SetTitleSize(font_size);
    hist_y_z_cut2->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_y_z_cut2_sim.pdf");

    //cut3

    hist_y_cut3->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_y_cut3->GetXaxis()->SetLabelSize(font_size);
    hist_y_cut3->GetYaxis()->SetLabelSize(font_size);
    hist_y_cut3->GetXaxis()->SetTitleSize(font_size);
    hist_y_cut3->GetYaxis()->SetTitleSize(font_size);
    hist_y_cut3->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_y_cut3_sim.pdf");

    hist_pt_cut3->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_pt_cut3->GetXaxis()->SetLabelSize(font_size);
    hist_pt_cut3->GetYaxis()->SetLabelSize(font_size);
    hist_pt_cut3->GetXaxis()->SetTitleSize(font_size);
    hist_pt_cut3->GetYaxis()->SetTitleSize(font_size);
    hist_pt_cut3->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_pt_cut3_sim.pdf");

    hist_ypt_cut3->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(big_margin);
    gPad->SetTopMargin(medium_margin);
    hist_ypt_cut3->GetXaxis()->SetLabelSize(font_size);
    hist_ypt_cut3->GetYaxis()->SetLabelSize(font_size);
    hist_ypt_cut3->GetXaxis()->SetTitleSize(font_size);
    hist_ypt_cut3->GetYaxis()->SetTitleSize(font_size);
    hist_ypt_cut3->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_ypt_cut3_sim.pdf");

    // TODO: FIT TO STRAIGHT LINE
    hist_phi_cut3->Fit("pol0f");
    hist_phi_cut3->Draw("hist");
    hist_phi_cut3->GetListOfFunctions()->FindObject("pol0f")->Draw("same");
    hist_phi_cut3->SetMinimum(0);
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_phi_cut3->GetXaxis()->SetLabelSize(font_size);
    hist_phi_cut3->GetYaxis()->SetLabelSize(font_size);
    hist_phi_cut3->GetXaxis()->SetTitleSize(font_size);
    hist_phi_cut3->GetYaxis()->SetTitleSize(font_size);
    hist_phi_cut3->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_phi_cut3_sim.pdf");

    hist_z_cut3->Draw("hist");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_z_cut3->GetXaxis()->SetLabelSize(font_size);
    hist_z_cut3->GetYaxis()->SetLabelSize(font_size);
    hist_z_cut3->GetXaxis()->SetTitleSize(font_size);
    hist_z_cut3->GetYaxis()->SetTitleSize(font_size);
    hist_z_cut3->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_delta_z_cut3_sim.pdf");

    hist_y_z_cut3->Draw("COLZ");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(medium_margin);
    hist_y_z_cut3->GetXaxis()->SetLabelSize(font_size);
    hist_y_z_cut3->GetYaxis()->SetLabelSize(font_size);
    hist_y_z_cut3->GetXaxis()->SetTitleSize(font_size);
    hist_y_z_cut3->GetYaxis()->SetTitleSize(font_size);
    hist_y_z_cut3->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_y_z_cut3_sim.pdf");
    /*
    hist_invmass->Draw("hist");
    c0->SaveAs("hist_invmass.pdf");

    ap->Draw("COLZ");
    c0->SaveAs("hist_ap_plot.pdf");
    */
    TF1 *func = new TF1("fit", fitf, -1, 1, 1);

    // total cosines distributions
    cout << "HISTOGRAMS TOTAL :  " << hist_cos_x->GetEntries() << endl;
    hist_cos_x->Scale(hist_cos_x->GetNbinsX() / (2 * hist_cos_x->GetEntries()));
    auto fitresultx = hist_cos_x->Fit("fit");
    TF1 *fitx = hist_cos_x->GetFunction("fit");
    Double_t x_p1 = fitx->GetParameter(0);
    Double_t x_e1 = fitx->GetParError(0);

    hist_cos_y->Scale(hist_cos_y->GetNbinsX() / (2 * hist_cos_y->GetEntries()));
    auto fitresulty = hist_cos_y->Fit("fit");
    TF1 *fity = hist_cos_y->GetFunction("fit");
    Double_t y_p1 = fity->GetParameter(0);
    Double_t y_e1 = fity->GetParError(0);

    hist_cos_z->Scale(hist_cos_z->GetNbinsX() / (2 * hist_cos_z->GetEntries()));
    auto fitresultz = hist_cos_z->Fit("fit");
    TF1 *fitz = hist_cos_z->GetFunction("fit");
    Double_t z_p1 = fitz->GetParameter(0);
    Double_t z_e1 = fitz->GetParError(0);

    maxes = {hist_cos_x->GetMaximum(), hist_cos_y->GetMaximum(), hist_cos_z->GetMaximum()};
    max_hist = 1.05 * (*max_element(maxes.begin(), maxes.end()));

    hist_cos_x->SetMaximum(max_hist);
    hist_cos_x->SetMinimum(0);
    hist_cos_x->Draw("hist");
    hist_cos_x->GetListOfFunctions()->FindObject("fit")->Draw("same");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(small_margin);
    hist_cos_x->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_x->GetXaxis()->SetLabelSize(font_size);
    hist_cos_x->GetYaxis()->SetLabelSize(font_size);
    hist_cos_x->GetXaxis()->SetTitleSize(font_size);
    hist_cos_x->GetYaxis()->SetTitleSize(font_size);
    hist_cos_x->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_cos_x_sim.pdf");

    hist_cos_y->SetMaximum(max_hist);
    hist_cos_y->SetMinimum(0);
    hist_cos_y->Draw("hist");
    hist_cos_y->GetListOfFunctions()->FindObject("fit")->Draw("same");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(small_margin);
    hist_cos_y->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_y->GetXaxis()->SetLabelSize(font_size);
    hist_cos_y->GetYaxis()->SetLabelSize(font_size);
    hist_cos_y->GetXaxis()->SetTitleSize(font_size);
    hist_cos_y->GetYaxis()->SetTitleSize(font_size);
    hist_cos_y->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_cos_y_sim.pdf");

    hist_cos_z->SetMaximum(max_hist);
    hist_cos_z->SetMinimum(0);
    hist_cos_z->Draw("hist");
    hist_cos_z->GetListOfFunctions()->FindObject("fit")->Draw("same");
    gPad->SetLeftMargin(big_margin);
    gPad->SetBottomMargin(big_margin);
    gPad->SetRightMargin(small_margin);
    gPad->SetTopMargin(small_margin);
    hist_cos_z->GetXaxis()->SetNdivisions(n_divisions);
    hist_cos_z->GetXaxis()->SetLabelSize(font_size);
    hist_cos_z->GetYaxis()->SetLabelSize(font_size);
    hist_cos_z->GetXaxis()->SetTitleSize(font_size);
    hist_cos_z->GetYaxis()->SetTitleSize(font_size);
    hist_cos_z->GetYaxis()->SetTitleOffset(title_offset);
    c0->SaveAs("hist_cos_z_sim.pdf");
}