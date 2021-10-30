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
void treeplot_rec() {
	TCanvas *c0 = new TCanvas("c0", "canvas", 0, 0, 1000, 800);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(11);

	TFile f("y-selection-rec.pp158.1307.root");
	double max_hist = 0.8, min_hist = 0.3;
	const Float_t font_size = 0.05, title_offset = 1.2, small_margin = 0.05, big_margin = 0.12, medium_margin = 0.09;
	const Int_t n_divisions = 805;
	const int Nbins = 100;
	const double min_mass = 1.08, max_mass = 1.25;
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

	TH1D* hist_p_pi_mass = new TH1D("hist_p_pi_mass", ";m_{p#pi^{-}} (GeV/c^{2});Entries", 3 * Nbins, min_mass, max_mass);
	TH1D* hist_p_pi_mass_cut = new TH1D("hist_p_pi_mass_cut", ";m_{p#pi^{-}} (GeV/c^{2});Entries", 3 * Nbins, min_mass, max_mass);
	TH1D* hist_p_pi_mass_cut1 = new TH1D("hist_p_pi_mass_cut1", ";m_{p#pi^{-}} (GeV/c^{2});Entries", 3 * Nbins, min_mass, max_mass);
	TH1D* hist_p_pi_mass_cut2 = new TH1D("hist_p_pi_mass_cut2", ";m_{p#pi^{-}} (GeV/c^{2});Entries", 3 * Nbins, min_mass, max_mass);
	TH1D* hist_p_pi_mass_cut3 = new TH1D("hist_p_pi_mass_cut3", ";m_{p#pi^{-}} (GeV/c^{2});Entries", 3 * Nbins, min_mass, max_mass);

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

	Double_t cos_x, cos_y, cos_z, DeltaZ, y;
	//double max_hist, min_hist;
	//vector<double> maxes;
	UInt_t nHits1, nHits2, nClusters1, nClusters2;
	Double_t v0Mass, p_pi_mass, cosPhi;
	UInt_t nClustersPass, XYTargetPass;
	UInt_t nS4;
	Double_t alpha, p_t;
	Double_t p_full_lambda, p_t_lambda, phi_lambda;

	TTree* t = (TTree*) f.Get("Lambdas");
	t->Print();
	t->SetBranchAddress("cos_x", &cos_x);
	t->SetBranchAddress("cos_y", &cos_y);
	t->SetBranchAddress("cos_z", &cos_z);
	//t->SetBranchAddress("nHits1", &nHits1);
	//t->SetBranchAddress("nHits2", &nHits2);
	//t->SetBranchAddress("nClusters1", &nClusters1);
	//t->SetBranchAddress("nClusters2", &nClusters2);
	t->SetBranchAddress("DeltaZ", &DeltaZ);
	t->SetBranchAddress("y", &y);
	t->SetBranchAddress("cosPhi", &cosPhi);
	t->SetBranchAddress("p_pi_mass", &p_pi_mass);
	t->SetBranchAddress("nS4", &nS4);
	t->SetBranchAddress("nClustersPass", &nClustersPass);
	t->SetBranchAddress("XYTargetPass", &XYTargetPass);
	t->SetBranchAddress("alpha", &alpha);
	t->SetBranchAddress("p_t", &p_t);
	t->SetBranchAddress("p_t_lambda", &p_t_lambda);
	t->SetBranchAddress("p_full_lambda", &p_full_lambda);
	t->SetBranchAddress("phi_lambda", &phi_lambda);

	const Long64_t nentries = t->GetEntries();
	for (Long64_t i = 0; i < nentries; ++i) {
		if (i % 100000 == 0)
			cout << "Entry: " << i << endl;
		t->GetEntry(i);
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
		hist_cos_x->Fill(cos_x);
		hist_cos_y->Fill(cos_y);
		hist_cos_z->Fill(cos_z);
		hist_p_pi_mass->Fill(p_pi_mass);
		hist_y->Fill(y);
		hist_ypt->Fill(y, p_t_lambda);
		hist_pt->Fill(p_t_lambda);
		hist_phi->Fill(phi_lambda);
		hist_z->Fill(DeltaZ);
		hist_y_z->Fill(y, DeltaZ);
		if (nClustersPass == 1)
		{
			hist_cos_x_cut->Fill(cos_x);
			hist_cos_y_cut->Fill(cos_y);
			hist_cos_z_cut->Fill(cos_z);
			hist_p_pi_mass_cut->Fill(p_pi_mass);
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
				hist_p_pi_mass_cut1->Fill(p_pi_mass);
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
				hist_p_pi_mass_cut2->Fill(p_pi_mass);
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
				hist_p_pi_mass_cut3->Fill(p_pi_mass);
				hist_y_cut3->Fill(y);
				hist_ypt_cut3->Fill(y, p_t_lambda);
				hist_pt_cut3->Fill(p_t_lambda);
				hist_phi_cut3->Fill(phi_lambda);
				hist_z_cut3->Fill(DeltaZ);
				hist_y_z_cut3->Fill(y, DeltaZ);
			}
		}
	}
	TF1 *func = new TF1("fit", fitf, -1, 1, 1);

	// total cosines distributions
	//cout << hist_cos_x->GetEntries() << endl;
	func->SetParNames("P_{x}:");
	hist_cos_x->Scale(hist_cos_x->GetNbinsX() / (2 * hist_cos_x->GetEntries()));
	auto fitresultx = hist_cos_x->Fit("fit");
	TF1 *fitx = hist_cos_x->GetFunction("fit");
	//Double_t x_p1 = fitx->GetParameter(0);
	//Double_t x_e1 = fitx->GetParError(0);
	func->SetParNames("P_{y}:");
	hist_cos_y->Scale(hist_cos_y->GetNbinsX() / (2 * hist_cos_y->GetEntries()));
	auto fitresulty = hist_cos_y->Fit("fit");
	TF1 *fity = hist_cos_y->GetFunction("fit");
	//Double_t y_p1 = fity->GetParameter(0);
	//Double_t y_e1 = fity->GetParError(0);
	func->SetParNames("P_{z}:");
	hist_cos_z->Scale(hist_cos_z->GetNbinsX() / (2 * hist_cos_z->GetEntries()));
	auto fitresultz = hist_cos_z->Fit("fit");
	TF1 *fitz = hist_cos_z->GetFunction("fit");
	//Double_t z_p1 = fitz->GetParameter(0);
	//Double_t z_e1 = fitz->GetParError(0);

	hist_cos_x->SetMaximum(max_hist);
	hist_cos_x->SetMinimum(min_hist);
	hist_cos_x->Draw("hist");
	hist_cos_x->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_x->GetXaxis()->SetLabelSize(font_size);
	hist_cos_x->GetYaxis()->SetLabelSize(font_size);
	hist_cos_x->GetXaxis()->SetTitleSize(font_size);
	hist_cos_x->GetYaxis()->SetTitleSize(font_size);
	hist_cos_x->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_x->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_x->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_x_rec.pdf");

	hist_cos_y->SetMaximum(max_hist);
	hist_cos_y->SetMinimum(min_hist);
	hist_cos_y->Draw("hist");
	hist_cos_y->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_y->GetXaxis()->SetLabelSize(font_size);
	hist_cos_y->GetYaxis()->SetLabelSize(font_size);
	hist_cos_y->GetXaxis()->SetTitleSize(font_size);
	hist_cos_y->GetYaxis()->SetTitleSize(font_size);
	hist_cos_y->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_y->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_y->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_y_rec.pdf");

	hist_cos_z->SetMaximum(max_hist);
	hist_cos_z->SetMinimum(min_hist);
	hist_cos_z->Draw("hist");
	hist_cos_z->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_z->GetXaxis()->SetLabelSize(font_size);
	hist_cos_z->GetYaxis()->SetLabelSize(font_size);
	hist_cos_z->GetXaxis()->SetTitleSize(font_size);
	hist_cos_z->GetYaxis()->SetTitleSize(font_size);
	hist_cos_z->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_z->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_z->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_z_rec.pdf");

	gStyle->SetOptStat(0);
	TF1 *func2 = new TF1("lorentzianPeak", lorentzianPeak, min_mass, max_mass, 3);
	func2->SetParNames ("Norm.Const", "m_{0}", "#Gamma");
	func2->SetParameters(hist_p_pi_mass->GetEntries() / hist_p_pi_mass->GetNbinsX(), 1.11566, 2.59380e-03);
	func2->SetNpx(2000);
	func2->SetLineWidth(1);
	func2->SetLineColor(kBlue);

	auto fitresult_v0Mass = hist_p_pi_mass->Fit("lorentzianPeak");
	TF1 *fit_v0Mass = hist_p_pi_mass->GetFunction("lorentzianPeak");
	hist_p_pi_mass->SetLineStyle(2);
	hist_p_pi_mass->SetLineWidth(2);
	hist_p_pi_mass->Draw("P");
	hist_p_pi_mass->SetMarkerSize(3);
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
	hist_p_pi_mass->GetYaxis()->SetNdivisions(508);
	c0->SaveAs("hist_p_pi_mass_rec.pdf");

	func->SetParNames("P_{x}:");
	hist_cos_x_cut->Scale(hist_cos_x_cut->GetNbinsX() / (2 * hist_cos_x_cut->GetEntries()));
	auto fitresultx_cut = hist_cos_x_cut->Fit("fit");
	TF1 *fitx_cut = hist_cos_x_cut->GetFunction("fit");
	//Double_t x_p1 = fitx->GetParameter(0);
	//Double_t x_e1 = fitx->GetParError(0);
	func->SetParNames("P_{y}:");
	hist_cos_y_cut->Scale(hist_cos_y_cut->GetNbinsX() / (2 * hist_cos_y_cut->GetEntries()));
	auto fitresulty_cut = hist_cos_y_cut->Fit("fit");
	TF1 *fity_cut = hist_cos_y_cut->GetFunction("fit");
	//Double_t y_p1 = fity->GetParameter(0);
	//Double_t y_e1 = fity->GetParError(0);
	func->SetParNames("P_{z}:");
	hist_cos_z_cut->Scale(hist_cos_z_cut->GetNbinsX() / (2 * hist_cos_z_cut->GetEntries()));
	auto fitresultz_cut = hist_cos_z_cut->Fit("fit");
	TF1 *fitz_cut = hist_cos_z_cut->GetFunction("fit");
	//Double_t z_p1 = fitz->GetParameter(0);
	//Double_t z_e1 = fitz->GetParError(0);

	hist_cos_x_cut->SetMaximum(max_hist);
	hist_cos_x_cut->SetMinimum(min_hist);
	hist_cos_x_cut->Draw("hist");
	hist_cos_x_cut->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_x_cut->GetXaxis()->SetLabelSize(font_size);
	hist_cos_x_cut->GetYaxis()->SetLabelSize(font_size);
	hist_cos_x_cut->GetXaxis()->SetTitleSize(font_size);
	hist_cos_x_cut->GetYaxis()->SetTitleSize(font_size);
	hist_cos_x_cut->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_x_cut->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_x_cut->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_x_cut_rec.pdf");

	hist_cos_y_cut->SetMaximum(max_hist);
	hist_cos_y_cut->SetMinimum(min_hist);
	hist_cos_y_cut->Draw("hist");
	hist_cos_y_cut->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_y_cut->GetXaxis()->SetLabelSize(font_size);
	hist_cos_y_cut->GetYaxis()->SetLabelSize(font_size);
	hist_cos_y_cut->GetXaxis()->SetTitleSize(font_size);
	hist_cos_y_cut->GetYaxis()->SetTitleSize(font_size);
	hist_cos_y_cut->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_y_cut->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_y_cut->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_y_cut_rec.pdf");

	hist_cos_z_cut->SetMaximum(max_hist);
	hist_cos_z_cut->SetMinimum(min_hist);
	hist_cos_z_cut->Draw("hist");
	hist_cos_z_cut->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_z_cut->GetXaxis()->SetLabelSize(font_size);
	hist_cos_z_cut->GetYaxis()->SetLabelSize(font_size);
	hist_cos_z_cut->GetXaxis()->SetTitleSize(font_size);
	hist_cos_z_cut->GetYaxis()->SetTitleSize(font_size);
	hist_cos_z_cut->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_z_cut->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_z_cut->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_z_cut_rec.pdf");

	gStyle->SetOptStat(0);
	TF1 *func_cut = new TF1("lorentzianPeak_cut", lorentzianPeak, min_mass, max_mass, 3);
	func_cut->SetParNames ("Norm.Const", "m_{0}", "#Gamma");
	func_cut->SetParameters(hist_p_pi_mass_cut->GetEntries() / hist_p_pi_mass_cut->GetNbinsX(), 1.11566, 2.59380e-03);
	func_cut->SetNpx(2000);
	func_cut->SetLineWidth(1);
	func_cut->SetLineColor(kBlue);

	auto fitresult_v0Mass_cut = hist_p_pi_mass_cut->Fit("lorentzianPeak_cut");
	TF1 *fit_v0Mass_cut = hist_p_pi_mass_cut->GetFunction("lorentzianPeak_cut");
	hist_p_pi_mass_cut->SetLineStyle(2);
	hist_p_pi_mass_cut->SetLineWidth(2);
	hist_p_pi_mass_cut->Draw("P");
	hist_p_pi_mass_cut->SetMarkerSize(3);
	hist_p_pi_mass_cut->GetListOfFunctions()->FindObject("lorentzianPeak_cut")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_p_pi_mass_cut->GetXaxis()->SetLabelSize(font_size);
	hist_p_pi_mass_cut->GetYaxis()->SetLabelSize(font_size);
	hist_p_pi_mass_cut->GetXaxis()->SetTitleSize(font_size);
	hist_p_pi_mass_cut->GetYaxis()->SetTitleSize(font_size);
	hist_p_pi_mass_cut->GetYaxis()->SetTitleOffset(title_offset);
	hist_p_pi_mass_cut->GetYaxis()->SetNdivisions(508);
	c0->SaveAs("hist_p_pi_mass_cut_rec.pdf");

	func->SetParNames("P_{x}:");
	hist_cos_x_cut1->Scale(hist_cos_x_cut1->GetNbinsX() / (2 * hist_cos_x_cut1->GetEntries()));
	auto fitresultx_cut1 = hist_cos_x_cut1->Fit("fit");
	TF1 *fitx_cut1 = hist_cos_x_cut1->GetFunction("fit");
	//Double_t x_p1 = fitx->GetParameter(0);
	//Double_t x_e1 = fitx->GetParError(0);
	func->SetParNames("P_{y}:");
	hist_cos_y_cut1->Scale(hist_cos_y_cut1->GetNbinsX() / (2 * hist_cos_y_cut1->GetEntries()));
	auto fitresulty_cut1 = hist_cos_y_cut1->Fit("fit");
	TF1 *fity_cut1 = hist_cos_y_cut1->GetFunction("fit");
	//Double_t y_p1 = fity->GetParameter(0);
	//Double_t y_e1 = fity->GetParError(0);
	func->SetParNames("P_{z}:");
	hist_cos_z_cut1->Scale(hist_cos_z_cut1->GetNbinsX() / (2 * hist_cos_z_cut1->GetEntries()));
	auto fitresultz_cut1 = hist_cos_z_cut1->Fit("fit");
	TF1 *fitz_cut1 = hist_cos_z_cut1->GetFunction("fit");
	//Double_t z_p1 = fitz->GetParameter(0);
	//Double_t z_e1 = fitz->GetParError(0);

	hist_cos_x_cut1->SetMaximum(max_hist);
	hist_cos_x_cut1->SetMinimum(min_hist);
	hist_cos_x_cut1->Draw("hist");
	hist_cos_x_cut1->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_x_cut1->GetXaxis()->SetLabelSize(font_size);
	hist_cos_x_cut1->GetYaxis()->SetLabelSize(font_size);
	hist_cos_x_cut1->GetXaxis()->SetTitleSize(font_size);
	hist_cos_x_cut1->GetYaxis()->SetTitleSize(font_size);
	hist_cos_x_cut1->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_x_cut1->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_x_cut1->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_x_cut1_rec.pdf");

	hist_cos_y_cut1->SetMaximum(max_hist);
	hist_cos_y_cut1->SetMinimum(min_hist);
	hist_cos_y_cut1->Draw("hist");
	hist_cos_y_cut1->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_y_cut1->GetXaxis()->SetLabelSize(font_size);
	hist_cos_y_cut1->GetYaxis()->SetLabelSize(font_size);
	hist_cos_y_cut1->GetXaxis()->SetTitleSize(font_size);
	hist_cos_y_cut1->GetYaxis()->SetTitleSize(font_size);
	hist_cos_y_cut1->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_y_cut1->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_y_cut1->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_y_cut1_rec.pdf");

	hist_cos_z_cut1->SetMaximum(max_hist);
	hist_cos_z_cut1->SetMinimum(min_hist);
	hist_cos_z_cut1->Draw("hist");
	hist_cos_z_cut1->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_z_cut1->GetXaxis()->SetLabelSize(font_size);
	hist_cos_z_cut1->GetYaxis()->SetLabelSize(font_size);
	hist_cos_z_cut1->GetXaxis()->SetTitleSize(font_size);
	hist_cos_z_cut1->GetYaxis()->SetTitleSize(font_size);
	hist_cos_z_cut1->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_z_cut1->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_z_cut1->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_z_cut1_rec.pdf");

	gStyle->SetOptStat(0);
	TF1 *func_cut1 = new TF1("lorentzianPeak_cut1", lorentzianPeak, min_mass, max_mass, 3);
	func_cut1->SetParNames ("Norm.Const", "m_{0}", "#Gamma");
	func_cut1->SetParameters(hist_p_pi_mass_cut1->GetEntries() / hist_p_pi_mass_cut1->GetNbinsX(), 1.11566, 2.59380e-03);
	func_cut1->SetNpx(2000);
	func_cut1->SetLineWidth(1);
	func_cut1->SetLineColor(kBlue);

	auto fitresult_v0Mass_cut1 = hist_p_pi_mass_cut1->Fit("lorentzianPeak_cut1");
	TF1 *fit_v0Mass_cut1 = hist_p_pi_mass_cut1->GetFunction("lorentzianPeak_cut1");
	hist_p_pi_mass_cut1->SetLineStyle(2);
	hist_p_pi_mass_cut1->SetLineWidth(2);
	hist_p_pi_mass_cut1->Draw("P");
	hist_p_pi_mass_cut1->SetMarkerSize(3);
	hist_p_pi_mass_cut1->GetListOfFunctions()->FindObject("lorentzianPeak_cut1")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_p_pi_mass_cut1->GetXaxis()->SetLabelSize(font_size);
	hist_p_pi_mass_cut1->GetYaxis()->SetLabelSize(font_size);
	hist_p_pi_mass_cut1->GetXaxis()->SetTitleSize(font_size);
	hist_p_pi_mass_cut1->GetYaxis()->SetTitleSize(font_size);
	hist_p_pi_mass_cut1->GetYaxis()->SetTitleOffset(title_offset);
	hist_p_pi_mass_cut1->GetYaxis()->SetNdivisions(508);
	c0->SaveAs("hist_p_pi_mass_cut1_rec.pdf");

	func->SetParNames("P_{x}:");
	hist_cos_x_cut2->Scale(hist_cos_x_cut2->GetNbinsX() / (2 * hist_cos_x_cut2->GetEntries()));
	auto fitresultx_cut2 = hist_cos_x_cut2->Fit("fit");
	TF1 *fitx_cut2 = hist_cos_x_cut2->GetFunction("fit");
	//Double_t x_p1 = fitx->GetParameter(0);
	//Double_t x_e1 = fitx->GetParError(0);
	func->SetParNames("P_{y}:");
	hist_cos_y_cut2->Scale(hist_cos_y_cut2->GetNbinsX() / (2 * hist_cos_y_cut2->GetEntries()));
	auto fitresulty_cut2 = hist_cos_y_cut2->Fit("fit");
	TF1 *fity_cut2 = hist_cos_y_cut2->GetFunction("fit");
	//Double_t y_p1 = fity->GetParameter(0);
	//Double_t y_e1 = fity->GetParError(0);
	func->SetParNames("P_{z}:");
	hist_cos_z_cut2->Scale(hist_cos_z_cut2->GetNbinsX() / (2 * hist_cos_z_cut2->GetEntries()));
	auto fitresultz_cut2 = hist_cos_z_cut2->Fit("fit");
	TF1 *fitz_cut2 = hist_cos_z_cut2->GetFunction("fit");
	//Double_t z_p1 = fitz->GetParameter(0);
	//Double_t z_e1 = fitz->GetParError(0);

	hist_cos_x_cut2->SetMaximum(max_hist);
	hist_cos_x_cut2->SetMinimum(min_hist);
	hist_cos_x_cut2->Draw("hist");
	hist_cos_x_cut2->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_x_cut2->GetXaxis()->SetLabelSize(font_size);
	hist_cos_x_cut2->GetYaxis()->SetLabelSize(font_size);
	hist_cos_x_cut2->GetXaxis()->SetTitleSize(font_size);
	hist_cos_x_cut2->GetYaxis()->SetTitleSize(font_size);
	hist_cos_x_cut2->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_x_cut2->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_x_cut2->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_x_cut2_rec.pdf");

	hist_cos_y_cut2->SetMaximum(max_hist);
	hist_cos_y_cut2->SetMinimum(min_hist);
	hist_cos_y_cut2->Draw("hist");
	hist_cos_y_cut2->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_y_cut2->GetXaxis()->SetLabelSize(font_size);
	hist_cos_y_cut2->GetYaxis()->SetLabelSize(font_size);
	hist_cos_y_cut2->GetXaxis()->SetTitleSize(font_size);
	hist_cos_y_cut2->GetYaxis()->SetTitleSize(font_size);
	hist_cos_y_cut2->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_y_cut2->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_y_cut2->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_y_cut2_rec.pdf");

	hist_cos_z_cut2->SetMaximum(max_hist);
	hist_cos_z_cut2->SetMinimum(min_hist);
	hist_cos_z_cut2->Draw("hist");
	hist_cos_z_cut2->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_z_cut2->GetXaxis()->SetLabelSize(font_size);
	hist_cos_z_cut2->GetYaxis()->SetLabelSize(font_size);
	hist_cos_z_cut2->GetXaxis()->SetTitleSize(font_size);
	hist_cos_z_cut2->GetYaxis()->SetTitleSize(font_size);
	hist_cos_z_cut2->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_z_cut2->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_z_cut2->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_z_cut2_rec.pdf");

	gStyle->SetOptStat(0);
	TF1 *func_cut2 = new TF1("lorentzianPeak_cut2", lorentzianPeak, min_mass, max_mass, 3);
	func_cut2->SetParNames ("Norm.Const", "m_{0}", "#Gamma");
	func_cut2->SetParameters(hist_p_pi_mass_cut2->GetEntries() / hist_p_pi_mass_cut2->GetNbinsX(), 1.11566, 2.59380e-03);
	func_cut2->SetNpx(2000);
	func_cut2->SetLineWidth(1);
	func_cut2->SetLineColor(kBlue);

	auto fitresult_v0Mass_cut2 = hist_p_pi_mass_cut2->Fit("lorentzianPeak_cut2");
	TF1 *fit_v0Mass_cut2 = hist_p_pi_mass_cut2->GetFunction("lorentzianPeak_cut2");
	hist_p_pi_mass_cut2->SetLineStyle(2);
	hist_p_pi_mass_cut2->SetLineWidth(2);
	hist_p_pi_mass_cut2->Draw("P");
	hist_p_pi_mass_cut2->SetMarkerSize(3);
	hist_p_pi_mass_cut2->GetListOfFunctions()->FindObject("lorentzianPeak_cut2")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_p_pi_mass_cut2->GetXaxis()->SetLabelSize(font_size);
	hist_p_pi_mass_cut2->GetYaxis()->SetLabelSize(font_size);
	hist_p_pi_mass_cut2->GetXaxis()->SetTitleSize(font_size);
	hist_p_pi_mass_cut2->GetYaxis()->SetTitleSize(font_size);
	hist_p_pi_mass_cut2->GetYaxis()->SetTitleOffset(title_offset);
	hist_p_pi_mass_cut2->GetYaxis()->SetNdivisions(508);
	c0->SaveAs("hist_p_pi_mass_cut2_rec.pdf");

	func->SetParNames("P_{x}:");
	hist_cos_x_cut3->Scale(hist_cos_x_cut3->GetNbinsX() / (2 * hist_cos_x_cut3->GetEntries()));
	auto fitresultx_cut3 = hist_cos_x_cut3->Fit("fit");
	TF1 *fitx_cut3 = hist_cos_x_cut3->GetFunction("fit");
	//Double_t x_p1 = fitx->GetParameter(0);
	//Double_t x_e1 = fitx->GetParError(0);
	func->SetParNames("P_{y}:");
	hist_cos_y_cut3->Scale(hist_cos_y_cut3->GetNbinsX() / (2 * hist_cos_y_cut3->GetEntries()));
	auto fitresulty_cut3 = hist_cos_y_cut3->Fit("fit");
	TF1 *fity_cut3 = hist_cos_y_cut3->GetFunction("fit");
	//Double_t y_p1 = fity->GetParameter(0);
	//Double_t y_e1 = fity->GetParError(0);
	func->SetParNames("P_{z}:");
	hist_cos_z_cut3->Scale(hist_cos_z_cut3->GetNbinsX() / (2 * hist_cos_z_cut3->GetEntries()));
	auto fitresultz_cut3 = hist_cos_z_cut3->Fit("fit");
	TF1 *fitz_cut3 = hist_cos_z_cut3->GetFunction("fit");
	//Double_t z_p1 = fitz->GetParameter(0);
	//Double_t z_e1 = fitz->GetParError(0);

	hist_cos_x_cut3->SetMaximum(max_hist);
	hist_cos_x_cut3->SetMinimum(min_hist);
	hist_cos_x_cut3->Draw("hist");
	hist_cos_x_cut3->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_x_cut3->GetXaxis()->SetLabelSize(font_size);
	hist_cos_x_cut3->GetYaxis()->SetLabelSize(font_size);
	hist_cos_x_cut3->GetXaxis()->SetTitleSize(font_size);
	hist_cos_x_cut3->GetYaxis()->SetTitleSize(font_size);
	hist_cos_x_cut3->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_x_cut3->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_x_cut3->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_x_cut3_rec.pdf");

	hist_cos_y_cut3->SetMaximum(max_hist);
	hist_cos_y_cut3->SetMinimum(min_hist);
	hist_cos_y_cut3->Draw("hist");
	hist_cos_y_cut3->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_y_cut3->GetXaxis()->SetLabelSize(font_size);
	hist_cos_y_cut3->GetYaxis()->SetLabelSize(font_size);
	hist_cos_y_cut3->GetXaxis()->SetTitleSize(font_size);
	hist_cos_y_cut3->GetYaxis()->SetTitleSize(font_size);
	hist_cos_y_cut3->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_y_cut3->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_y_cut3->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_y_cut3_rec.pdf");

	hist_cos_z_cut3->SetMaximum(max_hist);
	hist_cos_z_cut3->SetMinimum(min_hist);
	hist_cos_z_cut3->Draw("hist");
	hist_cos_z_cut3->GetListOfFunctions()->FindObject("fit")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_cos_z_cut3->GetXaxis()->SetLabelSize(font_size);
	hist_cos_z_cut3->GetYaxis()->SetLabelSize(font_size);
	hist_cos_z_cut3->GetXaxis()->SetTitleSize(font_size);
	hist_cos_z_cut3->GetYaxis()->SetTitleSize(font_size);
	hist_cos_z_cut3->GetYaxis()->SetTitleOffset(title_offset);
	hist_cos_z_cut3->GetXaxis()->SetNdivisions(n_divisions);
	hist_cos_z_cut3->GetYaxis()->SetNdivisions(n_divisions);
	c0->SaveAs("hist_cos_z_cut3_rec.pdf");

	gStyle->SetOptStat(0);
	TF1 *func_cut3 = new TF1("lorentzianPeak_cut3", lorentzianPeak, min_mass, max_mass, 3);
	func_cut3->SetParNames ("Norm.Const", "m_{0}", "#Gamma");
	func_cut3->SetParameters(hist_p_pi_mass_cut3->GetEntries() / hist_p_pi_mass_cut3->GetNbinsX(), 1.11566, 2.59380e-03);
	func_cut3->SetNpx(2000);
	func_cut3->SetLineWidth(1);
	func_cut3->SetLineColor(kBlue);

	auto fitresult_v0Mass_cut3 = hist_p_pi_mass_cut3->Fit("lorentzianPeak_cut3");
	TF1 *fit_v0Mass_cut3 = hist_p_pi_mass_cut3->GetFunction("lorentzianPeak_cut3");
	hist_p_pi_mass_cut3->SetLineStyle(2);
	hist_p_pi_mass_cut3->SetLineWidth(2);
	hist_p_pi_mass_cut3->Draw("P");
	hist_p_pi_mass_cut3->SetMarkerSize(3);
	hist_p_pi_mass_cut3->GetListOfFunctions()->FindObject("lorentzianPeak_cut3")->Draw("same");
	gPad->SetLeftMargin(big_margin);
	gPad->SetBottomMargin(big_margin);
	gPad->SetRightMargin(small_margin);
	gPad->SetTopMargin(medium_margin);
	hist_p_pi_mass_cut3->GetXaxis()->SetLabelSize(font_size);
	hist_p_pi_mass_cut3->GetYaxis()->SetLabelSize(font_size);
	hist_p_pi_mass_cut3->GetXaxis()->SetTitleSize(font_size);
	hist_p_pi_mass_cut3->GetYaxis()->SetTitleSize(font_size);
	hist_p_pi_mass_cut3->GetYaxis()->SetTitleOffset(title_offset);
	hist_p_pi_mass_cut3->GetYaxis()->SetNdivisions(508);
	c0->SaveAs("hist_p_pi_mass_cut3_rec.pdf");

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
	c0->SaveAs("hist_y_rec.pdf");

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
	c0->SaveAs("hist_pt_rec.pdf");

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
	c0->SaveAs("hist_phi_rec.pdf");

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
	c0->SaveAs("hist_y_z_rec.pdf");

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
	c0->SaveAs("hist_y_cut_rec.pdf");

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
	c0->SaveAs("hist_pt_cut_rec.pdf");

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
	c0->SaveAs("hist_phi_cut_rec.pdf");

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
	c0->SaveAs("hist_y_cut1_rec.pdf");

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
	c0->SaveAs("hist_ypt_cut1_rec.pdf");

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
	c0->SaveAs("hist_phi_cut1_rec.pdf");

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
	c0->SaveAs("hist_delta_z_cut1_rec.pdf");

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
	c0->SaveAs("hist_y_z_cut1_rec.pdf");

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
	c0->SaveAs("hist_y_cut2_rec.pdf");

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
	c0->SaveAs("hist_pt_cut2_rec.pdf");

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
	c0->SaveAs("hist_ypt_cut2_rec.pdf");

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
	c0->SaveAs("hist_phi_cut2_rec.pdf");

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
	c0->SaveAs("hist_delta_z_cut2_rec.pdf");

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
	c0->SaveAs("hist_y_z_cut2_rec.pdf");

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
	c0->SaveAs("hist_y_cut3_rec.pdf");

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
	c0->SaveAs("hist_pt_cut3_rec.pdf");

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
	c0->SaveAs("hist_ypt_cut3_rec.pdf");

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
	c0->SaveAs("hist_phi_cut3_rec.pdf");

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
	c0->SaveAs("hist_delta_z_cut3_rec.pdf");

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
	c0->SaveAs("hist_y_z_cut3_rec.pdf");

}
