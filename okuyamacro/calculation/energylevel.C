//K. Okuyama (Sep. 6, 2021)
void energylevel(){
gROOT->Reset();
gROOT->SetStyle("Plain");

TCanvas* c1 = new TCanvas("c1", "Energy Level", 800, 600);
c1->SetLeftMargin(0.15);
c1->SetBottomMargin(0.15);
c1->SetTopMargin(0.15);
c1->SetRightMargin(0.15);
c1->Draw();


TH1F* h1 = c1 -> DrawFrame(0.,-20.,9.,5.);
//h1->GetXaxis()->SetTitleFont(12);
h1->GetXaxis()->SetTitleSize(0.05);
h1->GetXaxis()->SetTitleOffset(1.05);
h1->GetXaxis()->SetTitle("");
h1->GetXaxis()->SetLabelOffset(0.02);
h1->GetXaxis()->SetNdivisions(0);
//h1->GetYaxis()->SetTitleFont(22);
h1->GetYaxis()->SetTitleOffset(1.1);
h1->GetYaxis()->SetTitleSize(0.05);
h1->GetYaxis()->SetTitle("Binding Energy [MeV]");


//Core States: Al27
double Al27_gs = -17.57;//MeV
double Al27_core[5] = {0.00, 0.84, 1.01, 2.21, 2.73};//MeV
TLine* lcore[5];
	for(int i=0;i<5;i++){
		lcore[i] = new TLine(3.5,Al27_core[i]+Al27_gs,5.5,Al27_core[i]+Al27_gs);
		lcore[i]->SetLineColor(kBlack);
		lcore[i]->SetLineWidth(2);
		lcore[i]->Draw("same");
	}
//Al28L: E01-011
double Al28L_exp[3] = {-17.57, -6.84, 1.82};//MeV
TLine* lhyp_exp[3];
	for(int i=0;i<3;i++){
		lhyp_exp[i] = new TLine(0.5,Al28L_exp[i],2.5,Al28L_exp[i]);
		lhyp_exp[i]->SetLineColor(kAzure);
		lhyp_exp[i]->SetLineWidth(4);
		lhyp_exp[i]->Draw("same");
	}
//Al28L: Shell model calc.
double Al28L_theo[5] = {-16.92, -8.60, -8.00, -0.29, 1.29};//MeV
TLine* lhyp_theo[5];
	for(int i=0;i<5;i++){
		lhyp_theo[i] = new TLine(6.5,Al28L_theo[i],8.5,Al28L_theo[i]);
		lhyp_theo[i]->SetLineColor(kRed);
		lhyp_theo[i]->SetLineWidth(4);
		lhyp_theo[i]->Draw("same");
	}
//TPave* p = new TPave(0.,-20.,9.,5.);
//p->SetFillColor(kGray);
//p->SetFillStyle(3003);
////p->SetBorderSize(4);
//p->Draw("same");


//TF1 *func = new TF1("func", "sqrt(x**2+0.511**2)/x", 200, 1400);
//func->Draw();
//TF1 *func2 = new TF1("func2", "sqrt(x**2+139**2)/x", 200, 1400);
//func2->Draw();
//TF1 *f1 = new TF1("f1", "[0]*pow(cos(x),[1])",0,90);

//f1->SetParameter(0,500);
//f1->SetParameter(1,2);

//f1->SetLineColor(2);
//f1->Draw("same");

//c1->Print("montecarlo_data0826_1kai.png");

}
