double sqrt_func(double *x, double *par){
	double x0 = x[0];
	double mass = par[0];

	return x0/TMath::Sqrt(x0*x0+mass*mass);
}

void cherenkov2(){
gROOT->Reset();
gROOT->SetStyle("Plain");

TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
c1->SetLeftMargin(0.15);
c1->SetBottomMargin(0.15);
c1->SetTopMargin(0.15);
c1->SetRightMargin(0.15);
c1->Draw();


TH1F* h1 = c1 -> DrawFrame(700.,0.9,2300.,1.02);
//h1->GetXaxis()->SetTitleFont(12);
h1->GetXaxis()->SetTitleSize(0.05);
h1->GetXaxis()->SetTitleOffset(1.05);
h1->GetXaxis()->SetTitle("p [MeV]");
h1->GetXaxis()->SetLabelOffset(0.02);
//h1->GetYaxis()->SetTitleFont(22);
h1->GetYaxis()->SetTitleOffset(1.1);
h1->GetYaxis()->SetTitleSize(0.05);
h1->GetYaxis()->SetTitle("#beta");


cout<<"aa"<<endl;
TF1* func = new TF1("func",sqrt_func, 700., 2300.,1);
func->SetNpx(600);
func->SetParameter(0,139);
func->SetLineColor(kOrange);
func->SetLineWidth(4);
func->Draw("same");

TF1* func2 = new TF1("func2",sqrt_func, 700., 2300.,1);
func2->SetNpx(600);
func2->SetParameter(0,494);
func2->SetLineColor(kGreen);
func2->SetLineWidth(4);
func2->Draw("same");

TF1* func3 = new TF1("func3",sqrt_func, 700., 2300.,1);
func3->SetNpx(600);
func3->SetParameter(0,979);
func3->SetLineColor(kRed);
func3->SetLineWidth(4);
func3->Draw("same");

TLine* n1 = new TLine(700.,1/1.015,2300.,1/1.015);
n1->SetLineColor(kAzure);
n1->Draw("same");
TLine* n2 = new TLine(700.,1/1.055,2300.,1/1.055);
n2->SetLineColor(kAzure);
n2->Draw("same");
TLine* l1 = new TLine(2000.,0.9,2000.,1.02);
//l1->SetLineColor(kRed);
l1->Draw("same");
TLine* l2 = new TLine(2200.,0.9,2200.,1.02);
//l2->SetLineColor(kRed);
l2->Draw("same");
TPave* p = new TPave(2000.,0.9,2200.,1.02);
p->SetFillColor(kGray);
p->SetFillStyle(3003);
//p->SetBorderSize(4);
p->Draw("same");


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
