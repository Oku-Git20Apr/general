double sqrt_func(double *x, double *par){
	double x0 = x[0];
	double mass = par[0];

	return TMath::Sqrt(TMath::Power(x0,2.)+TMath::Power(mass,2.))/x0;
}

void cherenkov_threshold(){
gROOT->Reset();
gROOT->SetStyle("Plain");

TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
c1->SetLeftMargin(0.15);
c1->SetBottomMargin(0.15);
c1->SetTopMargin(0.15);
c1->SetRightMargin(0.15);
c1->Draw();


TH1F* h1 = c1 -> DrawFrame(200.,0.99,1400.,1.02);
//h1->GetXaxis()->SetTitleFont(12);
h1->GetXaxis()->SetTitleSize(0.05);
h1->GetXaxis()->SetTitleOffset(1.05);
h1->GetXaxis()->SetTitle("p [MeV]");
h1->GetXaxis()->SetLabelOffset(0.02);
//h1->GetYaxis()->SetTitleFont(22);
h1->GetYaxis()->SetTitleOffset(1.1);
h1->GetYaxis()->SetTitleSize(0.05);
h1->GetYaxis()->SetTitle("1/#beta");


cout<<"aa"<<endl;
TF1* func = new TF1("func",sqrt_func, 200., 1400.,1);
func->SetNpx(600);
func->SetParameter(0,0.511);
func->SetLineColor(kOrange);
func->SetLineWidth(4);
func->Draw("same");

TF1* func2 = new TF1("func2",sqrt_func, 200., 1400.,1);
func2->SetNpx(600);
func2->SetParameter(0,139);
func2->SetLineColor(kBlue);
func2->SetLineWidth(4);
func2->Draw("same");

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
