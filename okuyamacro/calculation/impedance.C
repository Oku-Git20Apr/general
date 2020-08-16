double Impedance(double *x, double *par){
	double omega = 2.*TMath::Pi()*x[0];
	double Resistance = par[0];
	double Inductance = par[1];
	double Capacitance = par[2];
	double temp_func = omega*Inductance-1/(omega*Capacitance);
	return TMath::Sqrt(TMath::Power(Resistance,2.)+TMath::Power(temp_func,2.));
}

void impedance(){
gROOT->Reset();
gROOT->SetStyle("Plain");
gStyle->SetOptLogx(1);
gStyle->SetOptLogy(1);

TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
c1->SetLeftMargin(0.15);
c1->SetBottomMargin(0.15);
c1->SetTopMargin(0.15);
c1->SetRightMargin(0.15);
c1->Draw();

double R = 0.01;
double L = 2.*TMath::Power(10.,-9.);

double C1 = TMath::Power(10.,-8.);
double C2 = TMath::Power(10.,-9.);
double C3 = TMath::Power(10.,-10.);


TH1F* h1 = c1 -> DrawFrame(10000000.,0.01,1000000000.,1000.);
//h1->GetXaxis()->SetTitleFont(12);
h1->GetXaxis()->SetTitleSize(0.05);
h1->GetXaxis()->SetTitleOffset(1.1);
h1->GetXaxis()->SetTitle("frequency [Hz]");
h1->GetXaxis()->SetLabelOffset(0.02);
//h1->GetYaxis()->SetTitleFont(22);
h1->GetYaxis()->SetTitleOffset(1.1);
h1->GetYaxis()->SetTitleSize(0.05);
h1->GetYaxis()->SetTitle("Z [#Omega]");


cout<<"aa"<<endl;
TF1* func = new TF1("func",Impedance, 10000000, 1000000000.,3);
func->SetNpx(600);
func->SetParameter(0,R);
func->SetParameter(1,L);
func->SetParameter(2,C1);
func->SetLineColor(kGreen);
func->SetLineWidth(4);
func->Draw("same");

TF1* func2 = new TF1("func2",Impedance, 10000000, 1000000000.,3);
func2->SetNpx(600);
func2->SetParameter(0,R);
func2->SetParameter(1,L);
func2->SetParameter(2,C2);
func2->SetLineColor(kBlue);
func2->SetLineWidth(4);
func2->Draw("same");

TF1* func3 = new TF1("func3",Impedance, 10000000, 1000000000.,3);
func3->SetNpx(600);
func3->SetParameter(0,R);
func3->SetParameter(1,L);
func3->SetParameter(2,C3);
func3->SetLineColor(kRed);
func3->SetLineWidth(4);
func3->Draw("same");

 func->SetMarkerStyle(8);
 func->SetMarkerSize(1);
 func->SetMarkerColor(kGreen);
 func2->SetMarkerStyle(8);
 func2->SetMarkerSize(1); 
 func2->SetMarkerColor(kBlue);
 func3->SetMarkerStyle(8);
 func3->SetMarkerSize(1);
 func3->SetMarkerColor(kRed);



 TLegend *tl = new TLegend(0.65,0.75,0.85,0.85);
 tl->AddEntry(func,"10[nF]","p");
 tl->AddEntry(func2,"1[nF]","p");
 tl->AddEntry(func3,"0.1[nF]","p");;
 tl->Draw();

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
