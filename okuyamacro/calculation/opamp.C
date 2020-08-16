double Impedance(double *x, double *par){
	double omega = 2.*TMath::Pi()*x[0];
	double Resistance = par[0];
	double Inductance = par[1];
	double Capacitance = par[2];
	double temp_func = omega*Inductance-1/(omega*Capacitance);
	return TMath::Sqrt(TMath::Power(Resistance,2.)+TMath::Power(temp_func,2.));
}
//Inverse Differential Amplification
double Gain(double *x, double *par){
	double omega = 2.*TMath::Pi()*x[0];
	double R1 = par[0];
	double R2 = par[1];
	double C  = par[2];
	double gain = R2/sqrt((1/(omega*omega*C*C))+R1*R1);
	return 20.*TMath::Log10(gain);
}
double theta(double *x, double *par){//Phase shift
	double omega = 2.*TMath::Pi()*x[0];
	double R1 = par[0];
	double R2 = par[1];
	double C  = par[2];
	double theta = atan(-1./-(omega*R1*C))-TMath::Pi();
	return theta*180./TMath::Pi();
}

void opamp(){
gROOT->Reset();
gROOT->SetStyle("Plain");
gStyle->SetOptLogx(1);
gStyle->SetOptLogy(0);
gStyle->SetPadGridX(1);
gStyle->SetPadGridY(1);

TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
c1->SetLeftMargin(0.15);
c1->SetBottomMargin(0.15);
c1->SetTopMargin(0.15);
c1->SetRightMargin(0.15);
c2->SetLeftMargin(0.15);
c2->SetBottomMargin(0.15);
c2->SetTopMargin(0.15);
c2->SetRightMargin(0.15);

double R1 = 20;
double R2 = 1000;
double L = 2.*TMath::Power(10.,-9.);

double C1 = TMath::Power(10.,-8.);
double C2 = TMath::Power(10.,-9.);
double C3 = TMath::Power(10.,-10.);//100pF

c1->cd();
TH1F* h1 = c1 -> DrawFrame(1.,-160,TMath::Power(10.,11.),80.);
//h1->GetXaxis()->SetTitleFont(12);
h1->GetXaxis()->SetTitleSize(0.05);
h1->GetXaxis()->SetTitleOffset(1.1);
h1->GetXaxis()->SetTitle("#omega");
h1->GetXaxis()->SetLabelOffset(0.02);
//h1->GetYaxis()->SetTitleFont(22);
h1->GetYaxis()->SetTitleOffset(1.1);
h1->GetYaxis()->SetTitleSize(0.05);
h1->GetYaxis()->SetTitle("|G| [dB]");


cout<<"aa"<<endl;
TF1* func = new TF1("func",Gain, 1, TMath::Power(10.,11.),3);
func->SetNpx(600);
func->SetParameter(0,R1);
func->SetParameter(1,R2);
func->SetParameter(2,C3);
func->SetLineColor(kGreen);
func->SetLineWidth(4);
func->Draw("same");


 func->SetMarkerStyle(8);
 func->SetMarkerSize(1);
 func->SetMarkerColor(kGreen);

c2->cd();
TH1F* h2 = c2 -> DrawFrame(1.,-450,TMath::Power(10.,11.),-90.);
h2->GetXaxis()->SetTitleSize(0.05);
h2->GetXaxis()->SetTitleOffset(1.1);
h2->GetXaxis()->SetTitle("#omega");
h2->GetXaxis()->SetLabelOffset(0.02);
h2->GetYaxis()->SetTitleOffset(1.1);
h2->GetYaxis()->SetTitleSize(0.05);
h2->GetYaxis()->SetTitle("#theta [deg]");

TF1* func2 = new TF1("func2",theta, 1, TMath::Power(10.,11.),3);
func2->SetNpx(600);
func2->SetParameter(0,R1);
func2->SetParameter(1,R2);
func2->SetParameter(2,C3);
func2->SetLineColor(kGreen);
func2->SetLineWidth(4);
func2->SetLineStyle(5);
func2->Draw("same");


 func2->SetMarkerStyle(8);
 func2->SetMarkerSize(1);
 func2->SetMarkerColor(kGreen);

// TLegend *tl = new TLegend(0.65,0.75,0.85,0.85);
// tl->AddEntry(func,"10[nF]","p");
// tl->AddEntry(func2,"1[nF]","p");
// tl->AddEntry(func3,"0.1[nF]","p");;
// tl->Draw();

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
