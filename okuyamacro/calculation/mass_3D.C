//double sqrt_func(double *x, double *par){
//	double x0 = x[0];
//	double mass = par[0];
//
//	return TMath::Sqrt(TMath::Power(x0,2.)+TMath::Power(mass,2.))/x0;
//}
//
double sqrt(double *x){
double x0=x[0];
double f = TMath::Sqrt(x0);
return f;
}

double energy_part(double *x){
double pk=x[0];
double pe=x[1];
double M=940;
double f = TMath::Power((1851-TMath::Sqrt(TMath::Power(pe,2.)+TMath::Power(0.511,2.))+M-TMath::Sqrt(TMath::Power(pk,2.)+TMath::Power(494.,2.))),2.);
return f;
}

double momentum_part(double *x){
double pk=x[0];
double pe=x[1];
return TMath::Power(1851-pe-pk,2.);
}

double mass_func(double *x){
double pk=x[0];
double pe=x[1];
double M=940;
double f = TMath::Sqrt(TMath::Power((1851-TMath::Sqrt(TMath::Power(pe,2.)+TMath::Power(0.511,2.))+M-TMath::Sqrt(TMath::Power(pk,2.)+TMath::Power(494.,2.))),2.)-TMath::Power(1851-pe-pk,2.));
return f;
}

void mass_3D(){
gROOT->Reset();
gROOT->SetStyle("Plain");

//TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
//c1->SetLeftMargin(0.15);
//c1->SetBottomMargin(0.15);
//c1->SetTopMargin(0.15);
//c1->SetRightMargin(0.15);
//c1->Draw();


//TH2F* h1 = c1 -> DrawFrame(500.,1000.,1000.,2000.);
//h1->GetXaxis()->SetTitleFont(12);
//h1->GetXaxis()->SetTitleSize(0.05);
//h1->GetXaxis()->SetTitleOffset(1.05);
//h1->GetXaxis()->SetTitle("p [MeV]");
//h1->GetXaxis()->SetLabelOffset(0.02);
//h1->GetYaxis()->SetTitleFont(22);
//h1->GetYaxis()->SetTitleOffset(1.1);
//h1->GetYaxis()->SetTitleSize(0.05);
//h1->GetYaxis()->SetTitle("1/#beta");

double x,y;
double M = 940;
cout<<"aa"<<endl;
//TF1 *func = new TF1("func","sqrt(x)",1,5);
TF2* func = new TF2("func","TMath::Sqrt(TMath::Power((1851-TMath::Sqrt(TMath::Power(y,2.)+TMath::Power(0.511,2.))+940-TMath::Sqrt(TMath::Power(x,2.)+TMath::Power(494.,2.))),2.)-TMath::Power(1851-y-x,2.))",1100,1300,250,400);
//TF2* func = new TF2("func",")", 1,1000,1,2000);
func->Draw("surf");
TF2* f_Lambda = new TF2("f_Lambda","1115-0.00001*(x-y)",1100,1300,250,400);
//f_Lambda->Draw("same,surf");
//TF3* func2 = new TF3("func2",sqrt_func, 200., 1400.,1);
//func2->SetNpx(600);
//func2->SetParameter(0,139);
//func2->SetLineColor(kBlue);
//func2->SetLineWidth(4);
//func2->Draw("same");

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
