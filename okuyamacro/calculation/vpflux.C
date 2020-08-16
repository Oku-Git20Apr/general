double sqrt_func(double *x, double *par){
	double x0 = x[0];
	double mass = par[0];

	return TMath::Sqrt(TMath::Power(x0,2.)+TMath::Power(mass,2.))/x0;
}
double PI=4.*atan(1.);

double vpflux_doi(double *x, double *par){
	double Einc  = 4.3;//[GeV]
	//double Escat = 2.1;//[GeV]
	double Escat = par[0];//[GeV]
	double theta = x[0];	

	double Me=pow(511,-6.);//[GeV/c^2]
	double Mp=0.9382720;//[GeV/c^2]
	double Qsq=2*Einc*Escat*(1-cos(theta));
	double omega = Einc - Escat;
	double Eg = Einc - Escat - (Qsq/(2*Mp));
	double A=Me*Me*omega*omega/(4*Einc*Einc*Escat*Escat);
	double sinterm=sin(theta/2)*sin(theta/2);
	double a1=((Einc*Einc+Escat*Escat)/(2*Einc*Einc))/(A+sinterm);
	double a2=(Escat/Einc)*A/((A+sinterm)*(A+sinterm));
	//double a3=((Einc+Escat)*(Einc+Escat)/(4*Einc*Einc))/(omega*omega/(4*Einc*Escat)+sinterm);
	double a3=((Einc+Escat)*(Einc+Escat)/(4*Einc*Einc))/(A+sinterm);

	//double vpflux=Eg*(a1-a2-a3)/(137*4*PI*PI*omega*omega);
	double vpflux=(a1-a2-a3)/(137*4*PI*PI*omega);
	return vpflux;
}

double vpflux_kawama(double *x, double *par){
	double Einc  = 4.3;//[GeV]
	//double Escat = 2.1;//[GeV]
	double Escat = par[0];//[GeV]
	double theta = x[0];	

	double Me=pow(511,-6);//[GeV/c^2]
	double Mp=0.9382720;//[GeV/c^2]
	double Qsq=2*Einc*Escat*(1-cos(theta));
	double omega = Einc - Escat;
	double Eg = Einc - Escat - (Qsq/(2*Mp));
	double A=Me*Me*omega*omega/(4*Einc*Einc*Escat*Escat);
	double sinterm=sin(theta/2)*sin(theta/2);
	double a1=((Einc*Einc+Escat*Escat)/(2*Einc*Einc))/(A+sinterm);
	double a2=(Einc/Escat)*A/((A+sinterm)*(A+sinterm));
	double a3=((Einc+Escat)*(Einc+Escat)/(4*Einc*Einc))/(omega*omega/(4*Einc*Escat)+sinterm);

	double vpflux=Eg*(a1-a2-a3)/(137*4*PI*PI*omega*omega);
	return vpflux;
}
double vpflux_lab(double *x, double *par){
	double Einc  = 4.3;//[GeV]
	//double Escat = 2.1;//[GeV]
	double Escat = par[0];//[GeV]
	double theta = x[0];	

	double Me=pow(511,-6.);//[GeV/c^2]
	double Mp=0.9382720;//[GeV/c^2]
	double Qsq=2*Einc*Escat*(1-cos(theta));
	double omega = Einc - Escat;
	double q2=Qsq+omega*omega;
	double kg=omega-Qsq/(2*Mp);
	double eps=1/(1+2*(q2/Qsq)*tan(theta/2)*tan(theta/2));
//if(theta>0.230&&theta<0.231){
//cout<<"Qsq="<<Qsq<<endl;
//cout<<"q2="<<q2<<endl;
//cout<<"kg="<<kg<<endl;
//cout<<"eps="<<eps<<endl;
//}

	double vpflux=Escat*kg/(137*2*PI*PI*Einc*Qsq*(1-eps));
	return vpflux;
}

//double Escat = sqrt(pL*pL+Me*Me);
//double Einc = 4.5424;
//double omega = Einc - Escat;
//double A=Me*Me*omega*omega/(4*Einc*Einc*Escat*Escat);
//double sinterm=sin(theta_ee/2)*sin(theta_ee/2);
//double a1=((Einc*Einc+Escat*Escat)/(2*Einc*Einc))/(A+sinterm);
//double a2=(Escat/Einc)*A/((A+sinterm)*(A+sinterm));
//double a3=((Einc+Escat)*(Einc+Escat)/(4*Einc*Einc))/(omega*omega/(4*Einc*Escat)+sinterm);
//double vpflux=(a1-a2-a3)/(137*4*PI*PI*omega);
void vpflux(){
gROOT->Reset();
gROOT->SetStyle("Plain");

TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
c1->SetLeftMargin(0.15);
c1->SetBottomMargin(0.15);
c1->SetTopMargin(0.15);
c1->SetRightMargin(0.15);
c1->Draw();


TH1F* h1 = c1 -> DrawFrame(0.15,0.001,0.3,0.01);
//h1->GetXaxis()->SetTitleFont(12);
h1->GetXaxis()->SetTitleSize(0.05);
h1->GetXaxis()->SetTitleOffset(1.05);
h1->GetXaxis()->SetTitle("theta [rad]");
h1->GetXaxis()->SetLabelOffset(0.02);
//h1->GetYaxis()->SetTitleFont(22);
h1->GetYaxis()->SetTitleOffset(1.1);
h1->GetYaxis()->SetTitleSize(0.05);
h1->GetYaxis()->SetTitle("vpflux [/GeV/sr]");


TF1* func = new TF1("func",vpflux_doi, 0.001, 0.3,1);
func->SetNpx(600);
func->SetParameter(0,2.1);
func->SetLineColor(kAzure);
func->SetLineWidth(4);
func->Draw("same");
cout<<"val_doi"<<func->Eval(13.2*PI/180)<<endl;

TF1* func2 = new TF1("func2",vpflux_kawama, 0.001, 0.3,1);
func2->SetNpx(600);
func2->SetParameter(0,2.1);
func2->SetLineColor(kBlack);
func2->SetLineWidth(4);
func2->Draw("same");
cout<<"val_kawama"<<func2->Eval(13.2*PI/180)<<endl;

TF1* func3 = new TF1("func3",vpflux_lab, 0.001, 0.3,1);
func3->SetNpx(600);
func3->SetParameter(0,2.1);
func3->SetLineColor(kGreen);
func3->SetLineWidth(4);
func3->Draw("same");
cout<<"val_lab"<<func3->Eval(13.2*PI/180)<<endl;
//TF1* func2 = new TF1("func2",sqrt_func, 200., 1400.,1);
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
