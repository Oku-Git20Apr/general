//K. Okuyama (Sep. 10, 2021)
double mom_acpt(double *x, double *par){
	const double PI=4.*atan(1.);
	double x0 = x[0];
	double MT = par[0];//Target Mass
	double MM = par[1];//Missing Mass
	double MK = 494.;//MeV
	double MP = 938.272081;
	//double Ee = 4312.;//MeV
	double Ee = 4240.;//MeV
	//double cosine = TMath::Cos(0.05);
	//double cosine = TMath::Cos(13.2*PI/180.);
	//double cosine2 = TMath::Cos(2.*13.2*PI/180.);
	//double cosine = TMath::Cos(0.23);
	//double cosine2 = TMath::Cos(0.46);
	//double cosine = TMath::Cos(0.25);
	//double cosine2 = TMath::Cos(0.50);
	double cosine = TMath::Cos(0.20);//ee'
	double cosine2 = TMath::Cos(0.31);//ee'+eK
	double Ex=TMath::Sqrt(x0*x0+MP*MP);
	//double Qsq=2.*4318.*2200.*(1-cosine);
	//double Qsq=0.476*1000000.;//nnL
	double Qsq=0.15*1000000.;
	double Me=0.511;//MeV
	double pe=TMath::Sqrt(Ee*Ee-Me*Me);
	//double Eep=TMath::Sqrt(x0*x0+Me*Me);//Ee'
//	double Qsq=-(Ee*Ee-2*Ee*Eep+Eep*Eep)+pe*pe+x0*x0-2*pe*x0*TMath::Cos(13.2*PI/180.);
	
	//return 4300-(2*(MT+TMath::Sqrt(x0*x0+MK*MK))-TMath::Sqrt((MT+TMath::Sqrt(x0*x0+MK*MK))*((MT+TMath::Sqrt(x0*x0+MK*MK)))+2*x0*TMath::Cos(0.02)*(MM*MM-MT*MT-MK*MK+2*MT*TMath::Sqrt(x0*x0+MK*MK))))/2/x0/TMath::Cos(0.02);
//	return 4300.-(MM*MM-MT*MT-Ex*Ex+2*MT*Ex+x0*x0)/(2*MT-2*Ex+2*x0*cosine);
	//return Ee-(MM*MM-MT*MT-Ex*Ex+2*MT*Ex+x0*x0)/(2*MT-2*Ex+2*x0*cosine-Qsq);
	//cout<<"MM="<<MM<<endl;
	//cout<<"MT="<<MT<<endl;
	//cout<<"Ex="<<Ex<<endl;
	//cout<<"x0="<<x0<<endl;
	//cout<<"Ee="<<Ee<<endl;
	//cout<<"pe="<<pe<<endl;
	//cout<<"Qsq="<<Qsq<<endl;
	//return (-1.*MM*MM+MT*MT+Ex*Ex+2.*(MT-Ex)*Ee-2.*MT*Ex-x0*x0+x0*(pe-2.1)*cos(13.2*PI/180.))/(2.*(MT-Ex));
	//return Ee+(-1.*MM*MM+MT*MT+MK*MK-2*MT*Ex-Qsq)/(2*MT-2*Ex+2*x0*cos(3.*PI/180.));
	//double temp = (-1.*MM*MM+2.*Me*Me+MT*MT+Ex*Ex+2.*(MT-Ex)*Ee-2.*MT*Ex+2.*x0*pe*cos(13.2*PI/180.))/(2.*Ee+2.*(MT-Ex)-2.*pe*cos(13.2*PI/180.)+2.*x0*cos(26.4*PI/180.));
	//cout<<"pe'="<<temp<<endl;
	//return temp;
	 double a = 2.*(Ee+MT-Ex)-2.*pe*cosine+2.*x0*cosine2;
     double b = (Ee+MT-Ex)*(Ee+MT-Ex)-pe*pe-x0*x0+2.*x0*pe*cosine-MM*MM;
    // cout<<"a="<<a<<endl;
    // cout<<"b="<<b<<endl;
     return b/a;

	}

void next_acceptance(){
gROOT->Reset();
gROOT->SetStyle("Plain");

TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
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
double M = 938.272081;
double MK=493.677;
double Metap = 957.78;
double ML=1115.683;
double MS=1192.642;
//c1->cd()->DrawFrame(1730,1930,1950,2250);
//c1->cd()->DrawFrame(1720,1940,1980,2220);
//c1->cd()->DrawFrame(1730,1943,1950,2250);
//c1->cd()->DrawFrame(900,2500,1500,3200);
c1->cd()->DrawFrame(500,3000,500,3000);
TF1* func_acpt = new TF1("func_acpt",mom_acpt,900,1500,2);
func_acpt->SetNpx(600);
func_acpt->SetParameter(0,M);//target mass
func_acpt->SetParameter(1,ML);//hyp mass
func_acpt->SetLineColor(kAzure);
func_acpt->SetLineWidth(4);
func_acpt->Draw("same");
TF1* func_acpt2 = new TF1("func_acpt2",mom_acpt,900,1500,2);
func_acpt2->SetNpx(600);
func_acpt2->SetParameter(0,M);
func_acpt2->SetParameter(1,MS);
func_acpt2->SetLineColor(kCyan);
func_acpt2->SetLineWidth(4);
func_acpt2->Draw("same");

TLine *tl1, *tl2, *tl3, *tl4;
tl1 = new TLine(1760.,2010.,1760.,2160.);
tl2 = new TLine(1900.,2010.,1900.,2160.);
tl3 = new TLine(1760.,2010.,1900.,2010.);
tl4 = new TLine(1760.,2160.,1900.,2160.);
//tl1 = new TLine(1750.,1950.,1750.,2250.);
//tl2 = new TLine(1920.,1950.,1920.,2250.);
//tl3 = new TLine(1740.,2000.,1930.,2000.);
//tl4 = new TLine(1740.,2200.,1930.,2200.);
	//tl1->SetLineWidth(1);
	//tl1->SetLineColor(kBlack);
	//tl1->Draw("same");
	//tl2->SetLineWidth(1);
	//tl2->SetLineColor(kBlack);
	//tl2->Draw("same");
	//tl3->SetLineWidth(1);
	//tl3->SetLineColor(kBlack);
	//tl3->Draw("same");
	//tl4->SetLineWidth(1);
	//tl4->SetLineColor(kBlack);
	//tl4->Draw("same");

}
