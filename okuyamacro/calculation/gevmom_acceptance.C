double mom_acpt(double *x, double *par){
	double x0 = x[0];
	double MT = par[0];//Target Mass
	double MM = par[1];//Missing Mass
	double MK = 0.494;
	double Ee = 4.3;
	double cosine = TMath::Cos(0.05);
	double Ex=TMath::Sqrt(x0*x0+MK*MK);
	double Qsq=2*4.3*(1-cosine);
	double Me=0.000511;
	double pe=TMath::Sqrt(Ee*Ee-Me*Me);
	double Eep=TMath::Sqrt(x0*x0+Me*Me);//Ee'
//	double Qsq=-(Ee*Ee-2*Ee*Eep+Eep*Eep)+pe*pe+x0*x0-2*pe*x0*TMath::Cos(13.2*PI/180.);
	
	//return 4300-(2*(MT+TMath::Sqrt(x0*x0+MK*MK))-TMath::Sqrt((MT+TMath::Sqrt(x0*x0+MK*MK))*((MT+TMath::Sqrt(x0*x0+MK*MK)))+2*x0*TMath::Cos(0.02)*(MM*MM-MT*MT-MK*MK+2*MT*TMath::Sqrt(x0*x0+MK*MK))))/2/x0/TMath::Cos(0.02);
//	return 4300.-(MM*MM-MT*MT-Ex*Ex+2*MT*Ex+x0*x0)/(2*MT-2*Ex+2*x0*cosine);
	return Ee-(MM*MM-MT*MT-Ex*Ex+2*MT*Ex+x0*x0)/(2*MT-2*Ex+2*x0*cosine-Qsq);
	}

void gevmom_acceptance(){
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
double M  = 0.940;
double MK = 0.454;
double ML = 1.115;
double MS = 1.190;
c1->cd()->DrawFrame(1.730,1.930,1.950,2.250);
TF1* func_acpt = new TF1("func_acpt",mom_acpt,1.700,1.950,2);
func_acpt->SetNpx(600);
func_acpt->SetParameter(0,M);
func_acpt->SetParameter(1,ML);
func_acpt->SetLineColor(kAzure);
func_acpt->SetLineWidth(4);
func_acpt->Draw("same");
TF1* func_acpt2 = new TF1("func_acpt2",mom_acpt,1.700,1.950,2);
func_acpt2->SetNpx(600);
func_acpt2->SetParameter(0,M);
func_acpt2->SetParameter(1,MS);
func_acpt2->SetLineColor(kCyan);
func_acpt2->SetLineWidth(4);
func_acpt2->Draw("same");

TLine *tl1, *tl2, *tl3, *tl4;
tl1 = new TLine(1.750,1.950,1.750,2.250);
tl2 = new TLine(1.920,1.950,1.920,2.250);
tl3 = new TLine(1.740,2.000,1.930,2.000);
tl4 = new TLine(1.740,2.200,1.930,2.200);
	tl1->SetLineWidth(1);
	tl1->SetLineColor(kBlack);
	tl1->Draw("same");
	tl2->SetLineWidth(1);
	tl2->SetLineColor(kBlack);
	tl2->Draw("same");
	tl3->SetLineWidth(1);
	tl3->SetLineColor(kBlack);
	tl3->Draw("same");
	tl4->SetLineWidth(1);
	tl4->SetLineColor(kBlack);
	tl4->Draw("same");

}
