double Eg_func(double *x, double *par){
	double W = x[0];
	double MT = par[0];//Target Mass
	double Qsq = par[1];
	double temp = (W*W-MT*MT+Qsq)/(2.*MT);
	if(temp>=0.)return temp;
	else return -1000.;
	}

void W_Eg(){
gROOT->Reset();
gROOT->SetStyle("Plain");

TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
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
double M = 0.9382720;
TH1F* h1 = c1 -> DrawFrame(0.,0.,4.,10.);
//h1->GetXaxis()->SetTitleFont(12);
h1->GetXaxis()->SetTitleSize(0.05);
h1->GetXaxis()->SetTitleOffset(1.05);
h1->GetXaxis()->SetTitle("W [GeV]");
h1->GetXaxis()->SetLabelOffset(0.02);
//h1->GetYaxis()->SetTitleFont(22);
h1->GetYaxis()->SetTitleOffset(0.9);
h1->GetYaxis()->SetTitleSize(0.05);
h1->GetYaxis()->SetTitle("E_{#gamma} [GeV]");
TF1* func = new TF1("func",Eg_func,0.,4.0,2);
func->SetNpx(600);
func->SetParameter(0,M);
func->SetParameter(1,0.);//(GeV)^2
func->SetLineColor(kAzure);
func->SetLineWidth(4);
func->Draw("same");
TF1* func2 = new TF1("func2",Eg_func,0.,4.0,2);
func2->SetNpx(600);
func2->SetParameter(0,M);
func2->SetParameter(1,0.5);//(GeV)^2
func2->SetLineColor(kRed);
func2->SetLineWidth(4);
func2->Draw("same");
cout<<func->Eval(2.13)<<endl;
cout<<func2->Eval(2.13)<<endl;

}
