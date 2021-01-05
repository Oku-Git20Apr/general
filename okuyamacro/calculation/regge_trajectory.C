void regge_trajectory(){
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
double k_family_x[5];
double k_family_y[5];
double k_family_xe[5];
double k_family_ye[5];
double ks_family_x[5];
double ks_family_y[5];
double ks_family_xe[5];
double ks_family_ye[5];
k_family_x[0]=0.494*0.494;//K^+(494)
k_family_x[1]=1.270*1.270;//K_1(1270)
k_family_x[2]=1.770*1.770;//K_2(1770)
k_family_x[3]=2.320*2.320;//K_4(2320)
k_family_x[4]=2.500*2.500;//K_5(2500)
k_family_y[0]=0.;
k_family_y[1]=1.;
k_family_y[2]=2.;
k_family_y[3]=3.;
k_family_y[4]=4.;
k_family_xe[0]=0.;
k_family_xe[1]=0.;
k_family_xe[2]=0.;
k_family_xe[3]=0.;
k_family_xe[4]=0.;
k_family_ye[0]=0.;
k_family_ye[1]=0.;
k_family_ye[2]=0.;
k_family_ye[3]=0.;
k_family_ye[4]=0.;
//%%%%%%%%%%%%%%%%%%%%%%%%//
ks_family_x[0]=0.892*0.892;//K^*(892)
ks_family_x[1]=1.430*1.430;//K_2^*(1430)
ks_family_x[2]=1.780*1.780;//K_3^*(1780)
ks_family_x[3]=2.045*2.045;//K_4^*(2045)
ks_family_x[4]=2.380*2.380;//K_5^*(2380)
ks_family_y[0]=1.;
ks_family_y[1]=2.;
ks_family_y[2]=3.;
ks_family_y[3]=4.;
ks_family_y[4]=5.;
ks_family_xe[0]=0.;
ks_family_xe[1]=0.;
ks_family_xe[2]=0.;
ks_family_xe[3]=0.;
ks_family_xe[4]=0.;
ks_family_ye[0]=0.;
ks_family_ye[1]=0.;
ks_family_ye[2]=0.;
ks_family_ye[3]=0.;
ks_family_ye[4]=0.;

  TGraphErrors *gr_kfamily = new TGraphErrors(5, k_family_x, k_family_y, k_family_xe, k_family_ye);
  gr_kfamily->SetMarkerSize(1.5);
  gr_kfamily->SetMarkerStyle(21);
  gr_kfamily->SetFillStyle(3004);
  //gr_kfamily->SetLineColor(2);
  gr_kfamily->SetFillColor(kRed);
  gr_kfamily->SetMarkerColor(2);
  TGraphErrors *gr_ksfamily = new TGraphErrors(5, ks_family_x, ks_family_y, ks_family_xe, ks_family_ye);
  gr_ksfamily->SetMarkerSize(1.5);
  gr_ksfamily->SetMarkerStyle(21);
  gr_ksfamily->SetFillStyle(3004);
  //gr_ksfamily->SetLineColor(2);
  gr_ksfamily->SetFillColor(kAzure);
  gr_ksfamily->SetMarkerColor(kAzure);



TH1F* h1 = c1 -> DrawFrame(0.,0.,7.,7.);
//c1->SetGrid(0);
//h1->GetXaxis()->SetTitleFont(12);
h1->GetXaxis()->SetTitleSize(0.05);
h1->GetXaxis()->SetTitleOffset(1.05);
h1->GetXaxis()->SetTitle("t [GeV^{2}]");
h1->GetXaxis()->SetLabelOffset(0.02);
//h1->GetYaxis()->SetTitleFont(22);
h1->GetYaxis()->SetTitleOffset(0.9);
h1->GetYaxis()->SetTitleSize(0.05);
h1->GetYaxis()->SetTitle("#alpha(t)");
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.15);
	gPad->SetTopMargin(0.15);
	gPad->SetBottomMargin(0.15);
c1->Update();

TF1* f_k = new TF1("f_k" ,"[0]*(x-[1])",0.,10.);
TF1* f_ks= new TF1("f_ks","[0]*(x-[1])+1.",0.,10.);
f_k->SetParameter(0,0.70);
f_k->SetParameter(1,0.494*0.494);
f_k->SetLineColor(kRed);
f_ks->SetParameter(0,0.85);
f_ks->SetParameter(1,0.892*0.892);
f_ks->SetLineColor(kAzure);

gr_kfamily->Draw("Psame");
gr_ksfamily->Draw("Psame");
f_k->Draw("same");
f_ks->Draw("same");

}
