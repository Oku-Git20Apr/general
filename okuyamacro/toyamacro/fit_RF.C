//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gaus_pol0bg(double *x, double *par) {
  //par[0]=area (gaus)
  //par[1]=mean (gaus)
  //par[2]=sigma (gaus)
  //par[3]=pol const
  //par[4]=pol tilt
  double val;
  double bg;
  double ga;

ga=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
bg =  par[3];
val = ga + bg;
  return val;
}

//____________________________________________________________________________________________

void fit_RF(){
gStyle->SetOptStat(0);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.06);
gStyle->SetPadLeftMargin(0.16);
gStyle->SetPadBottomMargin(0.15);
TString ifname1=Form("./root/10246_rfana.root");//
  TFile *ifp1  = new TFile(ifname1  );
      cout<<"input filename : "<<ifname1<<endl;
TH1F *h_ctime_tdll_RF[10];

h_ctime_tdll_RF[0]= (TH1F*)ifp1->Get("h_ctime_tdll_RF2");



double min =-1.;
double max = 1.;
double bg;
TF1 *f_bg = new TF1("f_bg","pol0",min,max);
    f_bg ->SetLineWidth(2);  f_bg ->SetLineColor(2);f_bg ->SetLineStyle(3);

h_ctime_tdll_RF[0]->Fit(f_bg,"0QR","",min,-0.5);
bg=f_bg->GetParameter(0);

TF1 *fadd = new TF1("fadd",gaus_pol0bg,min,max,4);
    fadd ->SetLineWidth(2);  fadd ->SetLineColor(2);fadd ->SetLineStyle(1);fadd ->SetNpx(1000);
fadd -> SetParameter(0,150.);
fadd -> SetParameter(1,0.);
fadd -> SetParameter(2,0.2);
fadd -> SetParameter(3,bg);

h_ctime_tdll_RF[0]->Fit(fadd,"R","",min,max);
}
