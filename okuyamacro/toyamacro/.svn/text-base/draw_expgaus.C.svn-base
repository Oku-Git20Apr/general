//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetTF1(TF1 *f, int LColor=2, int LWidth=0, int LStyle=1, int Npx=1000){
  f->SetLineColor(LColor);
  f->SetLineWidth(LWidth);
  f->SetLineStyle(LStyle);
  f->SetNpx(Npx);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expgaus(double *x, double *par) {
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;

// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
  }
//  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));

  return val;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expgaus_tes(double *x, double *par) {
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;

// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);

  }
//  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));

  return val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expgaus_gaus(double *x, double *par){
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  return expgaus(x,&par[0]) + par[4]*TMath::Gaus(x[0],par[3],par[2],1);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expexpgaus_gaus(double *x, double *par){
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  return expgaus(x,&par[0]) + expgaus(x,&par[4]) + par[8]*TMath::Gaus(x[0],par[10],par[9],1);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void draw_expgaus(){
gStyle->SetPadTopMargin(0.04);
gStyle->SetPadBottomMargin(0.20);
gStyle->SetPadRightMargin(0.04);
gStyle->SetPadLeftMargin(0.07);
TString ifname1=Form("./root/ana10226_10350.root");//
//TString ifname1=Form("./root/test.root");//
TFile *ifp1  = new TFile(ifname1  );

  double mean = 0.0;
 // double Ltau = 270.;  // Lambda lifetime
  double Ltau = 263.2;  // Lambda lifetime PDG
  double Tau, Sigma, Mean;
  double eTau, eSigma, eMean;

  double Sigma_ga[2], Mean_ga[2];
  double eSigma_ga[2], eMean_ga[2];

  double Sigma_bg, Mean_bg;
  double eSigma_bg, eMean_bg;
  vector<double> vTau, veTau;

  double sigmaBG ;
  double sigmaH3L;
  double sigmaL  ;

 //  sigma = sigmaH3L;

  double fluc = 0.0;//fluctuation of num of event of Lambda and H3Lambda (+/- hoge[ratio])
  int nn = 0, NN = 0;  // No. of Requested events to get statistics
#if 0
  gROOT->SetBatch(1);
  if(N<400)        NN = 4000;
  else if(N<1000)  NN = 10000;
  else if(N<4000)  NN = 20000;
  else if(N<10000) NN = 40000;
  else NN = N;
#endif

  TF1 *f, *f_tes, *f_Lam, *f_bg;
  TF1 *f_ga[2];

//initial parameters
double area  = 100;
double tau   = 0.2;
double sigma = 0.2;
double mean  = 0.;

 
// Fitting //
    double par[11],par_bg[11];
    f = new TF1("f",expgaus,-1,3,4);
    f->SetParNames("area","tau","#sigma","mean");
    f->SetParameter(0,area); //20 bin num.
    f->SetParameter(1,tau);  f->SetParameter(2,sigma);  f->SetParameter(3,mean);

    f_tes = new TF1("f_tes",expgaus_tes,-1,3,4);
    f_tes->SetParNames("area","tau","#sigma","mean");
    f_tes->SetParameter(0,area); //20 bin num.
    f_tes->SetParameter(1,tau);  f_tes->SetParameter(2,sigma);  f_tes->SetParameter(3,mean);

// Drawing //
  SetTF1(f    ,1,3,1,2000);
  SetTF1(f_tes,2,3,2,2000);

  TCanvas *c[1];
  for(int i=0;i<1;i++){
    c[i] =new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1400,700);
    c[i]->Divide(1,1,1E-5,1E-5);
  }
 // TCanvas *c1 =new TCanvas("c1","c1",900,700);
  c[0]->cd(1)->SetMargin(0.10,0.05,0.18,0.07);
  gPad->SetLogy(1);
  f     ->Draw();
  f_tes ->Draw("same");


 // c[0]->Print("pdf/draw_expgaus.pdf[");
 // c[0]->Print("pdf/draw_expgaus.pdf" );
 // c[1]->Print("pdf/draw_expgaus.pdf" );
 // c[2]->Print("pdf/draw_expgaus.pdf" );
 // c[3]->Print("pdf/draw_expgaus.pdf" );
 // c[3]->Print("pdf/draw_expgaus.pdf]");
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
