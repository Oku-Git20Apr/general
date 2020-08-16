//2018.5.15
//Y. Toyama
//fitting pdf obtained from simulation
#include "Settings.cc"
const int color[]={2,3,4,6,7,8,9};
const int lifetime[21]={150, 160, 170, 180, 190 ,
                        200, 210, 220, 230, 240 ,
                        250, 260, 270, 280, 290 ,
                        300, 310, 320, 330, 340 ,350};

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
double expgaus_gaus(double *x, double *par){
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  //par[4]=area of additional gaus
  //par[5]=mean of additional gaus
  //par[6]=width of additional gaus
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  return expgaus(x,&par[0]) + par[4]*TMath::Gaus(x[0],par[5],par[6],1);
}
//****************//
void fit_pdf(){

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.12);
  Settings *set=new Settings();
  TString ifname2=Form("./root/pdf.root");
  TFile *ifp2  = new TFile(ifname2  );
  cout<<"input filename2 : "<<ifname2<<endl;

  TH1F*h_pdf[21];
  for(int i=0;i<21;i++){
    //cout<<i<<endl;
    h_pdf[i]       = (TH1F*)ifp2->Get(Form("h_pdf_L%dps_R%dps",lifetime[i],180));
    double nevent =h_pdf[i]->Integral();
    //cout<<"NBin : "<<h_dtimeBG->GetNbinsX()<<endl;
    //h_pdf[i] ->SetMarkerStyle( 20+i  ); 
    //h_pdf[i] ->SetMarkerColor( color[i%7]  ); 
  }

  TF1 *f_expgaus = new TF1("f_expgaus",expgaus,-2,2,4);
  TF1 *f_expgaus_gaus = new TF1("f_expgaus_gaus",expgaus_gaus,-2,2,7);
  double param[7];
  param[0] = 0.025*h_pdf[15]->Integral();
  param[1] = 0.20;
  param[2] = 0.180;
  param[3] = 0;
  f_expgaus->SetParameters(&param[0]);
  h_pdf[15]->Fit(f_expgaus,"0LL","",-1,2);

  
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  gPad->SetLogy();
  h_pdf[15]->Draw("PE");
  f_expgaus->Draw("same");
}
