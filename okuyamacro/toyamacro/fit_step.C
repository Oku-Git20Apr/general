#include <iostream>
#include <sstream>
#include <iomanip>
#include <csignal>
#include <stdlib.h>
#include <climits>
#include <fstream>
#include <math.h>
#include <string>
#include <unistd.h>
using namespace std;

#include "TApplication.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include <TLegend.h>
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TString.h"
#include "TPaveText.h"
#include "TSystem.h"

#include "TRandom.h"

void SetTH1(TH1F *h1, TString hname, TString xname, TString yname, int LColor, int FStyle, int FColor){
  h1->SetTitle(hname);
  h1->GetXaxis()->SetTitle(xname);
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->SetTitle(yname);
  h1->GetYaxis()->CenterTitle();
  h1->SetMinimum(0.8);
  h1->SetLineColor(LColor);
  h1->SetFillStyle(FStyle);
  h1->SetFillColor(FColor);
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->SetLineWidth(0);
  h1->SetTitleSize(0.05,"");
  h1->SetTitleSize(0.05,"x");
  h1->SetTitleSize(0.05,"y");
  ((TGaxis*)h1->GetYaxis())->SetMaxDigits(3);
}
////////////////
void SetTH1(TH1D *h1, TString hname, TString xname, TString yname, int LColor, int FStyle, int FColor){
  h1->SetTitle(hname);
  h1->GetXaxis()->SetTitle(xname);
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->SetTitle(yname);
  h1->GetYaxis()->CenterTitle();
  h1->SetMinimum(0.8);
  h1->SetLineColor(LColor);
  h1->SetFillStyle(FStyle);
  h1->SetFillColor(FColor);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->SetLineWidth(0);
  h1->SetTitleSize(0.06,"x");
  h1->SetTitleSize(0.05,"y");
  h1->SetTitleSize(0.05,"");
  ((TGaxis*)h1->GetYaxis())->SetMaxDigits(3);
  //h1->GetYaxis()->SetNdivisions(505);
}
////////////////
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double pol2bg(double *x, double *par) {
  //par[0]=amp
  //par[1]=pol
  double val;

if(x[0]>par[1])  val =  par[0]* (x[0] - par[1])*(x[0] - par[1]) ;
else  val = 0. ;

  return val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gaus_pol2bg(double *x, double *par) {
  //par[0]=area (gaus)
  //par[1]=mean (gaus)
  //par[2]=sigma (gaus)
  //par[3]=amp
  //par[4]=pol
  double val;
  double bg;
  double ga;

ga=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
if(x[0]>par[4])  bg =  par[3]* (x[0] - par[4])*(x[0] - par[4]) ;
else  bg = 0. ;
val = ga + bg;
  return val;
}

//____________________________________________________________________________________________

void fit_step(){
gStyle->SetOptStat(0);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.06);
gStyle->SetPadLeftMargin(0.16);
gStyle->SetPadBottomMargin(0.15);
int runnum1 =1000;
TString ifname1=Form("./root/H3L_650_QFL_65000.root");//H3L
  TFile *ifp1  = new TFile(ifname1  );
      cout<<"input filename : "<<ifname1<<endl;

TH1F *h1_BL_greso_all;//H3L mom
h1_BL_greso_all = (TH1F*)ifp1->Get("h1_BL_greso_all");
//SetTH1(h1_BL_greso_all,"", "-B_{#Lambda} [MeV]", "counts/(1MeV)"      ,4,1000,1);
h1_BL_greso_all -> SetFillColor(0);
h1_BL_greso_all -> SetLineWidth(3);
double param_bg[3];
double param[6], para_error[6];
double unbound_min =10;
double unbound_max =25;
TF1 *f5 = new TF1("f5",pol2bg,-20,25,2);
f5->SetParameter(0,100);
f5->SetParameter(1,-10);
    f5 ->SetLineWidth(3);  f5 ->SetLineColor(2);f5 ->SetLineStyle(3);
h1_BL_greso_all ->Fit(f5,"0QR","",unbound_min,unbound_max);
//h1_BL_greso_lam ->Fit(f5,"0QR","",-20,20);
f5 -> GetParameters(&param_bg[0]);
TLine *tline_ubmin = new TLine(unbound_min,0,unbound_min,900);
TLine *tline_ubmax = new TLine(unbound_max,0,unbound_max,900);
tline_ubmin->SetLineWidth(2.5);
tline_ubmin->SetLineColor(1);
tline_ubmin->SetLineStyle(1);
tline_ubmax->SetLineWidth(2.5);
tline_ubmax->SetLineColor(1);
tline_ubmax->SetLineStyle(1);
TF1 *fadd = new TF1("fadd",gaus_pol2bg,-20,25,5);
fadd -> SetParameter(0,1000.);
fadd -> SetParameter(1,-0.13);
fadd -> SetParameter(2,3.1);
fadd -> SetParameter(3,param_bg[0]);
fadd -> SetParameter(4,param_bg[1]);
//fadd -> FixParameter(3,param_bg[0]);
//fadd -> FixParameter(4,param_bg[1]);
    fadd ->SetLineWidth(2);  fadd ->SetLineColor(1);fadd ->SetLineStyle(1);
h1_BL_greso_all->Fit(fadd,"0QR","",-15,25);
fadd ->GetParameters(&param[0]);
double g_area  = param[0];
double g_mean  = param[1];
double g_sigma = param[2];
double bg_p0   = param[3];
double bg_p1   = param[4];
//cout<<"BG param(2nd step) ="<<bg_p0<<" "<<bg_p1<<endl;

double erg_area =fadd ->GetParError(0);
double erg_mean =fadd ->GetParError(1);
double erg_sigma=fadd ->GetParError(2);
double erbg_p0  =fadd ->GetParError(3);
double erbg_p1  =fadd ->GetParError(4);

TF1 *ga2 = new TF1("ga2","gausn",-20,20);
    ga2 ->SetLineWidth(2);  ga2 ->SetLineColor(2);ga2 ->SetLineStyle(2);
ga2 -> SetParameter(0, g_area);
ga2 -> SetParameter(1, g_mean);
ga2 -> SetParameter(2, g_sigma);

TF1 *f6 = new TF1("f6",pol2bg,-20,25,2);//to show fitting result
    f6 ->SetLineWidth(3);  f6 ->SetLineColor(4);f6 ->SetLineStyle(2);
f6->SetParameter(0,bg_p0);
f6->SetParameter(1,bg_p1);

int p_flag=0;
cout<<"Do you want to save canvas?"<<endl;
cout<<"yes :1, no :0"<<endl;
cin>>p_flag;

TCanvas *c[2];
for(int i=0;i<2;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800/2,1500/2);}
c[0]->Clear();
gPad->SetLogy(0);h1_BL_greso_all ->Draw("");f5->Draw("same");
tline_ubmin->Draw("same");tline_ubmax->Draw("same");
c[1]->Clear();
gPad->SetLogy(0);h1_BL_greso_all ->Draw("");fadd->Draw("same");ga2->Draw("same");f6->Draw("same");

if(p_flag==1){
//  c[0]->Print(Form("../pdf/kaondist_run%04d.pdf[",runnum1)  );
//  c[0]->Print(Form("../pdf/kaondist_run%04d.pdf" ,runnum1)  );
//  c[1]->Print(Form("../pdf/kaondist_run%04d.pdf" ,runnum1)  );
//  c[2]->Print(Form("../pdf/kaondist_run%04d.pdf" ,runnum1)  );
//  c[3]->Print(Form("../pdf/kaondist_run%04d.pdf" ,runnum1)  );
//  c[3]->Print(Form("../pdf/kaondist_run%04d.pdf]",runnum1)  );
//
  c[0]->Print(Form("../pdf/fit_example.pdf[",runnum1)  );
  c[0]->Print(Form("../pdf/fit_example.pdf" ,runnum1)  );
  c[1]->Print(Form("../pdf/fit_example.pdf" ,runnum1)  );
  c[1]->Print(Form("../pdf/fit_example.pdf]",runnum1)  );
 }

} 
