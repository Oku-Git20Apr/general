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


const int NIH   = 10;  // No. of IH
const int NTB   = 40;  // No. of TagB

////////////////
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
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->SetLineWidth(0);
  h1->SetTitleSize(0.05,"");
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
double gaus_pol0bg(double *x, double *par) {
  //par[0]=area (gaus)
  //par[1]=mean (gaus)
  //par[2]=sigma (gaus)
  //par[3]=bg (const)
  double val;
  double bg;
  double ga;

ga=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
bg =  par[3];
val = ga + bg;
  return val;
}

//____________________________________________________________________________________________

void TagB25_TDLR2SN(string root_file = "test.root"){
gStyle->SetOptStat("iem");
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.10);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadBottomMargin(0.15);
TString ifname1=;//H3L
  TFile *ifp1  = new TFile( Form("./root/%s",root_file.c_str()) );
      cout<<"input filename : "<<root_file<<endl;
TH1F *h_TagSumTDLL[NIH], *h_TagSumTDLR[NIH];
  //R
  h_TagSumTDLR[1]  = (TH1F*)ifp1    -> Get(Form("h_TagSumTDLR%d",1+1));
  SetTH1(h_TagSumTDLR[1], Form("ToF TagB Sum - TDL R%d",1+1), "coin. time[ns]", "count/0.1ns", 1, 0, 0);

  h_TagSumTDLR[1] ->  GetXaxis()-> SetRangeUser(-200,0);
int p_flag=1;
//cout<<"Do you want to save canvas?"<<endl;
//cout<<"yes :1, no :0"<<endl;
//cin>>p_flag;

TCanvas *c[1];
for(int i=0;i<1;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),900,1200);
                     }
#if 1
double peakpos =  h_TagSumTDLR[1] ->GetXaxis()->GetBinCenter(h_TagSumTDLR[1] ->GetMaximumBin() );
cout<<peakpos<<endl;
double xmin = peakpos - 60.;
double xmax = peakpos + 60.;
h_TagSumTDLR[1] ->GetXaxis()->SetRangeUser(xmin,xmax);
TF1 *f_bg = new TF1("f_bg","[0]",xmin+30,xmax-30);
double param_bg[1];
f_bg->SetParameter(0,100);
f_bg ->SetLineWidth(3);  f_bg ->SetLineColor(2);f_bg ->SetLineStyle(3);
h_TagSumTDLR[1] ->Fit(f_bg,"0QR","",peakpos-20,peakpos-15);
f_bg->GetParameters(&param_bg[0]);


TF1 *fadd = new TF1("fadd",gaus_pol0bg,xmin,xmax,4);
double param[4], para_error[4];
fadd -> SetParameter(0,1000.);//area
fadd -> SetParameter(1,peakpos);//peak pos
fadd -> SetParameter(2,7);//width
fadd -> SetParameter(3,param_bg[0]);
    fadd ->SetLineWidth(2);  fadd ->SetLineColor(2);fadd ->SetLineStyle(1);
h_TagSumTDLR[1] ->Fit(fadd,"0QR","",peakpos-15,peakpos+15);
fadd ->GetParameters(&param[0]);
fadd -> SetNpx(2000);
double g_area  = param[0];
double g_mean  = param[1];
double g_sigma = param[2];
double bg_p0   = param[3];

double S  = g_area*0.999;//3sigma
double N  = bg_p0*6*g_sigma;// +/- 3sigma
double N2 = bg_p0*40.;//B.G. level * coin. time gate 
TLine *tline_bmin = new TLine(g_mean-3*g_sigma,0,g_mean-3*g_sigma,bg_p0);
TLine *tline_bmax = new TLine(g_mean+3*g_sigma,0,g_mean+3*g_sigma,bg_p0);
tline_bmin->SetLineWidth(2.);
tline_bmin->SetLineColor(6);
tline_bmin->SetLineStyle(1);
tline_bmax->SetLineWidth(2.);
tline_bmax->SetLineColor(6);
tline_bmax->SetLineStyle(1);

TLine *tline_bmin2 = new TLine(g_mean-15.,0,g_mean-15.,bg_p0);
TLine *tline_bmax2 = new TLine(g_mean+15.,0,g_mean+15.,bg_p0);
tline_bmin2->SetLineWidth(2.);
tline_bmin2->SetLineColor(4);
tline_bmin2->SetLineStyle(2);
tline_bmax2->SetLineWidth(2.);
tline_bmax2->SetLineColor(4);
tline_bmax2->SetLineStyle(2);
//functions to show fitting result
TF1 *ga2 = new TF1("ga2","gausn",-600,200);
    ga2 ->SetLineWidth(2);  ga2 ->SetLineColor(2);ga2 ->SetLineStyle(2);
ga2 -> SetParameter(0, g_area);
ga2 -> SetParameter(1, g_mean);
ga2 -> SetParameter(2, g_sigma);
ga2 -> SetNpx(2000);
TF1 *f1 = new TF1("f1","[0]",xmin,xmax);//to show fitting result
    f1 ->SetLineWidth(2);  f1 ->SetLineColor(6);f1 ->SetLineStyle(1);
f1->SetParameter(0,bg_p0);
#endif
////
c[0]->Clear();
c[0]->Divide(1,3);
c[0]->cd(1);
gPad->SetLogy(1);
h_TagSumTDLR[1] ->Draw("");
//f_bg->Draw("same");
c[0]->cd(2);
gPad->SetLogy(0);
h_TagSumTDLR[1] ->Draw("");
fadd ->Draw("same");
f1   ->Draw("same");
ga2  ->Draw("same");
tline_bmin->Draw("same");tline_bmax->Draw("same");
tline_bmin2->Draw("same");tline_bmax2->Draw("same");

TH2F *h_frame = new TH2F("h_frame","h_frame",10,0,1,10,0,1);
  SetTH2(h_frame,"","","");
  h_frame->GetXaxis()->SetNdivisions(000);
  h_frame->GetYaxis()->SetNdivisions(000);
TLatex *tex_sigma  = new TLatex(0.5,0.7,Form("#sigma_{coin}=%.01lf ps",1000.*g_sigma));
TLatex *tex_signal = new TLatex(0.5,0.5,Form("Gaus area=%.01lf",S) );
TLatex *tex_noise  = new TLatex(0.5,0.3,Form("B.G. area=%.01lf",N) );
TLatex *tex_noise2 = new TLatex(0.5,0.1,Form("B.G.(wide) area=%.01lf",N2) );
string runname = root_file;
runname.erase(runname.size()-10);
runname.insert(0,"run");
TLatex *tex_runname = new TLatex(0.5,0.9,runname.c_str() );
tex_runname  -> SetTextSize(0.080);
tex_runname  -> SetTextAlign(22);

tex_sigma  -> SetTextSize(0.080);
tex_sigma  -> SetTextAlign(22);
tex_signal -> SetTextSize(0.080);
tex_signal -> SetTextAlign(22);
tex_noise  -> SetTextSize(0.080);
tex_noise  -> SetTextAlign(22);
tex_noise2 -> SetTextSize(0.080);
tex_noise2 -> SetTextAlign(22);
c[0]->cd(3);
//gStyle->SetOptStat("");
h_frame->SetStats(0);
h_frame->Draw("");
tex_runname ->Draw("same");
tex_sigma   ->Draw("same");
tex_signal  ->Draw("same");
tex_noise   ->Draw("same");
tex_noise2  ->Draw("same");

if(p_flag==1){
  string ofname_pdf = root_file;
  ofname_pdf.erase(ofname_pdf.size()-5);
  ofname_pdf.append("_SN.pdf");
  c[0]->Print(Form("./pdf/%s[",ofname_pdf.c_str())  );
  c[0]->Print(Form("./pdf/%s" ,ofname_pdf.c_str())  );
  c[0]->Print(Form("./pdf/%s]",ofname_pdf.c_str())  );
 }

} 

#if 0
c[4]->Clear();
c[4]->Divide(2,5);
for(int i=0;i<10;i++){
double xmin = h_TagSumTDLRU[24][i] ->GetXaxis()->GetBinCenter(h_TagSumTDLRU[24][i] ->GetMaximumBin() ) - 100;
double xmax = h_TagSumTDLRU[24][i] ->GetXaxis()->GetBinCenter(h_TagSumTDLRU[24][i] ->GetMaximumBin() ) + 100;
h_TagSumTDLRU[24][i] ->GetXaxis()->SetRangeUser(xmin,xmax);
c[4]->cd(i+1);
gPad->SetLogy(1);
h_TagSumTDLRU[24][i] ->Draw("");
}
c[5]->Clear();
c[5]->Divide(2,5);
for(int i=0;i<10;i++){
double xmin = h_TagSumTDLRD[24][i] ->GetXaxis()->GetBinCenter(h_TagSumTDLRD[24][i] ->GetMaximumBin() ) - 100;
double xmax = h_TagSumTDLRD[24][i] ->GetXaxis()->GetBinCenter(h_TagSumTDLRD[24][i] ->GetMaximumBin() ) + 100;
h_TagSumTDLRD[24][i] ->GetXaxis()->SetRangeUser(xmin,xmax);
c[5]->cd(i+1);
gPad->SetLogy(1);
h_TagSumTDLRD[24][i] ->Draw("");
c[1]->Clear();
c[1]->Divide(2,5);
for(int i=0;i<10;i++){
double xmin = h_TagSumTDLLU[24][i] ->GetXaxis()->GetBinCenter(h_TagSumTDLLU[24][i] ->GetMaximumBin() ) - 100;
double xmax = h_TagSumTDLLU[24][i] ->GetXaxis()->GetBinCenter(h_TagSumTDLLU[24][i] ->GetMaximumBin() ) + 100;
h_TagSumTDLLU[24][i] ->GetXaxis()->SetRangeUser(xmin,xmax);
c[1]->cd(i+1);
gPad->SetLogy(1);
h_TagSumTDLLU[24][i] ->Draw("");
}
c[2]->Clear();
c[2]->Divide(2,5);
for(int i=0;i<10;i++){
double xmin = h_TagSumTDLLD[24][i] ->GetXaxis()->GetBinCenter(h_TagSumTDLLD[24][i] ->GetMaximumBin() ) - 100;
double xmax = h_TagSumTDLLD[24][i] ->GetXaxis()->GetBinCenter(h_TagSumTDLLD[24][i] ->GetMaximumBin() ) + 100;
h_TagSumTDLLD[24][i] ->GetXaxis()->SetRangeUser(xmin,xmax);
c[2]->cd(i+1);
gPad->SetLogy(1);
h_TagSumTDLLD[24][i] ->Draw("");
}
}
#endif
