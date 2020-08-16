#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
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
#include "TBox.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGaxis.h"

#include "TRandom.h"

#include "HodoParamMan.hh"//aaa

#define Calibration

static const double PI = 4.0*atan(1.);
static const double mrad_to_deg = 1./1000*180./PI;
const double Mp = 938.272046;          // proton       mass (MeV/c2)
const double c = 0.299792458;          // speed of light in vacuum (m/ns)
const int NTB   = 40;  // No. of TagB
const int NTF   =160;  // No. of TagB
const int NMRPC =  6;  // No. of MRPC 

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetTH1(TH1 *h, TString name, TString xname, TString yname, int LColor=1, int FStyle=0, int FColor=0){
  h->SetTitle(name);
  h->SetLineColor(LColor);
  h->SetLineWidth(0);
  h->SetFillStyle(FStyle);
  h->SetFillColor(FColor);
  h->SetMinimum(0.8);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(3);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetTH2(TH2 *h, TString name, TString xname, TString yname, double min=0.8, double MStyle=1, double MSize=1.0){
  h->SetTitle(name);
  h->SetMinimum(min);
  h->SetLineWidth(0);
  h->SetTitleSize(0.05,"");
  h->SetMarkerStyle(1);
  h->SetMarkerSize(0.1);
  h->SetMarkerColor(1);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(0.7);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(3);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SetGrErr(TGraphErrors *gr, TString name, TString xname, TString yname, int LColor, int MColor, int MStyle, double Yoffset, double min, double max){
  gr->SetTitle(name);
  gr->SetName(name);
  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleOffset(1.0);
  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->SetTitleOffset(1.0);
  gr->GetYaxis()->CenterTitle();
  gr->SetLineColor(LColor);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(0.5);
  gr->GetYaxis()->SetTitleOffset(Yoffset);
//  gr->GetYaxis()->SetRangeUser(min,max);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetGr(TGraph *gr, TString name, TString xname, TString yname, int LColor, int MColor, int MStyle, double Yoffset){
  gr->SetTitle(name);
  gr->SetName(name);
  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleOffset(1.0);
  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->SetTitleOffset(1.0);
  gr->GetYaxis()->CenterTitle();
  gr->SetLineColor(LColor);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(0.5);
  gr->GetYaxis()->SetTitleOffset(Yoffset);
//  gr->GetYaxis()->SetRangeUser(min,max);
}
////////////////////////
void SetTF1(TF1 *f, int LColor, int LStyle,double LWidth){
  f->SetLineColor(LColor);
  f->SetLineStyle(LStyle);
  f->SetLineWidth(LWidth);
}
////////////////////////
void SetTLatex(TLatex *tex, int TColor, double TSize,int Align){
tex -> SetTextSize(TSize); 
tex -> SetTextColor(TColor);
tex -> SetTextAlign(Align);
tex -> SetTextFont(42);
}
////////////////////////
void SetTBox(TBox *box, int FColor, int FStyle){
box -> SetFillColor(FColor);
box -> SetFillStyle(FStyle);
}
////////////////////////

class checkparam
{
 public:
         checkparam();
        ~checkparam();
  void analysis();
  void fit();
  void draw(); 
  void savecanvas(); 
  void SetMaxRun( int N )  { RNumMax = N; }
  void SetMinRun( int N )  { RNumMin = N; }
  bool set_filename(string RunFile);
  void set_graph();

  private:
    FILE *fp;
    int RNumMax,RNumMin;
    int ENum;

  HodoParamMan *HodoParamManager;

  TGraph *gr_ihlu_adc_gain[nIHL],   *gr_ihru_adc_gain[nIHL];
  TGraph *gr_ihld_adc_gain[nIHL],   *gr_ihrd_adc_gain[nIHL];
  TGraph *gr_ohvlu_adc_gain[nOHVL], *gr_ohvru_adc_gain[nOHVL];
  TGraph *gr_ohvld_adc_gain[nOHVL], *gr_ohvrd_adc_gain[nOHVL];
  TGraph *gr_ohhlu_adc_gain[nOHHL], *gr_ohhru_adc_gain[nOHHL];
  TGraph *gr_ohhld_adc_gain[nOHHL], *gr_ohhrd_adc_gain[nOHHL];

  TGraph *gr_ihlu_adc_pede[nIHL]  , *gr_ihru_adc_pede[nIHL];
  TGraph *gr_ihld_adc_pede[nIHL]  , *gr_ihrd_adc_pede[nIHL];
  TGraph *gr_ohvlu_adc_pede[nOHVL], *gr_ohvru_adc_pede[nOHVL];
  TGraph *gr_ohvld_adc_pede[nOHVL], *gr_ohvrd_adc_pede[nOHVL];
  TGraph *gr_ohhlu_adc_pede[nOHHL], *gr_ohhru_adc_pede[nOHHL];
  TGraph *gr_ohhld_adc_pede[nOHHL], *gr_ohhrd_adc_pede[nOHHL];

  TGraph *gr_ihlu_tdc_offs[nIHL]  , *gr_ihru_tdc_offs[nIHL];
  TGraph *gr_ihld_tdc_offs[nIHL]  , *gr_ihrd_tdc_offs[nIHL];
  TGraph *gr_ohvlu_tdc_offs[nOHVL], *gr_ohvru_tdc_offs[nOHVL];
  TGraph *gr_ohvld_tdc_offs[nOHVL], *gr_ohvrd_tdc_offs[nOHVL];
  TGraph *gr_ohhlu_tdc_offs[nOHHL], *gr_ohhru_tdc_offs[nOHHL];
  TGraph *gr_ohhld_tdc_offs[nOHHL], *gr_ohhrd_tdc_offs[nOHHL];

  TLatex *tex;
  TBox *box;

  double ihlu_adc_pede[nIHL][600],ihru_adc_pede[nIHL][600];
  double ihld_adc_pede[nIHL][600],ihrd_adc_pede[nIHL][600];
  double ohvlu_adc_pede[nOHVL][600],ohvru_adc_pede[nOHVL][600];
  double ohvld_adc_pede[nOHVL][600],ohvrd_adc_pede[nOHVL][600];
  double ohhlu_adc_pede[nOHHL][600],ohhru_adc_pede[nOHHL][600];
  double ohhld_adc_pede[nOHHL][600],ohhrd_adc_pede[nOHHL][600];

  double ihlu_adc_gain[nIHL][600],ihru_adc_gain[nIHL][600];
  double ihld_adc_gain[nIHL][600],ihrd_adc_gain[nIHL][600];
  double ohvlu_adc_gain[nOHVL][600],ohvru_adc_gain[nOHVL][600];
  double ohvld_adc_gain[nOHVL][600],ohvrd_adc_gain[nOHVL][600];
  double ohhlu_adc_gain[nOHHL][600],ohhru_adc_gain[nOHHL][600];
  double ohhld_adc_gain[nOHHL][600],ohhrd_adc_gain[nOHHL][600];

  double ihlu_tdc_offs[nIHL][600],ihru_tdc_offs[nIHL][600];
  double ihld_tdc_offs[nIHL][600],ihrd_tdc_offs[nIHL][600];
  double ohvlu_tdc_offs[nOHVL][600],ohvru_tdc_offs[nOHVL][600];
  double ohvld_tdc_offs[nOHVL][600],ohvrd_tdc_offs[nOHVL][600];
  double ohhlu_tdc_offs[nOHHL][600],ohhru_tdc_offs[nOHHL][600];
  double ohhld_tdc_offs[nOHHL][600],ohhrd_tdc_offs[nOHHL][600];

  double ihlu_adc_pede_max[nIHL]  ,ihru_adc_pede_max[nIHL];
  double ihld_adc_pede_max[nIHL]  ,ihrd_adc_pede_max[nIHL];
  double ohvlu_adc_pede_max[nOHVL],ohvru_adc_pede_max[nOHVL];
  double ohvld_adc_pede_max[nOHVL],ohvrd_adc_pede_max[nOHVL];
  double ohhlu_adc_pede_max[nOHHL],ohhru_adc_pede_max[nOHHL];
  double ohhld_adc_pede_max[nOHHL],ohhrd_adc_pede_max[nOHHL];
  double ihlu_adc_pede_min[nIHL]  ,ihru_adc_pede_min[nIHL];
  double ihld_adc_pede_min[nIHL]  ,ihrd_adc_pede_min[nIHL];
  double ohvlu_adc_pede_min[nOHVL],ohvru_adc_pede_min[nOHVL];
  double ohvld_adc_pede_min[nOHVL],ohvrd_adc_pede_min[nOHVL];
  double ohhlu_adc_pede_min[nOHHL],ohhru_adc_pede_min[nOHHL];
  double ohhld_adc_pede_min[nOHHL],ohhrd_adc_pede_min[nOHHL];

  double ihlu_adc_gain_max[nIHL]  ,ihru_adc_gain_max[nIHL];
  double ihld_adc_gain_max[nIHL]  ,ihrd_adc_gain_max[nIHL];
  double ohvlu_adc_gain_max[nOHVL],ohvru_adc_gain_max[nOHVL];
  double ohvld_adc_gain_max[nOHVL],ohvrd_adc_gain_max[nOHVL];
  double ohhlu_adc_gain_max[nOHHL],ohhru_adc_gain_max[nOHHL];
  double ohhld_adc_gain_max[nOHHL],ohhrd_adc_gain_max[nOHHL];
  double ihlu_adc_gain_min[nIHL]  ,ihru_adc_gain_min[nIHL];
  double ihld_adc_gain_min[nIHL]  ,ihrd_adc_gain_min[nIHL];
  double ohvlu_adc_gain_min[nOHVL],ohvru_adc_gain_min[nOHVL];
  double ohvld_adc_gain_min[nOHVL],ohvrd_adc_gain_min[nOHVL];
  double ohhlu_adc_gain_min[nOHHL],ohhru_adc_gain_min[nOHHL];
  double ohhld_adc_gain_min[nOHHL],ohhrd_adc_gain_min[nOHHL];

  double ihlu_tdc_offs_max[nIHL]  ,ihru_tdc_offs_max[nIHL];
  double ihld_tdc_offs_max[nIHL]  ,ihrd_tdc_offs_max[nIHL];
  double ohvlu_tdc_offs_max[nOHVL],ohvru_tdc_offs_max[nOHVL];
  double ohvld_tdc_offs_max[nOHVL],ohvrd_tdc_offs_max[nOHVL];
  double ohhlu_tdc_offs_max[nOHHL],ohhru_tdc_offs_max[nOHHL];
  double ohhld_tdc_offs_max[nOHHL],ohhrd_tdc_offs_max[nOHHL];
  double ihlu_tdc_offs_min[nIHL]  ,ihru_tdc_offs_min[nIHL];
  double ihld_tdc_offs_min[nIHL]  ,ihrd_tdc_offs_min[nIHL];
  double ohvlu_tdc_offs_min[nOHVL],ohvru_tdc_offs_min[nOHVL];
  double ohvld_tdc_offs_min[nOHVL],ohvrd_tdc_offs_min[nOHVL];
  double ohhlu_tdc_offs_min[nOHHL],ohhru_tdc_offs_min[nOHHL];
  double ohhld_tdc_offs_min[nOHHL],ohhrd_tdc_offs_min[nOHHL];

  int run_ihlu_adc_pede_max[nIHL]  ,run_ihru_adc_pede_max[nIHL];
  int run_ihld_adc_pede_max[nIHL]  ,run_ihrd_adc_pede_max[nIHL];
  int run_ohvlu_adc_pede_max[nOHVL],run_ohvru_adc_pede_max[nOHVL];
  int run_ohvld_adc_pede_max[nOHVL],run_ohvrd_adc_pede_max[nOHVL];
  int run_ohhlu_adc_pede_max[nOHHL],run_ohhru_adc_pede_max[nOHHL];
  int run_ohhld_adc_pede_max[nOHHL],run_ohhrd_adc_pede_max[nOHHL];
  int run_ihlu_adc_pede_min[nIHL]  ,run_ihru_adc_pede_min[nIHL];
  int run_ihld_adc_pede_min[nIHL]  ,run_ihrd_adc_pede_min[nIHL];
  int run_ohvlu_adc_pede_min[nOHVL],run_ohvru_adc_pede_min[nOHVL];
  int run_ohvld_adc_pede_min[nOHVL],run_ohvrd_adc_pede_min[nOHVL];
  int run_ohhlu_adc_pede_min[nOHHL],run_ohhru_adc_pede_min[nOHHL];
  int run_ohhld_adc_pede_min[nOHHL],run_ohhrd_adc_pede_min[nOHHL];

  int run_ihlu_adc_gain_max[nIHL]  ,run_ihru_adc_gain_max[nIHL];
  int run_ihld_adc_gain_max[nIHL]  ,run_ihrd_adc_gain_max[nIHL];
  int run_ohvlu_adc_gain_max[nOHVL],run_ohvru_adc_gain_max[nOHVL];
  int run_ohvld_adc_gain_max[nOHVL],run_ohvrd_adc_gain_max[nOHVL];
  int run_ohhlu_adc_gain_max[nOHHL],run_ohhru_adc_gain_max[nOHHL];
  int run_ohhld_adc_gain_max[nOHHL],run_ohhrd_adc_gain_max[nOHHL];
  int run_ihlu_adc_gain_min[nIHL]  ,run_ihru_adc_gain_min[nIHL];
  int run_ihld_adc_gain_min[nIHL]  ,run_ihrd_adc_gain_min[nIHL];
  int run_ohvlu_adc_gain_min[nOHVL],run_ohvru_adc_gain_min[nOHVL];
  int run_ohvld_adc_gain_min[nOHVL],run_ohvrd_adc_gain_min[nOHVL];
  int run_ohhlu_adc_gain_min[nOHHL],run_ohhru_adc_gain_min[nOHHL];
  int run_ohhld_adc_gain_min[nOHHL],run_ohhrd_adc_gain_min[nOHHL];

  int run_ihlu_tdc_offs_max[nIHL]  ,run_ihru_tdc_offs_max[nIHL];
  int run_ihld_tdc_offs_max[nIHL]  ,run_ihrd_tdc_offs_max[nIHL];
  int run_ohvlu_tdc_offs_max[nOHVL],run_ohvru_tdc_offs_max[nOHVL];
  int run_ohvld_tdc_offs_max[nOHVL],run_ohvrd_tdc_offs_max[nOHVL];
  int run_ohhlu_tdc_offs_max[nOHHL],run_ohhru_tdc_offs_max[nOHHL];
  int run_ohhld_tdc_offs_max[nOHHL],run_ohhrd_tdc_offs_max[nOHHL];
  int run_ihlu_tdc_offs_min[nIHL]  ,run_ihru_tdc_offs_min[nIHL];
  int run_ihld_tdc_offs_min[nIHL]  ,run_ihrd_tdc_offs_min[nIHL];
  int run_ohvlu_tdc_offs_min[nOHVL],run_ohvru_tdc_offs_min[nOHVL];
  int run_ohvld_tdc_offs_min[nOHVL],run_ohvrd_tdc_offs_min[nOHVL];
  int run_ohhlu_tdc_offs_min[nOHHL],run_ohhru_tdc_offs_min[nOHHL];
  int run_ohhld_tdc_offs_min[nOHHL],run_ohhrd_tdc_offs_min[nOHHL];

  double rnum[600];
  int  nn;
 
    int run_num;
    int content;// number of filled event
    TCanvas *c[13];

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
checkparam::checkparam()
{

  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");
  gROOT->SetBatch(1);

  gStyle->SetOptDate(0);
  gStyle->SetOptFit(1);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetLineWidth(0);
  gStyle->SetOptDate(0);
//  gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleTextColor(1);
  gStyle->SetGridWidth(0);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetLineWidth(0);
  gStyle->SetNdivisions(510); // tertiary*10000 + secondary*100 + first
  gStyle->SetOptStat("iMen");
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);

  //const Int_t NRGBs = 5;
  //const Int_t NCont = 255;
  //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  //Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  //Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  //Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  //TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  //gStyle->SetNumberContours(NCont);
      
  for(int j=0; j<600; j++){
    for(int i=0; i<nIHL; i++){
      ihlu_adc_gain[i][j]=ihld_adc_gain[i][j]=0.;
      ihru_adc_gain[i][j]=ihrd_adc_gain[i][j]=0.;
    }
    for(int i=0; i<nOHVL; i++){
      ohvlu_adc_gain[i][j]=ohvld_adc_gain[i][j]=0.;
      ohvru_adc_gain[i][j]=ohvrd_adc_gain[i][j]=0.;
    }
    for(int i=0; i<nOHHL; i++){
      ohhlu_adc_gain[i][j]=ohhld_adc_gain[i][j]=0.;
      ohhru_adc_gain[i][j]=ohhrd_adc_gain[i][j]=0.;
    }
    
    rnum[j]=0.;
  }

 for(int i=0; i<nIHL; i++){
  ihlu_adc_pede_max[i] = ihru_adc_pede_max[i] =   -1.;
  ihld_adc_pede_max[i] = ihrd_adc_pede_max[i] =   -1.;
  ihlu_adc_pede_min[i] = ihru_adc_pede_min[i] = 9999.;
  ihld_adc_pede_min[i] = ihrd_adc_pede_min[i] = 9999.;

  ihlu_adc_gain_max[i] = ihru_adc_gain_max[i] =   -1.;
  ihld_adc_gain_max[i] = ihrd_adc_gain_max[i] =   -1.;
  ihlu_adc_gain_min[i] = ihru_adc_gain_min[i] = 9999.;
  ihld_adc_gain_min[i] = ihrd_adc_gain_min[i] = 9999.;

  ihlu_tdc_offs_max[i] = ihru_tdc_offs_max[i] =    0.;
  ihld_tdc_offs_max[i] = ihrd_tdc_offs_max[i] =    0.;
  ihlu_tdc_offs_min[i] = ihru_tdc_offs_min[i] = 9999.;
  ihld_tdc_offs_min[i] = ihrd_tdc_offs_min[i] = 9999.;
 }

 for(int i=0; i<nOHVL; i++){
  ohvlu_adc_pede_max[i] = ohvru_adc_pede_max[i] =   -1.;
  ohvld_adc_pede_max[i] = ohvrd_adc_pede_max[i] =   -1.;
  ohvlu_adc_pede_min[i] = ohvru_adc_pede_min[i] = 9999.;
  ohvld_adc_pede_min[i] = ohvrd_adc_pede_min[i] = 9999.;

  ohvlu_adc_gain_max[i] = ohvru_adc_gain_max[i] =   -1.;
  ohvld_adc_gain_max[i] = ohvrd_adc_gain_max[i] =   -1.;
  ohvlu_adc_gain_min[i] = ohvru_adc_gain_min[i] = 9999.;
  ohvld_adc_gain_min[i] = ohvrd_adc_gain_min[i] = 9999.;

  ohvlu_tdc_offs_max[i] = ohvru_tdc_offs_max[i] =   -1.;
  ohvld_tdc_offs_max[i] = ohvrd_tdc_offs_max[i] =   -1.;
  ohvlu_tdc_offs_min[i] = ohvru_tdc_offs_min[i] = 9999.;
  ohvld_tdc_offs_min[i] = ohvrd_tdc_offs_min[i] = 9999.;
 }

 for(int i=0; i<nOHHL; i++){
  ohhlu_adc_pede_max[i] = ohhru_adc_pede_max[i]=   -1.;
  ohhld_adc_pede_max[i] = ohhrd_adc_pede_max[i]=   -1.;
  ohhlu_adc_pede_min[i] = ohhru_adc_pede_min[i]= 9999.;
  ohhld_adc_pede_min[i] = ohhrd_adc_pede_min[i]= 9999.;

  ohhlu_adc_gain_max[i] = ohhru_adc_gain_max[i]=   -1.;
  ohhld_adc_gain_max[i] = ohhrd_adc_gain_max[i]=   -1.;
  ohhlu_adc_gain_min[i] = ohhru_adc_gain_min[i]= 9999.;
  ohhld_adc_gain_min[i] = ohhrd_adc_gain_min[i]= 9999.;

  ohhlu_tdc_offs_max[i] = ohhru_tdc_offs_max[i]=   -1.;
  ohhld_tdc_offs_max[i] = ohhrd_tdc_offs_max[i]=   -1.;
  ohhlu_tdc_offs_min[i] = ohhru_tdc_offs_min[i]= 9999.;
  ohhld_tdc_offs_min[i] = ohhrd_tdc_offs_min[i]= 9999.;
 }
  nn=0;

  tex = new TLatex(1,1,"");
  SetTLatex(tex,1,0.05,11);

  for(int i=0;i<13;i++){
    c[i]= new TCanvas(Form("c%d",i+1),Form("c%d",i+1),900,800 );
  }

}
////////////////////////////////////////////////////////////////////////////
checkparam::~checkparam(){
}
////////////////////////////////////////////////////////////////////////////
void checkparam::analysis(){

  int runnum;
  char str[256];
  while(fgets(str,144,fp)!=0){
    if(str[0]=='#') continue;
    if(sscanf(str,"%d",&runnum)==1){
      if(runnum==0);
      else {
	  cout << "runnum: " << runnum << endl;
      string paramname = Form("./param/hodo_combine/%d.param",runnum);
	  HodoParamManager = new HodoParamMan( paramname.c_str() );
	  HodoParamManager->Initialize();
      for(int i=0;i<nOHH;i++){
      //pedestal
	  ihlu_adc_pede[i][nn] = HodoParamManager->GetAdcOffset(0,CID_IH,i+1,0);//IHLU
	  ihld_adc_pede[i][nn] = HodoParamManager->GetAdcOffset(0,CID_IH,i+1,1);//IHLD
	  ihru_adc_pede[i][nn] = HodoParamManager->GetAdcOffset(1,CID_IH,i+1,0);//IHRU
	  ihrd_adc_pede[i][nn] = HodoParamManager->GetAdcOffset(1,CID_IH,i+1,1);//IHRD
         if(ihlu_adc_pede[i][nn]>ihlu_adc_pede_max[i]){ihlu_adc_pede_max[i]=ihlu_adc_pede[i][nn];run_ihlu_adc_pede_max[i]=runnum;}
         if(ihlu_adc_pede[i][nn]<ihlu_adc_pede_min[i]){ihlu_adc_pede_min[i]=ihlu_adc_pede[i][nn];run_ihlu_adc_pede_min[i]=runnum;}
         if(ihld_adc_pede[i][nn]>ihld_adc_pede_max[i]){ihld_adc_pede_max[i]=ihld_adc_pede[i][nn];run_ihld_adc_pede_max[i]=runnum;}
         if(ihld_adc_pede[i][nn]<ihld_adc_pede_min[i]){ihld_adc_pede_min[i]=ihld_adc_pede[i][nn];run_ihld_adc_pede_min[i]=runnum;}
         if(ihru_adc_pede[i][nn]>ihru_adc_pede_max[i]){ihru_adc_pede_max[i]=ihru_adc_pede[i][nn];run_ihru_adc_pede_max[i]=runnum;}
         if(ihru_adc_pede[i][nn]<ihru_adc_pede_min[i]){ihru_adc_pede_min[i]=ihru_adc_pede[i][nn];run_ihru_adc_pede_min[i]=runnum;}
         if(ihrd_adc_pede[i][nn]>ihrd_adc_pede_max[i]){ihrd_adc_pede_max[i]=ihrd_adc_pede[i][nn];run_ihrd_adc_pede_max[i]=runnum;}
         if(ihrd_adc_pede[i][nn]<ihrd_adc_pede_min[i]){ihrd_adc_pede_min[i]=ihrd_adc_pede[i][nn];run_ihrd_adc_pede_min[i]=runnum;}
      //Gain
	  ihlu_adc_gain[i][nn] = HodoParamManager->GetAdcGain(0,CID_IH,i+1,0);//IHLU
	  ihld_adc_gain[i][nn] = HodoParamManager->GetAdcGain(0,CID_IH,i+1,1);//IHLD
	  ihru_adc_gain[i][nn] = HodoParamManager->GetAdcGain(1,CID_IH,i+1,0);//IHRU
	  ihrd_adc_gain[i][nn] = HodoParamManager->GetAdcGain(1,CID_IH,i+1,1);//IHRD
         if(ihlu_adc_gain[i][nn]>ihlu_adc_gain_max[i]){ihlu_adc_gain_max[i]=ihlu_adc_gain[i][nn];run_ihlu_adc_gain_max[i]=runnum;}
         if(ihlu_adc_gain[i][nn]<ihlu_adc_gain_min[i]){ihlu_adc_gain_min[i]=ihlu_adc_gain[i][nn];run_ihlu_adc_gain_min[i]=runnum;}
         if(ihld_adc_gain[i][nn]>ihld_adc_gain_max[i]){ihld_adc_gain_max[i]=ihld_adc_gain[i][nn];run_ihld_adc_gain_max[i]=runnum;}
         if(ihld_adc_gain[i][nn]<ihld_adc_gain_min[i]){ihld_adc_gain_min[i]=ihld_adc_gain[i][nn];run_ihld_adc_gain_min[i]=runnum;}
         if(ihru_adc_gain[i][nn]>ihru_adc_gain_max[i]){ihru_adc_gain_max[i]=ihru_adc_gain[i][nn];run_ihru_adc_gain_max[i]=runnum;}
         if(ihru_adc_gain[i][nn]<ihru_adc_gain_min[i]){ihru_adc_gain_min[i]=ihru_adc_gain[i][nn];run_ihru_adc_gain_min[i]=runnum;}
         if(ihrd_adc_gain[i][nn]>ihrd_adc_gain_max[i]){ihrd_adc_gain_max[i]=ihrd_adc_gain[i][nn];run_ihrd_adc_gain_max[i]=runnum;}
         if(ihrd_adc_gain[i][nn]<ihrd_adc_gain_min[i]){ihrd_adc_gain_min[i]=ihrd_adc_gain[i][nn];run_ihrd_adc_gain_min[i]=runnum;}
      //t0
	  ihlu_tdc_offs[i][nn] = HodoParamManager->GetTdcOffset(0,CID_IH,i+1,0);//IHLU
	  ihld_tdc_offs[i][nn] = HodoParamManager->GetTdcOffset(0,CID_IH,i+1,1);//IHLD
	  ihru_tdc_offs[i][nn] = HodoParamManager->GetTdcOffset(1,CID_IH,i+1,0);//IHRU
	  ihrd_tdc_offs[i][nn] = HodoParamManager->GetTdcOffset(1,CID_IH,i+1,1);//IHRD
         if(ihlu_tdc_offs[i][nn]>ihlu_tdc_offs_max[i]){ihlu_tdc_offs_max[i]=ihlu_tdc_offs[i][nn];run_ihlu_tdc_offs_max[i]=runnum;}
         if(ihlu_tdc_offs[i][nn]<ihlu_tdc_offs_min[i]){ihlu_tdc_offs_min[i]=ihlu_tdc_offs[i][nn];run_ihlu_tdc_offs_min[i]=runnum;}
         if(ihld_tdc_offs[i][nn]>ihld_tdc_offs_max[i]){ihld_tdc_offs_max[i]=ihld_tdc_offs[i][nn];run_ihld_tdc_offs_max[i]=runnum;}
         if(ihld_tdc_offs[i][nn]<ihld_tdc_offs_min[i]){ihld_tdc_offs_min[i]=ihld_tdc_offs[i][nn];run_ihld_tdc_offs_min[i]=runnum;}
         if(ihru_tdc_offs[i][nn]>ihru_tdc_offs_max[i]){ihru_tdc_offs_max[i]=ihru_tdc_offs[i][nn];run_ihru_tdc_offs_max[i]=runnum;}
         if(ihru_tdc_offs[i][nn]<ihru_tdc_offs_min[i]){ihru_tdc_offs_min[i]=ihru_tdc_offs[i][nn];run_ihru_tdc_offs_min[i]=runnum;}
         if(ihrd_tdc_offs[i][nn]>ihrd_tdc_offs_max[i]){ihrd_tdc_offs_max[i]=ihrd_tdc_offs[i][nn];run_ihrd_tdc_offs_max[i]=runnum;}
         if(ihrd_tdc_offs[i][nn]<ihrd_tdc_offs_min[i]){ihrd_tdc_offs_min[i]=ihrd_tdc_offs[i][nn];run_ihrd_tdc_offs_min[i]=runnum;}
      }
      ///////
      //OHV//
      ///////
      for(int i=0;i<nOHV;i++){
      //pedestal
	  ohvlu_adc_pede[i][nn] = HodoParamManager->GetAdcOffset(2,CID_OH,i+1,0);//OHVLU
	  ohvld_adc_pede[i][nn] = HodoParamManager->GetAdcOffset(2,CID_OH,i+1,1);//OHVLD
	  ohvru_adc_pede[i][nn] = HodoParamManager->GetAdcOffset(3,CID_OH,i+1,0);//OHVRU
	  ohvrd_adc_pede[i][nn] = HodoParamManager->GetAdcOffset(3,CID_OH,i+1,1);//OHVRD
         if(ohvlu_adc_pede[i][nn]>ohvlu_adc_pede_max[i]){ohvlu_adc_pede_max[i]=ohvlu_adc_pede[i][nn];run_ohvlu_adc_pede_max[i]=runnum;}
         if(ohvlu_adc_pede[i][nn]<ohvlu_adc_pede_min[i]){ohvlu_adc_pede_min[i]=ohvlu_adc_pede[i][nn];run_ohvlu_adc_pede_min[i]=runnum;}
         if(ohvld_adc_pede[i][nn]>ohvld_adc_pede_max[i]){ohvld_adc_pede_max[i]=ohvld_adc_pede[i][nn];run_ohvld_adc_pede_max[i]=runnum;}
         if(ohvld_adc_pede[i][nn]<ohvld_adc_pede_min[i]){ohvld_adc_pede_min[i]=ohvld_adc_pede[i][nn];run_ohvld_adc_pede_min[i]=runnum;}
         if(ohvru_adc_pede[i][nn]>ohvru_adc_pede_max[i]){ohvru_adc_pede_max[i]=ohvru_adc_pede[i][nn];run_ohvru_adc_pede_max[i]=runnum;}
         if(ohvru_adc_pede[i][nn]<ohvru_adc_pede_min[i]){ohvru_adc_pede_min[i]=ohvru_adc_pede[i][nn];run_ohvru_adc_pede_min[i]=runnum;}
         if(ohvrd_adc_pede[i][nn]>ohvrd_adc_pede_max[i]){ohvrd_adc_pede_max[i]=ohvrd_adc_pede[i][nn];run_ohvrd_adc_pede_max[i]=runnum;}
         if(ohvrd_adc_pede[i][nn]<ohvrd_adc_pede_min[i]){ohvrd_adc_pede_min[i]=ohvrd_adc_pede[i][nn];run_ohvrd_adc_pede_min[i]=runnum;}
      //Gain
	  ohvlu_adc_gain[i][nn] = HodoParamManager->GetAdcGain(2,CID_OH,i+1,0);//OHVLU
	  ohvld_adc_gain[i][nn] = HodoParamManager->GetAdcGain(2,CID_OH,i+1,1);//OHVLD
	  ohvru_adc_gain[i][nn] = HodoParamManager->GetAdcGain(3,CID_OH,i+1,0);//OHVRU
	  ohvrd_adc_gain[i][nn] = HodoParamManager->GetAdcGain(3,CID_OH,i+1,1);//OHVRD
         if(ohvlu_adc_gain[i][nn]>ohvlu_adc_gain_max[i]){ohvlu_adc_gain_max[i]=ohvlu_adc_gain[i][nn];run_ohvlu_adc_gain_max[i]=runnum;}
         if(ohvlu_adc_gain[i][nn]<ohvlu_adc_gain_min[i]){ohvlu_adc_gain_min[i]=ohvlu_adc_gain[i][nn];run_ohvlu_adc_gain_min[i]=runnum;}
         if(ohvld_adc_gain[i][nn]>ohvld_adc_gain_max[i]){ohvld_adc_gain_max[i]=ohvld_adc_gain[i][nn];run_ohvld_adc_gain_max[i]=runnum;}
         if(ohvld_adc_gain[i][nn]<ohvld_adc_gain_min[i]){ohvld_adc_gain_min[i]=ohvld_adc_gain[i][nn];run_ohvld_adc_gain_min[i]=runnum;}
         if(ohvru_adc_gain[i][nn]>ohvru_adc_gain_max[i]){ohvru_adc_gain_max[i]=ohvru_adc_gain[i][nn];run_ohvru_adc_gain_max[i]=runnum;}
         if(ohvru_adc_gain[i][nn]<ohvru_adc_gain_min[i]){ohvru_adc_gain_min[i]=ohvru_adc_gain[i][nn];run_ohvru_adc_gain_min[i]=runnum;}
         if(ohvrd_adc_gain[i][nn]>ohvrd_adc_gain_max[i]){ohvrd_adc_gain_max[i]=ohvrd_adc_gain[i][nn];run_ohvrd_adc_gain_max[i]=runnum;}
         if(ohvrd_adc_gain[i][nn]<ohvrd_adc_gain_min[i]){ohvrd_adc_gain_min[i]=ohvrd_adc_gain[i][nn];run_ohvrd_adc_gain_min[i]=runnum;}
      //t0
	  ohvlu_tdc_offs[i][nn] = HodoParamManager->GetTdcOffset(2,CID_OH,i+1,0);//OHVLU
	  ohvld_tdc_offs[i][nn] = HodoParamManager->GetTdcOffset(2,CID_OH,i+1,1);//OHVLD
	  ohvru_tdc_offs[i][nn] = HodoParamManager->GetTdcOffset(3,CID_OH,i+1,0);//OHVRU
	  ohvrd_tdc_offs[i][nn] = HodoParamManager->GetTdcOffset(3,CID_OH,i+1,1);//OHVRD
         if(ohvlu_tdc_offs[i][nn]>ohvlu_tdc_offs_max[i]){ohvlu_tdc_offs_max[i]=ohvlu_tdc_offs[i][nn];run_ohvlu_tdc_offs_max[i]=runnum;}
         if(ohvlu_tdc_offs[i][nn]<ohvlu_tdc_offs_min[i]){ohvlu_tdc_offs_min[i]=ohvlu_tdc_offs[i][nn];run_ohvlu_tdc_offs_min[i]=runnum;}
         if(ohvld_tdc_offs[i][nn]>ohvld_tdc_offs_max[i]){ohvld_tdc_offs_max[i]=ohvld_tdc_offs[i][nn];run_ohvld_tdc_offs_max[i]=runnum;}
         if(ohvld_tdc_offs[i][nn]<ohvld_tdc_offs_min[i]){ohvld_tdc_offs_min[i]=ohvld_tdc_offs[i][nn];run_ohvld_tdc_offs_min[i]=runnum;}
         if(ohvru_tdc_offs[i][nn]>ohvru_tdc_offs_max[i]){ohvru_tdc_offs_max[i]=ohvru_tdc_offs[i][nn];run_ohvru_tdc_offs_max[i]=runnum;}
         if(ohvru_tdc_offs[i][nn]<ohvru_tdc_offs_min[i]){ohvru_tdc_offs_min[i]=ohvru_tdc_offs[i][nn];run_ohvru_tdc_offs_min[i]=runnum;}
         if(ohvrd_tdc_offs[i][nn]>ohvrd_tdc_offs_max[i]){ohvrd_tdc_offs_max[i]=ohvrd_tdc_offs[i][nn];run_ohvrd_tdc_offs_max[i]=runnum;}
         if(ohvrd_tdc_offs[i][nn]<ohvrd_tdc_offs_min[i]){ohvrd_tdc_offs_min[i]=ohvrd_tdc_offs[i][nn];run_ohvrd_tdc_offs_min[i]=runnum;}
      }
      ///////
      //OHH//
      ///////
      for(int i=0;i<nOHH;i++){
      //pedestal
	  ohhlu_adc_pede[i][nn] = HodoParamManager->GetAdcOffset(0,CID_OH,i+1,0);//OHHLU
	  ohhld_adc_pede[i][nn] = HodoParamManager->GetAdcOffset(0,CID_OH,i+1,1);//OHHLD
	  ohhru_adc_pede[i][nn] = HodoParamManager->GetAdcOffset(1,CID_OH,i+1,0);//OHHRU
	  ohhrd_adc_pede[i][nn] = HodoParamManager->GetAdcOffset(1,CID_OH,i+1,1);//OHHRD
         if(ohhlu_adc_pede[i][nn]>ohhlu_adc_pede_max[i]){ohhlu_adc_pede_max[i]=ohhlu_adc_pede[i][nn];run_ohhlu_adc_pede_max[i]=runnum;}
         if(ohhlu_adc_pede[i][nn]<ohhlu_adc_pede_min[i]){ohhlu_adc_pede_min[i]=ohhlu_adc_pede[i][nn];run_ohhlu_adc_pede_min[i]=runnum;}
         if(ohhld_adc_pede[i][nn]>ohhld_adc_pede_max[i]){ohhld_adc_pede_max[i]=ohhld_adc_pede[i][nn];run_ohhld_adc_pede_max[i]=runnum;}
         if(ohhld_adc_pede[i][nn]<ohhld_adc_pede_min[i]){ohhld_adc_pede_min[i]=ohhld_adc_pede[i][nn];run_ohhld_adc_pede_min[i]=runnum;}
         if(ohhru_adc_pede[i][nn]>ohhru_adc_pede_max[i]){ohhru_adc_pede_max[i]=ohhru_adc_pede[i][nn];run_ohhru_adc_pede_max[i]=runnum;}
         if(ohhru_adc_pede[i][nn]<ohhru_adc_pede_min[i]){ohhru_adc_pede_min[i]=ohhru_adc_pede[i][nn];run_ohhru_adc_pede_min[i]=runnum;}
         if(ohhrd_adc_pede[i][nn]>ohhrd_adc_pede_max[i]){ohhrd_adc_pede_max[i]=ohhrd_adc_pede[i][nn];run_ohhrd_adc_pede_max[i]=runnum;}
         if(ohhrd_adc_pede[i][nn]<ohhrd_adc_pede_min[i]){ohhrd_adc_pede_min[i]=ohhrd_adc_pede[i][nn];run_ohhrd_adc_pede_min[i]=runnum;}
      //Gain
	  ohhlu_adc_gain[i][nn] = HodoParamManager->GetAdcGain(0,CID_OH,i+1,0);//OHHLU
	  ohhld_adc_gain[i][nn] = HodoParamManager->GetAdcGain(0,CID_OH,i+1,1);//OHHLD
	  ohhru_adc_gain[i][nn] = HodoParamManager->GetAdcGain(1,CID_OH,i+1,0);//OHHRU
	  ohhrd_adc_gain[i][nn] = HodoParamManager->GetAdcGain(1,CID_OH,i+1,1);//OHHRD
         if(ohhlu_adc_gain[i][nn]>ohhlu_adc_gain_max[i]){ohhlu_adc_gain_max[i]=ohhlu_adc_gain[i][nn];run_ohhlu_adc_gain_max[i]=runnum;}
         if(ohhlu_adc_gain[i][nn]<ohhlu_adc_gain_min[i]){ohhlu_adc_gain_min[i]=ohhlu_adc_gain[i][nn];run_ohhlu_adc_gain_min[i]=runnum;}
         if(ohhld_adc_gain[i][nn]>ohhld_adc_gain_max[i]){ohhld_adc_gain_max[i]=ohhld_adc_gain[i][nn];run_ohhld_adc_gain_max[i]=runnum;}
         if(ohhld_adc_gain[i][nn]<ohhld_adc_gain_min[i]){ohhld_adc_gain_min[i]=ohhld_adc_gain[i][nn];run_ohhld_adc_gain_min[i]=runnum;}
         if(ohhru_adc_gain[i][nn]>ohhru_adc_gain_max[i]){ohhru_adc_gain_max[i]=ohhru_adc_gain[i][nn];run_ohhru_adc_gain_max[i]=runnum;}
         if(ohhru_adc_gain[i][nn]<ohhru_adc_gain_min[i]){ohhru_adc_gain_min[i]=ohhru_adc_gain[i][nn];run_ohhru_adc_gain_min[i]=runnum;}
         if(ohhrd_adc_gain[i][nn]>ohhrd_adc_gain_max[i]){ohhrd_adc_gain_max[i]=ohhrd_adc_gain[i][nn];run_ohhrd_adc_gain_max[i]=runnum;}
         if(ohhrd_adc_gain[i][nn]<ohhrd_adc_gain_min[i]){ohhrd_adc_gain_min[i]=ohhrd_adc_gain[i][nn];run_ohhrd_adc_gain_min[i]=runnum;}
      //t0
	  ohhlu_tdc_offs[i][nn] = HodoParamManager->GetTdcOffset(0,CID_OH,i+1,0);//OHHLU
	  ohhld_tdc_offs[i][nn] = HodoParamManager->GetTdcOffset(0,CID_OH,i+1,1);//OHHLD
	  ohhru_tdc_offs[i][nn] = HodoParamManager->GetTdcOffset(1,CID_OH,i+1,0);//OHHRU
	  ohhrd_tdc_offs[i][nn] = HodoParamManager->GetTdcOffset(1,CID_OH,i+1,1);//OHHRD
         if(ohhlu_tdc_offs[i][nn]>ohhlu_tdc_offs_max[i]){ohhlu_tdc_offs_max[i]=ohhlu_tdc_offs[i][nn];run_ohhlu_tdc_offs_max[i]=runnum;}
         if(ohhlu_tdc_offs[i][nn]<ohhlu_tdc_offs_min[i]){ohhlu_tdc_offs_min[i]=ohhlu_tdc_offs[i][nn];run_ohhlu_tdc_offs_min[i]=runnum;}
         if(ohhld_tdc_offs[i][nn]>ohhld_tdc_offs_max[i]){ohhld_tdc_offs_max[i]=ohhld_tdc_offs[i][nn];run_ohhld_tdc_offs_max[i]=runnum;}
         if(ohhld_tdc_offs[i][nn]<ohhld_tdc_offs_min[i]){ohhld_tdc_offs_min[i]=ohhld_tdc_offs[i][nn];run_ohhld_tdc_offs_min[i]=runnum;}
         if(ohhru_tdc_offs[i][nn]>ohhru_tdc_offs_max[i]){ohhru_tdc_offs_max[i]=ohhru_tdc_offs[i][nn];run_ohhru_tdc_offs_max[i]=runnum;}
         if(ohhru_tdc_offs[i][nn]<ohhru_tdc_offs_min[i]){ohhru_tdc_offs_min[i]=ohhru_tdc_offs[i][nn];run_ohhru_tdc_offs_min[i]=runnum;}
         if(ohhrd_tdc_offs[i][nn]>ohhrd_tdc_offs_max[i]){ohhrd_tdc_offs_max[i]=ohhrd_tdc_offs[i][nn];run_ohhrd_tdc_offs_max[i]=runnum;}
         if(ohhrd_tdc_offs[i][nn]<ohhrd_tdc_offs_min[i]){ohhrd_tdc_offs_min[i]=ohhrd_tdc_offs[i][nn];run_ohhrd_tdc_offs_min[i]=runnum;}
      }
     
      rnum[nn]=runnum;
      delete HodoParamManager;
      nn++;
      }
    }
  }
RNumMin=(int)rnum[0];
RNumMax=(int)rnum[nn-1];

}
////////////////////////////////////////////////////////////////////////////
bool checkparam::set_filename(string RunFile){
  if((fp=fopen(RunFile.c_str() ,"r"))==0){
    std::cerr << "file open fail" << std::endl;
    return false;
  }
  else return true;
}
////////////////////////////////////////////////////////////////////////////
void checkparam::set_graph(){
 /*IH*/
 for(int i=0;i<nIH;i++){
    //pedestal
    gr_ihlu_adc_pede[i]  = new TGraph( nn, rnum, ihlu_adc_pede[i]);//IHLU
    gr_ihld_adc_pede[i]  = new TGraph( nn, rnum, ihld_adc_pede[i]);//IHLD
    gr_ihru_adc_pede[i]  = new TGraph( nn, rnum, ihru_adc_pede[i]);//IHRU
    gr_ihrd_adc_pede[i]  = new TGraph( nn, rnum, ihrd_adc_pede[i]);//IHRD
    //Gain
    gr_ihlu_adc_gain[i]  = new TGraph( nn, rnum, ihlu_adc_gain[i]);//IHLU
    gr_ihld_adc_gain[i]  = new TGraph( nn, rnum, ihld_adc_gain[i]);//IHLD
    gr_ihru_adc_gain[i]  = new TGraph( nn, rnum, ihru_adc_gain[i]);//IHRU
    gr_ihrd_adc_gain[i]  = new TGraph( nn, rnum, ihrd_adc_gain[i]);//IHRD
    //t0
    gr_ihlu_tdc_offs[i]  = new TGraph( nn, rnum, ihlu_tdc_offs[i]);//IHLU
    gr_ihld_tdc_offs[i]  = new TGraph( nn, rnum, ihld_tdc_offs[i]);//IHLD
    gr_ihru_tdc_offs[i]  = new TGraph( nn, rnum, ihru_tdc_offs[i]);//IHRU
    gr_ihrd_tdc_offs[i]  = new TGraph( nn, rnum, ihrd_tdc_offs[i]);//IHRD
    //SetGr(TGraph *gr, TString name, TString xname, TString yname, int LColor, int MColor, int MStyle, double Yoffset)
    SetGr(gr_ihlu_adc_pede[i], Form("pedestal IHL%dU",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
    SetGr(gr_ihld_adc_pede[i], Form("pedestal IHL%dD",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
    SetGr(gr_ihru_adc_pede[i], Form("pedestal IHR%dU",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
    SetGr(gr_ihrd_adc_pede[i], Form("pedestal IHR%dD",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
    SetGr(gr_ihlu_adc_gain[i], Form("Gain IHL%dU"    ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
    SetGr(gr_ihld_adc_gain[i], Form("Gain IHL%dD"    ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
    SetGr(gr_ihru_adc_gain[i], Form("Gain IHR%dU"    ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
    SetGr(gr_ihrd_adc_gain[i], Form("Gain IHR%dD"    ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
    SetGr(gr_ihlu_tdc_offs[i], Form("t0 IHL%dU"      ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
    SetGr(gr_ihld_tdc_offs[i], Form("t0 IHL%dD"      ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
    SetGr(gr_ihru_tdc_offs[i], Form("t0 IHR%dU"      ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
    SetGr(gr_ihrd_tdc_offs[i], Form("t0 IHR%dD"      ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
 }
 /*OHV*/
 for(int i=0;i<nOHV;i++){
    //pedestal
    gr_ohvlu_adc_pede[i]  = new TGraph( nn, rnum, ohvlu_adc_pede[i]);//OHVLU
    gr_ohvld_adc_pede[i]  = new TGraph( nn, rnum, ohvld_adc_pede[i]);//OHVLD
    gr_ohvru_adc_pede[i]  = new TGraph( nn, rnum, ohvru_adc_pede[i]);//OHVRU
    gr_ohvrd_adc_pede[i]  = new TGraph( nn, rnum, ohvrd_adc_pede[i]);//OHVRD
    //Gain
    gr_ohvlu_adc_gain[i]  = new TGraph( nn, rnum, ohvlu_adc_gain[i]);//OHVLU
    gr_ohvld_adc_gain[i]  = new TGraph( nn, rnum, ohvld_adc_gain[i]);//OHVLD
    gr_ohvru_adc_gain[i]  = new TGraph( nn, rnum, ohvru_adc_gain[i]);//OHVRU
    gr_ohvrd_adc_gain[i]  = new TGraph( nn, rnum, ohvrd_adc_gain[i]);//OHVRD
    //t0
    gr_ohvlu_tdc_offs[i]  = new TGraph( nn, rnum, ohvlu_tdc_offs[i]);//OHVLU
    gr_ohvld_tdc_offs[i]  = new TGraph( nn, rnum, ohvld_tdc_offs[i]);//OHVLD
    gr_ohvru_tdc_offs[i]  = new TGraph( nn, rnum, ohvru_tdc_offs[i]);//OHVRU
    gr_ohvrd_tdc_offs[i]  = new TGraph( nn, rnum, ohvrd_tdc_offs[i]);//OHVRD
    SetGr(gr_ohvlu_adc_pede[i], Form("pedestal OHVL%dU",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohvld_adc_pede[i], Form("pedestal OHVL%dD",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohvru_adc_pede[i], Form("pedestal OHVR%dU",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohvrd_adc_pede[i], Form("pedestal OHVR%dD",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohvlu_adc_gain[i], Form("Gain OHVL%dU"    ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohvld_adc_gain[i], Form("Gain OHVL%dD"    ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohvru_adc_gain[i], Form("Gain OHVR%dU"    ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohvrd_adc_gain[i], Form("Gain OHVR%dD"    ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohvlu_tdc_offs[i], Form("t0 OHVL%dU"      ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohvld_tdc_offs[i], Form("t0 OHVL%dD"      ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohvru_tdc_offs[i], Form("t0 OHVR%dU"      ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohvrd_tdc_offs[i], Form("t0 OHVR%dD"      ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
 }
 /*OHH*/
 for(int i=0;i<nOHH;i++){
    //pedestal
    gr_ohhlu_adc_pede[i]  = new TGraph( nn, rnum, ohhlu_adc_pede[i]);//OHHLU
    gr_ohhld_adc_pede[i]  = new TGraph( nn, rnum, ohhld_adc_pede[i]);//OHHLD
    gr_ohhru_adc_pede[i]  = new TGraph( nn, rnum, ohhru_adc_pede[i]);//OHHRU
    gr_ohhrd_adc_pede[i]  = new TGraph( nn, rnum, ohhrd_adc_pede[i]);//OHHRD
    //Gain
    gr_ohhlu_adc_gain[i]  = new TGraph( nn, rnum, ohhlu_adc_gain[i]);//OHHLU
    gr_ohhld_adc_gain[i]  = new TGraph( nn, rnum, ohhld_adc_gain[i]);//OHHLD
    gr_ohhru_adc_gain[i]  = new TGraph( nn, rnum, ohhru_adc_gain[i]);//OHHRU
    gr_ohhrd_adc_gain[i]  = new TGraph( nn, rnum, ohhrd_adc_gain[i]);//OHHRD
    //t0
    gr_ohhlu_tdc_offs[i]  = new TGraph( nn, rnum, ohhlu_tdc_offs[i]);//OHHLU
    gr_ohhld_tdc_offs[i]  = new TGraph( nn, rnum, ohhld_tdc_offs[i]);//OHHLD
    gr_ohhru_tdc_offs[i]  = new TGraph( nn, rnum, ohhru_tdc_offs[i]);//OHHRU
    gr_ohhrd_tdc_offs[i]  = new TGraph( nn, rnum, ohhrd_tdc_offs[i]);//OHHRD
    SetGr(gr_ohhlu_adc_pede[i], Form("pedestal OHHL%dU",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohhld_adc_pede[i], Form("pedestal OHHL%dD",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohhru_adc_pede[i], Form("pedestal OHHR%dU",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohhrd_adc_pede[i], Form("pedestal OHHR%dD",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohhlu_adc_gain[i], Form("Gain OHHL%dU"    ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohhld_adc_gain[i], Form("Gain OHHL%dD"    ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohhru_adc_gain[i], Form("Gain OHHR%dU"    ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohhrd_adc_gain[i], Form("Gain OHHR%dD"    ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohhlu_tdc_offs[i], Form("t0 OHHL%dU"      ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohhld_tdc_offs[i], Form("t0 OHHL%dD"      ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohhru_tdc_offs[i], Form("t0 OHHR%dU"      ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
    SetGr(gr_ohhrd_tdc_offs[i], Form("t0 OHHR%dD"      ,i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.8);
 }
}
////////////////////////////////////////////////////////////////////////////
void checkparam::fit(){
}
////////////////////////////////////////////////////////////////////////////
void checkparam::draw(){

////////////
//Pedestal//
////////////
c[0]->Clear();
c[0]->Divide(4,5);
 for(int i=0;i<nIH;i++){
 c[0]->cd(i+1); gr_ihlu_adc_pede[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ihlu_adc_pede_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ihlu_adc_pede_max[i],run_ihlu_adc_pede_max[i],ihlu_adc_pede_min[i],run_ihlu_adc_pede_min[i]));
 c[0]->cd(i+11);gr_ihld_adc_pede[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ihld_adc_pede_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ihld_adc_pede_max[i],run_ihld_adc_pede_max[i],ihld_adc_pede_min[i],run_ihld_adc_pede_min[i]));
 }

c[1]->Clear();
c[1]->Divide(4,5);
 for(int i=0;i<nIH;i++){
 c[1]->cd(i+1); gr_ihru_adc_pede[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ihru_adc_pede_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ihru_adc_pede_max[i],run_ihru_adc_pede_max[i],ihru_adc_pede_min[i],run_ihru_adc_pede_min[i]));
 c[1]->cd(i+11);gr_ihrd_adc_pede[i]-> Draw("AP");                                                                                                                            
                tex-> DrawLatex(rnum[0], ihrd_adc_pede_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ihrd_adc_pede_max[i],run_ihrd_adc_pede_max[i],ihrd_adc_pede_min[i],run_ihrd_adc_pede_min[i]));
 }

c[2]->Clear();
c[2]->Divide(4,6);
 for(int i=0;i<nOHV;i++){
 c[2]->cd(i+1); gr_ohvlu_adc_pede[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ohvlu_adc_pede_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohvlu_adc_pede_max[i],run_ohvlu_adc_pede_max[i],ohvlu_adc_pede_min[i],run_ohvlu_adc_pede_min[i]));
 c[2]->cd(i+13);gr_ohvld_adc_pede[i]-> Draw("AP");                                                                                                                               
                tex-> DrawLatex(rnum[0], ohvld_adc_pede_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohvld_adc_pede_max[i],run_ohvld_adc_pede_max[i],ohvld_adc_pede_min[i],run_ohvld_adc_pede_min[i]));
 }

c[3]->Clear();
c[3]->Divide(4,6);
 for(int i=0;i<nOHV;i++){
 c[3]->cd(i+1); gr_ohvru_adc_pede[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ohvru_adc_pede_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohvru_adc_pede_max[i],run_ohvru_adc_pede_max[i],ohvru_adc_pede_min[i],run_ohvru_adc_pede_min[i]));
 c[3]->cd(i+13);gr_ohvrd_adc_pede[i]-> Draw("AP");                                                                                                                                                     
                tex-> DrawLatex(rnum[0], ohvrd_adc_pede_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohvrd_adc_pede_max[i],run_ohvrd_adc_pede_max[i],ohvrd_adc_pede_min[i],run_ohvrd_adc_pede_min[i]));
 }

c[4]->Clear();
c[4]->Divide(4,5);
 for(int i=0;i<nOHH;i++){
 c[4]->cd(i+1); gr_ohhlu_adc_pede[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ohhlu_adc_pede_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohhlu_adc_pede_max[i],run_ohhlu_adc_pede_max[i],ohhlu_adc_pede_min[i],run_ohhlu_adc_pede_min[i]));
 c[4]->cd(i+11); gr_ohhld_adc_pede[i]-> Draw("AP");                                                                                                                              
                tex-> DrawLatex(rnum[0], ohhld_adc_pede_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohhld_adc_pede_max[i],run_ohhld_adc_pede_max[i],ohhld_adc_pede_min[i],run_ohhld_adc_pede_min[i]));
 }

c[5]->Clear();
c[5]->Divide(4,5);
 for(int i=0;i<nOHH;i++){
 c[5]->cd(i+1); gr_ohhru_adc_pede[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ohhru_adc_pede_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohhru_adc_pede_max[i],run_ohhru_adc_pede_max[i],ohhru_adc_pede_min[i],run_ohhru_adc_pede_min[i]));
 c[5]->cd(i+11); gr_ohhrd_adc_pede[i]-> Draw("AP");                                                                                                                              
                tex-> DrawLatex(rnum[0], ohhrd_adc_pede_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohhrd_adc_pede_max[i],run_ohhrd_adc_pede_max[i],ohhrd_adc_pede_min[i],run_ohhrd_adc_pede_min[i]));
 }

////////
//Gain//
////////
/*
c[6]->Clear();
c[6]->Divide(4,5);
 for(int i=0;i<nIH;i++){
 c[6]->cd(i+1); gr_ihlu_adc_gain[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ihlu_adc_gain_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ihlu_adc_gain_max[i],run_ihlu_adc_gain_max[i],ihlu_adc_gain_min[i],run_ihlu_adc_gain_min[i]));
 c[6]->cd(i+11);gr_ihld_adc_gain[i]-> Draw("AP");                                                                                                                            
                tex-> DrawLatex(rnum[0], ihld_adc_gain_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ihld_adc_gain_max[i],run_ihld_adc_gain_max[i],ihld_adc_gain_min[i],run_ihld_adc_gain_min[i]));
 }

c[7]->Clear();
c[7]->Divide(4,5);
 for(int i=0;i<nIH;i++){
 c[7]->cd(i+1); gr_ihru_adc_gain[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ihru_adc_gain_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ihru_adc_gain_max[i],run_ihru_adc_gain_max[i],ihru_adc_gain_min[i],run_ihru_adc_gain_min[i]));
 c[7]->cd(i+11);gr_ihrd_adc_gain[i]-> Draw("AP");                                                                                                                            
                tex-> DrawLatex(rnum[0], ihrd_adc_gain_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ihrd_adc_gain_max[i],run_ihrd_adc_gain_max[i],ihrd_adc_gain_min[i],run_ihrd_adc_gain_min[i]));
 }

c[8]->Clear();
c[8]->Divide(4,6);
 for(int i=0;i<nOHV;i++){
 c[8]->cd(i+1); gr_ohvlu_adc_gain[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ohvlu_adc_gain_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohvlu_adc_gain_max[i],run_ohvlu_adc_gain_max[i],ohvlu_adc_gain_min[i],run_ohvlu_adc_gain_min[i]));
 c[8]->cd(i+13);gr_ohvld_adc_gain[i]-> Draw("AP");                                                                                                                               
                tex-> DrawLatex(rnum[0], ohvld_adc_gain_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohvld_adc_gain_max[i],run_ohvld_adc_gain_max[i],ohvld_adc_gain_min[i],run_ohvld_adc_gain_min[i]));
 }

c[9]->Clear();
c[9]->Divide(4,6);
 for(int i=0;i<nOHV;i++){
 c[9]->cd(i+1); gr_ohvru_adc_gain[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ohvru_adc_gain_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohvru_adc_gain_max[i],run_ohvru_adc_gain_max[i],ohvru_adc_gain_min[i],run_ohvru_adc_gain_min[i]));
 c[9]->cd(i+13);gr_ohvrd_adc_gain[i]-> Draw("AP");                                                                                                                               
                tex-> DrawLatex(rnum[0], ohvrd_adc_gain_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohvrd_adc_gain_max[i],run_ohvrd_adc_gain_max[i],ohvrd_adc_gain_min[i],run_ohvrd_adc_gain_min[i]));
 }

c[10]->Clear();
c[10]->Divide(4,5);
 for(int i=0;i<nOHH;i++){
 c[10]->cd(i+1); gr_ohhlu_adc_gain[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ohhlu_adc_gain_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohhlu_adc_gain_max[i],run_ohhlu_adc_gain_max[i],ohhlu_adc_gain_min[i],run_ohhlu_adc_gain_min[i]));
 c[10]->cd(i+13);gr_ohhld_adc_gain[i]-> Draw("AP");                                                                                                                              
                tex-> DrawLatex(rnum[0], ohhld_adc_gain_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohhld_adc_gain_max[i],run_ohhld_adc_gain_max[i],ohhld_adc_gain_min[i],run_ohhld_adc_gain_min[i]));
 }

c[11]->Clear();
c[11]->Divide(4,5);
 for(int i=0;i<nOHH;i++){
 c[11]->cd(i+1); gr_ohhru_adc_gain[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ohhru_adc_gain_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohhru_adc_gain_max[i],run_ohhru_adc_gain_max[i],ohhru_adc_gain_min[i],run_ohhru_adc_gain_min[i]));
 c[11]->cd(i+11);gr_ohhrd_adc_gain[i]-> Draw("AP");                                                                                                                              
                tex-> DrawLatex(rnum[0], ohhrd_adc_gain_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohhrd_adc_gain_max[i],run_ohhrd_adc_gain_max[i],ohhrd_adc_gain_min[i],run_ohhrd_adc_gain_min[i]));
 }
*/

//////
//t0//
//////
c[6]->Clear();
c[6]->Divide(4,5);
 for(int i=0;i<nIH;i++){
 c[6]->cd(i+1); gr_ihlu_tdc_offs[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ihlu_tdc_offs_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ihlu_tdc_offs_max[i],run_ihlu_tdc_offs_max[i],ihlu_tdc_offs_min[i],run_ihlu_tdc_offs_min[i]));
 c[6]->cd(i+11);gr_ihld_tdc_offs[i]-> Draw("AP");                                                                                                                            
                tex-> DrawLatex(rnum[0], ihld_tdc_offs_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ihld_tdc_offs_max[i],run_ihld_tdc_offs_max[i],ihld_tdc_offs_min[i],run_ihld_tdc_offs_min[i]));
 }

c[7]->Clear();
c[7]->Divide(4,5);
 for(int i=0;i<nIH;i++){
 c[7]->cd(i+1); gr_ihru_tdc_offs[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ihru_tdc_offs_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ihru_tdc_offs_max[i],run_ihru_tdc_offs_max[i],ihru_tdc_offs_min[i],run_ihru_tdc_offs_min[i]));
 c[7]->cd(i+11);gr_ihrd_tdc_offs[i]-> Draw("AP");                                                                                                                            
                tex-> DrawLatex(rnum[0], ihrd_tdc_offs_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ihrd_tdc_offs_max[i],run_ihrd_tdc_offs_max[i],ihrd_tdc_offs_min[i],run_ihrd_tdc_offs_min[i]));
 }

c[8]->Clear();
c[8]->Divide(4,6);
 for(int i=0;i<nOHV;i++){
 c[8]->cd(i+1); gr_ohvlu_tdc_offs[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ohvlu_tdc_offs_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohvlu_tdc_offs_max[i],run_ohvlu_tdc_offs_max[i],ohvlu_tdc_offs_min[i],run_ohvlu_tdc_offs_min[i]));
 c[8]->cd(i+13);gr_ohvld_tdc_offs[i]-> Draw("AP");                                                                                                                               
                tex-> DrawLatex(rnum[0], ohvld_tdc_offs_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohvld_tdc_offs_max[i],run_ohvld_tdc_offs_max[i],ohvld_tdc_offs_min[i],run_ohvld_tdc_offs_min[i]));
 }

c[9]->Clear();
c[9]->Divide(4,6);
 for(int i=0;i<nOHV;i++){
 c[9]->cd(i+1); gr_ohvru_tdc_offs[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ohvru_tdc_offs_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohvru_tdc_offs_max[i],run_ohvru_tdc_offs_max[i],ohvru_tdc_offs_min[i],run_ohvru_tdc_offs_min[i]));
 c[9]->cd(i+13);gr_ohvrd_tdc_offs[i]-> Draw("AP");                                                                                                                               
                tex-> DrawLatex(rnum[0], ohvrd_tdc_offs_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohvrd_tdc_offs_max[i],run_ohvrd_tdc_offs_max[i],ohvrd_tdc_offs_min[i],run_ohvrd_tdc_offs_min[i]));
 }

c[10]->Clear();
c[10]->Divide(4,5);
 for(int i=0;i<nOHH;i++){
 c[10]->cd(i+1); gr_ohhlu_tdc_offs[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ohhlu_tdc_offs_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohhlu_tdc_offs_max[i],run_ohhlu_tdc_offs_max[i],ohhlu_tdc_offs_min[i],run_ohhlu_tdc_offs_min[i]));
 c[10]->cd(i+13);gr_ohhld_tdc_offs[i]-> Draw("AP");                                                                                                                              
                tex-> DrawLatex(rnum[0], ohhld_tdc_offs_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohhld_tdc_offs_max[i],run_ohhld_tdc_offs_max[i],ohhld_tdc_offs_min[i],run_ohhld_tdc_offs_min[i]));
 }

c[11]->Clear();
c[11]->Divide(4,5);
 for(int i=0;i<nOHH;i++){
 c[11]->cd(i+1); gr_ohhru_tdc_offs[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], ohhru_tdc_offs_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohhru_tdc_offs_max[i],run_ohhru_tdc_offs_max[i],ohhru_tdc_offs_min[i],run_ohhru_tdc_offs_min[i]));
 c[11]->cd(i+11);gr_ohhrd_tdc_offs[i]-> Draw("AP");                                                                                                                              
                tex-> DrawLatex(rnum[0], ohhrd_tdc_offs_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",ohhrd_tdc_offs_max[i],run_ohhrd_tdc_offs_max[i],ohhrd_tdc_offs_min[i],run_ohhrd_tdc_offs_min[i]));
 }



}
////////////////////////////////////////////////////////////////////////////
void checkparam::savecanvas(){
cout<<RNumMin<<" - "<<RNumMax<<endl;
string pdf_name =Form("pdf/checkparam/run%d_%d.pdf",RNumMin, RNumMax  ); 
  c[0] ->Print(Form("%s[",pdf_name.c_str()  ) );
  c[0] ->Print(Form("%s" ,pdf_name.c_str()  ) );
  c[1] ->Print(Form("%s" ,pdf_name.c_str()  ) );
  c[2] ->Print(Form("%s" ,pdf_name.c_str()  ) );
  c[3] ->Print(Form("%s" ,pdf_name.c_str()  ) );
  c[4] ->Print(Form("%s" ,pdf_name.c_str()  ) );
  c[5] ->Print(Form("%s" ,pdf_name.c_str()  ) );
  c[6] ->Print(Form("%s" ,pdf_name.c_str()  ) );
  c[7] ->Print(Form("%s" ,pdf_name.c_str()  ) );
  c[8] ->Print(Form("%s" ,pdf_name.c_str()  ) );
  c[9] ->Print(Form("%s" ,pdf_name.c_str()  ) );
  c[10]->Print(Form("%s" ,pdf_name.c_str()  ) );
  c[11]->Print(Form("%s" ,pdf_name.c_str()  ) );
  c[11]->Print(Form("%s]",pdf_name.c_str()  ) );
  cout<<"saved : "<<pdf_name<<endl;
}
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "output.root";
  int ch;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool itr_flag = false;//iteration of PHC flag
  bool draw_flag = true;
  bool skip_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:bc:p:"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
    case 'h':
      cout<<"-f : input run list filename"<<endl;
      return 0;
      break;
    case '?':
      cout<<"unknown option...."<<endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;
    }
  }

  TApplication *theApp = new TApplication("App", &argc, argv);
  checkparam *ana = new checkparam();
  ana->set_filename(ifname);
  ana->analysis();
cout<<"after analysis"<<endl;
  ana->set_graph();
cout<<"after set graph"<<endl;
  ana->draw();
cout<<"after draw"<<endl;
  ana->savecanvas();
cout<<"after savecanvas"<<endl;
  delete ana;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}

