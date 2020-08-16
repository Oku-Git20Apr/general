#include <iostream>
#include <fstream>
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

#define Calibration

static const double PI = 4.0*atan(1.);
static const double mrad_to_deg = 1./1000*180./PI;
const double Mp = 938.272046;          // proton       mass (MeV/c2)
const double c = 0.299792458;          // speed of light in vacuum (m/ns)

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
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleSize(0.04);
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
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.5);
  h->SetMarkerColor(1);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(3);
}

void SetGrErr(TGraphErrors *gr, TString hname, TString xname, TString yname, int LColor, int MColor, int MStyle, double Yoffset, double min, double max){
  gr->SetTitle(hname);
  gr->SetName(hname);
  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->SetLineColor(LColor);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(0.5);
  gr->GetYaxis()->SetTitleOffset(Yoffset);
//  gr->GetYaxis()->SetRangeUser(min,max);
}

////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "output0000.dat";
  string ofname = "protonlike_extract.txt";
  int ch;
  int MaxNum = 0;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:w:n:bco"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;
    case 'w':
      output_flag = true;
      draw_flag = false;
      ofname = optarg;
      cout<<"output filename : "<<ofname<<endl;
      break;
    case 'n':
      MaxNum = atoi(optarg);
      break;
    case 'c':
      coin_flag = true;
      break;
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output txt filename"<<endl;
      cout<<"-n : maximum number of analysed events"<<endl;
      cout<<"-p : print png file"<<endl;
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

  TApplication theApp("App", &argc, argv);

  time_t start, end;
  start = time(NULL);
  time(&start);

  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");
  if(!draw_flag) gROOT->SetBatch(1);

  gStyle->SetOptDate(0);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetLineWidth(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptStat("ei");
//  gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleTextColor(1);

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
      
  gStyle->SetGridWidth(0);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetLineWidth(0);
  gStyle->SetNdivisions(510); // tertiary*10000 + secondary*100 + first

TFile *ifp = new TFile(Form("%s",ifname.c_str() ) );
TTree *tree = (TTree*)ifp->Get("tree");


gStyle->SetOptStat("iMen");
TH1F *h_IHRU[10];
TH1F *h_IHRD[10];
TH1F *h_IHLU[10];
TH1F *h_IHLD[10];
TH1F *h_OHVLU[12];
TH1F *h_OHVLD[12];
TH1F *h_OHVRU[12];
TH1F *h_OHVRD[12];
TH1F *h_OHHLU[9];
TH1F *h_OHHLD[9];
TH1F *h_OHHRU[9];
TH1F *h_OHHRD[9];
TH1F *h_EVSL[4];
TH1F *h_EVSR[4];
TH1F *h_EVLGL[2];
TH1F *h_EVLGR[2];

TH1F *h_IHRU_wt[10];
TH1F *h_IHRD_wt[10];
TH1F *h_IHLU_wt[10];
TH1F *h_IHLD_wt[10];
TH1F *h_OHVLU_wt[12];
TH1F *h_OHVLD_wt[12];
TH1F *h_OHVRU_wt[12];
TH1F *h_OHVRD_wt[12];
TH1F *h_OHHLU_wt[9];
TH1F *h_OHHLD_wt[9];
TH1F *h_OHHRU_wt[9];
TH1F *h_OHHRD_wt[9];
TH1F *h_EVSL_wt[4];
TH1F *h_EVSR_wt[4];
TH1F *h_EVLGL_wt[2];
TH1F *h_EVLGR_wt[2];

int IHRUmean[10],IHRDmean[10],IHLUmean[10],IHLDmean[10];
int OHVRUmean[12],OHVRDmean[12],OHVLUmean[12],OHVLDmean[12];
int OHHRUmean[9],OHHRDmean[9],OHHLUmean[9],OHHLDmean[9];
int EVSLmean[4],EVSRmean[4];
int EVLGLmean[2],EVLGRmean[2];

int IHRUrms[10],IHRDrms[10],IHLUrms[10],IHLDrms[10];
int OHVRUrms[12],OHVRDrms[12],OHVLUrms[12],OHVLDrms[12];
int OHHRUrms[9],OHHRDrms[9],OHHLUrms[9],OHHLDrms[9];
int EVSLrms[4],EVSRrms[4];
int EVLGLrms[4],EVLGRrms[4];
TLine *l_ihlu[10];
TLine *l_ihld[10];
TLine *l_ihru[10];
TLine *l_ihrd[10];
TLine *l_ohvlu[12];
TLine *l_ohvld[12];
TLine *l_ohvru[12];
TLine *l_ohvrd[12];
TLine *l_ohhlu[10];
TLine *l_ohhld[9];
TLine *l_ohhru[9];
TLine *l_ohhrd[9];
TLine *l_evl[4];
TLine *l_evr[4];
TLine *l_evlgl[2];
TLine *l_evlgr[2];

TText *t_ihlu[10];
TText *t_ihld[10];
TText *t_ihru[10];
TText *t_ihrd[10];
TText *t_ohvlu[12];
TText *t_ohvld[12];
TText *t_ohvru[12];
TText *t_ohvrd[12];
TText *t_ohhlu[10];
TText *t_ohhld[9];
TText *t_ohhru[9];
TText *t_ohhrd[9];
TText *t_evl[4];
TText *t_evr[4];
TText *t_evlgl[2];
TText *t_evlgr[2];
//////////////
//histograms//
//////////////
for(int i = 0; i<10 ; i++){//IH
h_IHLU[i]=new TH1F(Form("h_IHLU%d",i+1),Form("Pedestal IHLU %d",i+1),500,0,500);
h_IHLD[i]=new TH1F(Form("h_IHLD%d",i+1),Form("Pedestal IHLD %d",i+1),500,0,500);
h_IHRU[i]=new TH1F(Form("h_IHRU%d",i+1),Form("Pedestal IHRU %d",i+1),500,0,500);
h_IHRD[i]=new TH1F(Form("h_IHRD%d",i+1),Form("Pedestal IHRD %d",i+1),500,0,500);
SetTH1(h_IHLU[i], Form("Pedestal IHLU %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(h_IHLD[i], Form("Pedestal IHLD %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(h_IHRU[i], Form("Pedestal IHRU %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(h_IHRD[i], Form("Pedestal IHRD %d",i+1), "ADC[ch]", "count", 1, 0, 0);

h_IHLU_wt[i]=new TH1F(Form("h_IHLU_wt%d",i+1),Form("IHLU %d w/ TDC",i+1),500,0,500);
h_IHLD_wt[i]=new TH1F(Form("h_IHLD_wt%d",i+1),Form("IHLD %d w/ TDC",i+1),500,0,500);
h_IHRU_wt[i]=new TH1F(Form("h_IHRU_wt%d",i+1),Form("IHRU %d w/ TDC",i+1),500,0,500);
h_IHRD_wt[i]=new TH1F(Form("h_IHRD_wt%d",i+1),Form("IHRD %d w/ TDC",i+1),500,0,500);
SetTH1(h_IHLU_wt[i], Form("IHLU %d w/ TDC",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(h_IHLD_wt[i], Form("IHLD %d w/ TDC",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(h_IHRU_wt[i], Form("IHRU %d w/ TDC",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(h_IHRD_wt[i], Form("IHRD %d w/ TDC",i+1), "ADC[ch]", "count", 2, 0, 0);

tree->Project(Form("h_IHLU%d",i+1),Form("ihluadc[%d]",i),Form("ihlutdc[%d]<=0",i));
tree->Project(Form("h_IHLD%d",i+1),Form("ihldadc[%d]",i),Form("ihldtdc[%d]<=0",i));
tree->Project(Form("h_IHRU%d",i+1),Form("ihruadc[%d]",i),Form("ihrutdc[%d]<=0",i));
tree->Project(Form("h_IHRD%d",i+1),Form("ihrdadc[%d]",i),Form("ihrdtdc[%d]<=0",i));

tree->Project(Form("h_IHLU_wt%d",i+1),Form("ihluadc[%d]",i),Form("ihlutdc[%d]>0",i));
tree->Project(Form("h_IHLD_wt%d",i+1),Form("ihldadc[%d]",i),Form("ihldtdc[%d]>0",i));
tree->Project(Form("h_IHRU_wt%d",i+1),Form("ihruadc[%d]",i),Form("ihrutdc[%d]>0",i));
tree->Project(Form("h_IHRD_wt%d",i+1),Form("ihrdadc[%d]",i),Form("ihrdtdc[%d]>0",i));
 }
cout<<"IH filled"<<endl;
for(int i = 0; i<12 ; i++){//OHV
h_OHVLU[i]=new TH1F(Form("h_OHVLU%d",i+1),Form("Pedestal OHVLU %d",i+1),1000,0,1000);
h_OHVLD[i]=new TH1F(Form("h_OHVLD%d",i+1),Form("Pedestal OHVLD %d",i+1),1000,0,1000);
h_OHVRU[i]=new TH1F(Form("h_OHVRU%d",i+1),Form("Pedestal OHVRU %d",i+1),1000,0,1000);
h_OHVRD[i]=new TH1F(Form("h_OHVRD%d",i+1),Form("Pedestal OHVRD %d",i+1),1000,0,1000);
SetTH1(h_OHVLU[i], Form("Pedestal OHVLU %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(h_OHVLD[i], Form("Pedestal OHVLD %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(h_OHVRU[i], Form("Pedestal OHVRU %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(h_OHVRD[i], Form("Pedestal OHVRD %d",i+1), "ADC[ch]", "count", 1, 0, 0);

h_OHVLU_wt[i]=new TH1F(Form("h_OHVLU_wt%d",i+1),Form("OHVLU %d w/ TDC",i+1),1000,0,1000);
h_OHVLD_wt[i]=new TH1F(Form("h_OHVLD_wt%d",i+1),Form("OHVLD %d w/ TDC",i+1),1000,0,1000);
h_OHVRU_wt[i]=new TH1F(Form("h_OHVRU_wt%d",i+1),Form("OHVRU %d w/ TDC",i+1),1000,0,1000);
h_OHVRD_wt[i]=new TH1F(Form("h_OHVRD_wt%d",i+1),Form("OHVRD %d w/ TDC",i+1),1000,0,1000);
SetTH1(h_OHVLU_wt[i], Form("OHVLU %d",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(h_OHVLD_wt[i], Form("OHVLD %d",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(h_OHVRU_wt[i], Form("OHVRU %d",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(h_OHVRD_wt[i], Form("OHVRD %d",i+1), "ADC[ch]", "count", 2, 0, 0);

tree->Project(Form("h_OHVLU%d",i+1),Form("ohvluadc[%d]",i),Form("ohvlutdc[%d]<=0",i));
tree->Project(Form("h_OHVLD%d",i+1),Form("ohvldadc[%d]",i),Form("ohvldtdc[%d]<=0",i));
tree->Project(Form("h_OHVRU%d",i+1),Form("ohvruadc[%d]",i),Form("ohvrutdc[%d]<=0",i));
tree->Project(Form("h_OHVRD%d",i+1),Form("ohvrdadc[%d]",i),Form("ohvrdtdc[%d]<=0",i));

tree->Project(Form("h_OHVLU_wt%d",i+1),Form("ohvluadc[%d]",i),Form("ohvlutdc[%d]>0",i));
tree->Project(Form("h_OHVLD_wt%d",i+1),Form("ohvldadc[%d]",i),Form("ohvldtdc[%d]>0",i));
tree->Project(Form("h_OHVRU_wt%d",i+1),Form("ohvruadc[%d]",i),Form("ohvrutdc[%d]>0",i));
tree->Project(Form("h_OHVRD_wt%d",i+1),Form("ohvrdadc[%d]",i),Form("ohvrdtdc[%d]>0",i));
}
cout<<"OHV filled"<<endl;
for(int i = 0; i<9 ; i++){//OHH
h_OHHLU[i]=new TH1F(Form("h_OHHLU%d",i+1),Form("Pedestal OHHLU %d",i+1),1000,0,1000);
h_OHHLD[i]=new TH1F(Form("h_OHHLD%d",i+1),Form("Pedestal OHHLD %d",i+1),1000,0,1000);
h_OHHRU[i]=new TH1F(Form("h_OHHRU%d",i+1),Form("Pedestal OHHRU %d",i+1),1000,0,1000);
h_OHHRD[i]=new TH1F(Form("h_OHHRD%d",i+1),Form("Pedestal OHHRD %d",i+1),1000,0,1000);
SetTH1(h_OHHLU[i], Form("OHHLU %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(h_OHHLD[i], Form("OHHLD %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(h_OHHRU[i], Form("OHHRU %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(h_OHHRD[i], Form("OHHRD %d",i+1), "ADC[ch]", "count", 1, 0, 0);

h_OHHLU_wt[i]=new TH1F(Form("h_OHHLU_wt%d",i+1),Form("OHHLU %d w/ TDC",i+1),1000,0,1000);
h_OHHLD_wt[i]=new TH1F(Form("h_OHHLD_wt%d",i+1),Form("OHHLD %d w/ TDC",i+1),1000,0,1000);
h_OHHRU_wt[i]=new TH1F(Form("h_OHHRU_wt%d",i+1),Form("OHHRU %d w/ TDC",i+1),1000,0,1000);
h_OHHRD_wt[i]=new TH1F(Form("h_OHHRD_wt%d",i+1),Form("OHHRD %d w/ TDC",i+1),1000,0,1000);
SetTH1(h_OHHLU_wt[i], Form("OHHLU %d",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(h_OHHLD_wt[i], Form("OHHLD %d",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(h_OHHRU_wt[i], Form("OHHRU %d",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(h_OHHRD_wt[i], Form("OHHRD %d",i+1), "ADC[ch]", "count", 2, 0, 0);

tree->Project(Form("h_OHHLU%d",i+1),Form("ohhluadc[%d]",i),Form("ohhlutdc[%d]<=0",i));
tree->Project(Form("h_OHHLD%d",i+1),Form("ohhldadc[%d]",i),Form("ohhldtdc[%d]<=0",i));
tree->Project(Form("h_OHHRU%d",i+1),Form("ohhruadc[%d]",i),Form("ohhrutdc[%d]<=0",i));
tree->Project(Form("h_OHHRD%d",i+1),Form("ohhrdadc[%d]",i),Form("ohhrdtdc[%d]<=0",i));

tree->Project(Form("h_OHHLU_wt%d",i+1),Form("ohhluadc[%d]",i),Form("ohhlutdc[%d]>0",i));
tree->Project(Form("h_OHHLD_wt%d",i+1),Form("ohhldadc[%d]",i),Form("ohhldtdc[%d]>0",i));
tree->Project(Form("h_OHHRU_wt%d",i+1),Form("ohhruadc[%d]",i),Form("ohhrutdc[%d]>0",i));
tree->Project(Form("h_OHHRD_wt%d",i+1),Form("ohhrdadc[%d]",i),Form("ohhrdtdc[%d]>0",i));
}
cout<<"OHH filled"<<endl;

for(int i = 0; i<4 ; i++){//EVS
h_EVSL[i]=new TH1F(Form("h_EVSL%d",i+1),Form("Pedestal EVSL %d",i+1),1000,0,1000);
h_EVSR[i]=new TH1F(Form("h_EVSR%d",i+1),Form("Pedestal EVSR %d",i+1),1000,0,1000);
SetTH1(h_EVSL[i], Form("EVSL %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(h_EVSR[i], Form("EVSR %d",i+1), "ADC[ch]", "count", 1, 0, 0);

h_EVSL_wt[i]=new TH1F(Form("h_EVSL_wt%d",i+1),Form("EVSL %d",i+1),1000,0,1000);
h_EVSR_wt[i]=new TH1F(Form("h_EVSR_wt%d",i+1),Form("EVSR %d",i+1),1000,0,1000);
SetTH1(h_EVSL_wt[i], Form("EVSL %d",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(h_EVSR_wt[i], Form("EVSR %d",i+1), "ADC[ch]", "count", 2, 0, 0);

tree->Project(Form("h_EVSL%d",i+1),Form("evladc[%d]",i),Form("evltdc[%d]<=0",i));
tree->Project(Form("h_EVSR%d",i+1),Form("evradc[%d]",i),Form("evrtdc[%d]<=0",i));
tree->Project(Form("h_EVSL_wt%d",i+1),Form("evladc[%d]",i),Form("evltdc[%d]>0",i));
tree->Project(Form("h_EVSR_wt%d",i+1),Form("evradc[%d]",i),Form("evrtdc[%d]>0",i));
}

cout<<"EV filled"<<endl;
for(int i = 0; i<2 ; i++){//EVLG
h_EVLGL[i]=new TH1F(Form("h_EVLGL%d",i+1),Form("Pedestal EVLGL %d",i+1),1000,0,1000);
h_EVLGR[i]=new TH1F(Form("h_EVLGR%d",i+1),Form("Pedestal EVLGR %d",i+1),1000,0,1000);
SetTH1(h_EVLGL[i], Form("EVLGL %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(h_EVLGR[i], Form("EVLGR %d",i+1), "ADC[ch]", "count", 1, 0, 0);

h_EVLGL_wt[i]=new TH1F(Form("h_EVLGL_wt%d",i+1),Form("EVLGL %d",i+1),1000,0,1000);
h_EVLGR_wt[i]=new TH1F(Form("h_EVLGR_wt%d",i+1),Form("EVLGR %d",i+1),1000,0,1000);
SetTH1(h_EVLGL_wt[i], Form("EVLGL %d",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(h_EVLGR_wt[i], Form("EVLGR %d",i+1), "ADC[ch]", "count", 2, 0, 0);

tree->Project(Form("h_EVLGL%d",i+1),Form("evlgladc[%d]",i),Form("evlgltdc[%d]<=0",i));
tree->Project(Form("h_EVLGR%d",i+1),Form("evlgradc[%d]",i),Form("evlgrtdc[%d]<=0",i));
tree->Project(Form("h_EVLGL_wt%d",i+1),Form("evlgladc[%d]",i),Form("evlgltdc[%d]>0",i));
tree->Project(Form("h_EVLGR_wt%d",i+1),Form("evlgradc[%d]",i),Form("evlgrtdc[%d]>0",i));
}
cout<<"EVLG filled"<<endl;
  std::ofstream ofs(ofname.c_str());
for(int i = 0; i<10 ; i++){//IH
IHLUmean[i] = h_IHLU[i] -> GetMaximumBin();IHLUrms[i] = h_IHLU[i] -> GetRMS();
//cout<<Form("IHLU%d   ",i+1)<<IHLUmean[i] + 2*IHLUrms[i]<<endl;
ofs<< Form("IHLU%d   ",i+1)<<IHLUmean[i] + 2*IHLUrms[i]<<endl;
     l_ihlu[i]= new TLine(IHLUmean[i] + 2*IHLUrms[i],0,IHLUmean[i] + 2*IHLUrms[i],1000);
     l_ihlu[i]->SetLineColor(4); l_ihlu[i]->SetLineWidth(1);

     t_ihlu[i]= new TText( IHLUmean[i] + 2*IHLUrms[i],1000,Form("%.0d ch",IHLUmean[i] + 2*IHLUrms[i]));
     t_ihlu[i]->SetTextColor(6); t_ihlu[i]->SetTextSize(0.11);
}
for(int i = 0; i<10 ; i++){//IH
IHLDmean[i] = h_IHLD[i] -> GetMaximumBin();IHLDrms[i] = h_IHLD[i] -> GetRMS();
ofs<<Form("IHLD%d   ",i+1)<<IHLDmean[i] + 2*IHLDrms[i]<<endl;
     l_ihld[i]= new TLine(IHLDmean[i] + 2*IHLDrms[i],0,IHLDmean[i] + 2*IHLDrms[i],1000);
     l_ihld[i]->SetLineColor(4); l_ihld[i]->SetLineWidth(1);

     t_ihld[i]= new TText( IHLDmean[i] + 2*IHLDrms[i],1000,Form("%.0d ch",IHLDmean[i] + 2*IHLDrms[i]));
     t_ihld[i]->SetTextColor(6); t_ihld[i]->SetTextSize(0.11);
}
for(int i = 0; i<10 ; i++){//IH
IHRUmean[i] = h_IHRU[i] -> GetMaximumBin();IHRUrms[i] = h_IHRU[i] -> GetRMS();
ofs<<Form("IHRU%d   ",i+1)<<IHRUmean[i] + 2*IHRUrms[i]<<endl;
     l_ihru[i]= new TLine(IHRUmean[i] + 2*IHRUrms[i],0,IHRUmean[i] + 2*IHRUrms[i],1000);
     l_ihru[i]->SetLineColor(4); l_ihru[i]->SetLineWidth(1);

     t_ihru[i]= new TText( IHRUmean[i] + 2*IHRUrms[i],1000,Form("%.0d ch",IHRUmean[i] + 2*IHRUrms[i]));
     t_ihru[i]->SetTextColor(6); t_ihru[i]->SetTextSize(0.11);
}
for(int i = 0; i<10 ; i++){//IH
IHRDmean[i] = h_IHRD[i] -> GetMaximumBin();IHRDrms[i] = h_IHRD[i] -> GetRMS();
ofs<<Form("IHRD%d   ",i+1)<<IHRDmean[i] + 2*IHRDrms[i]<<endl;
     l_ihrd[i]= new TLine(IHRDmean[i] + 2*IHRDrms[i],0,IHRDmean[i] + 2*IHRDrms[i],1000);
     l_ihrd[i]->SetLineColor(4); l_ihrd[i]->SetLineWidth(1);

     t_ihrd[i]= new TText( IHRDmean[i] + 2*IHRDrms[i],1000,Form("%.0d ch",IHRDmean[i] + 2*IHRDrms[i]));
     t_ihrd[i]->SetTextColor(6); t_ihrd[i]->SetTextSize(0.11);
}
for(int i = 0; i<12 ; i++){//OHV
OHVLUmean[i] = h_OHVLU[i] -> GetMaximumBin();OHVLUrms[i] = h_OHVLU[i] -> GetRMS();
ofs<<Form("OHVLU%d   ",i+1)<<OHVLUmean[i] + 2*OHVLUrms[i]<<endl;
     l_ohvlu[i]= new TLine(OHVLUmean[i] + 2*OHVLUrms[i],0,OHVLUmean[i] + 2*OHVLUrms[i],1000);
     l_ohvlu[i]->SetLineColor(4); l_ohvlu[i]->SetLineWidth(1);

     t_ohvlu[i]= new TText( OHVLUmean[i] + 2*OHVLUrms[i],1000,Form("%.0d ch",OHVLUmean[i] + 2*OHVLUrms[i]));
     t_ohvlu[i]->SetTextColor(6); t_ohvlu[i]->SetTextSize(0.11);
}
for(int i = 0; i<12 ; i++){//OHV
OHVLDmean[i] = h_OHVLD[i] -> GetMaximumBin();OHVLDrms[i] = h_OHVLD[i] -> GetRMS();
ofs<<Form("OHVLD%d   ",i+1)<<OHVLDmean[i] + 2*OHVLDrms[i]<<endl;
     l_ohvld[i]= new TLine(OHVLDmean[i] + 2*OHVLDrms[i],0,OHVLDmean[i] + 2*OHVLDrms[i],1000);
     l_ohvld[i]->SetLineColor(4); l_ohvld[i]->SetLineWidth(1);

     t_ohvld[i]= new TText( OHVLDmean[i] + 2*OHVLDrms[i],1000,Form("%.0d ch",OHVLDmean[i] + 2*OHVLDrms[i]));
     t_ohvld[i]->SetTextColor(6); t_ohvld[i]->SetTextSize(0.11);
}
for(int i = 0; i<12 ; i++){//OHV
OHVRUmean[i] = h_OHVRU[i] -> GetMaximumBin();OHVRUrms[i] = h_OHVRU[i] -> GetRMS();
ofs<<Form("OHVRU%d   ",i+1)<<OHVRUmean[i] + 2*OHVRUrms[i]<<endl;
     l_ohvru[i]= new TLine(OHVRUmean[i] + 2*OHVRUrms[i],0,OHVRUmean[i] + 2*OHVRUrms[i],1000);
     l_ohvru[i]->SetLineColor(4); l_ohvru[i]->SetLineWidth(1);

     t_ohvru[i]= new TText( OHVRUmean[i] + 2*OHVRUrms[i],1000,Form("%.0d ch",OHVLUmean[i] + 2*OHVLUrms[i]));
     t_ohvru[i]->SetTextColor(6); t_ohvru[i]->SetTextSize(0.11);
}
for(int i = 0; i<12 ; i++){//OHV
OHVRDmean[i] = h_OHVRD[i] -> GetMaximumBin();OHVRDrms[i] = h_OHVRD[i] -> GetRMS();
ofs<<Form("OHVRD%d   ",i+1)<<OHVRDmean[i] + 2*OHVRDrms[i]<<endl;
     l_ohvrd[i]= new TLine(OHVRDmean[i] + 2*OHVRDrms[i],0,OHVRDmean[i] + 2*OHVRDrms[i],1000);
     l_ohvrd[i]->SetLineColor(4); l_ohvrd[i]->SetLineWidth(1);

     t_ohvrd[i]= new TText( OHVRDmean[i] + 2*OHVRDrms[i],1000,Form("%.0d ch",OHVLDmean[i] + 2*OHVLDrms[i]));
     t_ohvrd[i]->SetTextColor(6); t_ohvrd[i]->SetTextSize(0.11);
}
for(int i = 0; i<9 ; i++){//OHH
OHHLUmean[i] = h_OHHLU[i] -> GetMaximumBin();OHHLUrms[i] = h_OHHLU[i] -> GetRMS();
ofs<<Form("OHHLU%d   ",i+1)<<OHHLUmean[i] + 2*OHHLUrms[i]<<endl;
     l_ohhlu[i]= new TLine(OHHLUmean[i] + 2*OHHLUrms[i],0,OHHLUmean[i] + 2*OHHLUrms[i],1000);
     l_ohhlu[i]->SetLineColor(4); l_ohhlu[i]->SetLineWidth(1);

     t_ohhlu[i]= new TText( OHHLUmean[i] + 2*OHHLUrms[i],1000,Form("%.0d ch",OHHLUmean[i] + 2*OHHLUrms[i]));
     t_ohhlu[i]->SetTextColor(6); t_ohhlu[i]->SetTextSize(0.11);
}
for(int i = 0; i<9 ; i++){//OHH
OHHLDmean[i] = h_OHHLD[i] -> GetMaximumBin();OHHLDrms[i] = h_OHHLD[i] -> GetRMS();
ofs<<Form("OHHLD%d   ",i+1)<<OHHLDmean[i] + 2*OHHLDrms[i]<<endl;
     l_ohhld[i]= new TLine(OHHLDmean[i] + 2*OHHLDrms[i],0,OHHLDmean[i] + 2*OHHLDrms[i],1000);
     l_ohhld[i]->SetLineColor(4); l_ohhld[i]->SetLineWidth(1);

     t_ohhld[i]= new TText( OHHLDmean[i] + 2*OHHLDrms[i],1000,Form("%.0d ch",OHHLDmean[i] + 2*OHHLDrms[i]));
     t_ohhld[i]->SetTextColor(6); t_ohhld[i]->SetTextSize(0.11);
}
for(int i = 0; i<9 ; i++){//OHH
OHHRUmean[i] = h_OHHRU[i] -> GetMaximumBin();OHHRUrms[i] = h_OHHRU[i] -> GetRMS();
ofs<<Form("OHHRU%d   ",i+1)<<OHHRUmean[i] + 2*OHHRUrms[i]<<endl;
     l_ohhru[i]= new TLine(OHHRUmean[i] + 2*OHHRUrms[i],0,OHHRUmean[i] + 2*OHHRUrms[i],1000);
     l_ohhru[i]->SetLineColor(4); l_ohhru[i]->SetLineWidth(1);

     t_ohhru[i]= new TText( OHHRUmean[i] + 2*OHHRUrms[i],1000,Form("%.0d ch",OHHRUmean[i] + 2*OHHRUrms[i]));
     t_ohhru[i]->SetTextColor(6); t_ohhru[i]->SetTextSize(0.11);
}
for(int i = 0; i<9 ; i++){//OHH
OHHRDmean[i] = h_OHHRD[i] -> GetMaximumBin();OHHRDrms[i] = h_OHHRD[i] -> GetRMS();
ofs<<Form("OHHRD%d   ",i+1)<<OHHRDmean[i] + 2*OHHRDrms[i]<<endl;
     l_ohhrd[i]= new TLine(OHHRDmean[i] + 2*OHHRDrms[i],0,OHHRDmean[i] + 2*OHHRDrms[i],1000);
     l_ohhrd[i]->SetLineColor(4); l_ohhrd[i]->SetLineWidth(1);

     t_ohhrd[i]= new TText( OHHRDmean[i] + 2*OHHRDrms[i],1000,Form("%.0d ch",OHHRDmean[i] + 2*OHHRDrms[i]));
     t_ohhrd[i]->SetTextColor(6); t_ohhrd[i]->SetTextSize(0.11);
}
for(int i = 0; i<4 ; i++){//EVS
EVSLmean[i] = h_EVSL[i] -> GetMaximumBin();EVSLrms[i] = h_EVSL[i] -> GetRMS();
ofs<<Form("EVSL%d   ",i+1)<<EVSLmean[i] + 2*EVSLrms[i]<<endl;
     l_evl[i]= new TLine(EVSLmean[i] + 2*EVSLrms[i],0,EVSLmean[i] + 2*EVSLrms[i],1000);
     l_evl[i]->SetLineColor(4); l_evl[i]->SetLineWidth(1);

     t_evl[i]= new TText( EVSLmean[i] + 2*EVSLrms[i],1000,Form("%.0d ch",EVSLmean[i] + 2*EVSLrms[i]));
     t_evl[i]->SetTextColor(6); t_evl[i]->SetTextSize(0.11);
}
for(int i = 0; i<4 ; i++){//EVS
EVSRmean[i] = h_EVSR[i] -> GetMaximumBin();EVSRrms[i] = h_EVSR[i] -> GetRMS();
ofs<<Form("EVSR%d   ",i+1)<<EVSRmean[i] + 2*EVSRrms[i]<<endl;
     l_evr[i]= new TLine(EVSRmean[i] + 2*EVSRrms[i],0,EVSRmean[i] + 2*EVSRrms[i],1000);
     l_evr[i]->SetLineColor(4); l_evr[i]->SetLineWidth(1);

     t_evr[i]= new TText( EVSRmean[i] + 2*EVSRrms[i],1000,Form("%.0d ch",EVSRmean[i] + 2*EVSRrms[i]));
     t_evr[i]->SetTextColor(6); t_evr[i]->SetTextSize(0.11);
}

for(int i = 0; i<2 ; i++){//EVLG
EVLGRmean[i] = h_EVLGR[i] -> GetMaximumBin();EVLGLrms[i] = h_EVLGL[i] -> GetRMS();
ofs<<Form("EVLGL%d   ",i+1)<<EVLGLmean[i] + 2*EVLGLrms[i]<<endl;
     l_evlgl[i]= new TLine(EVLGLmean[i] + 2*EVLGLrms[i],0,EVLGLmean[i] + 2*EVLGLrms[i],1000);
     l_evlgl[i]->SetLineColor(4); l_evlgl[i]->SetLineWidth(1);

     t_evlgl[i]= new TText( EVLGLmean[i] + 2*EVLGLrms[i],1000,Form("%.0d ch",EVLGLmean[i] + 2*EVLGLrms[i]));
     t_evlgl[i]->SetTextColor(6); t_evlgl[i]->SetTextSize(0.11);
}

for(int i = 0; i<2 ; i++){//EVLG
EVLGRmean[i] = h_EVLGR[i] -> GetMaximumBin();EVLGRrms[i] = h_EVLGR[i] -> GetRMS();
ofs<<Form("EVLGR%d   ",i+1)<<EVLGRmean[i] + 2*EVLGRrms[i]<<endl;
     l_evlgr[i]= new TLine(EVLGRmean[i] + 2*EVLGRrms[i],0,EVLGRmean[i] + 2*EVLGRrms[i],1000);
     l_evlgr[i]->SetLineColor(4); l_evlgr[i]->SetLineWidth(1);

     t_evlgr[i]= new TText( EVLGRmean[i] + 2*EVLGRrms[i],1000,Form("%.0d ch",EVLGRmean[i] + 2*EVLGRrms[i]));
     t_evlgr[i]->SetTextColor(6); t_evlgr[i]->SetTextSize(0.11);
}

cout<<"start drawing!"<<endl;
TCanvas *c[14];
for(int i=0;i<14;i++){
c[i] = new TCanvas(Form("c%d",i+1), Form("canvas%d",i+1) , 800,1500);}

c[0]->Divide(2,5);
for(int i=0;i<10;i++){
c[0]->cd(i+1); gPad->SetLogy(1); h_IHLU[i]  ->Draw();h_IHLU_wt[i]  ->Draw("same");l_ihlu[i]->Draw("same");t_ihlu[i]->Draw("same");
 }                                                                                                                                
c[1]->Divide(2,5);                                                                                                                
for(int i=0;i<10;i++){                                                                                                            
c[1]->cd(i+1); gPad->SetLogy(1); h_IHLD[i]  ->Draw();h_IHLD_wt[i]  ->Draw("same");l_ihld[i]->Draw("same");t_ihld[i]->Draw("same");
 }                                                                       
c[2]->Divide(2,5);                                                       
for(int i=0;i<10;i++){                                                   
c[2]->cd(i+1); gPad->SetLogy(1); h_IHRU[i]  ->Draw();h_IHRU_wt[i]  ->Draw("same");l_ihru[i]->Draw("same");t_ihru[i]->Draw("same");
 }                                                                                                                                
c[3]->Divide(2,5);                                                                                                                
for(int i=0;i<10;i++){                                                                                                            
c[3]->cd(i+1); gPad->SetLogy(1); h_IHRD[i]  ->Draw();h_IHRD_wt[i]  ->Draw("same");l_ihrd[i]->Draw("same");t_ihrd[i]->Draw("same");
 }
c[4]->Divide(2,6);
for(int i=0;i<12;i++){
c[4]->cd(i+1); gPad->SetLogy(1); h_OHVLU[i]  ->Draw();h_OHVLU_wt[i]  ->Draw("same");l_ohvlu[i]->Draw("same");t_ohvlu[i]->Draw("same");
 }                                                                                                                                    
c[5]->Divide(2,6);                                                                                                                    
for(int i=0;i<12;i++){                                                                                                                
c[5]->cd(i+1); gPad->SetLogy(1); h_OHVLD[i]  ->Draw();h_OHVLD_wt[i]  ->Draw("same");l_ohvld[i]->Draw("same");t_ohvld[i]->Draw("same");
 }                                                                                                                                    
c[6]->Divide(2,6);                                                                                                                    
for(int i=0;i<12;i++){                                                                                                                
c[6]->cd(i+1); gPad->SetLogy(1); h_OHVRU[i]  ->Draw();h_OHVRU_wt[i]  ->Draw("same");l_ohvru[i]->Draw("same");t_ohvru[i]->Draw("same");
 }                                                                                                                                    
cout<<"drawing"<<endl;                                                                                                                
c[7]->Divide(2,6);                                                                                                                    
for(int i=0;i<12;i++){                                                                                                                
c[7]->cd(i+1); gPad->SetLogy(1); h_OHVRD[i]  ->Draw();h_OHVRD_wt[i]  ->Draw("same");l_ohvrd[i]->Draw("same");t_ohvrd[i]->Draw("same");
 }                                                                                                                                    
c[8]->Divide(2,5);                                                                                                                    
for(int i=0;i<9;i++){                                                                                                                 
c[8]->cd(i+1); gPad->SetLogy(1); h_OHHLU[i]  ->Draw();h_OHHLU_wt[i]  ->Draw("same");l_ohhlu[i]->Draw("same");t_ohhlu[i]->Draw("same");
 }                                                                                                                                    
cout<<"drawing OHH"<<endl;                                                                                                            
c[9]->Divide(2,5);                                                                                                                    
for(int i=0;i<9;i++){                                                                                                                 
c[9]->cd(i+1); gPad->SetLogy(1); h_OHHLD[i]  ->Draw();h_OHHLD_wt[i]  ->Draw("same");l_ohhld[i]->Draw("same");t_ohhld[i]->Draw("same");
 }
c[10]->Divide(2,5);
for(int i=0;i<9;i++){
c[10]->cd(i+1); gPad->SetLogy(1); h_OHHRU[i]  ->Draw();h_OHHRU_wt[i]  ->Draw("same");l_ohhru[i]->Draw("same");t_ohhru[i]->Draw("same");
 }                                                                                                                                     
c[11]->Divide(2,5);                                                                                                                    
for(int i=0;i<9;i++){                                                                                                                  
c[11]->cd(i+1); gPad->SetLogy(1); h_OHHRD[i]  ->Draw();h_OHHRD_wt[i]  ->Draw("same");l_ohhrd[i]->Draw("same");t_ohhrd[i]->Draw("same");
 }
cout<<"drawing EVS"<<endl;
c[12]->Divide(2,4);
for(int i=0;i<4;i++){
c[12]->cd(i+1); gPad->SetLogy(1); h_EVSL[i]  ->Draw();h_EVSL_wt[i]  ->Draw("same");l_evl[i]->Draw("same");t_evl[i]->Draw("same");
c[12]->cd(i+5); gPad->SetLogy(1); h_EVSR[i]  ->Draw();h_EVSR_wt[i]  ->Draw("same");l_evr[i]->Draw("same");t_evr[i]->Draw("same");
 }
cout<<"drawing EVLG"<<endl;
c[13]->Divide(2,4);
for(int i=0;i<2;i++){
c[13]->cd(i+1); gPad->SetLogy(1); h_EVLGL[i]  ->Draw();h_EVLGL_wt[i]  ->Draw("same");l_evlgl[i]->Draw("same");t_evlgl[i]->Draw("same");
c[13]->cd(i+3); gPad->SetLogy(1); h_EVLGR[i]  ->Draw();h_EVLGR_wt[i]  ->Draw("same");l_evlgr[i]->Draw("same");t_evlgr[i]->Draw("same");
 }
                                                                                                                                 
  string ofname_pdf = ofname;
  ofname_pdf.erase(ofname_pdf.size()-4);
  ofname_pdf.append(".pdf");

  c[0] ->Print(Form("%s[",ofname_pdf.c_str()));
  c[0] ->Print(Form("%s" ,ofname_pdf.c_str()));
  c[1] ->Print(Form("%s" ,ofname_pdf.c_str()));
  c[2] ->Print(Form("%s" ,ofname_pdf.c_str()));
  c[3] ->Print(Form("%s" ,ofname_pdf.c_str()));
  c[4] ->Print(Form("%s" ,ofname_pdf.c_str()));
  c[5] ->Print(Form("%s" ,ofname_pdf.c_str()));
  c[6] ->Print(Form("%s" ,ofname_pdf.c_str()));
  c[7] ->Print(Form("%s" ,ofname_pdf.c_str()));
  c[8] ->Print(Form("%s" ,ofname_pdf.c_str()));
  c[9] ->Print(Form("%s" ,ofname_pdf.c_str()));
  c[10]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[11]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[12]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[13]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[13]->Print(Form("%s]",ofname_pdf.c_str()));


}

