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
  string ofname = "protonlike_extract.pdf";
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
TH1F *hQ_IHRU[10];
TH1F *hQ_IHRD[10];
TH1F *hQ_IHLU[10];
TH1F *hQ_IHLD[10];
TH1F *hQ_OHVLU[12];
TH1F *hQ_OHVLD[12];
TH1F *hQ_OHVRU[12];
TH1F *hQ_OHVRD[12];
TH1F *hQ_OHHLU[9];
TH1F *hQ_OHHLD[9];
TH1F *hQ_OHHRU[9];
TH1F *hQ_OHHRD[9];
TH1F *hQ_EVSL[4];
TH1F *hQ_EVSR[4];
TH1F *hQ_EVLGL[2];
TH1F *hQ_EVLGR[2];

TH1F *hT_IHRU[10];
TH1F *hT_IHRD[10];
TH1F *hT_IHLU[10];
TH1F *hT_IHLD[10];
TH1F *hT_OHVLU[12];
TH1F *hT_OHVLD[12];
TH1F *hT_OHVRU[12];
TH1F *hT_OHVRD[12];
TH1F *hT_OHHLU[9];
TH1F *hT_OHHLD[9];
TH1F *hT_OHHRU[9];
TH1F *hT_OHHRD[9];
TH1F *hT_EVSL[4];
TH1F *hT_EVSR[4];
TH1F *hT_EVLGL[2];
TH1F *hT_EVLGR[2];

TH1F *hQ_IHRU_wt[10];
TH1F *hQ_IHRD_wt[10];
TH1F *hQ_IHLU_wt[10];
TH1F *hQ_IHLD_wt[10];
TH1F *hQ_OHVLU_wt[12];
TH1F *hQ_OHVLD_wt[12];
TH1F *hQ_OHVRU_wt[12];
TH1F *hQ_OHVRD_wt[12];
TH1F *hQ_OHHLU_wt[9];
TH1F *hQ_OHHLD_wt[9];
TH1F *hQ_OHHRU_wt[9];
TH1F *hQ_OHHRD_wt[9];
TH1F *hQ_EVSL_wt[4];
TH1F *hQ_EVSR_wt[4];
TH1F *hQ_EVLGL_wt[2];
TH1F *hQ_EVLGR_wt[2];

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
hQ_IHLU[i]=new TH1F(Form("hQ_IHLU%d",i+1),Form("QDC IHLU %d",i+1),500,0,1500);
hQ_IHLD[i]=new TH1F(Form("hQ_IHLD%d",i+1),Form("QDC IHLD %d",i+1),500,0,1500);
hQ_IHRU[i]=new TH1F(Form("hQ_IHRU%d",i+1),Form("QDC IHRU %d",i+1),500,0,1500);
hQ_IHRD[i]=new TH1F(Form("hQ_IHRD%d",i+1),Form("QDC IHRD %d",i+1),500,0,1500);
SetTH1(hQ_IHLU[i], Form("QDC IHLU %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(hQ_IHLD[i], Form("QDC IHLD %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(hQ_IHRU[i], Form("QDC IHRU %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(hQ_IHRD[i], Form("QDC IHRD %d",i+1), "ADC[ch]", "count", 1, 0, 0);

hT_IHLU[i]=new TH1F(Form("hT_IHLU%d",i+1),Form("QDC IHLU %d",i+1),3500,0,3500);
hT_IHLD[i]=new TH1F(Form("hT_IHLD%d",i+1),Form("QDC IHLD %d",i+1),3500,0,3500);
hT_IHRU[i]=new TH1F(Form("hT_IHRU%d",i+1),Form("QDC IHRU %d",i+1),3500,0,3500);
hT_IHRD[i]=new TH1F(Form("hT_IHRD%d",i+1),Form("QDC IHRD %d",i+1),3500,0,3500);
SetTH1(hT_IHLU[i], Form("TDC IHLU %d",i+1), "TDC[ch]", "count", 4, 0, 0);
SetTH1(hT_IHLD[i], Form("TDC IHLD %d",i+1), "TDC[ch]", "count", 4, 0, 0);
SetTH1(hT_IHRU[i], Form("TDC IHRU %d",i+1), "TDC[ch]", "count", 4, 0, 0);
SetTH1(hT_IHRD[i], Form("TDC IHRD %d",i+1), "TDC[ch]", "count", 4, 0, 0);

hQ_IHLU_wt[i]=new TH1F(Form("hQ_IHLU_wt%d",i+1),Form("IHLU QDC %d w/ TDC",i+1),500,0,1500);
hQ_IHLD_wt[i]=new TH1F(Form("hQ_IHLD_wt%d",i+1),Form("IHLD QDC %d w/ TDC",i+1),500,0,1500);
hQ_IHRU_wt[i]=new TH1F(Form("hQ_IHRU_wt%d",i+1),Form("IHRU QDC %d w/ TDC",i+1),500,0,1500);
hQ_IHRD_wt[i]=new TH1F(Form("hQ_IHRD_wt%d",i+1),Form("IHRD QDC %d w/ TDC",i+1),500,0,1500);
SetTH1(hQ_IHLU_wt[i], Form("IHLU %d w/ TDC",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(hQ_IHLD_wt[i], Form("IHLD %d w/ TDC",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(hQ_IHRU_wt[i], Form("IHRU %d w/ TDC",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(hQ_IHRD_wt[i], Form("IHRD %d w/ TDC",i+1), "ADC[ch]", "count", 2, 0, 0);

tree->Project(Form("hQ_IHLU%d",i+1),Form("ihluadc[%d]",i),"");//Form("ihlutdc[%d]<=0",i)
tree->Project(Form("hQ_IHLD%d",i+1),Form("ihldadc[%d]",i),"");//Form("ihldtdc[%d]<=0",i)
tree->Project(Form("hQ_IHRU%d",i+1),Form("ihruadc[%d]",i),"");//Form("ihrutdc[%d]<=0",i)
tree->Project(Form("hQ_IHRD%d",i+1),Form("ihrdadc[%d]",i),"");//Form("ihrdtdc[%d]<=0",i)

tree->Project(Form("hT_IHLU%d",i+1),Form("ihlutdc[%d]",i),Form("ihlutdc[%d]>0",i));
tree->Project(Form("hT_IHLD%d",i+1),Form("ihldtdc[%d]",i),Form("ihldtdc[%d]>0",i));
tree->Project(Form("hT_IHRU%d",i+1),Form("ihrutdc[%d]",i),Form("ihrutdc[%d]>0",i));
tree->Project(Form("hT_IHRD%d",i+1),Form("ihrdtdc[%d]",i),Form("ihrdtdc[%d]>0",i));

tree->Project(Form("hQ_IHLU_wt%d",i+1),Form("ihluadc[%d]",i),Form("ihlutdc[%d]>0",i));
tree->Project(Form("hQ_IHLD_wt%d",i+1),Form("ihldadc[%d]",i),Form("ihldtdc[%d]>0",i));
tree->Project(Form("hQ_IHRU_wt%d",i+1),Form("ihruadc[%d]",i),Form("ihrutdc[%d]>0",i));
tree->Project(Form("hQ_IHRD_wt%d",i+1),Form("ihrdadc[%d]",i),Form("ihrdtdc[%d]>0",i));
 }
cout<<"IH filled"<<endl;
for(int i = 0; i<12 ; i++){//OHV
hQ_OHVLU[i]=new TH1F(Form("hQ_OHVLU%d",i+1),Form("Pedestal OHVLU %d",i+1),1000,0,1000);
hQ_OHVLD[i]=new TH1F(Form("hQ_OHVLD%d",i+1),Form("Pedestal OHVLD %d",i+1),1000,0,1000);
hQ_OHVRU[i]=new TH1F(Form("hQ_OHVRU%d",i+1),Form("Pedestal OHVRU %d",i+1),1000,0,1000);
hQ_OHVRD[i]=new TH1F(Form("hQ_OHVRD%d",i+1),Form("Pedestal OHVRD %d",i+1),1000,0,1000);
SetTH1(hQ_OHVLU[i], Form("Pedestal OHVLU %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(hQ_OHVLD[i], Form("Pedestal OHVLD %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(hQ_OHVRU[i], Form("Pedestal OHVRU %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(hQ_OHVRD[i], Form("Pedestal OHVRD %d",i+1), "ADC[ch]", "count", 1, 0, 0);

hT_OHVLU[i]=new TH1F(Form("hT_OHVLU%d",i+1),Form("OHVLU %d TDC",i+1),1000,0,3000);
hT_OHVLD[i]=new TH1F(Form("hT_OHVLD%d",i+1),Form("OHVLD %d TDC",i+1),1000,0,3000);
hT_OHVRU[i]=new TH1F(Form("hT_OHVRU%d",i+1),Form("OHVRU %d TDC",i+1),1000,0,3000);
hT_OHVRD[i]=new TH1F(Form("hT_OHVRD%d",i+1),Form("OHVRD %d TDC",i+1),1000,0,3000);
SetTH1(hT_OHVLU[i], Form("OHVLU %d",i+1), "TDC[ch]", "count", 4, 0, 0);
SetTH1(hT_OHVLD[i], Form("OHVLD %d",i+1), "TDC[ch]", "count", 4, 0, 0);
SetTH1(hT_OHVRU[i], Form("OHVRU %d",i+1), "TDC[ch]", "count", 4, 0, 0);
SetTH1(hT_OHVRD[i], Form("OHVRD %d",i+1), "TDC[ch]", "count", 4, 0, 0);

hQ_OHVLU_wt[i]=new TH1F(Form("hQ_OHVLU_wt%d",i+1),Form("OHVLU %d w/ TDC",i+1),1000,0,1000);
hQ_OHVLD_wt[i]=new TH1F(Form("hQ_OHVLD_wt%d",i+1),Form("OHVLD %d w/ TDC",i+1),1000,0,1000);
hQ_OHVRU_wt[i]=new TH1F(Form("hQ_OHVRU_wt%d",i+1),Form("OHVRU %d w/ TDC",i+1),1000,0,1000);
hQ_OHVRD_wt[i]=new TH1F(Form("hQ_OHVRD_wt%d",i+1),Form("OHVRD %d w/ TDC",i+1),1000,0,1000);
SetTH1(hQ_OHVLU_wt[i], Form("OHVLU %d",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(hQ_OHVLD_wt[i], Form("OHVLD %d",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(hQ_OHVRU_wt[i], Form("OHVRU %d",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(hQ_OHVRD_wt[i], Form("OHVRD %d",i+1), "ADC[ch]", "count", 2, 0, 0);

tree->Project(Form("hQ_OHVLU%d",i+1),Form("ohvluadc[%d]",i),"");//Form("ohvlutdc[%d]<=0",i)
tree->Project(Form("hQ_OHVLD%d",i+1),Form("ohvldadc[%d]",i),"");//Form("ohvldtdc[%d]<=0",i)
tree->Project(Form("hQ_OHVRU%d",i+1),Form("ohvruadc[%d]",i),"");//Form("ohvrutdc[%d]<=0",i)
tree->Project(Form("hQ_OHVRD%d",i+1),Form("ohvrdadc[%d]",i),"");//Form("ohvrdtdc[%d]<=0",i)

tree->Project(Form("hT_OHVLU%d",i+1),Form("ohvlutdc[%d]",i),Form("ohvlutdc[%d]>0",i));
tree->Project(Form("hT_OHVLD%d",i+1),Form("ohvldtdc[%d]",i),Form("ohvldtdc[%d]>0",i));
tree->Project(Form("hT_OHVRU%d",i+1),Form("ohvrutdc[%d]",i),Form("ohvrutdc[%d]>0",i));
tree->Project(Form("hT_OHVRD%d",i+1),Form("ohvrdtdc[%d]",i),Form("ohvrdtdc[%d]>0",i));

tree->Project(Form("hQ_OHVLU_wt%d",i+1),Form("ohvluadc[%d]",i),Form("ohvlutdc[%d]>0",i));
tree->Project(Form("hQ_OHVLD_wt%d",i+1),Form("ohvldadc[%d]",i),Form("ohvldtdc[%d]>0",i));
tree->Project(Form("hQ_OHVRU_wt%d",i+1),Form("ohvruadc[%d]",i),Form("ohvrutdc[%d]>0",i));
tree->Project(Form("hQ_OHVRD_wt%d",i+1),Form("ohvrdadc[%d]",i),Form("ohvrdtdc[%d]>0",i));
}
cout<<"OHV filled"<<endl;
for(int i = 0; i<9 ; i++){//OHH
hQ_OHHLU[i]=new TH1F(Form("hQ_OHHLU%d",i+1),Form("Pedestal OHHLU %d",i+1),1000,0,1000);
hQ_OHHLD[i]=new TH1F(Form("hQ_OHHLD%d",i+1),Form("Pedestal OHHLD %d",i+1),1000,0,1000);
hQ_OHHRU[i]=new TH1F(Form("hQ_OHHRU%d",i+1),Form("Pedestal OHHRU %d",i+1),1000,0,1000);
hQ_OHHRD[i]=new TH1F(Form("hQ_OHHRD%d",i+1),Form("Pedestal OHHRD %d",i+1),1000,0,1000);
SetTH1(hQ_OHHLU[i], Form("OHHLU %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(hQ_OHHLD[i], Form("OHHLD %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(hQ_OHHRU[i], Form("OHHRU %d",i+1), "ADC[ch]", "count", 1, 0, 0);
SetTH1(hQ_OHHRD[i], Form("OHHRD %d",i+1), "ADC[ch]", "count", 1, 0, 0);

hT_OHHLU[i]=new TH1F(Form("hT_OHHLU%d",i+1),Form("OHHLU %d TDC",i+1),1000,0,3000);
hT_OHHLD[i]=new TH1F(Form("hT_OHHLD%d",i+1),Form("OHHLD %d TDC",i+1),1000,0,3000);
hT_OHHRU[i]=new TH1F(Form("hT_OHHRU%d",i+1),Form("OHHRU %d TDC",i+1),1000,0,3000);
hT_OHHRD[i]=new TH1F(Form("hT_OHHRD%d",i+1),Form("OHHRD %d TDC",i+1),1000,0,3000);
SetTH1(hT_OHHLU[i], Form("OHHLU %d",i+1), "TDC[ch]", "count", 4, 0, 0);
SetTH1(hT_OHHLD[i], Form("OHHLD %d",i+1), "TDC[ch]", "count", 4, 0, 0);
SetTH1(hT_OHHRU[i], Form("OHHRU %d",i+1), "TDC[ch]", "count", 4, 0, 0);
SetTH1(hT_OHHRD[i], Form("OHHRD %d",i+1), "TDC[ch]", "count", 4, 0, 0);

hQ_OHHLU_wt[i]=new TH1F(Form("hQ_OHHLU_wt%d",i+1),Form("OHHLU %d w/ TDC",i+1),1000,0,1000);
hQ_OHHLD_wt[i]=new TH1F(Form("hQ_OHHLD_wt%d",i+1),Form("OHHLD %d w/ TDC",i+1),1000,0,1000);
hQ_OHHRU_wt[i]=new TH1F(Form("hQ_OHHRU_wt%d",i+1),Form("OHHRU %d w/ TDC",i+1),1000,0,1000);
hQ_OHHRD_wt[i]=new TH1F(Form("hQ_OHHRD_wt%d",i+1),Form("OHHRD %d w/ TDC",i+1),1000,0,1000);
SetTH1(hQ_OHHLU_wt[i], Form("OHHLU %d",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(hQ_OHHLD_wt[i], Form("OHHLD %d",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(hQ_OHHRU_wt[i], Form("OHHRU %d",i+1), "ADC[ch]", "count", 2, 0, 0);
SetTH1(hQ_OHHRD_wt[i], Form("OHHRD %d",i+1), "ADC[ch]", "count", 2, 0, 0);

tree->Project(Form("hQ_OHHLU%d",i+1),Form("ohhluadc[%d]",i),"");//Form("ohhlutdc[%d]<=0",i)
tree->Project(Form("hQ_OHHLD%d",i+1),Form("ohhldadc[%d]",i),"");//Form("ohhldtdc[%d]<=0",i)
tree->Project(Form("hQ_OHHRU%d",i+1),Form("ohhruadc[%d]",i),"");//Form("ohhrutdc[%d]<=0",i)
tree->Project(Form("hQ_OHHRD%d",i+1),Form("ohhrdadc[%d]",i),"");//Form("ohhrdtdc[%d]<=0",i)

tree->Project(Form("hT_OHHLU%d",i+1),Form("ohhlutdc[%d]",i),Form("ohhlutdc[%d]>0",i));
tree->Project(Form("hT_OHHLD%d",i+1),Form("ohhldtdc[%d]",i),Form("ohhldtdc[%d]>0",i));
tree->Project(Form("hT_OHHRU%d",i+1),Form("ohhrutdc[%d]",i),Form("ohhrutdc[%d]>0",i));
tree->Project(Form("hT_OHHRD%d",i+1),Form("ohhrdtdc[%d]",i),Form("ohhrdtdc[%d]>0",i));

tree->Project(Form("hQ_OHHLU_wt%d",i+1),Form("ohhluadc[%d]",i),Form("ohhlutdc[%d]>0",i));
tree->Project(Form("hQ_OHHLD_wt%d",i+1),Form("ohhldadc[%d]",i),Form("ohhldtdc[%d]>0",i));
tree->Project(Form("hQ_OHHRU_wt%d",i+1),Form("ohhruadc[%d]",i),Form("ohhrutdc[%d]>0",i));
tree->Project(Form("hQ_OHHRD_wt%d",i+1),Form("ohhrdadc[%d]",i),Form("ohhrdtdc[%d]>0",i));
}
cout<<"OHH filled"<<endl;

cout<<"start drawing!"<<endl;
TCanvas *c[14];
for(int i=0;i<14;i++){
c[i] = new TCanvas(Form("c%d",i+1), Form("canvas%d",i+1) , 800,1500);}

c[0]->Divide(2,10);
for(int i=0;i<10;i++){
c[0]->cd(2*i+2); gPad->SetLogy(1); hQ_IHLU[i]  ->Draw();hQ_IHLU_wt[i]  ->Draw("same");
c[0]->cd(2*i+1); gPad->SetLogy(1); hT_IHLU[i]  ->Draw();

 }                                                                                                        
c[1]->Divide(2,10);                                                                                       
for(int i=0;i<10;i++){                                                                                    
c[1]->cd(2*i+2); gPad->SetLogy(1); hQ_IHLD[i]  ->Draw();hQ_IHLD_wt[i]  ->Draw("same");
c[1]->cd(2*i+1); gPad->SetLogy(1); hT_IHLD[i]  ->Draw();
 }                                                                       
c[2]->Divide(2,10);                                                       
for(int i=0;i<10;i++){                                                   
c[2]->cd(2*i+2); gPad->SetLogy(1); hQ_IHRU[i]  ->Draw();hQ_IHRU_wt[i]  ->Draw("same");
c[2]->cd(2*i+1); gPad->SetLogy(1); hT_IHRU[i]  ->Draw();
 }                                                                                                        
c[3]->Divide(2,10);                                                                                     
for(int i=0;i<10;i++){                                                                                    
c[3]->cd(2*i+2); gPad->SetLogy(1); hQ_IHRD[i]  ->Draw();hQ_IHRD_wt[i]  ->Draw("same");
c[3]->cd(2*i+1); gPad->SetLogy(1); hT_IHRD[i]  ->Draw();
 }

c[4]->Divide(2,12);
for(int i=0;i<12;i++){
c[4]->cd(2*i+2); gPad->SetLogy(1); hQ_OHVLU[i]  ->Draw();hQ_OHVLU_wt[i]  ->Draw("same");
c[4]->cd(2*i+1); gPad->SetLogy(1); hT_OHVLU[i]  ->Draw();
 }                                                                                    
c[5]->Divide(2,12);                                                                    
for(int i=0;i<12;i++){                                                                
c[5]->cd(2*i+2); gPad->SetLogy(1); hQ_OHVLD[i]  ->Draw();hQ_OHVLD_wt[i]  ->Draw("same");
c[5]->cd(2*i+1); gPad->SetLogy(1); hT_OHVLD[i]  ->Draw();
 }                                                                                    
c[6]->Divide(2,12);                                                                    
for(int i=0;i<12;i++){                                                                
c[6]->cd(2*i+2); gPad->SetLogy(1); hQ_OHVRU[i]  ->Draw();hQ_OHVRU_wt[i]  ->Draw("same");
c[6]->cd(2*i+1); gPad->SetLogy(1); hT_OHVRU[i]  ->Draw();
 }                                                                                    
cout<<"drawing"<<endl;                                                                
c[7]->Divide(2,12);                                                                    
for(int i=0;i<12;i++){                                                                
c[7]->cd(2*i+2); gPad->SetLogy(1); hQ_OHVRD[i]  ->Draw();hQ_OHVRD_wt[i]  ->Draw("same");
c[7]->cd(2*i+1); gPad->SetLogy(1); hT_OHVRD[i]  ->Draw();
 }                                                                                    
c[8]->Divide(2,9);                                                                    
for(int i=0;i<9;i++){                                                                 
c[8]->cd(2*i+2); gPad->SetLogy(1); hQ_OHHLU[i]  ->Draw();hQ_OHHLU_wt[i]  ->Draw("same");
c[8]->cd(2*i+1); gPad->SetLogy(1); hT_OHHLU[i]  ->Draw();
 }                                                                                    
cout<<"drawing OHH"<<endl;                                                            
c[9]->Divide(2,9);                                                                    
for(int i=0;i<9;i++){                                                                 
c[9]->cd(2*i+2); gPad->SetLogy(1); hQ_OHHLD[i]  ->Draw();hQ_OHHLD_wt[i]  ->Draw("same");
c[9]->cd(2*i+1); gPad->SetLogy(1); hT_OHHLD[i]  ->Draw();
 }
c[10]->Divide(2,9);
for(int i=0;i<9;i++){
c[10]->cd(2*i+2); gPad->SetLogy(1); hQ_OHHRU[i]  ->Draw();hQ_OHHRU_wt[i]  ->Draw("same");
c[10]->cd(2*i+1); gPad->SetLogy(1); hT_OHHRU[i]  ->Draw();
 }                                                                                     
c[11]->Divide(2,9);                                                                    
for(int i=0;i<9;i++){                                                                  
c[11]->cd(2*i+2); gPad->SetLogy(1); hQ_OHHRD[i]  ->Draw();hQ_OHHRD_wt[i]  ->Draw("same");
c[11]->cd(2*i+1); gPad->SetLogy(1); hT_OHHRD[i]  ->Draw();
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
  c[11]->Print(Form("%s]",ofname_pdf.c_str()));


}

