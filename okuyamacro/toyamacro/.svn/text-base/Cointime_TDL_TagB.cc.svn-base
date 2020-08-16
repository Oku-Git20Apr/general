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
const int NIH   = 10;  // No. of IH
const int NTB   = 40;  // No. of TagB
const int NTF   =160;  // No. of TagB

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
  h->SetMarkerStyle(1);
  h->SetMarkerSize(0.1);
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
  string ofname = "root/hoge.root";
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
  if(!draw_flag) gROOT->SetBatch(0);

  gStyle->SetOptDate(0);
  gStyle->SetOptFit(1);
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

  TFile *ofp = new TFile();
  if(output_flag){ ofp = new TFile(Form("%s",ofname.c_str()),"recreate"); }

  int ENum;
  int runnum,evnum;
  int tdllutdc_l[NIH], tdlrutdc_l[NIH], tdlldtdc_l[NIH],   tdlrdtdc_l[NIH];   // TDC (leading)
  int tdllutdc_t[NIH], tdlrutdc_t[NIH], tdlldtdc_t[NIH],   tdlrdtdc_t[NIH];   // TDC (trailing)
  int tagbtdc_l[NTB], tagbtdc_t[NTB],tagftdc_l[NTF],tagftdc_t[NTF];
  double tdllutime_l[NIH], tdlrutime_l[NIH], tdlldtime_l[NIH],   tdlrdtime_l[NIH];   // time (leading)
  double tdllutime_t[NIH], tdlrutime_t[NIH], tdlldtime_t[NIH],   tdlrdtime_t[NIH];   // time (trailing)
  double tdllutime_w[NIH], tdlrutime_w[NIH], tdlldtime_w[NIH],   tdlrdtime_w[NIH];   // time (width)
  double tdllmtime[NIH], tdlrmtime[NIH];   // time (width)
  double tagbtime_l[NTB], tagbtime_t[NTB],tagftime_l[NTF],tagftime_t[NTF];
  double tagbtime_w[NTB], tagftime_w[NTF];
  double tagb_trig_time, tagf_trig_time1, tagf_trig_time2;//trigger time
  double trig_time1, trig_time2;//trigger time
  double t_offset = 0;
  ENum = tree->GetEntries();
  if(MaxNum>0 && MaxNum<ENum) ENum = MaxNum;
  tree->SetBranchAddress("runnum"       ,&runnum     );
  tree->SetBranchAddress("evnum"        ,&evnum      );
  tree->SetBranchAddress("tdllutdc_l"   , tdllutdc_l );
  tree->SetBranchAddress("tdlrutdc_l"   , tdlrutdc_l );
  tree->SetBranchAddress("tdlldtdc_l"   , tdlldtdc_l );
  tree->SetBranchAddress("tdlrdtdc_l"   , tdlrdtdc_l );
  tree->SetBranchAddress("tdllutdc_t"   , tdllutdc_t );
  tree->SetBranchAddress("tdlrutdc_t"   , tdlrutdc_t );
  tree->SetBranchAddress("tdlldtdc_t"   , tdlldtdc_t );
  tree->SetBranchAddress("tdlrdtdc_t"   , tdlrdtdc_t );
  tree->SetBranchAddress("tagbtdc_l"    , tagbtdc_l );
  tree->SetBranchAddress("tagbtdc_t"    , tagbtdc_t );
  tree->SetBranchAddress("tagftdc_l"    , tagftdc_l );
  tree->SetBranchAddress("tagftdc_t"    , tagftdc_t );

gStyle->SetOptStat("iMen");

TH1F *h_TagSumTDLL[NIH], *h_TagSumTDLR[NIH];
TH1F *h_TagTDLL[NTB][NIH], *h_TagTDLR[NTB][NIH];
TH1F *h_TagTDLLU[NTB][NIH], *h_TagTDLRU[NTB][NIH];
TH1F *h_TagTDLLD[NTB][NIH], *h_TagTDLRD[NTB][NIH];
TH1F *h1_luTDC_l[NIH], *h1_ldTDC_l[NIH], *h1_ruTDC_l[NIH], *h1_rdTDC_l[NIH];
TH1F *h1_luTDC_t[NIH], *h1_ldTDC_t[NIH], *h1_ruTDC_t[NIH], *h1_rdTDC_t[NIH];
TH1F *h1_lutime_w[NIH], *h1_ldtime_w[NIH], *h1_rutime_w[NIH], *h1_rdtime_w[NIH];
TH2D *h2_luTDCl_evnum[NIH], *h2_ldTDCl_evnum[NIH], *h2_ruTDCl_evnum[NIH], *h2_rdTDCl_evnum[NIH];
TH2D *h2_daq0trig_evnum;
TH2D *h2_daq1trig_evnum;
h2_daq0trig_evnum = new TH2D("h2_daq0trig_evnum","h2_daq0trig_evnum",1000,0,500000,1000,-110000,10000);
h2_daq1trig_evnum = new TH2D("h2_daq1trig_evnum","h2_daq1trig_evnum",1000,0,500000,1000,-11000,10000);


 for(int j=0;j<NIH;j++){
  for(int i=0;i<NTB;i++){
  //L
  h_TagTDLLU[i][j] = new TH1F(Form("h_Tag%dTDLL%dU",i+1,j+1),Form("h_Tag%dTDLL%dU",i+1,j+1),20000,-1000,1000);
  SetTH1(h_TagTDLLU[i][j], Form("ToF TagB %d- TDL L%dU",i+1,j+1), "coin. time[ns]", "count/0.1ns", 1, 0, 0);
  h_TagTDLLD[i][j] = new TH1F(Form("h_Tag%dTDLL%dD",i+1,j+1),Form("h_Tag%dTDLL%dD",i+1,j+1),20000,-1000,1000);
  SetTH1(h_TagTDLLD[i][j], Form("ToF TagB %d- TDL L%dD",i+1,j+1), "coin. time[ns]", "count/0.1ns", 1, 0, 0);
  h_TagTDLL[i][j] = new TH1F(Form("h_Tag%dTDLL%d",i+1,j+1),Form("h_Tag%dTDLL%d",i+1,j+1),20000,-1000,1000);
  SetTH1(h_TagTDLL[i][j], Form("ToF TagB %d- TDL L%d",i+1,j+1), "coin. time[ns]", "count/0.1ns", 1, 0, 0);
  //R
  h_TagTDLRU[i][j] = new TH1F(Form("h_Tag%dTDLR%dU",i+1,j+1),Form("h_Tag%dTDLR%dU",i+1,j+1),20000,-1000,1000);
  SetTH1(h_TagTDLRU[i][j], Form("ToF TagB %d- TDL R%dU",i+1,j+1), "coin. time[ns]", "count/0.1ns", 1, 0, 0);
  h_TagTDLRD[i][j] = new TH1F(Form("h_Tag%dTDLR%dD",i+1,j+1),Form("h_Tag%dTDLR%dD",i+1,j+1),20000,-1000,1000);
  SetTH1(h_TagTDLRD[i][j], Form("ToF TagB %d- TDL R%dD",i+1,j+1), "coin. time[ns]", "count/0.1ns", 1, 0, 0);
  h_TagTDLR[i][j] = new TH1F(Form("h_Tag%dTDLR%d",i+1,j+1),Form("h_Tag%dTDLR%d",i+1,j+1),20000,-1000,1000);
  SetTH1(h_TagTDLR[i][j], Form("ToF TagB %d- TDL R%d",i+1,j+1), "coin. time[ns]", "count/0.1ns", 1, 0, 0);
  }
 }

TH2F *h2R4R10ToF_R4UWidth  = new TH2F("h2R4R10ToF_R4UWidth" ,"h2R4R10ToF_R4UWidth" ,2000,-100,100,1000,0,700);
TH2F *h2R4R10ToF_R4DWidth  = new TH2F("h2R4R10ToF_R4DWidth" ,"h2R4R10ToF_R4DWidth" ,2000,-100,100,1000,0,700);
TH2F *h2R4R10ToF_R10UWidth = new TH2F("h2R4R10ToF_R10UWidth","h2R4R10ToF_R10UWidth",2000,-100,100,1000,0,700);
TH2F *h2R4R10ToF_R10DWidth = new TH2F("h2R4R10ToF_R10DWidth","h2R4R10ToF_R10DWidth",2000,-100,100,1000,0,700);

TH1F *hToF_wt = new TH1F("hToF_wt","hToF_wt",2000,-100,100);
SetTH1(hToF_wt, "ToF R3-R10(tekito counter) ", "ToF[ns]", "count", 1, 0, 0);
TH2F *h2_luWvsToF[NIH],*h2_luWvsToF_cut[NIH];
TH2F *h2_ldWvsToF[NIH],*h2_ldWvsToF_cut[NIH];
TH2F *h2_ruWvsToF[NIH],*h2_ruWvsToF_cut[NIH];
TH2F *h2_rdWvsToF[NIH],*h2_rdWvsToF_cut[NIH];

for(int i=0;i<NIH;i++){
h1_luTDC_l[i] = new TH1F(Form("h1_luTDC_l%d",i+1),Form("h1_luTDC_l%d",i+1),10000,0,100000);
h1_ldTDC_l[i] = new TH1F(Form("h1_ldTDC_l%d",i+1),Form("h1_ldTDC_l%d",i+1),10000,0,100000);
h1_ruTDC_l[i] = new TH1F(Form("h1_ruTDC_l%d",i+1),Form("h1_ruTDC_l%d",i+1),10000,0,100000);
h1_rdTDC_l[i] = new TH1F(Form("h1_rdTDC_l%d",i+1),Form("h1_rdTDC_l%d",i+1),10000,0,100000);

h1_luTDC_t[i] = new TH1F(Form("h1_luTDC_t%d",i+1),Form("h1_luTDC_t%d",i+1),10000,0,100000);
h1_ldTDC_t[i] = new TH1F(Form("h1_ldTDC_t%d",i+1),Form("h1_ldTDC_t%d",i+1),10000,0,100000);
h1_ruTDC_t[i] = new TH1F(Form("h1_ruTDC_t%d",i+1),Form("h1_ruTDC_t%d",i+1),10000,0,100000);
h1_rdTDC_t[i] = new TH1F(Form("h1_rdTDC_t%d",i+1),Form("h1_rdTDC_t%d",i+1),10000,0,100000);
//
//h2_luTDCl_evnum[i] = new TH2D(Form("h2_luTDCl%d_evnum",i+1),Form("h2_luTDCl%d_evnum",i+1),1000,0,500000,1000,-11000,100000);
//h2_ldTDCl_evnum[i] = new TH2D(Form("h2_ldTDCl%d_evnum",i+1),Form("h2_ldTDCl%d_evnum",i+1),1000,0,500000,1000,-11000,100000);
//h2_ruTDCl_evnum[i] = new TH2D(Form("h2_ruTDCl%d_evnum",i+1),Form("h2_ruTDCl%d_evnum",i+1),1000,0,500000,1000,-11000,100000);
//h2_rdTDCl_evnum[i] = new TH2D(Form("h2_rdTDCl%d_evnum",i+1),Form("h2_rdTDCl%d_evnum",i+1),1000,0,500000,1000,-11000,100000);

h1_lutime_w[i] = new TH1F(Form("h1_lutime_w%d",i+1),Form("h1_lutime_w%d",i+1),1000,1,60);
h1_ldtime_w[i] = new TH1F(Form("h1_ldtime_w%d",i+1),Form("h1_ldtime_w%d",i+1),1000,1,60);
h1_rutime_w[i] = new TH1F(Form("h1_rutime_w%d",i+1),Form("h1_rutime_w%d",i+1),1000,1,60);
h1_rdtime_w[i] = new TH1F(Form("h1_rdtime_w%d",i+1),Form("h1_rdtime_w%d",i+1),1000,1,60);
  h_TagSumTDLL[i] = new TH1F(Form("h_TagSumTDLL%d",i+1),Form("h_TagSumTDLL%d",i+1),20000,-1000,1000);
  h_TagSumTDLR[i] = new TH1F(Form("h_TagSumTDLR%d",i+1),Form("h_TagSumTDLR%d",i+1),20000,-1000,1000);
  SetTH1(h_TagSumTDLL[i], Form("ToF TagB Sum - TDL L%d",i+1), "coin. time[ns]", "count/0.1ns", 1, 0, 0);
  SetTH1(h_TagSumTDLR[i], Form("ToF TagB Sum - TDL R%d",i+1), "coin. time[ns]", "count/0.1ns", 1, 0, 0);
h2_luWvsToF[i] = new TH2F(Form("h2_luWvsToF%d",i+1),Form("h2_luWvsToF%d",i+1),1000,800,900,4000,0,100);
h2_ldWvsToF[i] = new TH2F(Form("h2_ldWvsToF%d",i+1),Form("h2_ldWvsToF%d",i+1),1000,800,900,4000,0,100);
h2_ruWvsToF[i] = new TH2F(Form("h2_ruWvsToF%d",i+1),Form("h2_ruWvsToF%d",i+1),1000,800,900,4000,0,100);
h2_rdWvsToF[i] = new TH2F(Form("h2_rdWvsToF%d",i+1),Form("h2_rdWvsToF%d",i+1),1000,800,900,4000,0,100);
SetTH2(h2_luWvsToF[i], Form("Width LUseg%d vs ToF",i+1), "ToF [ns]", Form("width LU%d [ns]",i+1));
SetTH2(h2_ldWvsToF[i], Form("Width LDseg%d vs ToF",i+1), "ToF [ns]", Form("width LD%d [ns]",i+1));
SetTH2(h2_ruWvsToF[i], Form("Width RUseg%d vs ToF",i+1), "ToF [ns]", Form("width RU%d [ns]",i+1));
SetTH2(h2_rdWvsToF[i], Form("Width RDseg%d vs ToF",i+1), "ToF [ns]", Form("width RD%d [ns]",i+1));
/////
h2_luWvsToF_cut[i] = new TH2F(Form("h2_luWvsToF_cut%d",i+1),Form("h2_luWvsToF_cut%d",i+1),100,-3,3,50,0,400);
h2_ldWvsToF_cut[i] = new TH2F(Form("h2_ldWvsToF_cut%d",i+1),Form("h2_ldWvsToF_cut%d",i+1),100,-3,3,50,0,400);
h2_ruWvsToF_cut[i] = new TH2F(Form("h2_ruWvsToF_cut%d",i+1),Form("h2_ruWvsToF_cut%d",i+1),100,-3,3,50,0,400);
h2_rdWvsToF_cut[i] = new TH2F(Form("h2_rdWvsToF_cut%d",i+1),Form("h2_rdWvsToF_cut%d",i+1),100,-3,3,50,0,400);
SetTH2(h2_luWvsToF_cut[i], Form("Width LUseg%d vs ToF cut",i+1), "ToF [ns]", Form("width LU%d",i+1));
SetTH2(h2_ldWvsToF_cut[i], Form("Width LDseg%d vs ToF cut",i+1), "ToF [ns]", Form("width LD%d",i+1));
SetTH2(h2_ruWvsToF_cut[i], Form("Width RUseg%d vs ToF cut",i+1), "ToF [ns]", Form("width RU%d",i+1));
SetTH2(h2_rdWvsToF_cut[i], Form("Width RDseg%d vs ToF cut",i+1), "ToF [ns]", Form("width RD%d",i+1));
 }
TCut trigger = " tdlrutdc_l[0]>0";

////////////////////////////////////////////////////////////////////
//////////////       Fill
////////////////////////////////////////////////////////////////////
// event for loop
  for(int n=0;n<ENum;n++){
    tree->GetEntry(n);

    if(n%10000==0){
     char clout[100];
     end = time(NULL);
     time(&end);
     sprintf(clout,"%.0f sec",difftime(end,start));
     cout<<n<<" / "<<ENum<<" : "<<clout<<endl;
    } // if n%100000==0 std::cout

///////////
// Param //
///////////
double Width_R4Umin,Width_R4Dmin,Width_R10Umin,Width_R10Dmin;
double Width_R4Umax,Width_R4Dmax,Width_R10Umax,Width_R10Dmax;
//run10077, Vb = 58V
if(runnum==10077){
Width_R4Umin = 155.;Width_R4Dmin =  80.;Width_R10Umin =  90.;Width_R10Dmin = 120.;
Width_R4Umax = 175.;Width_R4Dmax = 100.;Width_R10Umax = 110.;Width_R10Dmax = 140.;
}
//////run10078, Vb = 60V
else if(runnum==10078){
Width_R4Umin = 220.;Width_R4Dmin = 120.;Width_R10Umin = 135.;Width_R10Dmin = 190.;
Width_R4Umax = 240.;Width_R4Dmax = 140.;Width_R10Umax = 155.;Width_R10Dmax = 210.;
}
//////run10082, Vb = 62V
else if(runnum==10082){
Width_R4Umin = 280.;Width_R4Dmin = 180.;Width_R10Umin = 190.;Width_R10Dmin = 240.;
Width_R4Umax = 300.;Width_R4Dmax = 200.;Width_R10Umax = 210.;Width_R10Dmax = 260.;
}
//run10114, Vb = 58V
if(runnum==10114){
Width_R4Umin = 155.;Width_R4Dmin =  80.;Width_R10Umin =  90.;Width_R10Dmin = 120.;
Width_R4Umax = 175.;Width_R4Dmax = 100.;Width_R10Umax = 110.;Width_R10Dmax = 140.;
}
//run10115, Vb = 60V
if(runnum==10115){
Width_R4Umin = 155.;Width_R4Dmin =  80.;Width_R10Umin =  90.;Width_R10Dmin = 120.;
Width_R4Umax = 175.;Width_R4Dmax = 100.;Width_R10Umax = 110.;Width_R10Dmax = 140.;
}
//run10116, Vb = 62V
if(runnum==10115){
Width_R4Umin = 155.;Width_R4Dmin =  80.;Width_R10Umin =  90.;Width_R10Dmin = 120.;
Width_R4Umax = 175.;Width_R4Dmax = 100.;Width_R10Umax = 110.;Width_R10Dmax = 140.;
}
else {
Width_R4Umin = 18.;Width_R4Dmin = 18.;Width_R10Umin = 18.;Width_R10Dmin = 18.;
Width_R4Umax = 22.;Width_R4Dmax = 22.;Width_R10Umax = 22.;Width_R10Dmax = 22.;
}
//////////////
//initialize//
//////////////
  for(int i=0; i<NIH; i++){
//   tdllutdc_l[i] = tdlrutdc_l[i] = tdlldtdc_l[i] = tdlrdtdc_l[i] = -9999;
//   tdllutdc_t[i] = tdlrutdc_t[i] = tdlldtdc_t[i] = tdlrdtdc_t[i] = -9999;
  tdllutime_l[i] =  tdlrutime_l[i] =  tdlldtime_l[i] =    tdlrdtime_l[i] = -9999;   // time (leading)
  tdllutime_t[i] =  tdlrutime_t[i] =  tdlldtime_t[i] =    tdlrdtime_t[i] = -9999;   // time (trailing)
  tdllutime_w[i] =  tdlrutime_w[i] =  tdlldtime_w[i] =    tdlrdtime_w[i] = -9999;   // time (width)
  tdllmtime[i]   =  tdlrmtime[i]   = -9999;
  }

  for(int i=0; i<NTB; i++){
  tagbtime_l[i] = tagbtime_t[i] =tagbtime_w[i]=-9999;
  }
  for(int i=0; i<NTF; i++){
  tagftime_l[i] = tagftime_t[i]=-9999;
  }
  tagb_trig_time = tagf_trig_time1 = tagf_trig_time2=-9999;//trigger time
  trig_time1=-9999;
  trig_time2=-9999;//trigger time

if(tdllutdc_l[0]>0) trig_time1 = tdllutdc_l[0]*0.025;//QTC2
if(tdlrutdc_l[0]>0) trig_time2 = tdlrutdc_l[0]*0.025;//QTC3(MRPC datta yatsu)
if(tagbtdc_l[38]>0) tagb_trig_time = tagbtdc_l[38]*0.025;//QTC3(MRPC datta yatsu)
//cout<<"trigger timing ="<<trig_time1<<endl;
if(trig_time1<0 || trig_time2<0)continue;
//cout<<"trigger timing1 ="<<trig_time1<<"trigger timing2 ="<<trig_time2<<endl;

  for(int i=0; i<NIH; i++){
//LU
  tdllutime_l[i] = trig_time1 - tdllutdc_l[i]*0.025;   // time (leading)
  tdllutime_t[i] = trig_time1 - tdllutdc_t[i]*0.025;   // time (trailing)
  tdllutime_w[i] = - tdllutime_t[i] + tdllutime_l[i];   // time (width)
//RU
  tdlrutime_l[i] = trig_time1 - tdlrutdc_l[i]*0.025;
  tdlrutime_t[i] = trig_time1 - tdlrutdc_t[i]*0.025;
  tdlrutime_w[i] = - tdlrutime_t[i] + tdlrutime_l[i];
//LD
  tdlldtime_l[i] = trig_time1 - tdlldtdc_l[i]*0.025;
  tdlldtime_t[i] = trig_time1 - tdlldtdc_t[i]*0.025;
  tdlldtime_w[i] = - tdlldtime_t[i] + tdlldtime_l[i];
//RD
  tdlrdtime_l[i] = trig_time1 - tdlrdtdc_l[i]*0.025;
  tdlrdtime_t[i] = trig_time1 - tdlrdtdc_t[i]*0.025;
  tdlrdtime_w[i] = - tdlrdtime_t[i] + tdlrdtime_l[i];

   if(i==3){//L4U,L4D
   tdllutime_l[i] = trig_time2 - tdllutdc_l[i]*0.025;   // time (leading)
   tdllutime_t[i] = trig_time2 - tdllutdc_t[i]*0.025;   // time (trailing)
   tdllutime_w[i] = - tdllutime_t[i] + tdllutime_l[i];   // time (width)
   tdlldtime_l[i] = trig_time2 - tdlldtdc_l[i]*0.025;
   tdlldtime_t[i] = trig_time2 - tdlldtdc_t[i]*0.025;
   tdlldtime_w[i] = - tdlldtime_t[i] + tdlldtime_l[i];
   }
   if(i==5){//R6D
   tdlrdtime_l[i] = trig_time2 - tdlrdtdc_l[i]*0.025;
   tdlrdtime_t[i] = trig_time2 - tdlrdtdc_t[i]*0.025;
   tdlrdtime_w[i] = - tdlrdtime_t[i] + tdlrdtime_l[i];
   }
   if(i==7||i==8){//R8U,R8D,R9U,R9D
   tdlrutime_l[i] = trig_time2 - tdlrutdc_l[i]*0.025;
   tdlrutime_t[i] = trig_time2 - tdlrutdc_t[i]*0.025;
   tdlrutime_w[i] = - tdlrutime_t[i] + tdlrutime_l[i];
   tdlrdtime_l[i] = trig_time2 - tdlrdtdc_l[i]*0.025;
   tdlrdtime_t[i] = trig_time2 - tdlrdtdc_t[i]*0.025;
   tdlrdtime_w[i] = - tdlrdtime_t[i] + tdlrdtime_l[i];
   }
  //mean time
  tdllmtime[i] = (tdllutime_l[i] + tdlldtime_l[i])/2.;   // time (mean)
  tdlrmtime[i] = (tdlrutime_l[i] + tdlrdtime_l[i])/2.;   // time (mean)
  //if(tdllutdc_l[i]>0)cout<<" time = "<<tdllutime_l[i]<<", width = "<<tdllutime_w[i]<<endl;
  }

  for(int i=0;i<NTB;i++){
  tagbtime_l[i] = tagb_trig_time - tagbtdc_l[i]*0.025;//lead
  tagbtime_t[i] = tagb_trig_time - tagbtdc_t[i]*0.025;//trail
  tagbtime_w[i] = tagbtime_l - tagbtime_t;//width
  }

/////////
//Fill
////////
//TDL raw data
  for(int i=0; i<NIH; i++){
  //leading
  h1_luTDC_l[i] -> Fill(tdllutdc_l[i]);
  h1_ldTDC_l[i] -> Fill(tdlldtdc_l[i]);
  h1_ruTDC_l[i] -> Fill(tdlrutdc_l[i]);
  h1_rdTDC_l[i] -> Fill(tdlrdtdc_l[i]);
  //trail
  h1_luTDC_t[i] -> Fill(tdllutdc_t[i]);
  h1_ldTDC_t[i] -> Fill(tdlldtdc_t[i]);
  h1_ruTDC_t[i] -> Fill(tdlrutdc_t[i]);
  h1_rdTDC_t[i] -> Fill(tdlrdtdc_t[i]);

  h1_lutime_w[i] -> Fill(tdllutime_w[i]);
  h1_ldtime_w[i] -> Fill(tdlldtime_w[i]);
  h1_rutime_w[i] -> Fill(tdlrutime_w[i]);
  h1_rdtime_w[i] -> Fill(tdlrdtime_w[i]);
  //h2_luTDCl_evnum[i] -> Fill(tdllutdc_l[i], evnum);
  //h2_ldTDCl_evnum[i] -> Fill(tdlldtdc_l[i], evnum);
  //h2_ruTDCl_evnum[i] -> Fill(tdlrutdc_l[i], evnum);
  //h2_rdTDCl_evnum[i] -> Fill(tdlrdtdc_l[i], evnum);
  }

  h2_daq0trig_evnum -> Fill( evnum, tdllutdc_l[0]);
  h2_daq1trig_evnum -> Fill( evnum, tdlrutdc_l[0]);
  for(int i=0; i<NTB; i++){
   for(int j=0;j<NIH;j++){
    if(tagbtdc_l[i]>0&&tdllutdc_l[j]>0)h_TagTDLLU[i][j] -> Fill(tagbtime_l[i] - tdllutime_l[j]);
    if(tagbtdc_l[i]>0&&tdlldtdc_l[j]>0)h_TagTDLLD[i][j] -> Fill(tagbtime_l[i] - tdlldtime_l[j]);
    if(tagbtdc_l[i]>0&&tdllutdc_l[j]>0&&tdlldtdc_l[j]>0){
                                       h_TagTDLL[i][j] -> Fill(tagbtime_l[i] - tdllmtime[j]);
                                       h_TagSumTDLL[j] -> Fill(tagbtime_l[i] - tdllmtime[j]);}
    if(tagbtdc_l[i]>0&&tdlrutdc_l[j]>0)h_TagTDLRU[i][j] -> Fill(tagbtime_l[i] - tdlrutime_l[j]);
    if(tagbtdc_l[i]>0&&tdlrdtdc_l[j]>0)h_TagTDLRD[i][j] -> Fill(tagbtime_l[i] - tdlrdtime_l[j]);
    if(tagbtdc_l[i]>0&&tdlrutdc_l[j]>0&&tdlrdtdc_l[j]>0){
                                       h_TagTDLR[i][j] -> Fill(tagbtime_l[i] - tdlrmtime[j]);
                                       h_TagSumTDLR[j] -> Fill(tagbtime_l[i] - tdlrmtime[j]);}
   }//TDL
  }//TagB



  for(int i=0; i<NIH; i++){
	  if(trig_time1!=0){
		  //h2_luWvsToF[i] -> Fill(tdllutime_l[i] - 810,tdllutime_w[i]);
		  h2_luWvsToF[i] -> Fill(tdllmtime[i],tdllutime_w[i]);
		  h2_ldWvsToF[i] -> Fill(tdllmtime[i],tdlldtime_w[i]);
		  h2_ruWvsToF[i] -> Fill(tdlrmtime[i],tdlrutime_w[i]);
		  h2_rdWvsToF[i] -> Fill(tdlrmtime[i],tdlrdtime_w[i]);
	  }
	  }

	//  if(//trig_time1>0&&
	//		  tdllutime_w[1]>0&&tdllutime_w[1]<10000&&
	//		  tdlldtime_w[1]>0&&tdlldtime_w[1]<10000&&
	//		  tdlrutime_w[2]>0&&tdlrutime_w[2]<10000&&
	//		  tdlrdtime_w[2]>0&&tdlrdtime_w[2]<10000
	//	)hToF    -> Fill(tdllmtime[1] - tdlrmtime[2]);
	  if(//trig_time1>0&&  Tekito counter v.s. TDLR4
			  tdlrutime_w[3]>Width_R4Umin &&tdlrutime_w[3]<Width_R4Umax &&
			  tdlrdtime_w[3]>Width_R4Dmin &&tdlrdtime_w[3]<Width_R4Dmax &&
			  tdlrutime_w[9]>Width_R10Umin&&tdlrutime_w[9]<Width_R10Umax&&
			  tdlrdtime_w[9]>Width_R10Dmin&&tdlrdtime_w[9]<Width_R10Dmax
		){
          hToF_wt -> Fill(tdlrmtime[9] - tdlrmtime[3]);
h2R4R10ToF_R4UWidth  -> Fill(tdlrmtime[9] - tdlrmtime[3],tdlrutime_w[3]);
h2R4R10ToF_R4DWidth  -> Fill(tdlrmtime[9] - tdlrmtime[3],tdlrdtime_w[3]);
h2R4R10ToF_R10UWidth -> Fill(tdlrmtime[9] - tdlrmtime[3],tdlrutime_w[9]);
h2R4R10ToF_R10DWidth -> Fill(tdlrmtime[9] - tdlrmtime[3],tdlrdtime_w[9]);
//cout<<"good event :"<<evnum<<endl;
		  //h2_luWvsToF_cut[1] ->Fill(tdllmtime[1] - tdlrmtime[2],tdllutime_w[1]);
		  //h2_ldWvsToF_cut[1] ->Fill(tdllmtime[1] - tdlrmtime[2],tdlldtime_w[1]);
		  //h2_ruWvsToF_cut[2] ->Fill(tdllmtime[1] - tdlrmtime[2],tdlrutime_w[2]);
		  //h2_rdWvsToF_cut[2] ->Fill(tdllmtime[1] - tdlrmtime[2],tdlrdtime_w[2]);
	  }
  }//event for loop


    ofp->Write(); ofp->Close();

}

