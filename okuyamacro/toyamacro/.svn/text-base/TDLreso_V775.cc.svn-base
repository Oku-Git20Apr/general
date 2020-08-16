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


gStyle->SetOptStat("iMen");

TH1F *hToF_wt = new TH1F("hToF_wt","hToF_wt",100,-3,3);
SetTH1(hToF_wt, "ToF R4-R10 ", "ToF[ns]", "count", 1, 0, 0);
TH2F *hQ10uToF = new TH2F("hQ10uToF","hQ10uToF",1000,-3,3,100,0,400);
TH2F *hQ10dToF = new TH2F("hQ10dToF","hQ10dToF",1000,-3,3,100,0,400);
TH2F *hQ4uToF  = new TH2F("hQ4uToF" ,"hQ4uToF" ,1000,-3,3,100,0,400);
TH2F *hQ4dToF  = new TH2F("hQ4dToF" ,"hQ4dToF" ,1000,-3,3,100,0,400);
SetTH2(hQ10uToF, "QDC 10U vs ToF", "ToF [ns]", "QDC 10U");
SetTH2(hQ10dToF, "QDC 10D vs ToF", "ToF [ns]", "QDC 10D");
SetTH2(hQ4uToF , "QDC  4U vs ToF", "ToF [ns]", "QDC  4U");
SetTH2(hQ4dToF , "QDC  4D vs ToF", "ToF [ns]", "QDC  4D");
double t_offset = 24.5;
TCut qcut10u = " ihruadc[9]>0 && ihruadc[9]<4000 && ihrutdc[9]>0";
TCut qcut10d = " ihrdadc[9]>0 && ihrdadc[9]<4000 && ihrdtdc[9]>0";
TCut qcut4u  = " ihruadc[3]>0 && ihruadc[3]<4000 && ihrutdc[3]>0";
TCut qcut4d  = " ihrdadc[3]>0 && ihrdadc[3]<4000 && ihrdtdc[3]>0";
tree->Project("hToF_wt",Form("( (ihrutdc[9]+ihrdtdc[9])-(ihrutdc[2]+ihrdtdc[2]) )*0.035/2-%lf",t_offset)
              ,qcut10u && qcut10d && qcut4u && qcut4d);//Form("ihlutdc[%d]<=0",i)
tree->Project("hQ10uToF",Form("ihruadc[9]:( (ihrutdc[9]+ihrdtdc[9])-(ihrutdc[2]+ihrdtdc[2]) )*0.035/2-%lf",t_offset)
              ,qcut10u && qcut10d && qcut4u && qcut4d);//Form("ihlutdc[%d]<=0",i)
tree->Project("hQ10dToF",Form("ihrdadc[9]:( (ihrutdc[9]+ihrdtdc[9])-(ihrutdc[2]+ihrdtdc[2]) )*0.035/2-%lf",t_offset)
              ,qcut10u && qcut10d && qcut4u && qcut4d);//Form("ihlutdc[%d]<=0",i)
tree->Project("hQ4uToF" ,Form("ihruadc[2]:( (ihrutdc[9]+ihrdtdc[9])-(ihrutdc[2]+ihrdtdc[2]) )*0.035/2-%lf",t_offset)
              ,qcut10u && qcut10d && qcut4u && qcut4d);//Form("ihlutdc[%d]<=0",i)
tree->Project("hQ4dToF" ,Form("ihrdadc[2]:( (ihrutdc[9]+ihrdtdc[9])-(ihrutdc[2]+ihrdtdc[2]) )*0.035/2-%lf",t_offset)
              ,qcut10u && qcut10d && qcut4u && qcut4d);//Form("ihlutdc[%d]<=0",i)

double sigma_t=100;
TText *Text_sigma_t = new TText(0,0,Form("%lf ps",sigma_t));
//////////////
//histograms//
//////////////

cout<<"start drawing!"<<endl;
TCanvas *c[14];
for(int i=0;i<14;i++){
c[i] = new TCanvas(Form("c%d",i+1), Form("canvas%d",i+1) , 800,700);}

c[0]->Clear();
hToF_wt->Draw();
hToF_wt ->Fit("gaus");
c[1]->Divide(2,2);
c[1]->cd(1);hQ10uToF->Draw("colz");
c[1]->cd(2);hQ10dToF->Draw("colz");
c[1]->cd(3);hQ4uToF ->Draw("colz");
c[1]->cd(4);hQ4dToF ->Draw("colz");
                                                                                                                                 
  string ofname_pdf = ofname;
  ofname_pdf.erase(ofname_pdf.size()-4);
  ofname_pdf.append(".pdf");

  c[0] ->Print(Form("%s[",ofname_pdf.c_str()));
  c[0] ->Print(Form("%s" ,ofname_pdf.c_str()));
  c[1] ->Print(Form("%s" ,ofname_pdf.c_str()));
  c[1] ->Print(Form("%s]" ,ofname_pdf.c_str()));


}

