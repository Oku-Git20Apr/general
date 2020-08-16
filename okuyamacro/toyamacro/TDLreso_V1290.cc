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

  int ENum;
  int tdllutdc_l[NIH], tdlrutdc_l[NIH], tdlldtdc_l[NIH],   tdlrdtdc_l[NIH];   // TDC (leading)
  int tdllutdc_t[NIH], tdlrutdc_t[NIH], tdlldtdc_t[NIH],   tdlrdtdc_t[NIH];   // TDC (trailing)
  double tdllutime_l[NIH], tdlrutime_l[NIH], tdlldtime_l[NIH],   tdlrdtime_l[NIH];   // time (leading)
  double tdllutime_t[NIH], tdlrutime_t[NIH], tdlldtime_t[NIH],   tdlrdtime_t[NIH];   // time (trailing)
  double tdllutime_w[NIH], tdlrutime_w[NIH], tdlldtime_w[NIH],   tdlrdtime_w[NIH];   // time (width)
  double tdllmtime[NIH], tdlrmtime[NIH];   // time (width)
  double trig_time1, trig_time2;//trigger time
  double t_offset = 24.5;
  ENum = tree->GetEntries();
  if(MaxNum>0 && MaxNum<ENum) ENum = MaxNum;
  tree->SetBranchAddress("tdllutdc_l"   , tdllutdc_l );
  tree->SetBranchAddress("tdlrutdc_l"   , tdlrutdc_l );
  tree->SetBranchAddress("tdlldtdc_l"   , tdlldtdc_l );
  tree->SetBranchAddress("tdlrdtdc_l"   , tdlrdtdc_l );
  tree->SetBranchAddress("tdllutdc_t"   , tdllutdc_t );
  tree->SetBranchAddress("tdlrutdc_t"   , tdlrutdc_t );
  tree->SetBranchAddress("tdlldtdc_t"   , tdlldtdc_t );
  tree->SetBranchAddress("tdlrdtdc_t"   , tdlrdtdc_t );

gStyle->SetOptStat("iMen");

TH1F *hToF = new TH1F("hToF","hToF",25,-1,1);
SetTH1(hToF, "ToF L2-R3 ", "ToF[ns]", "count", 1, 0, 0);
TH1F *hToF_wt = new TH1F("hToF_wt","hToF_wt",25,-1,1);
SetTH1(hToF_wt, "ToF L2-R3 ", "ToF[ns]", "count", 1, 0, 0);
TH2F *h2_luWvsToF[NIH],*h2_luWvsToF_cut[NIH];
TH2F *h2_ldWvsToF[NIH],*h2_ldWvsToF_cut[NIH];
TH2F *h2_ruWvsToF[NIH],*h2_ruWvsToF_cut[NIH];
TH2F *h2_rdWvsToF[NIH],*h2_rdWvsToF_cut[NIH];
for(int i=0;i<NIH;i++){
h2_luWvsToF[i] = new TH2F(Form("h2_luWvsToF%d",i+1),Form("h2_luWvsToF%d",i+1),1000,780,850,50,0,200);
h2_ldWvsToF[i] = new TH2F(Form("h2_ldWvsToF%d",i+1),Form("h2_ldWvsToF%d",i+1),1000,780,850,50,0,200);
h2_ruWvsToF[i] = new TH2F(Form("h2_ruWvsToF%d",i+1),Form("h2_ruWvsToF%d",i+1),1000,780,850,50,0,200);
h2_rdWvsToF[i] = new TH2F(Form("h2_rdWvsToF%d",i+1),Form("h2_rdWvsToF%d",i+1),1000,780,850,50,0,200);
SetTH2(h2_luWvsToF[i], Form("Width LUseg%d vs ToF",i+1), "ToF [ns]", Form("width LU%d",i+1));
SetTH2(h2_ldWvsToF[i], Form("Width LDseg%d vs ToF",i+1), "ToF [ns]", Form("width LD%d",i+1));
SetTH2(h2_ruWvsToF[i], Form("Width RUseg%d vs ToF",i+1), "ToF [ns]", Form("width RU%d",i+1));
SetTH2(h2_rdWvsToF[i], Form("Width RDseg%d vs ToF",i+1), "ToF [ns]", Form("width RD%d",i+1));
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

  trig_time1=-9999;
  trig_time2=-9999;//trigger time

if(tdllutdc_l[0]>0) trig_time1 = tdllutdc_l[0]*0.025;
if(tdlrutdc_l[0]>0) trig_time2 = tdlrutdc_l[0]*0.025;//QTC3(MRPC datta yatsu)
//cout<<"trigger timing ="<<trig_time1<<endl;

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
  tdllmtime[i] = (tdllutime_l[i] + tdlldtime_l[i])/2.;   // time (width)
  tdlrmtime[i] = (tdlrutime_l[i] + tdlrdtime_l[i])/2.;   // time (width)
//cout<<" time = "<<tdllutime_l[i]<<", width = "<<tdllutime_w[i]<<endl;
 }


  for(int i=0; i<NIH; i++){
if(trig_time1>0){
//h2_luWvsToF[i] -> Fill(tdllutime_l[i] - 810,tdllutime_w[i]);
h2_luWvsToF[i] -> Fill(tdllmtime[i],tdllutime_w[i]);
h2_ldWvsToF[i] -> Fill(tdllmtime[i],tdlldtime_w[i]);
h2_ruWvsToF[i] -> Fill(tdlrmtime[i],tdlrutime_w[i]);
h2_rdWvsToF[i] -> Fill(tdlrmtime[i],tdlrdtime_w[i]);
if(i==6&&tdlrdtdc_l[i]>0&&tdlrutdc_l[i]>0)cout<<"mean time, width :"<<tdlrmtime[i]<<"  "<<tdlrdtime_w[i]<<endl;
  }
 }

if(//trig_time1>0&&
  tdllutime_w[1]>0&&tdllutime_w[1]<10000&&
  tdlldtime_w[1]>0&&tdlldtime_w[1]<10000&&
  tdlrutime_w[2]>0&&tdlrutime_w[2]<10000&&
  tdlrdtime_w[2]>0&&tdlrdtime_w[2]<10000
)hToF    -> Fill(tdllmtime[1] - tdlrmtime[2]);
if(//trig_time1>0&&
  tdllutime_w[1]>110&&tdllutime_w[1]<125&&
  tdlldtime_w[1]>130&&tdlldtime_w[1]<145&&
  tdlrutime_w[2]>100&&tdlrutime_w[2]<115&&
  tdlrdtime_w[2]>100&&tdlrdtime_w[2]<115
){hToF_wt -> Fill(tdllmtime[1] - tdlrmtime[2]);
h2_luWvsToF_cut[1] ->Fill(tdllmtime[1] - tdlrmtime[2],tdllutime_w[1]);
h2_ldWvsToF_cut[1] ->Fill(tdllmtime[1] - tdlrmtime[2],tdlldtime_w[1]);
h2_ruWvsToF_cut[2] ->Fill(tdllmtime[1] - tdlrmtime[2],tdlrutime_w[2]);
h2_rdWvsToF_cut[2] ->Fill(tdllmtime[1] - tdlrmtime[2],tdlrdtime_w[2]);
 }
}//event for loop

double sigma_t=100;
TText *Text_sigma_t = new TText(0,0,Form("%lf ps",sigma_t));
//////////////
//histograms//
//////////////

cout<<"start drawing!"<<endl;
TCanvas *c[14];
for(int i=0;i<14;i++){
c[i] = new TCanvas(Form("c%d",i+1), Form("canvas%d",i+1) , 800,1900);}

c[0]->Clear();
c[0]->Divide(1,2);
c[0]->cd(1);
hToF ->Draw();
hToF ->Fit("gaus");
c[0]->cd(2);
hToF_wt->Draw();
hToF_wt ->Fit("gaus");

c[1]->Divide(2,5);
c[1]->cd(1) ;gPad->SetLogz(1);h2_luWvsToF[1] ->Draw("colz");
c[1]->cd(2) ;gPad->SetLogz(1);h2_ldWvsToF[1] ->Draw("colz");
c[1]->cd(3) ;gPad->SetLogz(1);h2_ruWvsToF[2] ->Draw("colz");
c[1]->cd(4) ;gPad->SetLogz(1);h2_rdWvsToF[2] ->Draw("colz");

c[1]->cd(7) ;gPad->SetLogz(1);h2_luWvsToF_cut[1] ->Draw("colz");
c[1]->cd(8) ;gPad->SetLogz(1);h2_ldWvsToF_cut[1] ->Draw("colz");
c[1]->cd(9) ;gPad->SetLogz(1);h2_ruWvsToF_cut[2] ->Draw("colz");
c[1]->cd(10);gPad->SetLogz(1);h2_rdWvsToF_cut[2] ->Draw("colz");


c[2]->Divide(2,5);
for(int i=0;i<10;i++){
c[2]->cd(i+1) ;gPad->SetLogz(1);h2_luWvsToF[i] ->Draw("colz");
}

c[3]->Divide(2,5);
for(int i=0;i<10;i++){
c[3]->cd(i+1) ;gPad->SetLogz(1);h2_ldWvsToF[i] ->Draw("colz");
}

c[4]->Divide(2,5);
for(int i=0;i<10;i++){
c[4]->cd(i+1) ;gPad->SetLogz(1);h2_ruWvsToF[i] ->Draw("colz");
}

c[5]->Divide(2,5);
for(int i=0;i<10;i++){
c[5]->cd(i+1) ;gPad->SetLogz(1);h2_rdWvsToF[i] ->Draw("colz");
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
  c[5] ->Print(Form("%s]" ,ofname_pdf.c_str()));


}

