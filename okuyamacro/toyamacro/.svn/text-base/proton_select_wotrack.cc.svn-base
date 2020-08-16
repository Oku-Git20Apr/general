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

#include "Tree.h"
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



class proton : public Tree
{
 public:
         proton();
        ~proton();
  void makehist();
  void loop();
  void draw(); 
  void savecanvas(); 
  void SetRoot(string ifname); 
  void SetPdfFilename(string ifname); 
  void SetMaxEvent( int N )  { ENumMax = N; }

  private:

    string pdf_name;
    int ENumMax;
    int ENum;

    TH2F *h2_detof_ihl[10], *h2_detof_ihr[10];
    TH2F *h2_detof_ohvl[12], *h2_detof_ohvr[12], *h2_detof_ohhl[9], *h2_detof_ohhr[9];
    TH2F *h2_OHVLdEdx_IHL_OHVLToF[5][8], *h2_OHVRdEdx_IHR_OHVRToF[5][8];
    TH2F *h2_IHLdEdx_IHL_OHVLToF[5][8], *h2_IHRdEdx_IHR_OHVRToF[5][8];

    TH2F *h2_IHL2dEdx_IHL2_OHVL4ToF;
    TH2F *h2_IHL2dEdx_IHL2_OHVL6ToF;
    TH2F *h2_OHVL4dEdx_IHL2_OHVL4ToF;
    TH2F *h2_OHVL6dEdx_IHL2_OHVL6ToF;
    TCanvas *cc[12];

};
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proton::proton(){
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
      


}
//
proton::~proton(){
}
////////////////////////////////////////////////////////////////////////////
void proton::makehist(){

  for(int i=0;i<10;i++){
    h2_detof_ihl[i] = new TH2F(Form("h2_detof_ihl%d",i+1),Form("h2_detof_ihl%d",i+1),  100, 0.,  20.,  110,  0.,  10.);
    h2_detof_ihr[i] = new TH2F(Form("h2_detof_ihr%d",i+1),Form("h2_detof_ihr%d",i+1),  100, 0.,  20.,  110,  0.,  10.);
    SetTH2(h2_detof_ihl[i],Form("IHL%d dE vs ToF",i+1), "ToF[ns]","dE[MeV]");
    SetTH2(h2_detof_ihr[i],Form("IHR%d dE vs ToF",i+1), "ToF[ns]","dE[MeV]");
  }
  for(int i=0;i<12;i++){
    h2_detof_ohvl[i] = new TH2F(Form("h2_detof_ohvl%d",i+1),Form("h2_detof_ohvl%d",i+1),  100, 0.,  20.,  110,  0.,  10.);
    h2_detof_ohvr[i] = new TH2F(Form("h2_detof_ohvr%d",i+1),Form("h2_detof_ohvr%d",i+1),  100, 0.,  20.,  110,  0.,  10.);
    SetTH2(h2_detof_ohvl[i],Form("OHVL%d dE vs ToF",i+1), "ToF[ns]","dE[MeV]");
    SetTH2(h2_detof_ohvr[i],Form("OHVR%d dE vs ToF",i+1), "ToF[ns]","dE[MeV]");
  }
  for(int i=0;i<12;i++){
    h2_detof_ohhl[i] = new TH2F(Form("h2_detof_ohhl%d",i+1),Form("h2_detof_ohhl%d",i+1),  100, 0.,  20.,  110,  0.,  10.);
    h2_detof_ohhr[i] = new TH2F(Form("h2_detof_ohhr%d",i+1),Form("h2_detof_ohhr%d",i+1),  100, 0.,  20.,  110,  0.,  10.);
    SetTH2(h2_detof_ohhl[i],Form("OHHL%d dE vs ToF",i+1), "ToF[ns]","dE[MeV]");
    SetTH2(h2_detof_ohhr[i],Form("OHHR%d dE vs ToF",i+1), "ToF[ns]","dE[MeV]");
  }

for(int i=0;i<5;i++){//IHL seg
 for(int j=0;j<8;j++){//OHVL seg
   h2_OHVLdEdx_IHL_OHVLToF[i][j] = new TH2F(Form("h2_OHVL%ddEdx_IHL%d_OHVL%dToF",j+1,i+1,j+1),Form("h2_OHVL%ddEdx_IHL%d_OHVL%dToF",j+1,i+1,j+1), 340, 0., 14, 110, 0.,  22.);
   h2_IHLdEdx_IHL_OHVLToF[i][j]  = new TH2F(Form("h2_IHL%ddEdx_IHL%d_OHVL%dToF" ,i+1,i+1,j+1),Form("h2_IHL%ddEdx_IHL%d_OHVL%dToF" ,i+1,i+1,j+1), 340, 0., 14, 110, 0.,  22.);
  
   SetTH2(h2_OHVLdEdx_IHL_OHVLToF[i][j],Form("OHVL%d dEdx vs OHVL%d - IHL%d ToF"      ,j+1,j+1,i+1), Form("OHVL%d - IHL%d ToF [ns]",j+1,i+1), Form("OHVL%d dEdx [MeV/cm]",j+1) );
   SetTH2(h2_IHLdEdx_IHL_OHVLToF[i][j] ,Form("IHL%d  dE/dx vs OHVL%d - IHL%d ToF"     ,i+1,j+1,i+1), Form("OHVL%d - IHL%d ToF [ns]",j+1,i+1), Form("IHL%d dE/dx [MeV/cm]",i+1) );
  
   h2_OHVRdEdx_IHR_OHVRToF[i][j] = new TH2F(Form("h2_OHVR%ddEdx_IHR%d_OHVR%dToF",j+1,i+1,j+1),Form("h2_OHVR%ddEdx_IHR%d_OHVR%dToF",j+1,i+1,j+1), 340, 0., 14, 110, 0.,  22.);
   h2_IHRdEdx_IHR_OHVRToF[i][j]  = new TH2F(Form("h2_IHR%ddEdx_IHR%d_OHVR%dToF" ,i+1,i+1,j+1),Form("h2_IHR%ddEdx_IHR%d_OHVR%dToF" ,i+1,i+1,j+1), 340, 0., 14, 110, 0.,  22.);
  
   SetTH2(h2_OHVRdEdx_IHR_OHVRToF[i][j],Form("OHVR%d dEdx vs OHVR%d - IHR%d ToF"      ,j+1,j+1,i+1), Form("OHVR%d - IHR%d ToF [ns]",j+1,i+1), Form("OHVR%d dEdx [MeV/cm]",j+1) );
   SetTH2(h2_IHRdEdx_IHR_OHVRToF[i][j] ,Form("IHR%d  dE/dx vs OHVR%d - IHR%d ToF"     ,i+1,j+1,i+1), Form("OHVR%d - IHR%d ToF [ns]",j+1,i+1), Form("IHR%d dE/dx [MeV/cm]",i+1) );
   }
  }
    h2_IHL2dEdx_IHL2_OHVL4ToF  = new TH2F("h2_IHL2dEdx_IHL2_OHVL4ToF" ,"h2_IHL2dEdx_IHL2_OHVL4ToF" , 340, 0., 14, 110, 0.,  22.);
    h2_IHL2dEdx_IHL2_OHVL6ToF  = new TH2F("h2_IHL2dEdx_IHL2_OHVL6ToF" ,"h2_IHL2dEdx_IHL2_OHVL6ToF" , 340, 0., 14, 110, 0.,  22.);
    h2_OHVL4dEdx_IHL2_OHVL4ToF = new TH2F("h2_OHVL4dEdx_IHL2_OHVL4ToF","h2_OHVL4dEdx_IHL2_OHVL4ToF", 340, 0., 14, 110, 0.,  22.);
    h2_OHVL6dEdx_IHL2_OHVL6ToF = new TH2F("h2_OHVL6dEdx_IHL2_OHVL6ToF","h2_OHVL6dEdx_IHL2_OHVL6ToF", 340, 0., 14, 110, 0.,  22.);
}
////////////////////////////////////////////////////////////////////////////
void proton::loop(){

bool ihlseg_flag[10], ihrseg_flag[10], ohvlseg_flag[12], ohvrseg_flag[12];
  bool IHL_OHVLhit[5][5];
  bool IHR_OHVRhit[5][5];
bool IHL2_OHVL4_p, IHL2_OHVL6_p;

  if( ENumMax>0 && ENumMax<ENum) ENum = ENumMax;

  for(int n=0;n<ENum;n++){
    //tr->tree->GetEntry(n);
    tree->GetEntry(n);
    if(n%100000==0) cout<<n<<" / "<<ENum<<endl;

  for(int i=0;i<10;i++){
    ihlseg_flag[i]=ihrseg_flag[i]=false;
    if(ihlutdc[i]>0&&ihldtdc[i]>0)ihlseg_flag[i]=true;
    if(ihrutdc[i]>0&&ihrdtdc[i]>0)ihrseg_flag[i]=true;
  }

  for(int i=0;i<12;i++){
    ohvlseg_flag[i]=ohvrseg_flag[i]=false;
    if(ohvlutdc[i]>0&&ohvldtdc[i]>0)ohvlseg_flag[i]=true;
    if(ohvrutdc[i]>0&&ohvrdtdc[i]>0)ohvrseg_flag[i]=true;
  }

  IHL2_OHVL4_p = IHL2_OHVL6_p = false;
for(int i=0;i<5;i++){//IHL seg
 for(int j=0;j<8;j++){//OHVL seg
  IHL_OHVLhit[i][j]=false;
  IHR_OHVRhit[i][j]=false;
  }
 }
//IH, OH segment flags//
for(int i=0;i<5;i++){//IH seg
 for(int j=0;j<8;j++){//OHV seg
if(ihlutdc[i]>0 && ihldtdc[i]>0 && ohvlutdc[j]>0 && ohvldtdc[j]>0)IHL_OHVLhit[i][j]=true;
if(ihrutdc[i]>0 && ihrdtdc[i]>0 && ohvrutdc[j]>0 && ohvrdtdc[j]>0)IHR_OHVRhit[i][j]=true;
  }
 }


  for(int i=0;i<10;i++){
    for(int j=0;j<12;j++){
      if(ihlseg_flag[i]&&ohvlseg_flag[j]){//IHL&OHVL
       h2_detof_ihl[i]  -> Fill(ohvlmctime[j]-ihlmctime[i],ihlde[i]);
       h2_detof_ohvl[j] -> Fill(ohvlmctime[j]-ihlmctime[i],ohvlde[i]);
      }        

      if(ihlseg_flag[i]&&ohvrseg_flag[j]){//IHL&OHVR
       h2_detof_ihl[i]  -> Fill(ohvrmctime[j]-ihlmctime[i],ihlde[i]);
       h2_detof_ohvr[j] -> Fill(ohvrmctime[j]-ihlmctime[i],ohvrde[i]);
      }        

      if(ihrseg_flag[i]&&ohvlseg_flag[j]){//IHR&OHVL
       h2_detof_ihr[i]  -> Fill(ohvlmctime[j]-ihrmctime[i],ihrde[i]);
       h2_detof_ohvl[j] -> Fill(ohvlmctime[j]-ihrmctime[i],ohvlde[i]);
      }        

      if(ihrseg_flag[i]&&ohvrseg_flag[j]){//IHR&OHVR
       h2_detof_ihr[i]  -> Fill(ohvrmctime[j]-ihrmctime[i],ihrde[i]);
       h2_detof_ohvr[j] -> Fill(ohvrmctime[j]-ihrmctime[i],ohvrde[i]);
      }        
    }
  }

  for(int i=0;i<5;i++){//IHL seg
   for(int j=0;j<8;j++){//OHVL seg
     if(IHL_OHVLhit[i][j]){
        h2_OHVLdEdx_IHL_OHVLToF[i][j] ->Fill(ohvlmctime[j]-ihlmctime[i], ohvlde[j]);
        h2_IHLdEdx_IHL_OHVLToF[i][j]  ->Fill(ohvlmctime[j]-ihlmctime[i], ihlde[i] );
        if(i==1&&j==3&&ihlde[i]>3&&ohvlde[j]>7&&(ohvlmctime[j]-ihlmctime[i])>6){IHL2_OHVL4_p=true;h2_IHL2dEdx_IHL2_OHVL4ToF->Fill(ohvlmctime[j]-ihlmctime[i],ihlde[i]); h2_OHVL4dEdx_IHL2_OHVL4ToF->Fill(ohvlmctime[j]-ihlmctime[i],ohvlde[i]);}
        if(i==1&&j==5&&ihlde[i]>3&&ohvlde[j]>7&&(ohvlmctime[j]-ihlmctime[i])>6){IHL2_OHVL6_p=true;h2_IHL2dEdx_IHL2_OHVL6ToF->Fill(ohvlmctime[j]-ihlmctime[i],ihlde[i]); h2_OHVL6dEdx_IHL2_OHVL6ToF->Fill(ohvlmctime[j]-ihlmctime[i],ohvlde[i]);}
     }//if IHL_OHVL hit i,j
     if(IHR_OHVRhit[i][j]){
        h2_OHVRdEdx_IHR_OHVRToF[i][j] ->Fill(ohvrmctime[j]-ihrmctime[i], ohvrde[j]);
        h2_IHRdEdx_IHR_OHVRToF[i][j]  ->Fill(ohvrmctime[j]-ihrmctime[i], ihrde[i] );
     }//if IHR_OHVR hit i,j
    }
   }

   }//event loop
cout<<"end of loop"<<endl;
}
////////////////////////////////////////////////////////////////////////////
void proton::draw(){
  for(int i=0;i<12;i++){
  cc[i]= new TCanvas(Form("c_%d",i+1),Form("c_%d",i+1),1400,900 );
  }

cc[0]->Divide(4,3);
  for(int i=0;i<10;i++){
  cc[0]->cd(i+1);gPad->SetLogz(1);h2_detof_ihl[i]->Draw("colz");
  }
cout<<"ihl finish"<<endl;

cc[1]->Divide(3,4);
  for(int i=0;i<12;i++){
  cc[1]->cd(i+1);gPad->SetLogz(1);h2_detof_ihr[i]->Draw("colz");
  }
cout<<"ihr finish"<<endl;

cc[2]->Divide(4,3);
  for(int i=0;i<12;i++){
  cc[2]->cd(i+1);gPad->SetLogz(1);h2_detof_ohvl[i]->Draw("colz");
  }
cout<<"ohvl finish"<<endl;

cc[3]->Divide(4,3);
  for(int i=0;i<12;i++){
  cc[3]->cd(i+1);gPad->SetLogz(1);h2_detof_ohvr[i]->Draw("colz");
  }

cc[4]->Divide(5,8);
for(int i=0;i<5;i++){//IHL seg
 for(int j=0;j<8;j++){//OHVL seg
cc[4]->cd(5*i+j+1); gPad->SetLogz(1);       h2_IHLdEdx_IHL_OHVLToF[i][j]  ->Draw("colz");
  }
 }

cc[5]->Divide(5,8);
for(int i=0;i<5;i++){//IHL seg
 for(int j=0;j<8;j++){//OHVL seg
cc[5]->cd(5*i+j+1); gPad->SetLogz(1);       h2_OHVLdEdx_IHL_OHVLToF[i][j]  ->Draw("colz");
  }
 }

cc[6]->Divide(5,8);
for(int i=0;i<5;i++){//IHR seg
 for(int j=0;j<8;j++){//OHVR seg
cc[6]->cd(5*i+j+1); gPad->SetLogz(1);      h2_IHRdEdx_IHR_OHVRToF[i][j] ->Draw("colz");
  }
 }

cc[7]->Divide(5,8);
for(int i=0;i<5;i++){//IHR seg
 for(int j=0;j<8;j++){//OHVR seg
cc[7]->cd(5*i+j+1); gPad->SetLogz(1);      h2_OHVRdEdx_IHR_OHVRToF[i][j] ->Draw("colz");
  }
 }

cc[8]->Divide(2,2);
cc[8]->cd(1); gPad->SetLogz(1); h2_IHL2dEdx_IHL2_OHVL4ToF ->Draw("colz");
cc[8]->cd(2); gPad->SetLogz(1); h2_IHL2dEdx_IHL2_OHVL6ToF ->Draw("colz");
cc[8]->cd(3); gPad->SetLogz(1); h2_OHVL4dEdx_IHL2_OHVL4ToF->Draw("colz");
cc[8]->cd(4); gPad->SetLogz(1); h2_OHVL6dEdx_IHL2_OHVL6ToF->Draw("colz");
cout<<"ohvr finish"<<endl;
cout<<"finish drawing"<<endl;
}
////////////////////////////////////////////////////////////////////////////
void proton::savecanvas(){
  cc[0] ->Print(Form("%s[",pdf_name.c_str()) );
  cc[0] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[1] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[2] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[3] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[4] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[5] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[6] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[7] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[8] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[8] ->Print(Form("%s]",pdf_name.c_str()) );
cout<<pdf_name<<" saved!"<<endl;
}
////////////////////////////////////////////////////////////////////////////
void proton::SetPdfFilename(string ifname){
pdf_name = ifname;
} 
////////////////////////////////////////////////
void proton::SetRoot(string ifname)
{
cout<<"SetRoot"<<endl;
  add(ifname);
  readtree();
  ENum = GetEntries();
}
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "test.root";
  string pdfname = "protonlike_extract.pdf";
  int ch;
  int MaxNum = 0;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:p:n:bco"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;
    case 'p':
      output_flag = true;
      draw_flag = false;
      pdfname = optarg;
      cout<<"output pdf filename : "<<pdfname<<endl;
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
      cout<<"-f : input filename"<<endl;
      cout<<"-n : maximum number of analysed events"<<endl;
      cout<<"-p : output pdf file name"<<endl;
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
  proton *pro = new proton();
  pro->SetMaxEvent(MaxNum);
  pro->SetRoot(ifname);
  pro->makehist();
  pro->loop();
  pro->draw();
  pro->SetPdfFilename(pdfname);
  pro->savecanvas();
  delete pro;

  gSystem->Exit(1);
  theApp->Run();
  return 0;

} 
