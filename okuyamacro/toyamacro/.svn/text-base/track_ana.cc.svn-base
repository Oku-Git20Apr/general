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
#include "Track_tr.h"

#define Calibration

static const double PI = 4.0*atan(1.);
static const double mrad_to_deg = 1./1000*180./PI;
const double Mp = 938.272046;          // proton       mass (MeV/c2)
const double Mpi = 139.57018;          // charged pion mass (MeV/c2)
const double MK = 493.677;             // charged Kaon mass (MeV/c2)
const double c = 0.299792458;          // speed of light in vacuum (m/ns)
const int NIH   = 20;  // No. of IH
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
  h->GetYaxis()->SetTitleOffset(0.9);
  h->GetYaxis()->SetTitleSize(0.06);
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


class track_ana : public Track_tr
{
 public:
         track_ana();
        ~track_ana();
  void makehist(string ofname);
  void roop();
  void fit();
  void draw(); 
  void savecanvas(); 
  void SetRoot(string ifname); 
  void GetHist(string ifname); 
  void SetMaxEvent( int N )  { ENumMax = N; }
  void SetPdfFilename(string ifname); 

  private:
    //Track_tr *tr;
    TFile *ofp ;
    string pdf_name;
    int GetMaxEvent() { return ENumMax; }
    int ENumMax;
    int ENum;

    TH1F *h_ctime;
    TH1F *h_fl;
    TH2F *h_ctime_uwid;

    TH1F *h_msqr, *h_msqr_decut;
    TH2F *h2_pid_ohv;
    TH2F *h2_pid_ohh;
    TH2F *h2_pid_all;
    TH2F *h2_pid_decut;
    TH2F *h2_ohdedxbeta_all;
    TH2F *h2_ohdedxbeta_decut;
    TH2F *h2_ohdedxbeta_ohv;
    TH2F *h2_ohdedxbeta_ohh;
    TH2F *h2_ihdedxbeta_all;
    TH2F *h2_dedxbeta_ihl[10], *h2_dedxbeta_ihr[10];
    TH2F *h2_dedxbeta_ohvl[12], *h2_dedxbeta_ohvr[12], *h2_dedxbeta_ohhl[9], *h2_dedxbeta_ohhr[9];

    TH2F *h2_pid_ihl[10], *h2_pid_ihr[10];
    TH2F *h2_pid_ohvl[12], *h2_pid_ohvr[12], *h2_pid_ohhl[9], *h2_pid_ohhr[9];
    TH1F *h_fl_ohvl[12], *h_fl_ohvr[12], *h_fl_ohhl[9], *h_fl_ohhr[9];
    TH1F *h_ct_tb[40];
    TH2F *h2_ct_seg;

    TLine *line;
    int run_num;
TCanvas *c1,*c2,*c3,*c4,*c5;
TCanvas *cc[12];
TLatex *tex_fl[4][12];
TF1 *f_pi, *f_k, *f_p;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
track_ana::track_ana()
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
      

  for(int i=0;i<12;i++){
  cc[i]= new TCanvas(Form("c_%d",i+1),Form("c_%d",i+1),1400,800 );
  }

  line = new TLine(0,0,1,1);
  line ->SetLineColor(6);
  line ->SetLineWidth(1);
  line ->SetLineStyle(1);

f_pi = new TF1("f_pi","[0]/sqrt(x*x-1)", 1,7);
f_k  = new TF1("f_k" ,"[0]/sqrt(x*x-1)", 1,7);
f_p  = new TF1("f_p" ,"[0]/sqrt(x*x-1)", 1,7);
f_pi->SetParameter(0, 0.001*Mpi);
f_k ->SetParameter(0, 0.001*MK);
f_p ->SetParameter(0, 0.001*Mp);
f_pi->SetLineColor(6);f_pi->SetLineWidth(1.5);f_pi->SetNpx(1000);
f_k ->SetLineColor(3);f_k ->SetLineWidth(1.5);f_k ->SetNpx(1000);
f_p ->SetLineColor(2);f_p ->SetLineWidth(1.5);f_p ->SetNpx(1000);

}
////////////////////////////////////////////////////////////////////////////
track_ana::~track_ana(){
}
////////////////////////////////////////////////////////////////////////////
void track_ana::makehist(string ofname){
cout<<"makehist"<<endl;
  ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  ofp->mkdir("track"); ofp->cd("track");
  h_ctime    = new TH1F("h_ctime"   ,"h_ctime"   ,2400,-30,30);
  SetTH1(h_ctime  ,"Coincidence time","cointime (ns)","counts / 25ps");
  h_ctime_uwid  = new TH2F("h_ctime_uwid" ,"h_ctime_uwid" ,200,-10,10,200,0,50);
  SetTH2(h_ctime_uwid,"Cointime vs Width","cointime (ns)","TDLU width (ns)");

    h2_pid_all  = new TH2F("h2_pid_all"   ,"h2_pid_all"   , 500, 0, 7, 1000, -1, 1);
    h2_pid_decut= new TH2F("h2_pid_decut" ,"h2_pid_decut" , 500, 0, 7, 1000, -1, 1);
    SetTH2(h2_pid_all  ,"PID plot all"             ,"1/#beta","momentum[GeV/#it{c}]");
    SetTH2(h2_pid_decut,"PID plot all w/ dEdx cut" ,"1/#beta","momentum[GeV/#it{c}]");

  ofp->cd();
  ofp->mkdir("ih"); ofp->cd("ih");
  for(int i=0;i<10;i++){
    h2_pid_ihl[i] = new TH2F(Form("h2_pid_ihl%d",i+1) ,Form("h2_pid_ihl%d",i+1) , 500, 0, 9, 1000, -1, 1);
    h2_pid_ihr[i] = new TH2F(Form("h2_pid_ihr%d",i+1) ,Form("h2_pid_ihr%d",i+1) , 500, 0, 9, 1000, -1, 1);
    SetTH2(h2_pid_ihl[i],Form("PID plot ihL%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
    SetTH2(h2_pid_ihr[i],Form("PID plot ihR%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
  }

  ofp->cd();
  ofp->mkdir("oh"); ofp->cd("oh");

    h2_pid_ohv = new TH2F("h2_pid_ohv"    ,"h2_pid_ohv"   , 500, 0, 7, 1000, -1, 1);
    h2_pid_ohh = new TH2F("h2_pid_ohh"    ,"h2_pid_ohh"   , 500, 0, 7, 1000, -1, 1);
    h2_ohdedxbeta_all  = new TH2F("h2_ohdedxbeta_all"   ,"h2_ohdedxbeta_all"    , 500, 0, 2, 1000, 0, 9);
    h2_ohdedxbeta_decut= new TH2F("h2_ohdedxbeta_decut" ,"h2_ohdedxbeta_decut"  , 500, 0, 2, 1000, 0, 9);
    h2_ohdedxbeta_ohv  = new TH2F("h2_ohdedxbeta_ohv"   ,"h2_ohdedxbeta_ohv"    , 500, 0, 2, 1000, 0, 9);
    h2_ohdedxbeta_ohh  = new TH2F("h2_ohdedxbeta_ohh"   ,"h2_ohdedxbeta_ohh"    , 500, 0, 2, 1000, 0, 9);
    h2_ihdedxbeta_all  = new TH2F("h2_ihdedxbeta_all"   ,"h2_ihdedxbeta_all"    , 500, 0, 2, 1000, 0, 9);
    SetTH2(h2_pid_ohv  ,"PID plot OHV"             ,"1/#beta","momentum[GeV/#it{c}]");
    SetTH2(h2_pid_ohh  ,"PID plot OHH"             ,"1/#beta","momentum[GeV/#it{c}]");
    SetTH2(h2_ohdedxbeta_all  ,"de/dx(OH) vs #beta all"         ,"#beta","momentum[GeV/#it{c}]");
    SetTH2(h2_ohdedxbeta_decut,"de/dx(OH) vs #beta all(dE cut)" ,"#beta","momentum[GeV/#it{c}]");
    SetTH2(h2_ohdedxbeta_ohv  ,"de/dx(OH) vs #beta OHV1-8"      ,"#beta","momentum[GeV/#it{c}]");
    SetTH2(h2_ohdedxbeta_ohh  ,"de/dx(OH) vs #beta OHH4-6"      ,"#beta","momentum[GeV/#it{c}]");
    SetTH2(h2_ihdedxbeta_all  ,"de/dx(IH) vs #beta all"         ,"#beta","momentum[GeV/#it{c}]");


  for(int i=0;i<12;i++){
    h2_pid_ohvl[i] = new TH2F(Form("h2_pid_ohvl%d",i+1) ,Form("h2_pid_ohvl%d",i+1) , 500, 0, 9, 1000, -1, 1);
    h2_pid_ohvr[i] = new TH2F(Form("h2_pid_ohvr%d",i+1) ,Form("h2_pid_ohvr%d",i+1) , 500, 0, 9, 1000, -1, 1);
    h_fl_ohvl[i]  = new TH1F(Form("h_fl_ohvl%d",i+1) ,Form("h_fl_ohvl%d",i+1) , 1000, 0, 1200);
    h_fl_ohvr[i]  = new TH1F(Form("h_fl_ohvr%d",i+1) ,Form("h_fl_ohvr%d",i+1) , 1000, 0, 1200);
    SetTH2(h2_pid_ohvl[i],Form("PID plot OHVL%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
    SetTH2(h2_pid_ohvr[i],Form("PID plot OHVR%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
    SetTH1(h_fl_ohvl[i] ,Form("Flight length OHVL%d",i+1) ,"flight length[cm]","counts", 1, 3001, 3);
    SetTH1(h_fl_ohvr[i] ,Form("Flight length OHVR%d",i+1) ,"flight length[cm]","counts", 1, 3001, 4);
  }
  for(int i=0;i<9;i++){
    h2_pid_ohhl[i] = new TH2F(Form("h2_pid_ohhl%d",i+1) ,Form("h2_pid_ohhl%d",i+1) , 500, 0, 9, 1000, -1, 1);
    h2_pid_ohhr[i] = new TH2F(Form("h2_pid_ohhr%d",i+1) ,Form("h2_pid_ohhr%d",i+1) , 500, 0, 9, 1000, -1, 1);;
    h_fl_ohhl[i]  = new TH1F(Form("h_fl_ohhl%d",i+1) ,Form("h_fl_ohhl%d",i+1) , 1000, 0, 1200);
    h_fl_ohhr[i]  = new TH1F(Form("h_fl_ohhr%d",i+1) ,Form("h_fl_ohhr%d",i+1) , 1000, 0, 1200);
    SetTH2(h2_pid_ohhl[i],Form("PID plot OHHL%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
    SetTH2(h2_pid_ohhr[i],Form("PID plot OHHR%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
    SetTH1(h_fl_ohhl[i] ,Form("Flight length OHHL%d",i+1) ,"flight length[cm]","counts", 1, 3001, 8);
    SetTH1(h_fl_ohhr[i] ,Form("Flight length OHHR%d",i+1) ,"flight length[cm]","counts", 1, 3001, 7);
  }

}
////////////////////////////////////////////////////////////////////////////
void track_ana::roop(){
cout<<"start roop"<<endl;

bool dedx_flag, ohseg_flag;

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){
    //tr->tree->GetEntry(n);
    tree->GetEntry(n);
    dedx_flag=ohseg_flag=false;
    if(n%100000==0) cout<<n<<" / "<<ENum<<endl;
    //int ohlr =tr->ohlr;
    //int ohseg=tr->ohseg;
    //int ihlr =tr->ihlr;
    //int ihseg=tr->ihseg;
    //double mom   =tr->mom;
    //double beta  =tr->beta;
    //double ohdedx=tr->ohdedx;
    //double ihdedx=tr->ihdedx;

//make flags
    if(ohdedx> -6.*beta+6.0 &&ohdedx>3.636*beta-3.2727 
     &&ihdedx> -6.*beta+6.0 &&ihdedx>3.636*beta-3.2727)dedx_flag=true;

    if( ((ohlr==-2||ohlr==2)&&ohseg>3&&ohseg<7) || ((ohlr==-1||ohlr==1)&&ohseg>0&&ohseg<9))ohseg_flag=true;

      h2_pid_all           ->Fill(1./beta,mom);
      if(dedx_flag && ohseg_flag){
        h2_pid_decut           ->Fill(1./beta,mom);
        h2_ohdedxbeta_decut    ->Fill(beta,ohdedx);
      }

      h2_ohdedxbeta_all      ->Fill(beta,ohdedx);
      h2_ihdedxbeta_all      ->Fill(beta,ihdedx);
      if(ihlr==-1){//IHL
        h2_pid_ihl[ihseg-1] ->Fill(1./beta,mom);
      }
      if(ihlr== 1){//IHR
        h2_pid_ihr[ihseg-1] ->Fill(1./beta,mom);
      }

      if(ohlr==-2){//OHHL
        h2_pid_ohh           ->Fill(1./beta,mom);
        h2_pid_ohhl[ohseg-1] ->Fill(1./beta,mom);
        h_fl_ohhl[ohseg-1]  ->Fill(fl);
        if(ohseg>3&&ohseg<7)h2_ohdedxbeta_ohh    ->Fill(beta,ohdedx);
      }
      else if(ohlr==-1){//OHVL
        h2_pid_ohv           ->Fill(1./beta,mom);
        if(ohseg<9)h2_ohdedxbeta_ohv    ->Fill(beta,ohdedx);
        h2_pid_ohvl[ohseg-1] ->Fill(1./beta,mom);
        h_fl_ohvl[ohseg-1]  ->Fill(fl);
      }
      else if(ohlr==2){//OHHR
        h2_pid_ohh           ->Fill(1./beta,mom);
        h2_pid_ohhr[ohseg-1] ->Fill(1./beta,mom);
        h_fl_ohhr[ohseg-1]  ->Fill(fl);
        if(ohseg>3&&ohseg<7)h2_ohdedxbeta_ohh    ->Fill(beta,ohdedx);
      }
      else if(ohlr==1){//OHVR
        h2_pid_ohv           ->Fill(1./beta,mom);
        if(ohseg<9)h2_ohdedxbeta_ohv    ->Fill(beta,ohdedx);
        h2_pid_ohvr[ohseg-1] ->Fill(1./beta,mom);
        h_fl_ohvr[ohseg-1]  ->Fill(fl);
      }

    }
  ofp->Write();
}
////////////////////////////////////////////////////////////////////////////
void track_ana::fit(){
TF1 *f_tdl = new TF1("f_tdl","gaus",-1,1);
//h_tdl_c       ->Fit(f_tdl);
double peak[4][12];
double tof[4][12];
double max[4][12];
      for(int i=0;i<12;i++){
      peak[0][i]=h_fl_ohvl[i]->GetXaxis()->GetBinCenter(h_fl_ohvl[i]->GetMaximumBin() );
      peak[1][i]=h_fl_ohvr[i]->GetXaxis()->GetBinCenter(h_fl_ohvr[i]->GetMaximumBin() );
      max[0][i]=h_fl_ohvl[i]->GetMaximum() ;
      max[1][i]=h_fl_ohvr[i]->GetMaximum() ;
      tex_fl[0][i] = new TLatex(peak[0][i]+100,max[0][i]/2.,Form("FL=%.01lf cm. tof=%.02lf ns",peak[0][i],0.01*peak[0][i]/c) );
      tex_fl[1][i] = new TLatex(peak[1][i]+100,max[1][i]/2.,Form("FL=%.01lf cm. tof=%.02lf ns",peak[1][i],0.01*peak[1][i]/c) );
      }
      for(int i=0;i<8;i++){
      peak[2][i]=h_fl_ohhl[i]->GetXaxis()->GetBinCenter(h_fl_ohhl[i]->GetMaximumBin() );
      peak[3][i]=h_fl_ohhr[i]->GetXaxis()->GetBinCenter(h_fl_ohhr[i]->GetMaximumBin() );
      max[2][i]=h_fl_ohhl[i]->GetMaximum() ;
      max[3][i]=h_fl_ohhr[i]->GetMaximum() ;
      tex_fl[2][i] = new TLatex(peak[2][i]+100,max[2][i]/2.,Form("FL=%.01lf cm. tof=%.02lf ns",peak[2][i],0.01*peak[2][i]/c) );
      tex_fl[3][i] = new TLatex(peak[3][i]+100,max[3][i]/2.,Form("FL=%.01lf cm. tof=%.02lf ns",peak[3][i],0.01*peak[3][i]/c) );
      }

   //for(int i=0;i<4;i++){
   // cout<<"oh "<<endl;
   // for(int j=0;j<8;j++){
   // cout<<0.01*peak[i][j]/c<<endl;
   // }
   //}

}
////////////////////////////////////////////////////////////////////////////
void track_ana::draw(){

cc[0]->Divide(4,2);
cc[0]->cd(1);gPad->SetLogz(1); h2_pid_ohv  ->Draw("colz");
cc[0]->cd(2);gPad->SetLogz(1); h2_pid_ohh  ->Draw("colz");
cc[0]->cd(3);gPad->SetLogz(1); h2_pid_all  ->Draw("colz");
cc[0]->cd(4);gPad->SetLogz(1); h2_pid_decut->Draw("colz");
cc[0]->cd(5);gPad->SetLogz(1); h2_ohdedxbeta_decut->Draw("colz");line ->DrawLine(0,6,1,0);line ->DrawLine(0.9,0,2,4);
//cc[0]->cd(5);gPad->SetLogz(1); h_ohdedxbeta_all->Draw("colz");
cc[0]->cd(6);gPad->SetLogz(1); h2_ohdedxbeta_ohv->Draw("colz");line ->DrawLine(0,6,1,0);line ->DrawLine(0.9,0,2,4);
cc[0]->cd(7);gPad->SetLogz(1); h2_ohdedxbeta_ohh->Draw("colz");line ->DrawLine(0,6,1,0);line ->DrawLine(0.9,0,2,4);
cc[0]->cd(8);gPad->SetLogz(1); h2_ihdedxbeta_all->Draw("colz");line ->DrawLine(0,6,1,0);line ->DrawLine(0.9,0,2,4);

cc[1]->Divide(4,3);
for(int i=0;i<12;i++){
cc[1]->cd(i+1);gPad->SetLogz(1);h2_pid_ohvl[i]->Draw("colz");
}
cc[2]->Divide(4,3);
for(int i=0;i<12;i++){
cc[2]->cd(i+1);gPad->SetLogz(1);h2_pid_ohvr[i]->Draw("colz");
}

cc[3]->Divide(4,3);
for(int i=0;i<8;i++){
cc[3]->cd(i+1);gPad->SetLogz(1);h2_pid_ohhl[i]->Draw("colz");
}

cc[4]->Divide(4,3);
for(int i=0;i<8;i++){
cc[4]->cd(i+1);gPad->SetLogz(1);h2_pid_ohhr[i]->Draw("colz");
}
cc[5]->Divide(4,3);
for(int i=0;i<10;i++){
cc[5]->cd(i+1);gPad->SetLogz(1);h2_pid_ihl[i]->Draw("colz");
}
cc[6]->Divide(4,3);
for(int i=0;i<10;i++){
cc[6]->cd(i+1);gPad->SetLogz(1);h2_pid_ihr[i]->Draw("colz");
}


cc[7]->Divide(4,3);
for(int i=0;i<12;i++){
cc[7]->cd(i+1);gPad->SetLogy(1);h_fl_ohvl[i]->Draw("");
tex_fl[0][i]->Draw("same");
}

cc[8]->Divide(4,3);
for(int i=0;i<12;i++){
cc[8]->cd(i+1);gPad->SetLogy(1);h_fl_ohvr[i]->Draw("");
tex_fl[1][i]->Draw("same");
}

cc[8]->Divide(4,3);
for(int i=0;i<8;i++){
cc[8]->cd(i+1);gPad->SetLogy(1);h_fl_ohhl[i]->Draw("");
tex_fl[2][i]->Draw("same");
}

cc[9]->Divide(4,3);
for(int i=0;i<8;i++){
cc[9]->cd(i+1);gPad->SetLogy(1);h_fl_ohhr[i]->Draw("");
tex_fl[3][i]->Draw("same");
}

cc[10]->Clear();
gPad->SetLogz(1);  h2_pid_decut->SetStats(0);h2_pid_decut->Draw("colz");
f_pi->Draw("same");
f_k ->Draw("same");
f_p ->Draw("same");


}
////////////////////////////////////////////////////////////////////////////
void track_ana::SetPdfFilename(string ifname){
pdf_name = ifname;
} 
////////////////////////////////////////////////////////////////////////////
void track_ana::savecanvas(){
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
  cc[9] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[10]->Print(Form("%s" ,pdf_name.c_str()) );
  cc[10]->Print(Form("%s]",pdf_name.c_str()) );
cout<<pdf_name<<" saved!"<<endl;
}
////////////////////////////////////////////////
void track_ana::SetRoot(string ifname)
{
cout<<"SetRoot"<<endl;
  //tr = new Track_tr();
  //tr->add(ifname);

  //tr->readtree();
  //ENum = tr->GetEntries();
  add(ifname);
  readtree();
  ENum = GetEntries();
}
////////////////////////////////////////////////
void track_ana::GetHist(string ifname)
{
cout<<"Get Histogram from root file"<<endl;
TFile *file = new TFile(ifname.c_str() );
//h =(TH2F*)file->Get("h");
}
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "output0000.dat";
  string ofname = "root/hoge.root";
  string pdfname = "pdf/track/hoge.pdf";
  int ch;
  int MaxNum = 0;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:w:n:bcop:"))!=-1){
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
    case 'p':
      pdfname = optarg;
      break;
    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output root filename"<<endl;
      cout<<"-n : maximum number of analysed events"<<endl;
      cout<<"-p : print pdf file"<<endl;
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
  track_ana *ana = new track_ana();

  ana->SetMaxEvent(MaxNum);
  ana->SetRoot(ifname);
  ana->makehist(ofname);
  ana->roop();
  ana->fit();
  ana->draw();
  ana->SetPdfFilename(pdfname);
  ana->savecanvas();
  delete ana;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}

