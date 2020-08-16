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


class tdl_ana
{
 public:
         tdl_ana();
        ~tdl_ana();
  void makehist();
  void loop();
  void fit();
  void draw(); 
  void savecanvas(); 
    void SetMaxEvent( int N )  { ENumMax = N; }

  private:
    Tree *tr;
    int GetMaxEvent() { return ENumMax; }
    int ENumMax;
    int ENum;

    TH1F *h_ctime;
    TH2F *h_ctime_uwid , *h_ctime_dwid , *h_ctime_bwid ;

    TH1F *h_ctime_c;
    TH2F *h_ctime_uwid_c , *h_ctime_dwid_c , *h_ctime_bwid_c;

    TH1F *h_tdl;
    TH2F *h_tdl_uwid1 , *h_tdl_dwid1;
    TH2F *h_tdl_uwid2 , *h_tdl_dwid2;

    TH1F *h_tdl_c;
    TH2F *h_tdl_uwid1_c , *h_tdl_dwid1_c;
    TH2F *h_tdl_uwid2_c , *h_tdl_dwid2_c;

    TH1F *h_tb;
    TH2F *h_tb_bwid1 , *h_tb_bwid2;
    TH1F *h_tb_c;
    TH2F *h_tb_bwid1_c , *h_tb_bwid2_c;

    int run_num;
TCanvas *c1,*c2,*c3,*c4,*c5;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
tdl_ana::tdl_ana()
{

  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");
  gROOT->SetBatch(0);

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
      

  tr = new Tree();
  tr->add("../../root/all/10212.root");
  tr->add("../../root/all/10213.root");

 // tr->add("../../root/all/10229.root");
 // tr->add("../../root/all/10230.root");

  tr->readtree();
  ENum = tr->GetEntries();
  c1= new TCanvas("c1","c1",1400,800 );
  c2= new TCanvas("c2","c2",1400,800 );
  c3= new TCanvas("c3","c3",1400,800 );
  c4= new TCanvas("c4","c4",1400,800 );
  c5= new TCanvas("c5","c5",1400,800 );

}
////////////////////////////////////////////////////////////////////////////
tdl_ana::~tdl_ana(){
}
////////////////////////////////////////////////////////////////////////////
void tdl_ana::makehist(){
  h_ctime    = new TH1F("h_ctime"   ,"h_ctime"   ,2400,-30,30);
  h_ctime_c  = new TH1F("h_ctime_c" ,"h_ctime_c" ,2400,-30,30);
  SetTH1(h_ctime  ,"Coincidence time","cointime (ns)","counts / 25ps");
  SetTH1(h_ctime_c,"Coincidence time w/ cut","cointime (ns)","counts / 25ps");
  h_ctime_uwid  = new TH2F("h_ctime_uwid" ,"h_ctime_uwid" ,200,-10,10,200,0,50);
  h_ctime_dwid  = new TH2F("h_ctime_dwid" ,"h_ctime_dwid" ,200,-10,10,200,0,50);
  h_ctime_bwid  = new TH2F("h_ctime_bwid" ,"h_ctime_bwid" ,200,-10,10,200,0,50);
  SetTH2(h_ctime_uwid,"Cointime vs Width","cointime (ns)","TDLU width (ns)");
  SetTH2(h_ctime_dwid,"Cointime vs Width","cointime (ns)","TDLD width (ns)");
  SetTH2(h_ctime_bwid,"Cointime vs Width","cointime (ns)","TagB width (ns)");

  h_ctime_uwid_c  = new TH2F("h_ctime_uwid_c" ,"h_ctime_uwid_c" ,200,-10,10,200,0,50);
  h_ctime_dwid_c  = new TH2F("h_ctime_dwid_c" ,"h_ctime_dwid_c" ,200,-10,10,200,0,50);
  h_ctime_bwid_c  = new TH2F("h_ctime_bwid_c" ,"h_ctime_bwid_c" ,200,-10,10,200,0,50);
  SetTH2(h_ctime_uwid_c,"Cointime vs Width","cointime (ns)","TDLU width (ns)");
  SetTH2(h_ctime_dwid_c,"Cointime vs Width","cointime (ns)","TDLD width (ns)");
  SetTH2(h_ctime_bwid_c,"Cointime vs Width","cointime (ns)","TagB width (ns)");

  h_tdl    = new TH1F("h_tdl"   ,"h_tdl"   ,2400,-30,30);
  h_tdl_c  = new TH1F("h_tdl_c" ,"h_tdl_c" ,2400,-30,30);
  SetTH1(h_tdl  ,"ToF TDL-TDL"       ,"ToF(ns)","counts / 25ps");
  SetTH1(h_tdl_c,"ToF TDL-TDL w/ cut","ToF(ns)","counts / 25ps");
  h_tdl_uwid1 = new TH2F("h_tdl_uwid1" ,"h_tdl_uwid1" ,200,-10,10,200,0,50);
  h_tdl_dwid1 = new TH2F("h_tdl_dwid1" ,"h_tdl_dwid1" ,200,-10,10,200,0,50);
  h_tdl_uwid2 = new TH2F("h_tdl_uwid2" ,"h_tdl_uwid2" ,200,-10,10,200,0,50);
  h_tdl_dwid2 = new TH2F("h_tdl_dwid2" ,"h_tdl_dwid2" ,200,-10,10,200,0,50);
  SetTH2(h_tdl_uwid1,"ToF TDL-TDL vs Width","cointime (ns)","TDL4U width (ns)");
  SetTH2(h_tdl_dwid1,"ToF TDL-TDL vs Width","cointime (ns)","TDL4D width (ns)");
  SetTH2(h_tdl_uwid2,"ToF TDL-TDL vs Width","cointime (ns)","TDL10U width (ns)");
  SetTH2(h_tdl_dwid2,"ToF TDL-TDL vs Width","cointime (ns)","TDL10D width (ns)");
  h_tdl_uwid1_c = new TH2F("h_tdl_uwid1_c" ,"h_tdl_uwid1_c" ,200,-10,10,200,0,50);
  h_tdl_dwid1_c = new TH2F("h_tdl_dwid1_c" ,"h_tdl_dwid1_c" ,200,-10,10,200,0,50);
  h_tdl_uwid2_c = new TH2F("h_tdl_uwid2_c" ,"h_tdl_uwid2_c" ,200,-10,10,200,0,50);
  h_tdl_dwid2_c = new TH2F("h_tdl_dwid2_c" ,"h_tdl_dwid2_c" ,200,-10,10,200,0,50);
  SetTH2(h_tdl_uwid1_c,"ToF TDL-TDL vs Width","cointime (ns)","TDL4U width (ns)");
  SetTH2(h_tdl_dwid1_c,"ToF TDL-TDL vs Width","cointime (ns)","TDL4D width (ns)");
  SetTH2(h_tdl_uwid2_c,"ToF TDL-TDL vs Width","cointime (ns)","TDL10U width (ns)");
  SetTH2(h_tdl_dwid2_c,"ToF TDL-TDL vs Width","cointime (ns)","TDL10D width (ns)");

  h_tb    = new TH1F("h_tb"   ,"h_tb"   ,2400,-29.1,30.9);
  h_tb_c  = new TH1F("h_tb_c" ,"h_tb_c" ,2400,-29.1,30.9);
  SetTH1(h_tb  ,"ToF TagB-TagB"       ,"ToF(ns)","counts / 25ps");
  SetTH1(h_tb_c,"ToF TagB-TagB w/ cut","ToF(ns)","counts / 25ps");
  h_tb_bwid1 = new TH2F("h_tb_bwid1" ,"h_tb_bwid1" ,800,-9.1,10.9,200,0,50);
  h_tb_bwid2 = new TH2F("h_tb_bwid2" ,"h_tb_bwid2" ,800,-9.1,10.9,200,0,50);
  SetTH2(h_tb_bwid1,"ToF TagB-TagB vs Width","ToF(ns)","TagB width (ns)");
  SetTH2(h_tb_bwid2,"ToF TagB-TagB vs Width","ToF(ns)","TagB width (ns)");
  h_tb_bwid1_c = new TH2F("h_tb_bwid1_c" ,"h_tb_bwid1_c" ,800,-9.1,10.9,200,0,50);
  h_tb_bwid2_c = new TH2F("h_tb_bwid2_c" ,"h_tb_bwid2_c" ,800,-9.1,10.9,200,0,50);
  SetTH2(h_tb_bwid1_c,"ToF TagB-TagB vs Width","ToF(ns)","TagB width (ns)");
  SetTH2(h_tb_bwid2_c,"ToF TagB-TagB vs Width","ToF(ns)","TagB width (ns)");

}
////////////////////////////////////////////////////////////////////////////
void tdl_ana::loop(){
double tdl_offset=2.;
double tagb_offset=0.;
double coin_offset=37.4;
//double coin_offset=12.4;
double Width_R4Umin,Width_R4Dmin,Width_R10Umin,Width_R10Dmin;
double Width_R4Umax,Width_R4Dmax,Width_R10Umax,Width_R10Dmax;
Width_R4Umin = 25.;Width_R4Dmin = 20.;Width_R10Umin =25. ;Width_R10Dmin =20.;
Width_R4Umax = 30.;Width_R4Dmax = 25.;Width_R10Umax =30. ;Width_R10Dmax =25.;
double Width_B1max,Width_B2max,Width_B1min,Width_B2min;
Width_B1max = 27.;
Width_B2max = 33.;
Width_B1min = 22.;
Width_B2min = 28.;

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){
    tr->tree->GetEntry(n);
    run_num = tr->runnum;
    if(n%100000==0) cout<<n<<" / "<<ENum<<endl;
		double	tdl_uwid1 = tr->tdlrutime_w[3];
		double	tdl_dwid1 = tr->tdlrdtime_w[3];
		double	tdl_uwid2 = tr->tdlrutime_w[9];
		double	tdl_dwid2 = tr->tdlrdtime_w[9];
		double	tdl_mt1   = tr->tdlrmtime[3];
		double	tdl_mt2   = tr->tdlrmtime[9];

		double	tagb_wid1 = tr->tagbtime_w[22];
		double	tagb_wid2 = tr->tagbtime_w[23];
        double  tagb_mt1  = (tr->tagbtdc_lm[38][0] - tr->tagbtdc_lm[22][0]) * 0.025;
        double  tagb_mt2  = (tr->tagbtdc_lm[38][0] - tr->tagbtdc_lm[23][0]) * 0.025;
        double  tagb_tof  = (tr->tagbtdc_lm[22][0] - tr->tagbtdc_lm[23][0]) * 0.025;
		int	tagb_tdc1  = tr->tagbtdc_l[22];
		int	tagb_tdc2  = tr->tagbtdc_l[23];
		int	tagb_trigtdc  = tr->tagbtdc_l[38];

	  if(//Tekito counter v.s. TDLR4
			  tdl_uwid1>Width_R4Umin &&tdl_uwid1<Width_R4Umax &&
			  tdl_dwid1>Width_R4Dmin &&tdl_dwid1<Width_R4Dmax &&
			  tdl_uwid2>Width_R10Umin&&tdl_uwid2<Width_R10Umax&&
			  tdl_dwid2>Width_R10Dmin&&tdl_dwid2<Width_R10Dmax
		){
          h_tdl_c       -> Fill(tdl_mt1 - tdl_mt2 + tdl_offset);
          h_tdl_uwid1_c -> Fill(tdl_mt1 - tdl_mt2 + tdl_offset,tdl_uwid1);
          h_tdl_dwid1_c -> Fill(tdl_mt1 - tdl_mt2 + tdl_offset,tdl_dwid1);
          h_tdl_uwid2_c -> Fill(tdl_mt1 - tdl_mt2 + tdl_offset,tdl_uwid2);
          h_tdl_dwid2_c -> Fill(tdl_mt1 - tdl_mt2 + tdl_offset,tdl_dwid2);
	      }

	  if(//trig_time1>0&&
			  tdl_uwid1>0&&tdl_uwid1<10000&&
			  tdl_dwid1>0&&tdl_dwid1<10000&&
			  tdl_uwid2>0&&tdl_uwid2<10000&&
			  tdl_dwid2>0&&tdl_dwid2<10000
		){h_tdl       -> Fill(tdl_mt1 - tdl_mt2 + tdl_offset);
          h_tdl_uwid1 -> Fill(tdl_mt1 - tdl_mt2 + tdl_offset,tdl_uwid1);
          h_tdl_dwid1 -> Fill(tdl_mt1 - tdl_mt2 + tdl_offset,tdl_dwid1);
          h_tdl_uwid2 -> Fill(tdl_mt1 - tdl_mt2 + tdl_offset,tdl_uwid2);
          h_tdl_dwid2 -> Fill(tdl_mt1 - tdl_mt2 + tdl_offset,tdl_dwid2);
         }

	  if(//TagB 23 v.s. TagB24
             tagb_trigtdc>0&& tagb_tdc1>0          &&tagb_tdc2>0&&
			  tagb_wid1>Width_B1min&&tagb_wid1<Width_B1max&&
			  tagb_wid2>Width_B2min&&tagb_wid2<Width_B2max
		){
          h_tb_c       -> Fill(tagb_tof + tagb_offset);
          h_tb_bwid1_c -> Fill(tagb_tof + tagb_offset,tagb_wid1);
          h_tb_bwid2_c -> Fill(tagb_tof + tagb_offset,tagb_wid2);
	      }
	  if(//TagB 23 v.s. TagB24
             tagb_trigtdc>0&& tagb_tdc1>0          &&tagb_tdc2>0&&
			  tagb_wid1>0&&tagb_wid1<10000&&
			  tagb_wid2>0&&tagb_wid2<10000
		){
          h_tb       -> Fill(tagb_tof + tagb_offset);
          h_tb_bwid1 -> Fill(tagb_tof + tagb_offset,tagb_wid1);
          h_tb_bwid2 -> Fill(tagb_tof + tagb_offset,tagb_wid2);
	      }
	  if(//TagB 23 v.s. TDL R4
             tagb_trigtdc>0&& tagb_tdc1>0  &&
			  tagb_wid1>0&&tagb_wid1<10000&&
			  tdl_uwid1>0&&tdl_uwid1<10000&&
			  tdl_dwid1>0&&tdl_dwid1<10000
       ){
          h_ctime       -> Fill(tagb_mt1 - tdl_mt1 + coin_offset);
          h_ctime_uwid  -> Fill(tagb_mt1 - tdl_mt1 + coin_offset,tdl_uwid1);
          h_ctime_dwid  -> Fill(tagb_mt1 - tdl_mt1 + coin_offset,tdl_dwid1);
          h_ctime_bwid  -> Fill(tagb_mt1 - tdl_mt1 + coin_offset,tagb_wid1);
        }

	  if(//TagB 23 v.s. TDL R4
             tagb_trigtdc>0&& tagb_tdc1>0  &&
			  tagb_wid1>Width_B1min &&tagb_wid1<Width_B1max &&
			  tdl_uwid1>Width_R4Umin&&tdl_uwid1<Width_R4Umax&&
			  tdl_dwid1>Width_R4Dmin&&tdl_dwid1<Width_R4Dmax
       ){
          h_ctime_c       -> Fill(tagb_mt1 - tdl_mt1 + coin_offset);
          h_ctime_uwid_c  -> Fill(tagb_mt1 - tdl_mt1 + coin_offset,tdl_uwid1);
          h_ctime_dwid_c  -> Fill(tagb_mt1 - tdl_mt1 + coin_offset,tdl_dwid1);
          h_ctime_bwid_c  -> Fill(tagb_mt1 - tdl_mt1 + coin_offset,tagb_wid1);
        }
 }

}
////////////////////////////////////////////////////////////////////////////
void tdl_ana::fit(){
TF1 *f_tdl = new TF1("f_tdl","gaus",-1,1);
h_tdl_c       ->Fit(f_tdl);
TF1 *f_tb = new TF1("f_tb","gaus",-1,1);
h_tb_c       ->Fit(f_tb);
TF1 *f_ct = new TF1("f_ct","gaus",-2,2);
h_ctime_c       ->Fit(f_ct,"R");

}
////////////////////////////////////////////////////////////////////////////
void tdl_ana::draw(){

c1->Clear();
c2->Clear();
c3->Clear();
c4->Clear();
c5->Clear();

c1->Divide(3,2);
c1->cd(1);h_tdl       ->Draw("");
c1->cd(2);h_tdl_uwid1 ->Draw("colz");
c1->cd(3);h_tdl_dwid1 ->Draw("colz");
c1->cd(4);h_tdl_uwid2 ->Draw("colz");
c1->cd(5);h_tdl_dwid2 ->Draw("colz");

c2->Divide(3,2);
c2->cd(1);h_tdl_c       ->Draw("");
c2->cd(2);h_tdl_uwid1_c ->Draw("colz");
c2->cd(3);h_tdl_dwid1_c ->Draw("colz");
c2->cd(4);h_tdl_uwid2_c ->Draw("colz");
c2->cd(5);h_tdl_dwid2_c ->Draw("colz");

c3->Divide(3,2);
c3->cd(1);h_tb         ->Draw("");
c3->cd(2);h_tb_bwid1   ->Draw("colz");
c3->cd(3);h_tb_bwid2   ->Draw("colz");
c3->cd(4);h_tb_c       ->Draw("");
c3->cd(5);h_tb_bwid1_c ->Draw("colz");
c3->cd(6);h_tb_bwid2_c ->Draw("colz");

c4->Divide(3,2);
c4->cd(1); h_ctime      ->Draw("");
c4->cd(2); h_ctime_uwid ->Draw("colz");
c4->cd(3); h_ctime_dwid ->Draw("colz");
c4->cd(5); h_ctime_bwid ->Draw("colz");

c5->Divide(3,2);
c5->cd(1); h_ctime_c      ->Draw("");
c5->cd(2); h_ctime_uwid_c ->Draw("colz");
c5->cd(3); h_ctime_dwid_c ->Draw("colz");
c5->cd(5); h_ctime_bwid_c ->Draw("colz");
}
////////////////////////////////////////////////////////////////////////////
void tdl_ana::savecanvas(){
  c1->Print(Form("pdf/%d.pdf[",run_num) );
  c1->Print(Form("pdf/%d.pdf" ,run_num) );
  c2->Print(Form("pdf/%d.pdf" ,run_num) );
  c3->Print(Form("pdf/%d.pdf" ,run_num) );
  c4->Print(Form("pdf/%d.pdf" ,run_num) );
  c5->Print(Form("pdf/%d.pdf" ,run_num) );
  c5->Print(Form("pdf/%d.pdf]",run_num) );
cout<<"saved"<<endl;
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

  TApplication *theApp = new TApplication("App", &argc, argv);
  tdl_ana *ana = new tdl_ana();

  ana->SetMaxEvent(MaxNum);
  ana->makehist();
  ana->loop();
  ana->fit();
  ana->draw();
//  ana->savecanvas();
  delete ana;

  theApp->Run();
  return 0;
}

