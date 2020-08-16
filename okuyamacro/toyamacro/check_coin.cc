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
#include "Settings.h"

#define Calibration

static const double PI = 4.0*atan(1.);
static const double mrad_to_deg = 1./1000*180./PI;
const double Mp = 938.272046;          // proton       mass (MeV/c2)
const double c = 0.299792458;          // speed of light in vacuum (m/ns)
const int NTB   = 40;  // No. of TagB
const int NTF   =160;  // No. of TagB
const int NMRPC =  6;  // No. of MRPC 


class check_coin
{
 public:
         check_coin();
        ~check_coin();
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
  Settings *set;

  TGraph *gr_tdllu_t0[10], *gr_tdlru_t0[10];
  TGraph *gr_tdlld_t0[10], *gr_tdlrd_t0[10];
  TLatex *tex;
  TBox *box;

  double tdllu_t0_max[10]  ,tdlru_t0_max[10];
  double tdlld_t0_max[10]  ,tdlrd_t0_max[10];
  double tdllu_t0_min[10]  ,tdlru_t0_min[10];
  double tdlld_t0_min[10]  ,tdlrd_t0_min[10];
  double tdllu_t0[10][600]  ,tdlru_t0[10][600];
  double tdlld_t0[10][600]  ,tdlrd_t0[10][600];

  int run_tdllu_t0_max[10]  ,run_tdlru_t0_max[10];
  int run_tdlld_t0_max[10]  ,run_tdlrd_t0_max[10];
  int run_tdllu_t0_min[10]  ,run_tdlru_t0_min[10];
  int run_tdlld_t0_min[10]  ,run_tdlrd_t0_min[10];


  double rnum[600];
  int  nn;
 
    int run_num;
    int content;// number of filled event
    TCanvas *c[13];

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
check_coin::check_coin()
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
    for(int i=0; i<10; i++){
      tdllu_t0[i][j]=tdlld_t0[i][j]=0.;
      tdlru_t0[i][j]=tdlrd_t0[i][j]=0.;
    }
    rnum[j]=0.;
  }

 for(int i=0; i<10; i++){
  tdllu_t0_max[i] = tdlru_t0_max[i] =   -1.;
  tdlld_t0_max[i] = tdlrd_t0_max[i] =   -1.;
  tdllu_t0_min[i] = tdlru_t0_min[i] = 999999.;
  tdlld_t0_min[i] = tdlrd_t0_min[i] = 999999.;
 }

  nn=0;

  set = new Settings();
  tex = new TLatex(1,1,"");
  set->SetTLatex(tex,1,0.05,11);

  for(int i=0;i<13;i++){
    c[i]= new TCanvas(Form("c%d",i+1),Form("c%d",i+1),900,800 );
  }

}
////////////////////////////////////////////////////////////////////////////
check_coin::~check_coin(){
}
////////////////////////////////////////////////////////////////////////////
void check_coin::analysis(){

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
      for(int i=0;i<10;i++){
        //t0
	    tdllu_t0[i][nn] = HodoParamManager->GetTdcOffset(0,CID_TDL,i+1,0);//TDLLU
	    tdlld_t0[i][nn] = HodoParamManager->GetTdcOffset(0,CID_TDL,i+1,1);//TDLLD
	    tdlru_t0[i][nn] = HodoParamManager->GetTdcOffset(1,CID_TDL,i+1,0);//TDLRU
	    tdlrd_t0[i][nn] = HodoParamManager->GetTdcOffset(1,CID_TDL,i+1,1);//TDLRD
        if(tdllu_t0[i][nn]>tdllu_t0_max[i]){tdllu_t0_max[i]=tdllu_t0[i][nn];run_tdllu_t0_max[i]=runnum;}
        if(tdllu_t0[i][nn]<tdllu_t0_min[i]){tdllu_t0_min[i]=tdllu_t0[i][nn];run_tdllu_t0_min[i]=runnum;}
        if(tdlld_t0[i][nn]>tdlld_t0_max[i]){tdlld_t0_max[i]=tdlld_t0[i][nn];run_tdlld_t0_max[i]=runnum;}
        if(tdlld_t0[i][nn]<tdlld_t0_min[i]){tdlld_t0_min[i]=tdlld_t0[i][nn];run_tdlld_t0_min[i]=runnum;}
        if(tdlru_t0[i][nn]>tdlru_t0_max[i]){tdlru_t0_max[i]=tdlru_t0[i][nn];run_tdlru_t0_max[i]=runnum;}
        if(tdlru_t0[i][nn]<tdlru_t0_min[i]){tdlru_t0_min[i]=tdlru_t0[i][nn];run_tdlru_t0_min[i]=runnum;}
        if(tdlrd_t0[i][nn]>tdlrd_t0_max[i]){tdlrd_t0_max[i]=tdlrd_t0[i][nn];run_tdlrd_t0_max[i]=runnum;}
        if(tdlrd_t0[i][nn]<tdlrd_t0_min[i]){tdlrd_t0_min[i]=tdlrd_t0[i][nn];run_tdlrd_t0_min[i]=runnum;}
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
bool check_coin::set_filename(string RunFile){
  if((fp=fopen(RunFile.c_str() ,"r"))==0){
    std::cerr << "file open fail" << std::endl;
    return false;
  }
  else return true;
}
////////////////////////////////////////////////////////////////////////////
void check_coin::set_graph(){
 /*IH*/
 for(int i=0;i<10;i++){
    //pedestal
    gr_tdllu_t0[i]  = new TGraph( nn, rnum, tdllu_t0[i]);//IHLU
    gr_tdlld_t0[i]  = new TGraph( nn, rnum, tdlld_t0[i]);//IHLD
    gr_tdlru_t0[i]  = new TGraph( nn, rnum, tdlru_t0[i]);//IHRU
    gr_tdlrd_t0[i]  = new TGraph( nn, rnum, tdlrd_t0[i]);//IHRD
    //SetGr(TGraph *gr, TString name, TString xname, TString yname, int LColor, int MColor, int MStyle, double Yoffset)
    set->SetGr(gr_tdllu_t0[i], Form("T0 TDLL%dU",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
    set->SetGr(gr_tdlld_t0[i], Form("T0 TDLL%dD",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
    set->SetGr(gr_tdlru_t0[i], Form("T0 TDLR%dU",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
    set->SetGr(gr_tdlrd_t0[i], Form("T0 TDLR%dD",i+1), "run num","ch", 2+i, 2+i%8, 20+i, 0.9);
 }
}
////////////////////////////////////////////////////////////////////////////
void check_coin::fit(){
}
////////////////////////////////////////////////////////////////////////////
void check_coin::draw(){

////////////
//Pedestal//
////////////
c[0]->Clear();
c[0]->Divide(4,5);
 for(int i=0;i<10;i++){
 c[0]->cd(i+1); gr_tdllu_t0[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], tdllu_t0_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",tdllu_t0_max[i],run_tdllu_t0_max[i],tdllu_t0_min[i],run_tdllu_t0_min[i]));
 c[0]->cd(i+11);gr_tdlld_t0[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], tdlld_t0_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",tdlld_t0_max[i],run_tdlld_t0_max[i],tdlld_t0_min[i],run_tdlld_t0_min[i]));
 }

c[1]->Clear();
c[1]->Divide(4,5);
 for(int i=0;i<10;i++){
 c[1]->cd(i+1); gr_tdlru_t0[i]-> Draw("AP");
                tex-> DrawLatex(rnum[0], tdlru_t0_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",tdlru_t0_max[i],run_tdlru_t0_max[i],tdlru_t0_min[i],run_tdlru_t0_min[i]));
 c[1]->cd(i+11);gr_tdlrd_t0[i]-> Draw("AP");                                                                                                                            
                tex-> DrawLatex(rnum[0], tdlrd_t0_max[i], Form("max:%.02lf(run%d), min:%.02lf(run%d)",tdlrd_t0_max[i],run_tdlrd_t0_max[i],tdlrd_t0_min[i],run_tdlrd_t0_min[i]));
 }


}
////////////////////////////////////////////////////////////////////////////
void check_coin::savecanvas(){
cout<<RNumMin<<" - "<<RNumMax<<endl;
string pdf_name =Form("pdf/check_coin/run%d_%d.pdf",RNumMin, RNumMax  ); 
  c[0] ->Print(Form("%s[",pdf_name.c_str()  ) );
  c[0] ->Print(Form("%s" ,pdf_name.c_str()  ) );
  c[1] ->Print(Form("%s" ,pdf_name.c_str()  ) );
  c[1]->Print(Form("%s]",pdf_name.c_str()  ) );
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
  check_coin *ana = new check_coin();
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

