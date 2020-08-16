//draw_tdlparam.cc  2017.8.3 Y.Toyama
//combine parameters


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

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"


#include "Settings.h"
#include "RKV_ParamMan.h"

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
class draw_tdlparam
{
 public:
         draw_tdlparam();
        ~draw_tdlparam();
  void SetParamRK(string ifname); 
  bool SetFileName(string RunFile);
  void SetValue(); 
  void draw(); 
  RKV_ParamMan *ParamManRKV;
  Settings *set;

  private:
    int run_num;
    FILE *fp;
    TGraph *tg_tdll[10], *tg_tdlr[10];
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
draw_tdlparam::draw_tdlparam()
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

  set = new Settings();

  ParamManRKV  = new RKV_ParamMan();
  for(int i=0;i<10;i++){
    tg_tdll[i] = new TGraph();
    tg_tdlr[i] = new TGraph();
  }
}
////////////////////////////////////////////////////////////////////////////
draw_tdlparam::~draw_tdlparam(){
}
////////////////////////////////////////////////////////////////////////////
void draw_tdlparam::SetValue(){
  int lr,cid;
  int nn=0;
  int runnum;
  char str[256];
  while(fgets(str,144,fp)!=0){
    if(str[0]=='#') continue;
    if(sscanf(str,"%d",&runnum)==1){
      if(runnum==0)cout<<runnum<<endl;
      else {
        ParamManRKV->read_param(Form("param/%d.param",runnum));
        for(int i=0;i<10;i++){//TDL
          cid = 15;
          lr=0;  tg_tdll[i] ->SetPoint(nn, runnum, ParamManRKV->GetT0(-1,cid,i+1));//TDLL
          lr=1;  tg_tdlr[i] ->SetPoint(nn, runnum, ParamManRKV->GetT0( 1,cid,i+1));//TDLR
        }
        nn++;
      }
    }
  }


}
////////////////////////////////////////////////////////////////////////////
bool draw_tdlparam::SetFileName(string RunFile){
  if((fp=fopen(RunFile.c_str() ,"r"))==0){
    std::cerr << "file open fail" << std::endl;
    return false;
  }
  else return true;
}
////////////////////////////////////////////////
void draw_tdlparam::draw(){
  for(int i=0;i<10;i++){
    set->SetGr(tg_tdll[i],Form("TDLL%d",i+1),"run","t_{0}[ns]",1,1+i%8, 20+i%8);
    set->SetGr(tg_tdlr[i],Form("TDLR%d",i+1),"run","t_{0}[ns]",1,1+i%8, 20+i%8);
  }

  TCanvas*c[2];
  for(int i=0;i<2;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),840,1188 );
  }

  c[0]->Divide(2,5);
  for(int i=0;i<10;i++){
    c[0]->cd(i+1);tg_tdll[i]->Draw("AP");
  }

  c[1]->Divide(2,5);
  for(int i=0;i<10;i++){
    c[1]->cd(i+1);tg_tdlr[i]->Draw("AP");
  }

  c[0]->Print("pdf/tdl_param.pdf[");
  c[0]->Print("pdf/tdl_param.pdf");
  c[1]->Print("pdf/tdl_param.pdf");
  c[1]->Print("pdf/tdl_param.pdf]");
  cout<<"pdf/tdl_param.pdf : saved"<<endl;
}

////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "runlist.dat";
  int ch;
  extern char *optarg;

  draw_tdlparam *tdl = new draw_tdlparam();

  if(!tdl->SetFileName(ifname))return -1;
  tdl->SetValue();
  tdl->draw();
  delete tdl;

  return 0;
}

