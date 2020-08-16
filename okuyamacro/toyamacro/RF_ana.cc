#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
using namespace std;

#include "TApplication.h"
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
#include "Setting.h"

static const double PI = 4.0*atan(1.);
static const double mrad_to_deg = 1./1000*180./PI;
const double Mp = 938.272046;          // proton       mass (MeV/c2)
//const double c = 0.299792458;          // speed of light in vacuum (m/ns)
const int NCanvas = 14;

double tb_offs[24]={ 0.,  //TagB1 
                     0.,  //TagB2 
                     0.,  //TagB3 
                     0.,  //TagB4 
                     0.,  //TagB5 
                     0.,  //TagB6 
                     0.,  //TagB7 
                     0.,  //TagB8 
                     0.,  //TagB9 
                     0.,  //TagB10
                     0.,  //TagB11
                     0.,  //TagB12
                     0.,  //TagB13
                     0.,  //TagB14
                     0.,  //TagB15
                     0.,  //TagB16
                     0.,  //TagB17
                     0.,  //TagB18
                     0.,  //TagB19
                     0.,  //TagB20
                     0.,  //TagB21
                     0.,  //TagB22
                     0.,  //TagB23
                     0.}; //TagB24

double tb_w_min[24]={0., //23.,
                     0., //20.,
                     0., //20.,
                     0., //10.,
                     0., //15.,
                     0., //20.,
                     0., //15.,
                     0., //15.,
                     0., //15.,
                     0., //20.,
                     0., //15.,
                     0., //20.,
                     0., //13.,
                     0., //20.,
                     0., //20.,
                     0., //18.,
                     0., //22.,
                     0., //21.,
                     0., //22.,
                     0., //25.,
                     0., //16.,
                     0., //25.,
                     0., //21.,
                     0.}; //25.};
double tb_w_max[24]={50., //35.,
                     50., //30.,
                     50., //30.,
                     50., //15.,
                     50., //25.,
                     50., //30.,
                     50., //25.,
                     50., //30.,
                     50., //25.,
                     50., //25.,
                     50., //25.,
                     50., //35.,
                     50., //20.,
                     50., //30.,
                     50., //30.,
                     50., //30.,
                     50., //35.,
                     50., //30.,
                     50., //35.,
                     50., //35.,
                     50., //30.,
                     50., //40.,
                     50., //35.,
                     50.}; //40.};


double tdll_offs[NTDL]={0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                       -0.9};

double tdllu_w_min[NTDL]={0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.};

double tdllu_w_max[NTDL]={1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.};

double tdlld_w_min[NTDL]={0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.,
                          0.};

double tdlld_w_max[NTDL]={1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.,
                          1000.};
//////////////////////////////////////////+++++++++++////////////////
class RF_ana : public Tree
{
 public:
         RF_ana();
        ~RF_ana();
  void makehist();
  void loop();
  void fit();
  void draw(); 
  void savecanvas(string name); 
  void show_tfratio(); 
  void SetRoot(string ifname );
  void SetOfname(string name );
  void SetMaxEvent( int N )  { ENumMax = N; }
  Setting *set;

  private:
    int GetMaxEvent() { return ENumMax; }
    int ENumMax;
    int ENum;
    TFile *ofp;
    string ofname;

    TH1F *h_ctime_tb_RF[40], *h_ctime_tb_RF_cut[40], *h_ctime_tb_RF_f_or[40], *h_ctime_tb_RF_fcut[40][4], *h_time_tb_RF[40];
    TH1F *h_ctime_tb_RF_tbtb[40], *h_ctime_tb_RF_tbtb_bar[40], *h_ctime_tb_RF_fbor[40];
    TH1F *h_ctime_tdll_RF[NTDL], *h_ctime_tdll_RF_cut[NTDL], *h_time_tdll_RF[NTDL];
    TH1F *h_time_RF, *h_time_RFdiff[4];
    TH1F *h_ctof_tb[12], *h_ctof_tb_cut[12], *h_ctof_tb_fcut[12];

    TH1F *h_tagfhitpat_tb[40];
    
    TH1F *h_cointime_tb[40], *h_cointime_tdll[NTDL];
    TH2F *h2_cointime_tb_wid[40], *h2_cointime_tdll_uwid[NTDL], *h2_cointime_tdll_dwid[NTDL];
    
    TH2F *h2_ctime_tb_RF[40], *h2_ctime_tb_RF_cut[40], *h2_ctime_tb_RF_fcut[40][4], *h2_time_tb_RF[40];
    TH2F *h2_ctime_tdllu_RF[NTDL], *h2_ctime_tdllu_RF_cut[NTDL], *h2_time_tdllu_RF[NTDL];
    TH2F *h2_ctime_tdlld_RF[NTDL], *h2_ctime_tdlld_RF_cut[NTDL], *h2_time_tdlld_RF[NTDL];
    TH2F *h2_ctof_tb_w1[12], *h2_ctof_tb_w_cut1[12];
    TH2F *h2_ctof_tb_w2[12], *h2_ctof_tb_w_cut2[12];
    TH2F *h2_ctof_tb_RF1[12], *h2_ctof_tb_RF2[12];
    TH2F *h2_ctof_tb_tdc1[12], *h2_ctof_tb_tdc2[12];

    int run_num;
    int tb_total[40],tf_total[40];
    TCanvas *c[NCanvas];
    double peak_tb[40], max_tb[40];
    TLine *tl_tb;

    bool tdll_flag[NTDL];
    bool tb_flag[40], tb_rf_flag[40],tf_flag[40];
    bool tbtb_flag[20];//TagB odd seg and TagB even seg coincidence
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
RF_ana::RF_ana()
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

  set = new Setting();
  tl_tb = new TLine();

  for(int i=0;i<NCanvas;i++){
    c[i]= new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1400,800 );
  }

}
////////////////////////////////////////////////////////////////////////////
RF_ana::~RF_ana(){
}
////////////////////////////////////////////////////////////////////////////
void RF_ana::SetRoot(string ifname ){
  add(ifname.c_str());
  readtree();
  ENum = GetEntries();

}
////////////////////////////////////////////////////////////////////////////
void RF_ana::SetOfname(string name ){
  ofname = name;
}
////////////////////////////////////////////////////////////////////////////
void RF_ana::makehist(){
  ofp = new TFile(Form("%s",ofname.c_str()),"recreate");

  for(int i=0;i<40;i++){
    tb_total[i] = tf_total[i] = 0;
    h_ctime_tb_RF[i]     =new TH1F(Form("h_ctime_tb_RF%d",i+1)     ,Form("h_ctime_tb_RF%d",i+1)     ,1000,   -3, 3);
    h_ctime_tb_RF_cut[i] =new TH1F(Form("h_ctime_tb_RF_cut%d",i+1) ,Form("h_ctime_tb_RF_cut%d",i+1) ,1000,   -3, 3);
    h_ctime_tb_RF_f_or[i]=new TH1F(Form("h_ctime_tb_RF_f_or%d",i+1),Form("h_ctime_tb_RF_f_or%d",i+1),1000,   -3, 3);
    h_time_tb_RF[i]      =new TH1F(Form("h_time_tb_RF%d",i+1)      ,Form("h_time_tb_RF%d",i+1)      ,2000,  -10,10);
    h_ctime_tb_RF_tbtb[i]=new TH1F(Form("h_ctime_tb_RF_tbtb%d",i+1),Form("h_ctime_tb_RF_tbtb%d",i+1),1000,   -3, 3);
    h_ctime_tb_RF_fbor[i]=new TH1F(Form("h_ctime_tb_RF_fbor%d",i+1),Form("h_ctime_tb_RF_fbor%d",i+1),1000,   -3, 3);
    h_ctime_tb_RF_tbtb_bar[i]=new TH1F(Form("h_ctime_tb_RF_tbtb_bar%d",i+1),Form("h_ctime_tb_RF_tbtb_bar%d",i+1),1000,   -3, 3);
    h_tagfhitpat_tb[i]   =new TH1F(Form("h_tagfhitpat_tb%d",i+1)   ,Form("h_tagfhitpat_tb%d",i+1)   , 120,  0.5, 120.5);//TagF hit pattern with TagB hit
    set->SetTH1(h_ctime_tb_RF[i]     ,Form("TagB%d - RF(w/ TWC)"                              ,i+1),"TagB - RF[ns]"   ,"Counts"       , 1, 3000, 0);
    set->SetTH1(h_ctime_tb_RF_cut[i] ,Form("TagB%d - RF(w/ TWC&Width cut)"                    ,i+1),"TagB - RF[ns]"   ,"Counts"       , 1, 3000, 0);
    set->SetTH1(h_ctime_tb_RF_f_or[i],Form("TagB%d - RF(w/ TWC&Width&TagF cut)"               ,i+1),"TagB - RF[ns]"   ,"Counts"       , 4, 3000, 0);
    set->SetTH1(h_ctime_tb_RF_tbtb[i],Form("TagB%d - RF(w/ TWC&Width&TagBTagB cut)"           ,i+1),"TagB - RF[ns]"   ,"Counts"       , 3, 3000, 0);
    set->SetTH1(h_ctime_tb_RF_tbtb_bar[i],Form("TagB%d - RF(w/ TWC&Width&!TagBTagB cut)"      ,i+1),"TagB - RF[ns]"   ,"Counts"       , 2, 3000, 0);
    set->SetTH1(h_ctime_tb_RF_fbor[i],Form("TagB%d - RF(w/ TWC&Width&(TagBTagB or TagF) cut)" ,i+1),"TagB - RF[ns]"   ,"Counts"       , 6, 3000, 0);
    set->SetTH1(h_time_tb_RF[i]      ,Form("TagB%d - RF(w/o TWC)"                             ,i+1),"TagB - RF[ns]"   ,"Counts"       , 1, 3000, 0);
    int fcol =2;
    if(i%2==0)fcol = 4;
    set->SetTH1(h_tagfhitpat_tb[i]   ,Form("TagF hit pattern(w/ TagB hit)"          ,i+1),"TagF segment"      ,"Counts"       , 1, 3001, fcol);

    for(int j=0;j<4;j++){
      h_ctime_tb_RF_fcut[i][j]=new TH1F(Form("h_ctime_tb_RF_f%dcut%d",j+1,i+1),Form("h_ctime_tb_RF_f%dcut%d",j+1,i+1),1000, -3, 3);
      h2_ctime_tb_RF_fcut[i][j] =new TH2F(Form("h2_ctime_tb_RF_f%dcut%d",j+1,i+1),Form("h2_ctime_tb_RF_f%dcut%d",j+1,i+1),1000, -3, 3,1000,0,50);
      set->SetTH2(h2_ctime_tb_RF_fcut[i][j],Form("TagB%d - RF(w/ TWC&Width&F cut)" ,i+1),"TagB - RF[ns]","Width[ns]");
      set->SetTH1(h_ctime_tb_RF_fcut[i][j],Form("TagB%d - RF(w/ TWC&Width&F%d cut)" ,i+1,j+1),"TagB - RF[ns]"      ,"Counts"       , 2+2*j, 3000, 0);
    }


    h2_ctime_tb_RF[i]      =new TH2F(Form("h2_ctime_tb_RF%d",i+1)     ,Form("h2_ctime_tb_RF%d",i+1)     ,1000, -3, 3,1000,0,50);
    h2_ctime_tb_RF_cut[i]  =new TH2F(Form("h2_ctime_tb_RF_cut%d",i+1) ,Form("h2_ctime_tb_RF_cut%d",i+1) ,1000, -3, 3,1000,0,50);
    h2_time_tb_RF[i]       =new TH2F(Form("h2_time_tb_RF%d",i+1)      ,Form("h2_time_tb_RF%d",i+1)     ,2000,-10,10,1000,0,50);
    set->SetTH2(h2_ctime_tb_RF[i]     ,Form("TagB%d - RF(w/ TWC)"             ,i+1),"TagB - RF[ns]","Width[ns]");
    set->SetTH2(h2_ctime_tb_RF_cut[i] ,Form("TagB%d - RF(w/ TWC&Width cut)"   ,i+1),"TagB - RF[ns]","Width[ns]");
    set->SetTH2(h2_time_tb_RF[i]      ,Form("TagB%d - RF(w/o TWC)"            ,i+1),"TagB - RF[ns]","Width[ns]");
  }

  for(int i=0;i<NTDL;i++){
    h_ctime_tdll_RF[i]     =new TH1F(Form("h_ctime_tdll_RF%d",i+1)    ,Form("h_ctime_tdll_RF%d",i+1)    ,1000, -3, 3);
    h_ctime_tdll_RF_cut[i] =new TH1F(Form("h_ctime_tdll_RF_cut%d",i+1),Form("h_ctime_tdll_RF_cut%d",i+1),1000, -3, 3);
    h_time_tdll_RF[i]      =new TH1F(Form("h_time_tdll_RF%d",i+1)     ,Form("h_time_tdll_RF%d",i+1)     ,2000,-10,10);
    set->SetTH1(h_ctime_tdll_RF[i]    ,Form("TDLL%d - RF(w/ TWC)"             ,i+1),"TDLL - RF[ns]"      ,"Counts"       , 1, 3000, 0);
    set->SetTH1(h_ctime_tdll_RF_cut[i],Form("TDLL%d - RF(w/ TWC&Width cut)"   ,i+1),"TDLL - RF[ns]"      ,"Counts"       , 1, 3000, 0);
    set->SetTH1(h_time_tdll_RF[i]     ,Form("TDLL%d - RF(w/o TWC)"            ,i+1),"TDLL - RF[ns]"      ,"Counts"       , 1, 3000, 0);


    h2_ctime_tdllu_RF[i]      =new TH2F(Form("h2_ctime_tdllu_RF%d",i+1)    ,Form("h2_ctime_tdllu_RF%d",i+1)    ,1000, -3, 3,1000,0,50);
    h2_ctime_tdllu_RF_cut[i]  =new TH2F(Form("h2_ctime_tdllu_RF_cut%d",i+1),Form("h2_ctime_tdllu_RF_cut%d",i+1),1000, -3, 3,1000,0,50);
    h2_time_tdllu_RF[i]       =new TH2F(Form("h2_time_tdllu_RF%d",i+1)     ,Form("h2_time_tdllu_RF%d",i+1)     ,2000,-10,10,1000,0,50);
    h2_ctime_tdlld_RF[i]      =new TH2F(Form("h2_ctime_tdlld_RF%d",i+1)    ,Form("h2_ctime_tdlld_RF%d",i+1)    ,1000, -3, 3,1000,0,50);
    h2_ctime_tdlld_RF_cut[i]  =new TH2F(Form("h2_ctime_tdlld_RF_cut%d",i+1),Form("h2_ctime_tdlld_RF_cut%d",i+1),1000, -3, 3,1000,0,50);
    h2_time_tdlld_RF[i]       =new TH2F(Form("h2_time_tdlld_RF%d",i+1)     ,Form("h2_time_tdlld_RF%d",i+1)     ,2000,-10,10,1000,0,50);
    set->SetTH2(h2_ctime_tdllu_RF[i]    ,Form("TDLLU%d - RF(w/ TWC)"           ,i+1),"TDLL - RF[ns]","TDLLU Width[ns]");
    set->SetTH2(h2_ctime_tdllu_RF_cut[i],Form("TDLLU%d - RF(w/ TWC&Width cut)" ,i+1),"TDLL - RF[ns]","TDLLU Width[ns]");
    set->SetTH2(h2_time_tdllu_RF[i]     ,Form("TDLLU%d - RF(w/o TWC)"          ,i+1),"TDLL - RF[ns]","TDLLU Width[ns]");
    set->SetTH2(h2_ctime_tdlld_RF[i]    ,Form("TDLLD%d - RF(w/ TWC)"           ,i+1),"TDLL - RF[ns]","TDLLD Width[ns]");
    set->SetTH2(h2_ctime_tdlld_RF_cut[i],Form("TDLLD%d - RF(w/ TWC&Width cut)" ,i+1),"TDLL - RF[ns]","TDLLD Width[ns]");
    set->SetTH2(h2_time_tdlld_RF[i]     ,Form("TDLLD%d - RF(w/o TWC)"          ,i+1),"TDLL - RF[ns]","TDLLD Width[ns]");
  }

  h_time_RF         = new TH1F("h_time_RF"    ,"h_time_RF"    ,10000,2000,3200);

  for(int i=0;i<4;i++){
    h_time_RFdiff[i]     = new TH1F(Form("h_time_RFdiff%d",i+1),Form("h_time_RFdiff%d",i+1), 40, 165.95*(i+1)-5., 165.95*(i+1)+5.);
    set->SetTH1(h_time_RFdiff[i] ,Form("RF 1st - %d",i+2)            ,"diff. time[ns]"      ,"Counts/25ps"       , 1, 3000, 0);
  }

  for(int i=0;i<12;i++){
    h_ctof_tb[i]         = new TH1F(Form("h_ctof_tb_%d"     ,i+1), Form("h_ctof_tb_%d"     ,i+1), 80,-1,1);
    h_ctof_tb_cut[i]     = new TH1F(Form("h_ctof_tb_cut_%d" ,i+1), Form("h_ctof_tb_cut_%d" ,i+1), 80,-1,1);
    h_ctof_tb_fcut[i]    = new TH1F(Form("h_ctof_tb_fcut_%d",i+1), Form("h_ctof_tb_fcut_%d",i+1), 80,-1,1);
    set->SetTH1(h_ctof_tb[i]            ,Form("TagB%d - TagB%d(w/ TWC)"    ,2*i+1,2*i+2)      ,Form("TagB%d - TagB%d[ns]",2*i+1,2*i+2)    ,"Counts/25ps"       , 1, 3000, 0);
    set->SetTH1(h_ctof_tb_cut[i]        ,Form("TagB%d - TagB%d(w/ Qcut)"   ,2*i+1,2*i+2)      ,Form("TagB%d - TagB%d[ns]",2*i+1,2*i+2)    ,"Counts/25ps"       , 1, 3000, 0);
    set->SetTH1(h_ctof_tb_fcut[i]       ,Form("TagB%d - TagB%d(w/ Q&F cut)",2*i+1,2*i+2)      ,Form("TagB%d - TagB%d[ns]",2*i+1,2*i+2)    ,"Counts/25ps"       , 1, 3000, 0);

    h2_ctof_tb_w1[i]          = new TH2F(Form("h2_ctof_tb_w1_%d"      ,i+1), Form("h2_ctof_tb_w1_%d"     ,i+1),2000, -1, 1,1000,0,50);
    h2_ctof_tb_w2[i]          = new TH2F(Form("h2_ctof_tb_w2_%d"      ,i+1), Form("h2_ctof_tb_w2_%d"     ,i+1),2000, -1, 1,1000,0,50);
    h2_ctof_tb_w_cut1[i]      = new TH2F(Form("h2_ctof_tb_w_cut1_%d"  ,i+1), Form("h2_ctof_tb_w_cut1_%d" ,i+1),2000, -1, 1,1000,0,50);
    h2_ctof_tb_w_cut2[i]      = new TH2F(Form("h2_ctof_tb_w_cut2_%d"  ,i+1), Form("h2_ctof_tb_w_cut2_%d" ,i+1),2000, -1, 1,1000,0,50);
    h2_ctof_tb_RF1[i]         = new TH2F(Form("h2_ctof_tb_RF1_%d"     ,i+1), Form("h2_ctof_tb_RF1_%d"     ,i+1),2000, -1, 1,2000,-1,1);
    h2_ctof_tb_RF2[i]         = new TH2F(Form("h2_ctof_tb_RF2_%d"     ,i+1), Form("h2_ctof_tb_RF2_%d"     ,i+1),2000, -1, 1,2000,-1,1);
    h2_ctof_tb_tdc1[i]        = new TH2F(Form("h2_ctof_tb_tdc1_%d"    ,i+1), Form("h2_ctof_tb_tdc1_%d"    ,i+1),2000, -1, 1,400,3000,7000);
    h2_ctof_tb_tdc2[i]        = new TH2F(Form("h2_ctof_tb_tdc2_%d"    ,i+1), Form("h2_ctof_tb_tdc2_%d"    ,i+1),2000, -1, 1,400,3000,7000);
    set->SetTH2(h2_ctof_tb_w1[i]        ,Form("TagB%d - TagB%d" ,2*i+1,2*i+2) ,Form("TagB%d - TagB%d[ns]",2*i+1,2*i+2),Form("TagB%d Width[ns]",2*i+1));
    set->SetTH2(h2_ctof_tb_w2[i]        ,Form("TagB%d - TagB%d" ,2*i+1,2*i+2) ,Form("TagB%d - TagB%d[ns]",2*i+1,2*i+2),Form("TagB%d Width[ns]",2*i+2));
    set->SetTH2(h2_ctof_tb_w_cut1[i]    ,Form("TagB%d - TagB%d" ,2*i+1,2*i+2) ,Form("TagB%d - TagB%d[ns]",2*i+1,2*i+2),Form("TagB%d Width[ns]",2*i+1));
    set->SetTH2(h2_ctof_tb_w_cut2[i]    ,Form("TagB%d - TagB%d" ,2*i+1,2*i+2) ,Form("TagB%d - TagB%d[ns]",2*i+1,2*i+2),Form("TagB%d Width[ns]",2*i+2));
    set->SetTH2(h2_ctof_tb_RF1[i]       ,Form("TagB%d - TagB%d" ,2*i+1,2*i+2) ,Form("TagB%d - TagB%d[ns]",2*i+1,2*i+2),Form("TagB%d - RF[ns]",2*i+1));
    set->SetTH2(h2_ctof_tb_RF2[i]       ,Form("TagB%d - TagB%d" ,2*i+1,2*i+2) ,Form("TagB%d - TagB%d[ns]",2*i+1,2*i+2),Form("TagB%d - RF[ns]",2*i+2));
    set->SetTH2(h2_ctof_tb_tdc1[i]      ,Form("TagB%d - TagB%d" ,2*i+1,2*i+2) ,Form("TagB%d - TagB%d[ns]",2*i+1,2*i+2),Form("TagB%d TDC[ch]",2*i+1));
    set->SetTH2(h2_ctof_tb_tdc2[i]      ,Form("TagB%d - TagB%d" ,2*i+1,2*i+2) ,Form("TagB%d - TagB%d[ns]",2*i+1,2*i+2),Form("TagB%d TDC[ch]",2*i+2));
  }


  for(int i=0;i<NTDL;i++){//TDL seg
    h_cointime_tdll[i] = new TH1F(Form("h_cointime_tdll%d",i+1),Form("h_cointime_tdll%d",i+1) , 200,-3,3);
    h2_cointime_tdll_uwid[i] = new TH2F(Form("h2_cointime_tdll%d_uwid",i+1),Form("h2_cointime_tdll%d_uwid",i+1) , 200,-3,3,1000,0,70);
    h2_cointime_tdll_dwid[i] = new TH2F(Form("h2_cointime_tdll%d_dwid",i+1),Form("h2_cointime_tdll%d_dwid",i+1) , 200,-3,3,1000,0,70);
    set->SetTH1(h_cointime_tdll[i],Form("cointime TDLL%d",i+1)  ,"cointime[ns]"      ,"Counts"       , 1, 3001, 2);
  }

  set->SetTLine(tl_tb , 6,1,1);
}
////////////////////////////////////////////////////////////////////////////
void RF_ana::loop(){
  double RF_offset=0.;//[ns]
  double RFtime[5];
  int bunch[40];
  
  double tagb_offset=0.;
  double tmp_tb,val_tb[24];
  double tmp_tdll,val_tdll[NTDL];

  bool tagf73_flag,tagf74_flag,tagf75_flag,tagf76_flag;
  bool tagb19_flag,tagb20_flag;
  int ntf73,ntf74,ntf75,ntf76;
  int ntf73_tb19,ntf74_tb19,ntf75_tb19,ntf76_tb19;
  int ntf73_tb20,ntf74_tb20,ntf75_tb20,ntf76_tb20;
  int ntf73_tb1920,ntf74_tb1920,ntf75_tb1920,ntf76_tb1920;
  ntf73 = ntf74 = ntf75 = ntf76 = 0;
  ntf73_tb19 = ntf74_tb19 = ntf75_tb19 = ntf76_tb19 = 0;
  ntf73_tb20 = ntf74_tb20 = ntf75_tb20 = ntf76_tb20 = 0;
  ntf73_tb1920 = ntf74_tb1920 = ntf75_tb1920 = ntf76_tb1920 = 0;

  //cut param for TDL24
  double tdlu_min =  70;
  double tdlu_max = 317;
  double tdld_min =  80;
  double tdld_max = 380;

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){
    tree->GetEntry(n);
    if(n%10000==0) cout<<n<<" / "<<ENum<<endl;

    //////////////
    //Initialize//
    //////////////
    for(int i=0;i<NTDL;i++){
      tdll_flag[i]=false;
    }
    for(int i=0;i<40;i++){
      tb_flag[i] = tb_rf_flag[i]=tf_flag[i] = false;
    }
    for(int i=0;i<20;i++)tbtb_flag[i] = false;

    tagf73_flag=tagf74_flag=tagf75_flag=tagf76_flag = false;
    tagb19_flag = tagb20_flag = false;

    ////////
    //Fill//
    ////////
    for(int i=0;i<5;i++){
      RFtime[i] = tagbtime_l[39][i]+RF_offset;
      h_time_RF -> Fill(RFtime[i]);
    }

    for(int i=0;i<4;i++){
      h_time_RFdiff[i] -> Fill(RFtime[0]-RFtime[i+1]);
      //cout<<i<<"  "<<RFtime[0]-RFtime[i+1]<<endl;
    }

    //TagB - TagB ToF
    for(int i=0;i<12;i++){
      if(tagbctime[2*i]>-100&&tagbctime[2*i+1]>-100){
        h_ctof_tb[i]     ->Fill(tagbctime[2*i]-tagbctime[2*i+1]);
        h2_ctof_tb_w1[i] ->Fill(tagbctime[2*i]-tagbctime[2*i+1], tagbtime_w[2*i]);
        h2_ctof_tb_w2[i] ->Fill(tagbctime[2*i]-tagbctime[2*i+1], tagbtime_w[2*i+1]);
        if(tagbtime_w[2*i]>tb_w_min[2*i]&&tagbtime_w[2*i]<tb_w_max[2*i]&&tagbtime_w[2*i+1]>tb_w_min[2*i+1]&&tagbtime_w[2*i+1]<tb_w_max[2*i+1]){
          if(abs(tagbctime[2*i]-tagbctime[2*i+1]) < 1.)tbtb_flag[i]=true;
          h_ctof_tb_cut[i]   ->Fill(tagbctime[2*i]-tagbctime[2*i+1]);
          h2_ctof_tb_RF1[i]  ->Fill(tagbctime[2*i]-tagbctime[2*i+1], val_tb[2*i]);
          h2_ctof_tb_RF2[i]  ->Fill(tagbctime[2*i]-tagbctime[2*i+1], val_tb[2*i+1]);
          //for(int j=0;j<10;j++){//multi hit
          //  h2_ctof_tb_tdc1[i] ->Fill(tagbctime[2*i]-tagbctime[2*i+1], tagbtdc_lm[2*i][j]);
          //  h2_ctof_tb_tdc2[i] ->Fill(tagbctime[2*i]-tagbctime[2*i+1], tagbtdc_lm[2*i+1][j]);
          //}
          //if(tagftdc_lm[4*2*i+3][0]>0){
          //  h_ctof_tb_fcut[i]     ->Fill(tagbctime[2*i]-tagbctime[2*i+1]);
          //}
          //if((val_tb[2*i]>-1&&val_tb[2*i]<-0.7)||(val_tb[2*i]<1&&val_tb[2*i]>0.7)){
          ////  h_ctof_tb_fcut[i]     ->Fill(tagbctime[2*i]-tagbctime[2*i+1]);
          //}
          //if(abs(val_tb[2*i])<0.5&&abs(val_tb[2*i+1])<0.5){
          ////  h_ctof_tb_fcut[i]     ->Fill(tagbctime[2*i]-tagbctime[2*i+1]);
          //}
          h2_ctof_tb_w_cut1[i] ->Fill(tagbctime[2*i]-tagbctime[2*i+1], tagbtime_w[2*i]);
          h2_ctof_tb_w_cut2[i] ->Fill(tagbctime[2*i]-tagbctime[2*i+1], tagbtime_w[2*i+1]);
        }
      }
    }

    //TagB - RF
    for(int i=0;i<24;i++){
      bunch[i]  =0;
      val_tb[i] =1000.;
      if(tagbctime[i]>-100){
        for(int k=0;k<166;k++){
          tmp_tb = RFtime[1]+k*(RFtime[0]-RFtime[1])/83. - tagbctime[i]-tb_offs[i];
          if(abs(tmp_tb)<abs(val_tb[i])){val_tb[i]=tmp_tb;bunch[i]=k;}
        }
 
        if(tagbtime_w[i]>0.){
          h_ctime_tb_RF[i]  ->Fill(val_tb[i]);
          h2_ctime_tb_RF[i] ->Fill(val_tb[i], tagbtime_w[i]);
        }
        if(tagbtime_w[i]>tb_w_min[i]&&tagbtime_w[i]<tb_w_max[i]){
          tb_flag[i]=true;
          if(abs(val_tb[i])<0.5)tb_rf_flag[i]=true;
          h_ctime_tb_RF_cut[i]  ->Fill(val_tb[i]);
          h2_ctime_tb_RF_cut[i] ->Fill(val_tb[i], tagbtime_w[i]);

          
          for(int j=0;j<8;j++){
            int tf_j = 4*i+j;
            if(i%2==0)tf_j =4*i+j -4;
            for(int l=0;l<tagf_lsize[tf_j];l++){
              if(fabs(tagftime_l[tf_j][l])<50.){
                tf_flag[i] = true;
                h_ctime_tb_RF_fcut[i][j]  ->Fill(val_tb[i]);
                h2_ctime_tb_RF_fcut[i][j] ->Fill(val_tb[i], tagbtime_w[i]);
                if(tf_j == 72)tagf73_flag = true;
                if(tf_j == 73)tagf74_flag = true;
                if(tf_j == 74)tagf75_flag = true;
                if(tf_j == 75)tagf76_flag = true;
              }
            }
          }
          if(tf_flag[i])h_ctime_tb_RF_f_or[i]         ->Fill(val_tb[i]);
          if(tbtb_flag[i/2])h_ctime_tb_RF_tbtb[i]     ->Fill(val_tb[i]);
          if(!tbtb_flag[i/2])h_ctime_tb_RF_tbtb_bar[i]     ->Fill(val_tb[i]);
          if(tbtb_flag[i/2] || tf_flag[i])h_ctime_tb_RF_fbor[i]     ->Fill(val_tb[i]);
        }
      }
      for(int j=0;j<10;j++){
        if(tagbtime_l[i][j]>-100){
          h_time_tb_RF[i]   ->Fill(RFtime[1] - tagbtime_l[i][j] +(double)bunch[i]*(RFtime[0]-RFtime[1])/83.);
          h2_time_tb_RF[i]  ->Fill(RFtime[1] - tagbtime_l[i][j] +(double)bunch[i]*(RFtime[0]-RFtime[1])/83., tagbtime_w[i]);
        }
      }

      if(tb_rf_flag[i]){
        for(int k=0;k<96;k++){
          for(int l=0;l<tagf_lsize[k];l++){
            if(fabs(tagftime_l[k][l])<50.)h_tagfhitpat_tb[i] ->Fill(k+1);
          }
        }
      }

      if(tb_flag[i])tb_total[i]++;
      if(tf_flag[i])tf_total[i]++;
    }


    if(tagf73_flag)ntf73++;
    if(tagf74_flag)ntf74++;
    if(tagf75_flag)ntf75++;
    if(tagf76_flag)ntf76++;
    if(tagf73_flag&&tb_rf_flag[18])ntf73_tb19++; 
    if(tagf74_flag&&tb_rf_flag[18])ntf74_tb19++; 
    if(tagf75_flag&&tb_rf_flag[18])ntf75_tb19++; 
    if(tagf76_flag&&tb_rf_flag[18])ntf76_tb19++;
    if(tagf73_flag&&tb_rf_flag[19])ntf73_tb20++;
    if(tagf74_flag&&tb_rf_flag[19])ntf74_tb20++;
    if(tagf75_flag&&tb_rf_flag[19])ntf75_tb20++;
    if(tagf76_flag&&tb_rf_flag[19])ntf76_tb20++;
    if(tagf73_flag&&tb_rf_flag[18]&&tb_rf_flag[19])ntf73_tb1920++;
    if(tagf74_flag&&tb_rf_flag[18]&&tb_rf_flag[19])ntf74_tb1920++;
    if(tagf75_flag&&tb_rf_flag[18]&&tb_rf_flag[19])ntf75_tb1920++;
    if(tagf76_flag&&tb_rf_flag[18]&&tb_rf_flag[19])ntf76_tb1920++;

   //TDL - RF
    for(int i=0;i<NTDL;i++){
      bunch[i]  =0;
      val_tdll[i] =1000.;
      if(tdlluctime[i]>-100 && tdlldctime[i]>-100){
        for(int k=0;k<166;k++){
          tmp_tdll = RFtime[1]+k*(RFtime[0]-RFtime[1])/83. - 0.5*(tdlluctime[i]+tdlldctime[i])-tdll_offs[i];
          if(abs(tmp_tdll)<abs(val_tdll[i])){val_tdll[i]=tmp_tdll;bunch[i]=k;}
        }
 
        if(i==23){
          if(tdllutime_w[i]>tdlu_min && tdllutime_w[i]<tdlu_max && 
             tdlldtime_w[i]>tdld_min && tdlldtime_w[i]<tdld_max && 
             tdllutime_wn[23]>tdllutime_wn[20] &&
             tdllutime_wn[23]>tdllutime_wn[21] &&
             tdllutime_wn[23]>tdllutime_wn[22]
          ){
            h_ctime_tdll_RF[i]  ->Fill(val_tdll[i]);
            h2_ctime_tdllu_RF[i] ->Fill(val_tdll[i], tdllutime_w[i]);
            h2_ctime_tdlld_RF[i] ->Fill(val_tdll[i], tdlldtime_w[i]);
            //cout<<val_tdll[i]<<endl;
          }
        }
        else{
          h_ctime_tdll_RF[i]  ->Fill(val_tdll[i]);
          h2_ctime_tdllu_RF[i] ->Fill(val_tdll[i], tdllutime_w[i]);
          h2_ctime_tdlld_RF[i] ->Fill(val_tdll[i], tdlldtime_w[i]);
        }
      }
    }

/*
     //TDL
    for(int i=0;i<NTDL;i++){
      bunch[i]  =0;
      val_tdll =1000.;

      //TDLL
      if(tdllmctime[i]>-100){
        for(int k=0;k<166;k++){
          tmp_tdll = RFtime[1]+k*(RFtime[0]-RFtime[1])/83. - tdllmctime[i]-tdll_offs[i];
          if(abs(tmp_tdll)<abs(val_tdll)){val_tdll=tmp_tdll;bunch[i]=k;}
        }
 
        h_ctime_tdll_RF[i]  ->Fill(val_tdll);
        h2_ctime_tdllu_RF[i] ->Fill(val_tdll, tdllutime_w[i]);
        h2_ctime_tdlld_RF[i] ->Fill(val_tdll, tdlldtime_w[i]);
        if(tdllutime_w[i]>tdllu_w_min[i]&&tdllutime_w[i]<tdllu_w_max[i]&&tdlldtime_w[i]>tdlld_w_min[i]&&tdlldtime_w[i]<tdlld_w_max[i]){
          tdll_flag[i]=true;
          h_ctime_tdll_RF_cut[i]  ->Fill(val_tdll);
          h2_ctime_tdllu_RF_cut[i] ->Fill(val_tdll, tdllutime_w[i]);
          h2_ctime_tdlld_RF_cut[i] ->Fill(val_tdll, tdlldtime_w[i]);
        }
      }

      for(int j=0;j<10;j++){
        if(tdllutime_l[i][j]>-100&&tdlldtime_l[i][j]>-100){
          h_time_tdll_RF[i]    ->Fill(RFtime[1] - (tdllutime_l[i][j]+tdlldtime_l[i][j])/2. +(double)bunch[i]*(RFtime[0]-RFtime[1])/83.);
          h2_time_tdllu_RF[i]  ->Fill(RFtime[1] - (tdllutime_l[i][j]+tdlldtime_l[i][j])/2. +(double)bunch[i]*(RFtime[0]-RFtime[1])/83., tdllutime_w[i]);
          h2_time_tdlld_RF[i]  ->Fill(RFtime[1] - (tdllutime_l[i][j]+tdlldtime_l[i][j])/2. +(double)bunch[i]*(RFtime[0]-RFtime[1])/83., tdllutime_w[i]);
        }
      }

    }
   


    for(int i=0;i<24;i++){//TagB
      if(tdll_flag[1]){
        if(tb_flag[i]){
        }
      }

      for(int j=0;j<NTDL;j++){//TDL
        if(tdll_flag[j]){
          if(tb_flag[i]){
            h_cointime_tdll[j]         ->Fill(tagbctime[i]-tdllmctime[j]);
            h_cointime_tb[i]           ->Fill(tagbctime[i]-tdllmctime[j]);
            h2_cointime_tb_wid[i]      ->Fill(tagbctime[i]-tdllmctime[j],tagbtime_w[i]);
            h2_cointime_tdll_uwid[j]   ->Fill(tagbctime[i]-tdllmctime[j],tdllutime_w[j]);
            h2_cointime_tdll_dwid[j]   ->Fill(tagbctime[i]-tdllmctime[j],tdlldtime_w[j]);
          }
        }
      }
    }
*/
  }//event loop

  double eff_73 =(double)ntf73_tb1920/(double)ntf73_tb20;
  double eff_74 =(double)ntf74_tb1920/(double)ntf74_tb20;
  double eff_75 =(double)ntf75_tb1920/(double)ntf75_tb20;
  double eff_76 =(double)ntf76_tb1920/(double)ntf76_tb20;
  //error calculated by binary distribution
  double er_73  = sqrt((double)ntf73_tb20 * (1.-eff_73))/(double)ntf73_tb20;
  double er_74  = sqrt((double)ntf74_tb20 * (1.-eff_74))/(double)ntf74_tb20;
  double er_75  = sqrt((double)ntf75_tb20 * (1.-eff_75))/(double)ntf75_tb20;
  double er_76  = sqrt((double)ntf76_tb20 * (1.-eff_76))/(double)ntf76_tb20;

  //cout<<"TagF73(x)TagB "<<ntf73_tb1920 << " / "<<ntf73_tb20<<" = "<<eff_73<<" +/- "<<er_73<<endl;
  //cout<<"TagF74(x)TagB "<<ntf74_tb1920 << " / "<<ntf74_tb20<<" = "<<eff_74<<" +/- "<<er_74<<endl;
  //cout<<"TagF75(x)TagB "<<ntf75_tb1920 << " / "<<ntf75_tb20<<" = "<<eff_75<<" +/- "<<er_75<<endl;  
  //cout<<"TagF76(x)TagB "<<ntf76_tb1920 << " / "<<ntf76_tb20<<" = "<<eff_76<<" +/- "<<er_76<<endl;  

}
////////////////////////////////////////////////////////////////////////////
void RF_ana::fit(){
  TF1 *f_tdl = new TF1("f_tdl","gaus",-1,1);
  TF1 *f_tb[40];

  for(int i=0;i<24;i++){
    f_tb[i]= new TF1(Form("f_tb%d",i+1),"gaus",-2.5,2.5);
    set->SetTF1(f_tb[i],2,1,1.0 );
    h_ctime_tb_RF_cut[i]->Fit(f_tb[i],"0QR","",-2.5,2.5);

    peak_tb[i] = h_ctime_tb_RF[i]->GetBinCenter(h_ctime_tb_RF[i]->GetMaximumBin());
    max_tb[i]  = h_ctime_tb_RF[i]->GetBinContent(h_ctime_tb_RF[i]->GetMaximumBin());
    //cout<<"    "<<peak_tb[i]<<endl;  
  }
  TF1 *f_ct = new TF1("f_ct","gaus",-2,2);

}
////////////////////////////////////////////////////////////////////////////
void RF_ana::draw(){

  for(int i=0;i<14;i++){
    c[i]->Clear();
  }

  c[0]->Divide(6,4);
  for(int i=0;i<24;i++){
    c[0]->cd(i+1);gPad->SetLogy(1);h_ctime_tb_RF[i]->Draw("");
                                tl_tb->DrawLine(peak_tb[i],0.8,peak_tb[i],max_tb[i]);
  }

  c[1]->Divide(6,4);
  for(int i=0;i<24;i++){
    c[1]->cd(i+1);gPad->SetLogy(1);h_ctime_tb_RF_cut[i]->Draw("");
                                   //h_ctime_tb_RF_fbor[i]   ->Draw("same");
                                   h_ctime_tb_RF_tbtb[i]   ->Draw("same");
                                   h_ctime_tb_RF_tbtb_bar[i]   ->Draw("same");
                                   //h_ctime_tb_RF_f_or[i]   ->Draw("same");
                                   //h_ctime_tb_RF_fcut[i][0]->Draw("same");
                                   //h_ctime_tb_RF_fcut[i][1]->Draw("same");
                                   //h_ctime_tb_RF_fcut[i][2]->Draw("same");
                                   //h_ctime_tb_RF_fcut[i][3]->Draw("same");
  }

  c[2]->Divide(6,4);
  for(int i=0;i<24;i++){
    c[2]->cd(i+1);gPad->SetLogy(1);h_time_tb_RF[i]->Draw("");
  }

  c[3]->Divide(6,4);
  for(int i=0;i<24;i++){
    c[3]->cd(i+1);gPad->SetLogz(1);h2_ctime_tb_RF[i]->Draw("colz");
  }

  c[4]->Divide(6,4);
  for(int i=0;i<24;i++){
    c[4]->cd(i+1);gPad->SetLogz(1);//h2_ctime_tb_RF_fcut[i]->Draw("colz");
    //c[4]->cd(i+1);gPad->SetLogz(1);h2_ctime_tb_RF_cut[i]->Draw("colz");
  }

  c[5]->Divide(6,4);
  for(int i=0;i<24;i++){
    c[5]->cd(i+1);gPad->SetLogz(1);h2_time_tb_RF[i]->Draw("colz");
  }

  c[6]->Divide(5,4);
  for(int i=0;i<10;i++){
    c[6]->cd(i+1) ;gPad->SetLogz(1);h2_time_tdllu_RF[i]->Draw("colz");
    c[6]->cd(i+11);gPad->SetLogz(1);h2_time_tdlld_RF[i]->Draw("colz");
  }

  c[7]->Divide(5,4);
  for(int i=0;i<10;i++){
    c[7]->cd(i+1) ;gPad->SetLogz(1);h2_ctime_tdllu_RF[i]->Draw("colz");
    c[7]->cd(i+11);gPad->SetLogz(1);h2_ctime_tdlld_RF[i]->Draw("colz");
  }

  c[8]->Divide(5,5);
  for(int i=0;i<24;i++){
    c[8]->cd(i+1) ;gPad->SetLogy(1);h_ctime_tdll_RF[i]->Draw("");
    //c[8]->cd(i+11);gPad->SetLogy(1);h_time_tdll_RF[i] ->Draw("");
  }


  c[9]->Divide(6,4);
  c[9]->cd(1); gPad->SetLogy(1);h_time_RF            ->Draw("");
  c[9]->cd(2); gPad->SetLogy(1);h_time_RFdiff[0]     ->Draw("");

  c[9]->cd(12);gPad->SetLogy(1);h_ctof_tb[0]        ->Draw("");
  c[9]->cd(13);gPad->SetLogz(1);h2_ctof_tb_w1[0]    ->Draw("colz");
  c[9]->cd(14);gPad->SetLogz(1);h2_ctof_tb_w2[0]    ->Draw("colz");
  c[9]->cd(15);gPad->SetLogy(1);h_ctof_tb_cut[0]    ->Draw("");
  c[9]->cd(16);gPad->SetLogy(1);h_ctof_tb_fcut[0]   ->Draw("");
  c[9]->cd(17);gPad->SetLogz(1);h2_ctof_tb_w_cut1[0]->Draw("colz");
  c[9]->cd(18);gPad->SetLogz(1);h2_ctof_tb_w_cut2[0]->Draw("colz");

  c[10]->Divide(7,6);
  for(int i=0;i<6;i++){
    c[10]->cd(7*i+1); gPad->SetLogy(1);h_ctof_tb[i]        ->Draw("");
    c[10]->cd(7*i+2); gPad->SetLogz(1);h2_ctof_tb_w1[i]    ->Draw("colz");
    c[10]->cd(7*i+3); gPad->SetLogz(1);h2_ctof_tb_w2[i]    ->Draw("colz");
    c[10]->cd(7*i+4); gPad->SetLogy(1);h_ctof_tb_cut[i]    ->Draw("");
    c[10]->cd(7*i+5); gPad->SetLogy(1);h_ctof_tb_fcut[i]   ->Draw("");
    c[10]->cd(7*i+6); gPad->SetLogz(1);h2_ctof_tb_w_cut1[i]->Draw("colz");
    c[10]->cd(7*i+7); gPad->SetLogz(1);h2_ctof_tb_w_cut2[i]->Draw("colz");
  }

  c[11]->Divide(7,6);
  for(int i=0;i<6;i++){
    c[11]->cd(7*i+1); gPad->SetLogy(1);h_ctof_tb[i+6]        ->Draw("");
    c[11]->cd(7*i+2); gPad->SetLogz(1);h2_ctof_tb_w1[i+6]    ->Draw("colz");
    c[11]->cd(7*i+3); gPad->SetLogz(1);h2_ctof_tb_w2[i+6]    ->Draw("colz");
    c[11]->cd(7*i+4); gPad->SetLogy(1);h_ctof_tb_cut[i+6]    ->Draw("");
    c[11]->cd(7*i+5); gPad->SetLogy(1);h_ctof_tb_fcut[i+6]   ->Draw("");
    c[11]->cd(7*i+6); gPad->SetLogz(1);h2_ctof_tb_w_cut1[i+6]->Draw("colz");
    c[11]->cd(7*i+7); gPad->SetLogz(1);h2_ctof_tb_w_cut2[i+6]->Draw("colz");
  }


  c[12]->Divide(6,4);
  for(int i=0;i<12;i++){
    c[12]->cd(2*i+1); gPad->SetLogz(1);h2_ctof_tb_tdc1[i]    ->Draw("colz");
    c[12]->cd(2*i+2); gPad->SetLogz(1);h2_ctof_tb_tdc2[i]    ->Draw("colz");
  }

  c[13]->Divide(1,1);
  c[13]->cd(1);h_tagfhitpat_tb[23] ->Draw("");
  for(int i=0;i<23;i++){
    h_tagfhitpat_tb[i] ->Draw("same");
  }
} 
////////////////////////////////////////////////////////////////////////////
void RF_ana::savecanvas(string name){
  c[0] ->Print(Form("%s[",name.c_str()) );
  c[0] ->Print(Form("%s" ,name.c_str()) );
  c[1] ->Print(Form("%s" ,name.c_str()) );
  c[2] ->Print(Form("%s" ,name.c_str()) );
  c[3] ->Print(Form("%s" ,name.c_str()) );
  c[4] ->Print(Form("%s" ,name.c_str()) );
  c[5] ->Print(Form("%s" ,name.c_str()) );
  c[6] ->Print(Form("%s" ,name.c_str()) );
  c[7] ->Print(Form("%s" ,name.c_str()) );
  c[8] ->Print(Form("%s" ,name.c_str()) );
  c[9] ->Print(Form("%s" ,name.c_str()) );
  c[10]->Print(Form("%s" ,name.c_str()) );
  c[11]->Print(Form("%s" ,name.c_str()) );
  c[12]->Print(Form("%s" ,name.c_str()) );
  c[13]->Print(Form("%s" ,name.c_str()) );
  c[13]->Print(Form("%s]",name.c_str()) );
  cout<<name<<" saved"<<endl;
  ofp->Write(); ofp->Close();
}
////////////////////////////////////////////////////////////////////////////
void RF_ana::show_tfratio(){

  for(int i=0;i<24;i++){
    double tf_ratio = (double)tf_total[i]/(double)tb_total[i];
    cout<<"TagB"<<i+1<<" "<<tf_ratio<<" +/- "<<sqrt((double)tb_total[i]*tf_ratio*(1.-tf_ratio))/(double)tb_total[i]<<endl;
  }
}
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "output0000.dat";
  string ofname = "root/RF_ana.root";
  string pdfname = "pdf/RF_ana.pdf";
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
    case 'p':
      pdfname = optarg;
      break;
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output root filename"<<endl;
      cout<<"-n : maximum number of analysed events"<<endl;
      cout<<"-p : print pdf filename"<<endl;
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
  RF_ana *ana = new RF_ana();

  ana->SetMaxEvent(MaxNum);
  ana->SetRoot(ifname);
  ana->SetOfname(ofname);
  ana->makehist();
  ana->loop();
  ana->fit();
  ana->draw();
  ana->savecanvas(pdfname);
  ana->show_tfratio();
  delete ana;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}

