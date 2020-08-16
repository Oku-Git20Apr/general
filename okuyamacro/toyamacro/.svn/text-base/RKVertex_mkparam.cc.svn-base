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
#include "RKVertex_tr.h"
#include "RKV_ParamMan.h"

static const double PI = 4.0*atan(1.);
static const double mrad_to_deg = 1./1000*180./PI;
const double Mp = 938.272046;          // proton       mass (MeV/c2)
const double Mpi = 139.57018;          // charged pion mass (MeV/c2)
const double MK = 493.677;             // charged Kaon mass (MeV/c2)
const double c = 0.299792458;          // speed of light in vacuum (m/ns)

const double MaxCoinTime =  1.3;//IH-Tagger Coincidence Timing[ns]
const double MinCoinTime = -1.3;//IH-Tagger Coincidence Timing[ns]
const int min_stat_for_fit = 400;//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class RKVertex_mkparam : public RKVertex_tr
{
 public:
         RKVertex_mkparam();
        ~RKVertex_mkparam();
  void makehist(string ofname);
  void SetPIDfunc();
  void loop();
  void fit_TagB();
  void fit_TDL();
  void fit_TDL2_6();
  void fit_TDL7_9();
  void fit_OHVFW();
  void fit_OHVBW();
  void fit_OHHmid();
  void fit_OHHtil();
  void draw(); 
  void savecanvas(); 
  void SetRoot(string ifname); 
  void GetHist(string ifname); 
  void SetMaxEvent( int N )  { ENumMax = N; }
  void SetPdfFilename(string ifname); 
  void SetInputParam(string ifname); 
  void SetOutputParam(string ifname); 
  void WriteParam();
  double MissingMass(double Eg, TLorentzVector P1, TLorentzVector P2);
  int TagBHit(double *TagBTime, double IHTime, int *seg);
  bool SetEgamma(void);
  void FitGaus(TH1F *h, double &gamin, double &gamax, double range=2.5,int itr = 10);

  bool Proton(int charge, double momentum, double beta);
  bool PiMinus(int charge, double momentum, double beta);
  bool PiPlus(int charge, double momentum, double beta);
  Settings *set;
  RKV_ParamMan *ParamMan;

  private:
    TFile *ofp ;
    string pdf_name;
    string input_param;
    string output_param;
    int GetMaxEvent() { return ENumMax; }
    int ENumMax;
    int ENum;
    int runmin,runmax;
    //PID Mass Range
    double MassRange_proton_min;
    double MassRange_proton_max;
    double MassRange_pion_min;
    double MassRange_pion_max;
    double MassRange_K_min;
    double MassRange_K_max;
    double TPE[40]; //Tagged Photon Energy
    double oa_min;
    double oa_max;

    double gamin, gamax;
    double mean,e_mean,p0;//p0->default param
    

    TLorentzVector P_inv;    //lorentz vector of invariant mass
    TLorentzVector P_pi, P_p;//lorentz vector of pi^{-] and proton
    TVector3 vec[3];//mom(vector) of pi^{-] and proton

    TH1F *h_TDiffOHVL_IHLR2[12], *h_TDiffOHVR_IHLR2[12];
    TH1F *h_TDiffOHVL_IHL[12][10], *h_TDiffOHVL_IHR[12][10], *h_TDiffOHVR_IHL[12][10], *h_TDiffOHVR_IHR[12][10];
    TH1F *h_TDiffOHHL_IHL[9][10], *h_TDiffOHHR_IHR[9][10];
    TH1F *h_TDiffOHVL[12], *h_TDiffOHVR[12], *h_TDiffOHHR[9], *h_TDiffOHHL[9];

    TH1F *h_cointime_tb[40], *h_cointime_tdlr[10], *h_cointime_tdll[10];
    TH2F *h2_cointime_tdlseg, *h2_cointime_tagseg;

    TH1F *h_taggertime, *h_tdlmctime, *h_ihctime, *h_ohctime;
    TH2F *h2_tagger_tdl, *h2_tdl_ih, *h2_ih_tagger, *h2_ih_oh;
    TH2F *h2_tagger_oh_c1, *h2_tagger_oh, *h2_tagger_tdiff, *h2_oh_tdiff;


    TH2F *h2_dvol, *h2_dvol_oa, *h2_dvol_ppi;
    TH2F *h2_im_vs_mm;
    TH1F *h_inv_mass, *h_mis_mass;
    TH1F *h_invm_eg1, *h_invm_eg2, *h_invm_eg3;
    TH1F *h_oa_ee, *h_oa_ppi;

    TH1F *h_ctime_allpi, *h_ctime_ppipi;
    TH1F *h_fl;
    TH2F *h_ctime_uwid;
    TH1F *h_tag_nhit;
    TH1F *h_tag_seg; 

    TH1F *h_msqr, *h_msqr_decut;
    TH2F *h2_pid_ohv;
    TH2F *h2_pid_ohh;
    TH2F *h2_pid_all;
    TH2F *h2_pid_ppi;
    TH2F *h2_pid_tag;
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

    TH2F *h_frame;
    TLatex *tex;
    TLine *line;
    int run_num;
    double beta[3];
    TCanvas *c1,*c2,*c3,*c4,*c5;
    TCanvas *cc[12];
    TLatex *tex_fl[4][12];
    TF1 *f_pi, *f_k, *f_p;
    TF1 *f_pi_min, *f_p_min;
    TF1 *f_pi_max, *f_p_max;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
RKVertex_mkparam::RKVertex_mkparam()
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

  MassRange_proton_min = 1.0;
  MassRange_proton_max = 1.5;
  MassRange_pion_min   = 0.0;
  MassRange_pion_max   = 0.8;       
  MassRange_K_min = 0.45;
  MassRange_K_max = 0.55;
  oa_min = -0.95;
  oa_max =  0.95;

  for(int i=0;i<12;i++){
  cc[i]= new TCanvas(Form("c_%d",i+1),Form("c_%d",i+1),1400,800 );
  }


  tex  = new TLatex(0.5,0.7,"cut condition");
  tex  -> SetTextSize(0.050);
  tex  -> SetTextAlign(22);

  line = new TLine(0,0,1,1);
  line ->SetLineColor(6);
  line ->SetLineWidth(1);
  line ->SetLineStyle(1);

  
  set = new Settings();
  ParamMan = new RKV_ParamMan();
}
////////////////////////////////////////////////////////////////////////////
RKVertex_mkparam::~RKVertex_mkparam(){
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::SetPIDfunc(){
  f_pi = new TF1("f_pi","[0]/sqrt(x*x-1)", 1,7);
  f_k  = new TF1("f_k" ,"[0]/sqrt(x*x-1)", 1,7);
  f_p  = new TF1("f_p" ,"[0]/sqrt(x*x-1)", 1,7);
  f_pi->SetParameter(0, 0.001*Mpi);
  f_k ->SetParameter(0, 0.001*MK);
  f_p ->SetParameter(0, 0.001*Mp);
  set->SetTF1(f_pi ,6,1,1.5);
  set->SetTF1(f_k  ,3,1,1.5);
  set->SetTF1(f_p  ,2,1,1.5);

  f_pi_min = new TF1("f_pi_min","[0]/sqrt(x*x-1)", 1,7);
  f_p_min  = new TF1("f_p_min" ,"[0]/sqrt(x*x-1)", 1,7);
  f_pi_min->SetParameter(0, MassRange_pion_min  );
  f_p_min ->SetParameter(0, MassRange_proton_min);
  set->SetTF1(f_pi_min ,6,1,1.5);
  set->SetTF1(f_p_min  ,2,1,1.5);

  f_pi_max = new TF1("f_pi_max","[0]/sqrt(x*x-1)", 1,7);
  f_p_max  = new TF1("f_p_max" ,"[0]/sqrt(x*x-1)", 1,7);
  f_pi_max->SetParameter(0, MassRange_pion_max  );
  f_p_max ->SetParameter(0, MassRange_proton_max);
  set->SetTF1(f_pi_max ,6,1,1.5);
  set->SetTF1(f_p_max  ,2,1,1.5);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::makehist(string ofname){
cout<<"makehist"<<endl;
  ofp = new TFile(Form("%s",ofname.c_str()),"recreate");

  h_frame = new TH2F("h_frame","h_frame",10,0,1,10,0,1);
  set->SetTH2(h_frame,"","","");
  h_frame->GetXaxis()->SetNdivisions(000);
  h_frame->GetYaxis()->SetNdivisions(000);
  h_frame->SetStats(0);

  ofp->mkdir("TDiffOHV"); ofp->cd("TDiffOHV");
  for(int i=0; i<12; i++){
    h_TDiffOHVL_IHLR2[i] =  new TH1F(Form("TDiffOHVL%d_IHLR2", i+1), Form("TDiff DC-(OHVL%d-IHLR2)", i+1), 200, -10., 10.);
    h_TDiffOHVR_IHLR2[i] =  new TH1F(Form("TDiffOHVR%d_IHLR2", i+1), Form("TDiff DC-(OHVR%d-IHLR2)", i+1), 200, -10., 10.);
  }

  for(int i=0; i<12; i++){
    h_TDiffOHVL[i] = new TH1F(Form("TDiffOHVL%d" , i+1), Form("TDiff DC-(OHVL%d-TagB)" , i+1), 200, -10., 10.);
    h_TDiffOHVR[i] = new TH1F(Form("TDiffOHVR%d" , i+1), Form("TDiff DC-(OHVR%d-TagB)" , i+1), 200, -10., 10.);
    for(int j=0; j<10; j++){
      h_TDiffOHVL_IHL[i][j] = new TH1F(Form("TDiffOHVL%d_IHL%d" , i+1, j+1), Form("TDiff DC-(OHVL%d-IHL%d)" , i+1, j+1), 200, -10., 10.);
      h_TDiffOHVL_IHR[i][j] = new TH1F(Form("TDiffOHVL%d_IHR%d" , i+1, j+1), Form("TDiff DC-(OHVL%d-IHR%d)" , i+1, j+1), 200, -10., 10.);
      h_TDiffOHVR_IHL[i][j] = new TH1F(Form("TDiffOHVR%d_IHL%d" , i+1, j+1), Form("TDiff DC-(OHVR%d-IHL%d)" , i+1, j+1), 200, -10., 10.);
      h_TDiffOHVR_IHR[i][j] = new TH1F(Form("TDiffOHVR%d_IHR%d" , i+1, j+1), Form("TDiff DC-(OHVR%d-IHR%d)" , i+1, j+1), 200, -10., 10.);
    }
  }
  ofp->cd();

  ofp->mkdir("TDiffOHH"); ofp->cd("TDiffOHH");
  for(int i=0;i<9;i++){//ohh seg
    h_TDiffOHHL[i] = new TH1F(Form("TDiffOHHL%d" , i+1), Form("TDiff DC-(OHHL%d-TagB)" , i+1), 200, -10., 10.);
    h_TDiffOHHR[i] = new TH1F(Form("TDiffOHHR%d" , i+1), Form("TDiff DC-(OHHR%d-TagB)" , i+1), 200, -10., 10.);
    for(int j=0;j<10;j++){//ih seg
      h_TDiffOHHL_IHL[i][j] = new TH1F(Form("TDiffOHHL%d_IHL%d",i+1,j+1), Form("TDiffOHHL%d_IHL%d",i+1,j+1), 200,-10,10);
      h_TDiffOHHR_IHR[i][j] = new TH1F(Form("TDiffOHHR%d_IHR%d",i+1,j+1), Form("TDiffOHHR%d_IHR%d",i+1,j+1), 200,-10,10);
    }
  }
  ofp->cd();

  ofp->mkdir("cointime"); ofp->cd("cointime");
  for(int i=0;i<10;i++){//TDL seg
    h_cointime_tdll[i] = new TH1F(Form("h_cointime_tdll%d",i+1),Form("h_cointime_tdll%d",i+1) , 400,-5,5);
    h_cointime_tdlr[i] = new TH1F(Form("h_cointime_tdlr%d",i+1),Form("h_cointime_tdlr%d",i+1) , 400,-5,5);
    set->SetTH1(h_cointime_tdll[i],Form("cointime TDLL%d",i+1)  ,"cointime[ns]"      ,"Counts/0.025ns"       , 1, 3001, 2);
    set->SetTH1(h_cointime_tdlr[i],Form("cointime TDLR%d",i+1)  ,"cointime[ns]"      ,"Counts/0.025ns"       , 1, 3001, 2);
  }
  for(int i=0;i<40;i++){//TagB seg
    h_cointime_tb[i] = new TH1F(Form("h_cointime_tb%d",i+1),Form("h_cointime_tb%d",i+1) , 200,-10,10);
    set->SetTH1(h_cointime_tb[i],Form("cointime TagB%d",i+1)  ,"cointime[ns]"      ,"Counts"       , 1, 3001, 4);
  }
  h_taggertime = new TH1F("h_taggertime"  ,"h_taggertime" , 400, -20., 20.);
  h_tdlmctime  = new TH1F("h_tdlmctime"   ,"h_tdlmctime"  , 400, -20., 20.);
  h_ihctime    = new TH1F("h_ihctime"     ,"h_ihctime"    , 400, -20., 20.);
  h_ohctime    = new TH1F("h_ohctime"     ,"h_ohctime"    , 400, -20., 20.);
  h2_cointime_tdlseg   = new TH2F("h2_cointime_tdlseg"    ,"h2_cointime_tdlseg"    ,  300, -10, 10,  21,   -10.5,    10.5);
  h2_cointime_tagseg   = new TH2F("h2_cointime_tagseg"    ,"h2_cointime_tagseg"    ,  300, -10, 10,  21,     0.5,    27.5);
  h2_tagger_tdl        = new TH2F("h2_tagger_tdl"         ,"h2_tagger_tdl"         ,  400,   -20., 20.,  400,    -20.,  20.);
  h2_tdl_ih            = new TH2F("h2_tdl_ih"             ,"h2_tdl_ih"             ,  400,   -20., 20.,  400,    -20.,  20.);
  h2_ih_tagger         = new TH2F("h2_ih_tagger"          ,"h2_ih_tagger"          ,  400,   -20., 20.,  400,    -20.,  20.);
  h2_ih_oh             = new TH2F("h2_ih_oh"              ,"h2_ih_oh"              ,  400,   -20., 20.,  400,    -20.,  20.);

  h2_tagger_oh   = new TH2F("h2_tagger_oh"     ,"h2_tagger_oh"        ,  400,   -20., 20.,   400,    -20.,   20.);
  h2_tagger_oh_c1= new TH2F("h2_tagger_oh_c1"  ,"h2_tagger_oh_c1"     ,  400,   -20., 20.,   400,    -20.,   20.);
  h2_tagger_tdiff= new TH2F("h2_tagger_tdiff"  ,"h2_tagger_tdiff"     ,  400,   -20., 20.,  1200,    -20.,  100.);
  h2_oh_tdiff    = new TH2F("h2_oh_tdiff"      ,"h2_oh_tdiff"         ,  400,   -20., 20.,  1200,    -20.,  100.);
  set->SetTH1(h_taggertime ,"tagger time" ,"time[ns]"        ,"Counts"       , 1, 1000, 0);
  set->SetTH1(h_tdlmctime  ,"TDL time"    ,"time[ns]"        ,"Counts"       , 1, 1000, 0);
  set->SetTH1(h_ihctime    ,"IH time"     ,"time[ns]"        ,"Counts"       , 1, 1000, 0);
  set->SetTH1(h_ohctime    ,"OH time"     ,"time[ns]"        ,"Counts"       , 1, 1000, 0);
  set->SetTH2(h2_cointime_tdlseg    ,"coin. time vs TDL seg"     ,"cointime[ns]"     ,"TDL seg"      );
  set->SetTH2(h2_cointime_tagseg    ,"coin. time vs TagB seg"    ,"cointime[ns]"     ,"TagB seg"     );
  set->SetTH2(h2_tagger_tdl         ,"Tagger_TDL"                ,"Tagger[ns]"       ,"TDL[ns]"      );
  set->SetTH2(h2_tdl_ih             ,"TDL_IH"                    ,"TDL[ns]"          ,"IH[ns]"       );
  set->SetTH2(h2_ih_tagger          ,"IH_Tagger"                 ,"IH[ns]"           ,"Tagger[ns]"   );
  set->SetTH2(h2_ih_oh              ,"IH_OH"                     ,"IH[ns]"           ,"OH[ns]"       );
  set->SetTH2(h2_tagger_oh          ,"Tagger_OH   "              ,"Tagger[ns]"       ,"OH[ns]"       );
  set->SetTH2(h2_tagger_oh_c1       ,"Tagger_OH (TDiff>20)"      ,"Tagger[ns]"       ,"OH[ns]"       );
  set->SetTH2(h2_tagger_tdiff       ,"Tagger_TDiff"              ,"Tagger[ns]"       ,"TDiff[ns]"    );
  set->SetTH2(h2_oh_tdiff           ,"OH_TDiff    "              ,"OH[ns]"           ,"TDiff[ns]"    );

  ofp->cd();


  ofp->mkdir("vertex"); ofp->cd("vertex");
  h2_dvol     = new TH2F("h2_dvol"    ,"h2_dvol"     , 200, -10., 10., 106, -10.6, 10.6);
  h2_dvol_oa  = new TH2F("h2_dvol_oa" ,"h2_dvol_oa"  , 200, -10., 10., 106, -10.6, 10.6);
  h2_dvol_ppi = new TH2F("h2_dvol_ppi","h2_dvol_ppi" , 200, -10., 10., 106, -10.6, 10.6);
  h2_im_vs_mm = new TH2F("h2_im_vs_mm","h2_im_vs_mm" ,  50,   1., 1.5,  35,    0.,  0.7);
  h_inv_mass  = new TH1F("h_inv_mass" ,"h_inv_mass", 200, 1., 1.5);
  h_invm_eg1  = new TH1F("h_invm_eg1" ,"h_invm_eg1", 200, 1., 1.5);
  h_invm_eg2  = new TH1F("h_invm_eg2" ,"h_invm_eg2", 200, 1., 1.5);
  h_invm_eg3  = new TH1F("h_invm_eg3" ,"h_invm_eg3", 200, 1., 1.5);
  h_mis_mass  = new TH1F("h_mis_mass" ,"h_mis_mass", 200, 0., 1.);
  h_oa_ee     = new TH1F("h_oa_ee"  ,"h_oa_ee" , 200, -1., 1.);
  h_oa_ppi    = new TH1F("h_oa_ppi" ,"h_oa_ppi", 200, -1., 1.);
  set->SetTH2(h2_dvol    ,"Vertex Position"           ,"x [cm]","y [cm]");
  set->SetTH2(h2_dvol_oa ,"Vertex Position(OA cut)"   ,"x [cm]","y [cm]");
  set->SetTH2(h2_dvol_ppi,"Vertex Position(p-#pi^{-})","x [cm]","y [cm]");
  set->SetTH2(h2_im_vs_mm,"M_{inv} vs M_{x} p#pi^{-}" ,"Invariant Mass [GeV/#it{c}^{2}]","Missing Mass [GeV/#it{c}^{2}]");
  set->SetTH1(h_inv_mass,"Invariant Mass"                       ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.5MeV/#it{c}^{2}"       , 1, 3000, 1);
  set->SetTH1(h_invm_eg1,"Invariant Mass(E_{#gamma}<1.0)"       ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.5MeV/#it{c}^{2}"       , 1, 3000, 1);
  set->SetTH1(h_invm_eg2,"Invariant Mass(E_{#gamma}=1.0-1.1)"   ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.5MeV/#it{c}^{2}"       , 1, 3000, 1);
  set->SetTH1(h_invm_eg3,"Invariant Mass(E_{#gamma}>1.1)"       ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.5MeV/#it{c}^{2}"       , 1, 3000, 1);
  set->SetTH1(h_mis_mass,"Missing Mass"     ,"Missing Mass [GeV/#it{c}^{2}]"        ,"Counts"       , 1, 3000, 1);
  set->SetTH1(h_oa_ee   ,"Opening Angle(e^{+} e{-})"     ,"cos#theta"        ,"Counts"       , 1, 3000, 1);
  set->SetTH1(h_oa_ppi  ,"Opening Angle(p #pi^{-})"      ,"cos#theta"        ,"Counts"       , 1, 3000, 1);
  ofp->cd();

  ofp->mkdir("track"); ofp->cd("track");

  h2_pid_all  = new TH2F("h2_pid_all"   ,"h2_pid_all"   , 500, 0, 7, 1000, -1, 1);
  h2_pid_ppi  = new TH2F("h2_pid_ppi"   ,"h2_pid_ppi"   , 500, 0, 7, 1000, -1, 1);
  h2_pid_tag  = new TH2F("h2_pid_tag"   ,"h2_pid_tag"   , 500, 0, 7, 1000, -1, 1);
  set->SetTH2(h2_pid_all  ,"PID plot all(OH-TDL)"   ,"1/#beta","momentum[GeV/#it{c}]");
  set->SetTH2(h2_pid_ppi  ,"PID plot ppi"           ,"1/#beta","momentum[GeV/#it{c}]");
  set->SetTH2(h2_pid_tag  ,"PID plot all (OH-TagB)" ,"1/#beta","momentum[GeV/#it{c}]");
  ofp->cd();

  ofp->mkdir("tag"); ofp->cd("tag");
  h_ctime_allpi   = new TH1F("h_ctime_allpi"   ,"h_ctime_allpi"   ,800,-10,10);
  h_ctime_ppipi   = new TH1F("h_ctime_ppipi"   ,"h_ctime_ppipi"   ,800,-10,10);
  h_tag_nhit = new TH1F("h_tag_nhit", "h_tag_nhit", 15, 0., 15.);
  h_tag_seg  = new TH1F("h_tag_seg",  "h_tag_seg" , 45, 0., 45.);
  h_ctime_uwid  = new TH2F("h_ctime_uwid" ,"h_ctime_uwid" ,200,-10,10,200,0,50);
  set->SetTH1(h_ctime_allpi     ,"IH#otimes Tag coin. time(all:#pi^{-})"      ,"cointime (ns)","counts / 25ps", 1, 3001, 3);
  set->SetTH1(h_ctime_ppipi     ,"IH#otimes Tag coin. time(p&#pi^{-}:#pi^{-})","cointime (ns)","counts / 25ps", 1, 3000, 1);
  set->SetTH1(h_tag_nhit  ,"TagB Multiplicity"               ,"Num. of Hit"  ,"Counts"       , 1, 3000, 1);
  set->SetTH1(h_tag_seg   ,"TagB Segment"                    ,"Segment"      ,"Counts"       , 1, 3000, 1);
  set->SetTH2(h_ctime_uwid,"Cointime vs Width"               ,"cointime (ns)","TDLU width (ns)");
  ofp->cd();

  ofp->mkdir("ih"); ofp->cd("ih");
  for(int i=0;i<10;i++){
    h2_pid_ihl[i] = new TH2F(Form("h2_pid_ihl%d",i+1) ,Form("h2_pid_ihl%d",i+1) , 500, 0, 9, 1000, -1, 1);
    h2_pid_ihr[i] = new TH2F(Form("h2_pid_ihr%d",i+1) ,Form("h2_pid_ihr%d",i+1) , 500, 0, 9, 1000, -1, 1);
    set->SetTH2(h2_pid_ihl[i],Form("PID plot ihL%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
    set->SetTH2(h2_pid_ihr[i],Form("PID plot ihR%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
  }
    h2_ihdedxbeta_all  = new TH2F("h2_ihdedxbeta_all"   ,"h2_ihdedxbeta_all"    , 500, 0, 2, 1000, 0, 9);
    set->SetTH2(h2_ihdedxbeta_all  ,"de/dx(IH) vs #beta all"         ,"#beta","momentum[GeV/#it{c}]");

  ofp->cd();
  ofp->mkdir("oh"); ofp->cd("oh");

    h2_pid_ohv = new TH2F("h2_pid_ohv"    ,"h2_pid_ohv"   , 500, 0, 7, 1000, -1, 1);
    h2_pid_ohh = new TH2F("h2_pid_ohh"    ,"h2_pid_ohh"   , 500, 0, 7, 1000, -1, 1);
    h2_ohdedxbeta_all  = new TH2F("h2_ohdedxbeta_all"   ,"h2_ohdedxbeta_all"    , 500, 0, 2, 1000, 0, 9);
    h2_ohdedxbeta_decut= new TH2F("h2_ohdedxbeta_decut" ,"h2_ohdedxbeta_decut"  , 500, 0, 2, 1000, 0, 9);
    h2_ohdedxbeta_ohv  = new TH2F("h2_ohdedxbeta_ohv"   ,"h2_ohdedxbeta_ohv"    , 500, 0, 2, 1000, 0, 9);
    h2_ohdedxbeta_ohh  = new TH2F("h2_ohdedxbeta_ohh"   ,"h2_ohdedxbeta_ohh"    , 500, 0, 2, 1000, 0, 9);
    set->SetTH2(h2_pid_ohv  ,"PID plot OHV"             ,"1/#beta","momentum[GeV/#it{c}]");
    set->SetTH2(h2_pid_ohh  ,"PID plot OHH"             ,"1/#beta","momentum[GeV/#it{c}]");
    set->SetTH2(h2_ohdedxbeta_all  ,"de/dx(OH) vs #beta all"         ,"#beta","momentum[GeV/#it{c}]");
    set->SetTH2(h2_ohdedxbeta_decut,"de/dx(OH) vs #beta all(dE cut)" ,"#beta","momentum[GeV/#it{c}]");
    set->SetTH2(h2_ohdedxbeta_ohv  ,"de/dx(OH) vs #beta OHV1-8"      ,"#beta","momentum[GeV/#it{c}]");
    set->SetTH2(h2_ohdedxbeta_ohh  ,"de/dx(OH) vs #beta OHH4-6"      ,"#beta","momentum[GeV/#it{c}]");


  for(int i=0;i<12;i++){
    h2_pid_ohvl[i] = new TH2F(Form("h2_pid_ohvl%d",i+1) ,Form("h2_pid_ohvl%d",i+1) , 500, 0, 9, 1000, -1, 1);
    h2_pid_ohvr[i] = new TH2F(Form("h2_pid_ohvr%d",i+1) ,Form("h2_pid_ohvr%d",i+1) , 500, 0, 9, 1000, -1, 1);
    h_fl_ohvl[i]  = new TH1F(Form("h_fl_ohvl%d",i+1) ,Form("h_fl_ohvl%d",i+1) , 1000, 0, 1200);
    h_fl_ohvr[i]  = new TH1F(Form("h_fl_ohvr%d",i+1) ,Form("h_fl_ohvr%d",i+1) , 1000, 0, 1200);
    set->SetTH2(h2_pid_ohvl[i],Form("PID plot OHVL%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
    set->SetTH2(h2_pid_ohvr[i],Form("PID plot OHVR%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
    set->SetTH1(h_fl_ohvl[i] ,Form("Flight length OHVL%d",i+1) ,"flight length[cm]","counts", 1, 3001, 3);
    set->SetTH1(h_fl_ohvr[i] ,Form("Flight length OHVR%d",i+1) ,"flight length[cm]","counts", 1, 3001, 4);
  }
  for(int i=0;i<9;i++){
    h2_pid_ohhl[i] = new TH2F(Form("h2_pid_ohhl%d",i+1) ,Form("h2_pid_ohhl%d",i+1) , 500, 0, 9, 1000, -1, 1);
    h2_pid_ohhr[i] = new TH2F(Form("h2_pid_ohhr%d",i+1) ,Form("h2_pid_ohhr%d",i+1) , 500, 0, 9, 1000, -1, 1);;
    h_fl_ohhl[i]  = new TH1F(Form("h_fl_ohhl%d",i+1) ,Form("h_fl_ohhl%d",i+1) , 1000, 0, 1200);
    h_fl_ohhr[i]  = new TH1F(Form("h_fl_ohhr%d",i+1) ,Form("h_fl_ohhr%d",i+1) , 1000, 0, 1200);
    set->SetTH2(h2_pid_ohhl[i],Form("PID plot OHHL%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
    set->SetTH2(h2_pid_ohhr[i],Form("PID plot OHHR%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
    set->SetTH1(h_fl_ohhl[i] ,Form("Flight length OHHL%d",i+1) ,"flight length[cm]","counts", 1, 3001, 8);
    set->SetTH1(h_fl_ohhr[i] ,Form("Flight length OHHR%d",i+1) ,"flight length[cm]","counts", 1, 3001, 7);
  }

}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::loop(){
cout<<"start loop"<<endl;
runmin=99999;runmax=-99;

int ppicounter=0;
int pcounter=0;
int picounter=0;
bool ppi_flag;//proton & pi^{-} vertex
bool oa_flag, tag_flag, dedx_flag, ohseg_flag;
int pi_id, p_id;//track ID number for pi- and proton
double Egamma, mm, inv_mass;
double cos_oa;
double TDiff;

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){
    tree->GetEntry(n);
    if(runnum>runmax)runmax=runnum;
    if(runnum<runmin)runmin=runnum;

    //////////////
    //Initialize//
    //////////////
    ppi_flag=false;
    pi_id = p_id =-1;
    oa_flag = tag_flag = dedx_flag=ohseg_flag=false;

    if(n%100000==0) cout<<n<<" / "<<ENum<<endl;

    h2_dvol    ->Fill(vertex[0],vertex[1]);//x,y
    //h2_dvol_oa ->Fill(vertex[0],vertex[1]);

    

      ///////////////
      //search pion//
      ///////////////
    for(int i=0;i<ntr;i++){
      if(ass_noh[i]<1||ass_nih[i]<1)continue;//if no hit in hodoscopes
      beta[i] = (fl_oh[i] - fl_ih[i])/(ohct[i] - ihct[i])/(c*100.);
      h2_pid_all->Fill(1./beta[i], mom[i]*charge[i]);
      if(PiMinus(charge[i],mom[i],beta[i])) {pi_id=i;}
      if(PiPlus(charge[i],mom[i],beta[i]))  {pi_id=i;}
      vec[i].SetXYZ(momx[i], momy[i], momz[i]);
    }//for each track


    //////////////
    //TagB info.//
    //////////////
    double tdlmctime,taggertime=-999.;
    if(pi_id>-1){
      if(ihlr[pi_id]==-1)    tdlmctime=tdllmctime[ihseg[pi_id]-1];
      else if(ihlr[pi_id]==1)tdlmctime=tdlrmctime[ihseg[pi_id]-1];
      int seg[3];
      int NHit = TagBHit(tagbctime, tdlmctime, seg);
      h_tag_nhit ->Fill(NHit);
      taggertime=tagbctime[seg[0]-1];
      Egamma = -1.; 
      for(int i=0;i<24;i++){
        if(tagbctime[i]>-100){
          h_ctime_allpi ->Fill(tdlmctime-tagbctime[i]);
          //if(ihlr[pi_id]==-1)    h_cointime_tdll[ihseg[pi_id]-1] ->Fill(tdlmctime-tagbctime[i]);
          //else if(ihlr[pi_id]==1)h_cointime_tdlr[ihseg[pi_id]-1] ->Fill(tdlmctime-tagbctime[i]);
        }
      }

      if(NHit == 1){
        Egamma = TPE[seg[0]-1];
        //h_ctime_allpi ->Fill(tagbctime[seg[0]-1]-tdlmctime);
        h_cointime_tb[seg[0]-1] ->Fill(tdlmctime-tagbctime[seg[0]-1]);
        if(ihlr[pi_id]==-1)    h_cointime_tdll[ihseg[pi_id]-1] ->Fill(tdlmctime-taggertime);
        else if(ihlr[pi_id]==1)h_cointime_tdlr[ihseg[pi_id]-1] ->Fill(tdlmctime-taggertime);
        h2_cointime_tdlseg ->Fill(tdlmctime-taggertime , ihseg[pi_id]*ihlr[pi_id]);
        h2_cointime_tagseg ->Fill(tdlmctime-taggertime , seg[0]);
        //h_ctime_allpi ->Fill(tagbctime[seg[0]]-ihct[pi_id]);
        if(ppi_flag)h_ctime_ppipi ->Fill(tagbctime[seg[0]]-ihct[pi_id]);
        tag_flag = true;
      }
      else if(NHit == 2 && abs(seg[0]-seg[1])==1){
        Egamma = (TPE[seg[0]-1]+TPE[seg[1]-1])/2; 
        //h_ctime_allpi ->Fill(tagbctime[seg[0]]-ihct[pi_id]);
        //h_ctime_allpi ->Fill(tagbctime[seg[0]-1]-tdlmctime);
        if(ppi_flag)h_ctime_ppipi ->Fill(tdlmctime-tagbctime[seg[0]-1]);
        tag_flag = true;
      }
      if(NHit == 1){h_tag_seg -> Fill(seg[0]);}
      else if((NHit == 2) && (abs(seg[0]-seg[1])==1)){h_tag_seg -> Fill((seg[0]+seg[1])/2.);}
    }

      ////////////////
      //recalc. beta//
      ////////////////
    for(int i=0;i<ntr;i++){
      if(ass_noh[i]<1||taggertime<-990)continue;//if no hit in OH
      //if(ass_noh[i]<1||!tag_flag)continue;//if no hit in OH
      //beta[i] = (fl_oh[i]-fl_ih[i])/(ohct[i] - tdlmctime)/(c*100.);
      beta[i] = (fl_oh[i])/(ohct[i] + taggertime)/(c*100.);
      //beta[i] = fl_oh[i]/(ohct[i] - taggertime)/(c*100.);
      h2_pid_tag  ->Fill(1./beta[i], mom[i]*charge[i]);
      //TDiff =(ohct[i]-taggertime);//
      //TDiff =fl_oh[i]/(c*100.) - (ohct[i]-taggertime);//
      //if(ass_nih[i]>0)TDiff =(fl_oh[i]-fl_ih[i])/(c*100.) - (ohct[i]-ihct[i]);//
      //else            TDiff =(fl_oh[i]-fl_ih[i])/(c*100.) - (ohct[i]-taggertime);//

      //TDiff =(fl_oh[i]-fl_ih[i])/(c*100.) - (ohct[i]+taggertime);//
      TDiff =(fl_oh[i])/(c*100.) - (ohct[i]+taggertime);//

      h_taggertime    ->Fill(taggertime);
      h_tdlmctime     ->Fill(tdlmctime );
      h_ihctime       ->Fill(ihct[i]);
      h_ohctime       ->Fill(ohct[i]);
      h2_tagger_tdl   ->Fill(taggertime , tdlmctime  );
      h2_tdl_ih       ->Fill(tdlmctime  , ihct[i]    );
      h2_ih_tagger    ->Fill(ihct[i]    , taggertime );
      h2_ih_oh        ->Fill(ihct[i]    , ohct[i]    );
      h2_tagger_oh    ->Fill(taggertime , ohct[i]    );
      if(TDiff>20.)h2_tagger_oh_c1 ->Fill(taggertime , ohct[i]    );
      h2_tagger_tdiff ->Fill(taggertime , TDiff      );
      h2_oh_tdiff     ->Fill(ohct[i]    , TDiff      );


      if(Proton(charge[i],mom[i],beta[i]) ) {p_id =i; pcounter++; }
      if(PiMinus(charge[i],mom[i],beta[i])) {pi_id=i; picounter++; }

        if(ohlr[i]==-1){//OHVL
          h_TDiffOHVL[ohseg[i]-1]->Fill(TDiff);
          h2_pid_ohvl[ohseg[i]-1]->Fill(1./beta[i], mom[i]*charge[i]);
        }
        if(ohlr[i]== 1){//OHVR
          h_TDiffOHVR[ohseg[i]-1] ->Fill(TDiff);
          h2_pid_ohvr[ohseg[i]-1] ->Fill(1./beta[i], mom[i]*charge[i]);
        }


        if(ohlr[i]==-2){//OHHL
          h_TDiffOHHL_IHL[ohseg[i]-1][ihseg[i]-1]->Fill(TDiff);
          h_TDiffOHHL[ohseg[i]-1] ->Fill(TDiff);
          h2_pid_ohh              ->Fill(1./beta[i], mom[i]*charge[i]);
          h2_pid_ohhl[ohseg[i]-1] ->Fill(1./beta[i], mom[i]*charge[i]);
        }
        if(ohlr[i]== 2){//OHHR
          h_TDiffOHHR_IHR[ohseg[i]-1][ihseg[i]-1]->Fill(TDiff);
          h_TDiffOHHR[ohseg[i]-1] ->Fill(TDiff);
          h2_pid_ohh              ->Fill(1./beta[i], mom[i]*charge[i]);
          h2_pid_ohhr[ohseg[i]-1] ->Fill(1./beta[i], mom[i]*charge[i]);
        }
      vec[i].SetXYZ(momx[i], momy[i], momz[i]);
    }//recalc. beta

    if(pi_id>-1&&p_id>-1)ppi_flag=true;
    


    //anlysis of proton & pi^{-} event
    if(ppi_flag){
      ppicounter++;
      P_pi.SetPxPyPzE(momx[pi_id], momy[pi_id], momz[pi_id], sqrt(0.000001*Mpi*Mpi +mom[pi_id]*mom[pi_id] ));
      P_p.SetPxPyPzE( momx[p_id] , momy[p_id] , momz[p_id] , sqrt(0.000001*Mp*Mp   +mom[p_id] *mom[p_id]    ));
      P_inv = P_p + P_pi;
      inv_mass = P_inv.M();
      cos_oa = vec[pi_id]*vec[p_id]/(vec[pi_id].Mag()*vec[p_id].Mag());

      h_oa_ppi   ->Fill(cos_oa);
      if(cos_oa>oa_min&&cos_oa<oa_max)oa_flag=true;
      h2_dvol_ppi->Fill(vertex[0],vertex[1]);
      h2_pid_ppi ->Fill(1./beta[pi_id], mom[pi_id]*charge[pi_id]);
      h2_pid_ppi ->Fill(1./beta[p_id] , mom[p_id]*charge[p_id]  );
    }//ppi_flag


   ////////////////
   //Missing mass//
   ////////////////
   if(tag_flag&&ppi_flag&&oa_flag){
     mm = MissingMass(Egamma, P_pi, P_p);
       if(mm<MassRange_K_max&&mm>MassRange_K_min){
         h_inv_mass  -> Fill(inv_mass);
         if(Egamma<1.0)            h_invm_eg1  -> Fill(inv_mass);
         if(Egamma<1.1&&Egamma>1.0)h_invm_eg2  -> Fill(inv_mass);
         if(Egamma>1.1)            h_invm_eg3  -> Fill(inv_mass);
       }
     h_mis_mass  -> Fill(mm);
     h2_im_vs_mm -> Fill(inv_mass, mm);
   }



  }//event loop
cout<<"pi- event: "<<picounter<<endl;
  ofp->Write();
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::fit_TagB(){

  for(int i=0; i<40; i++){//Fitting TagB
    FitGaus(h_cointime_tb[i],gamin,gamax);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_cointime_tb[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0( 0,18,i+1);
    if(e_mean<1.)ParamMan->SetT0( 0,18,i+1,p0+mean);
  }


}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::fit_TDL(){
int lr,cid;
cid=15;
  for(int i=0; i<10; i++){//Fitting TDLL
    lr=-1;
    if(h_cointime_tdll[i]->Integral()<min_stat_for_fit)continue;
    FitGaus(h_cointime_tdll[i],gamin,gamax);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_cointime_tdll[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0( lr,cid,i+1);
    if(e_mean<1.)ParamMan->SetT0( lr,cid,i+1,p0-1.*mean);
  }

  for(int i=0; i<10; i++){//Fitting TDLR
    lr= 1;
    if(h_cointime_tdlr[i]->Integral()<min_stat_for_fit)continue;
    FitGaus(h_cointime_tdlr[i],gamin,gamax);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_cointime_tdlr[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0( lr,15,i+1);
    if(e_mean<1.)ParamMan->SetT0( lr,15,i+1,p0-1.*mean);
  }

}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::fit_TDL2_6(){
int lr,cid;
cid=15;
  for(int i=1; i<6; i++){//Fitting TDLL
    lr=-1;
    if(h_cointime_tdll[i]->Integral()<min_stat_for_fit)continue;
    if(h_cointime_tdll[i]->Integral()<1000)h_cointime_tdll[i]->Rebin(2);
    FitGaus(h_cointime_tdll[i],gamin,gamax,2.0);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_cointime_tdll[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0( lr,cid,i+1);
    if(e_mean<1.)ParamMan->SetT0( lr,cid,i+1,p0-1.*mean);
  }

  for(int i=1; i<6; i++){//Fitting TDLR
    lr= 1;
    if(h_cointime_tdlr[i]->Integral()<min_stat_for_fit)continue;
    if(h_cointime_tdlr[i]->Integral()<1000)h_cointime_tdlr[i]->Rebin(2);
    FitGaus(h_cointime_tdlr[i],gamin,gamax,2.0);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_cointime_tdlr[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0( lr,15,i+1);
    if(e_mean<1.)ParamMan->SetT0( lr,15,i+1,p0-1.*mean);
  }

}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::fit_TDL7_9(){
int lr,cid;
cid=15;
  for(int i=6; i<9; i++){//Fitting TDLL
    lr=-1;
    if(h_cointime_tdll[i]->Integral()<min_stat_for_fit)continue;
    FitGaus(h_cointime_tdll[i],gamin,gamax,1.5);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_cointime_tdll[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0( lr,cid,i+1);
    if(e_mean<1.)ParamMan->SetT0( lr,cid,i+1,p0-1.*mean);
  }

  for(int i=6; i<9; i++){//Fitting TDLR
    lr= 1;
    if(h_cointime_tdlr[i]->Integral()<min_stat_for_fit)continue;
    FitGaus(h_cointime_tdlr[i],gamin,gamax,1.5);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_cointime_tdlr[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0( lr,15,i+1);
    if(e_mean<1.)ParamMan->SetT0( lr,15,i+1,p0-1.*mean);
  }

}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::fit_OHVFW(){
  int lr,cid;
  cid=2;
  lr=-1;
  for(int i=0; i<8; i++){//Fitting OHVL1-8 - TagB
    FitGaus(h_TDiffOHVL[i],gamin,gamax);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_TDiffOHVL[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0( lr,cid,i+1);
    if(e_mean<1.)ParamMan->SetT0(lr,cid,i+1,p0+mean);
  }

  lr=1;
  for(int i=0; i<8; i++){//Fitting OHVR1-8 - TagB
    FitGaus(h_TDiffOHVR[i],gamin,gamax);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_TDiffOHVR[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0( lr,cid,i+1);
    if(e_mean<1.)ParamMan->SetT0( lr,cid,i+1,p0+mean);
  }
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::fit_OHVBW(){
  int lr,cid;
  cid=2;
  lr=-1;
  for(int i=8; i<12; i++){//Fitting OHVL8-12 - TagB
    FitGaus(h_TDiffOHVL[i],gamin,gamax);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_TDiffOHVL[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0(lr,cid,i+1);
    if(e_mean<1.)ParamMan->SetT0(lr,cid,i+1,p0+mean);
  }

  lr=1;
  for(int i=8; i<12; i++){//Fitting OHVR8-12 - TagB
    FitGaus(h_TDiffOHVR[i],gamin,gamax);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_TDiffOHVR[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0( lr,cid,i+1);
    if(e_mean<1.)ParamMan->SetT0( lr,cid,i+1,p0+mean);
  }
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::fit_OHHmid(){//mid plane
  int lr,cid;
  cid=2;
  lr=-2;
  for(int i=3; i<6; i++){//Fitting OHHL4-6 - TagB
    FitGaus(h_TDiffOHHL[i],gamin,gamax);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_TDiffOHHL[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0(lr,cid,i+1);
    if(e_mean<1.)ParamMan->SetT0(lr,cid,i+1,p0+mean);
  }

  lr=2;
  for(int i=3; i<6; i++){//Fitting OHHR4-6 - TagB
    FitGaus(h_TDiffOHHR[i],gamin,gamax);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_TDiffOHHR[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0( lr,cid,i+1);
    if(e_mean<1.)ParamMan->SetT0( lr,cid,i+1,p0+mean);
  }
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::fit_OHHtil(){//tilted counters
  int lr,cid;
  cid=2;
  for(int i=0; i<9; i++){//Fitting OHHL1-3,7-9 - TagB
    if(i>2&&i<6)continue;
    lr=-2;
    FitGaus(h_TDiffOHHL[i],gamin,gamax);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_TDiffOHHL[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0(-1,2,i+1);
    if(e_mean<1.)ParamMan->SetT0(lr,cid,i+1,p0+mean);
  }

  for(int i=0; i<9; i++){//Fitting OHHR4-6 - TagB
    if(i>2&&i<6)continue;
    lr=2;
    FitGaus(h_TDiffOHHR[i],gamin,gamax);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_TDiffOHHR[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    p0=ParamMan->GetT0( lr,cid,i+1);
    if(e_mean<1.)ParamMan->SetT0( lr,cid,i+1,p0+mean);
  }
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::draw(){
  cc[0]->Divide(6,5);
double maxy=h_ctime_allpi ->GetBinContent(h_ctime_allpi ->GetMaximumBin());
  cc[0]->cd(1);h_ctime_allpi ->Draw("");line->DrawLine(MaxCoinTime,0,MaxCoinTime,maxy);line->DrawLine(MinCoinTime,0,MinCoinTime,maxy);
  cc[0]->cd(2);gPad->SetLogz(1);h2_pid_all  ->Draw("colz");
  cc[0]->cd(3);gPad->SetLogz(1);h2_pid_tag  ->Draw("colz");
  for(int i=0;i<12;i++){
    cc[0]->cd(i+7);   h_TDiffOHVL[i] -> Draw("");
    cc[0]->cd(i+12+7);h_TDiffOHVR[i] -> Draw("");
  }

  cc[1]->Divide(5,4);
  for(int i=0;i<9;i++){
    cc[1]->cd(i+1); h_TDiffOHHL[i] -> Draw("");
    cc[1]->cd(i+11);h_TDiffOHHR[i] -> Draw("");
  }

  cc[2]->Divide(5,3);
  cc[2]->cd(1); gPad->SetLogy(0);h_taggertime    ->Draw("");
  cc[2]->cd(2); gPad->SetLogy(0);h_tdlmctime     ->Draw("");
  cc[2]->cd(3); gPad->SetLogy(0);h_ihctime       ->Draw("");
  cc[2]->cd(4); gPad->SetLogy(0);h_ohctime       ->Draw("");
  cc[2]->cd(5); gPad->SetLogz(1);h2_tagger_tdl   ->Draw("colz");
  cc[2]->cd(6); gPad->SetLogz(1);h2_tdl_ih       ->Draw("colz");
  cc[2]->cd(7); gPad->SetLogz(1);h2_ih_tagger    ->Draw("colz");
  cc[2]->cd(8); gPad->SetLogz(1);h2_ih_oh        ->Draw("colz");
  cc[2]->cd(9); gPad->SetLogz(1);h2_tagger_oh    ->Draw("colz");
  cc[2]->cd(10);gPad->SetLogz(1);h2_tagger_tdiff ->Draw("colz");
  cc[2]->cd(11);gPad->SetLogz(1);h2_oh_tdiff     ->Draw("colz");
  cc[2]->cd(12);gPad->SetLogz(1);h2_tagger_oh_c1 ->Draw("colz");
  cc[2]->cd(13);gPad->SetLogz(1);h2_cointime_tdlseg ->Draw("colz");
  cc[2]->cd(14);gPad->SetLogz(1);h2_cointime_tagseg ->Draw("colz");

  cc[3]->Divide(4,6);
  for(int i=0;i<24;i++){
    cc[3]->cd(i+1);h_cointime_tb[i]->Draw("");
  }
  cc[4]->Divide(4,5);
  for(int i=0;i<10;i++){
    h_cointime_tdll[i]->GetXaxis()->SetRangeUser(-3,3);
    h_cointime_tdlr[i]->GetXaxis()->SetRangeUser(-3,3);
    cc[4]->cd(i+1); h_cointime_tdll[i]->Draw("");
    cc[4]->cd(i+11);h_cointime_tdlr[i]->Draw("");
  }
  cc[5]->Divide(4,6);
  for(int i=0;i<12;i++){
    cc[5]->cd(i+1); gPad->SetLogz(1); h2_pid_ohvl[i]->Draw("colz");
    cc[5]->cd(i+12);gPad->SetLogz(1); h2_pid_ohvr[i]->Draw("colz");
  }
//  for(int i=0;i<10;i++){
//  cc[i]->Divide(3,4);
//    for(int j=2;j<12;j++){
//      cc[i]->cd(i+1);h_TDiffOHVL_IHL[i][j] -> Draw("");
//    }
//  }

}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::SetPdfFilename(string ifname){
  pdf_name = ifname;
} 
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::SetInputParam(string ifname){
  input_param = ifname;
  ParamMan->read_param(input_param.c_str());
} 
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::SetOutputParam(string ifname){
  output_param = ifname;
} 
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::WriteParam(){
  ParamMan->write_param(output_param.c_str());
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam::savecanvas(){
  cc[0] ->Print(Form("%s[",pdf_name.c_str()) );
  cc[0] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[1] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[2] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[3] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[4] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[5] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[5] ->Print(Form("%s]",pdf_name.c_str()) );
  cout<<pdf_name<<" saved!"<<endl;
}
////////////////////////////////////////////////
void RKVertex_mkparam::SetRoot(string ifname)
{
  cout<<"SetRoot"<<endl;
  add(ifname);
  readtree();
  ENum = GetEntries();
}
////////////////////////////////////////////////
void RKVertex_mkparam::GetHist(string ifname)
{
  cout<<"Get Histogram from root file"<<endl;
  TFile *file = new TFile(ifname.c_str() );
  for(int i=0; i<12; i++){
    h_TDiffOHVL_IHLR2[i] =  (TH1F*)file->Get(Form("TDiffOHV/TDiffOHVL%d_IHLR2", i+1));
    h_TDiffOHVR_IHLR2[i] =  (TH1F*)file->Get(Form("TDiffOHV/TDiffOHVR%d_IHLR2", i+1));
  }

  for(int i=0; i<12; i++){
    h_TDiffOHVL[i] = (TH1F*)file->Get(Form("TDiffOHV/TDiffOHVL%d" , i+1));
    h_TDiffOHVR[i] = (TH1F*)file->Get(Form("TDiffOHV/TDiffOHVR%d" , i+1));
    for(int j=0; j<10; j++){
      h_TDiffOHVL_IHL[i][j] = (TH1F*)file->Get(Form("TDiffOHV/TDiffOHVL%d_IHL%d" , i+1, j+1));
      h_TDiffOHVL_IHR[i][j] = (TH1F*)file->Get(Form("TDiffOHV/TDiffOHVL%d_IHR%d" , i+1, j+1));
      h_TDiffOHVR_IHL[i][j] = (TH1F*)file->Get(Form("TDiffOHV/TDiffOHVR%d_IHL%d" , i+1, j+1));
      h_TDiffOHVR_IHR[i][j] = (TH1F*)file->Get(Form("TDiffOHV/TDiffOHVR%d_IHR%d" , i+1, j+1));
    }
  }

  for(int i=0;i<9;i++){//ohh seg
    h_TDiffOHHL[i] = (TH1F*)file->Get(Form("TDiffOHH/TDiffOHHL%d" , i+1));
    h_TDiffOHHR[i] = (TH1F*)file->Get(Form("TDiffOHH/TDiffOHHR%d" , i+1));
    for(int j=0;j<10;j++){//ih seg
      h_TDiffOHHL_IHL[i][j] =  (TH1F*)file->Get(Form("TDiffOHH/TDiffOHHL%d_IHL%d",i+1,j+1));
      h_TDiffOHHR_IHR[i][j] =  (TH1F*)file->Get(Form("TDiffOHH/TDiffOHHR%d_IHR%d",i+1,j+1));
    }
  }

  for(int i=0;i<10;i++){//TDL seg
    h_cointime_tdll[i] = (TH1F*)file->Get(Form("cointime/h_cointime_tdll%d",i+1));
    h_cointime_tdlr[i] = (TH1F*)file->Get(Form("cointime/h_cointime_tdlr%d",i+1));
  }
  for(int i=0;i<40;i++){//TagB seg
    h_cointime_tb[i] = (TH1F*)file->Get(Form("cointime/h_cointime_tb%d",i+1));
  }
  h_taggertime = (TH1F*)file->Get("cointime/h_taggertime" );
  h_tdlmctime  = (TH1F*)file->Get("cointime/h_tdlmctime"  );
  h_ihctime    = (TH1F*)file->Get("cointime/h_ihctime"    );
  h_ohctime    = (TH1F*)file->Get("cointime/h_ohctime"    );
  h2_tagger_tdl= (TH2F*)file->Get("cointime/h2_tagger_tdl");
  h2_tdl_ih    = (TH2F*)file->Get("cointime/h2_tdl_ih"    );
  h2_ih_tagger = (TH2F*)file->Get("cointime/h2_ih_tagger" );
  h2_ih_oh     = (TH2F*)file->Get("cointime/h2_ih_oh"     );

  h2_tagger_oh   =(TH2F*)file->Get("cointime/h2_tagger_oh"   );
  h2_tagger_oh_c1=(TH2F*)file->Get("cointime/h2_tagger_oh_c1");
  h2_tagger_tdiff=(TH2F*)file->Get("cointime/h2_tagger_tdiff");
  h2_oh_tdiff    =(TH2F*)file->Get("cointime/h2_oh_tdiff"    );
  h2_cointime_tdlseg   = (TH2F*)file->Get("cointime/h2_cointime_tdlseg" );
  h2_cointime_tagseg   = (TH2F*)file->Get("cointime/h2_cointime_tagseg" );

  h2_dvol     = (TH2F*)file->Get("vertex/h2_dvol"    );
  h2_dvol_oa  = (TH2F*)file->Get("vertex/h2_dvol_oa" );
  h2_dvol_ppi = (TH2F*)file->Get("vertex/h2_dvol_ppi");
  h2_im_vs_mm = (TH2F*)file->Get("vertex/h2_im_vs_mm");
  h_inv_mass  = (TH1F*)file->Get("vertex/h_inv_mass" );
  h_invm_eg1  = (TH1F*)file->Get("vertex/h_invm_eg1" );
  h_invm_eg2  = (TH1F*)file->Get("vertex/h_invm_eg2" );
  h_invm_eg3  = (TH1F*)file->Get("vertex/h_invm_eg3" );
  h_mis_mass  = (TH1F*)file->Get("vertex/h_mis_mass" );
  h_oa_ee     = (TH1F*)file->Get("vertex/h_oa_ee"    );
  h_oa_ppi    = (TH1F*)file->Get("vertex/h_oa_ppi"   );

  h2_pid_all  = (TH2F*)file->Get("track/h2_pid_all");
  h2_pid_ppi  = (TH2F*)file->Get("track/h2_pid_ppi");
  h2_pid_tag  = (TH2F*)file->Get("track/h2_pid_tag");

  h_ctime_allpi   = (TH1F*)file->Get("tag/h_ctime_allpi");
  h_ctime_ppipi   = (TH1F*)file->Get("tag/h_ctime_ppipi");
  h_tag_nhit      = (TH1F*)file->Get("tag/h_tag_nhit"       );
  h_tag_seg       = (TH1F*)file->Get("tag/h_tag_seg"        );
  h_ctime_uwid    = (TH2F*)file->Get("tag/h_ctime_uwid"     );

  for(int i=0;i<10;i++){
    h2_pid_ihl[i] = (TH2F*)file->Get(Form("ih/h2_pid_ihl%d",i+1));
    h2_pid_ihr[i] = (TH2F*)file->Get(Form("ih/h2_pid_ihr%d",i+1));
  }
    h2_ihdedxbeta_all  =  (TH2F*)file->Get("ih/h2_ihdedxbeta_all"  );

    h2_pid_ohv         = (TH2F*)file->Get("oh/h2_pid_ohv"          );
    h2_pid_ohh         = (TH2F*)file->Get("oh/h2_pid_ohh"          );
    h2_ohdedxbeta_all  = (TH2F*)file->Get("oh/h2_ohdedxbeta_all"   );
    h2_ohdedxbeta_decut= (TH2F*)file->Get("oh/h2_ohdedxbeta_decut" );
    h2_ohdedxbeta_ohv  = (TH2F*)file->Get("oh/h2_ohdedxbeta_ohv"   );
    h2_ohdedxbeta_ohh  = (TH2F*)file->Get("oh/h2_ohdedxbeta_ohh"   );


  for(int i=0;i<12;i++){
    h2_pid_ohvl[i] = (TH2F*)file->Get(Form("oh/h2_pid_ohvl%d",i+1));
    h2_pid_ohvr[i] = (TH2F*)file->Get(Form("oh/h2_pid_ohvr%d",i+1));
    h_fl_ohvl[i]   = (TH1F*)file->Get(Form("oh/h_fl_ohvl%d",i+1) );
    h_fl_ohvr[i]   = (TH1F*)file->Get(Form("oh/h_fl_ohvr%d",i+1) );
  }
  for(int i=0;i<9;i++){
    h2_pid_ohhl[i]= (TH2F*)file->Get(Form("oh/h2_pid_ohhl%d",i+1));
    h2_pid_ohhr[i]= (TH2F*)file->Get(Form("oh/h2_pid_ohhr%d",i+1));
    h_fl_ohhl[i]  = (TH1F*)file->Get(Form("oh/h_fl_ohhl%d",i+1)  );
    h_fl_ohhr[i]  = (TH1F*)file->Get(Form("oh/h_fl_ohhr%d",i+1)  );
  }
cout<<"end of gethist"<<endl;
}

////////////////////////////////////////////////
int RKVertex_mkparam::TagBHit(double *TagBTime, double IHTime, int *seg){

  seg[0] = seg[1] =seg[2] = -1;
  int NHit = 0;

  for(int i=0; i<27; i++){
    if( (IHTime-TagBTime[i])>MinCoinTime && (IHTime-TagBTime[i])<MaxCoinTime){
      if(NHit<2)seg[NHit] = i+1;
      NHit++;
    }
  }

  return NHit;
}
////////////////////////////////////////////////

bool RKVertex_mkparam::SetEgamma(void){

  for(int i=0; i<40; i++){TPE[i] = 0.;}

  string str;
  int id; double val;
  ifstream File("PhotonEnergy_1310MeV.dat");
  if(!File){cout << "SetEgamma: File Open Error" << endl; return false;}
  while(!File.fail()){
    getline(File, str);
    if((str.c_str())[0] != '#'){
      sscanf(str.c_str(), "%d %lf", &id, &val);
      TPE[id-1] = val;
    }
  }

  return true;

}

////////////////////////////////////////////////

double RKVertex_mkparam::MissingMass(double Eg, TLorentzVector P1, TLorentzVector P2){
  TLorentzVector P_gamma(Eg, 0., 0., Eg);
  TLorentzVector P_target(0., 0., 0., 0.001*Mp);

  TLorentzVector P_X = P_gamma + P_target - P1 - P2;
  return P_X.M();

}

////////////////////////////////////////////////

bool RKVertex_mkparam::Proton(int charge, double momentum, double beta){

  if(charge < 0){return false;}
  if(((momentum*sqrt(1-pow(beta, 2.))/beta) > MassRange_proton_min) &&
     ((momentum*sqrt(1-pow(beta, 2.))/beta) < MassRange_proton_max)){
    return true;
  }
  else return false;

}

////////////////////////////////////////////////

bool RKVertex_mkparam::PiMinus(int charge, double momentum, double beta){

  //momentum>0
  if(charge > 0){return false;}
  if(beta>1){return true;}
  else{
    if(((momentum*sqrt(1-pow(beta, 2.))/beta) > MassRange_pion_min) &&
       ((momentum*sqrt(1-pow(beta, 2.))/beta) < MassRange_pion_max)){
      //  if((momentum*sqrt(1-pow(beta, 2.))/beta) > MassRange_pion[1]){
      return true;
    }
  }
  return false;

}
////////////////////////////////////////////////

bool RKVertex_mkparam::PiPlus(int charge, double momentum, double beta){

  //momentum>0
  if(charge < 0){return false;}
  if(beta>1){return true;}
  else{
    if(((momentum*sqrt(1-pow(beta, 2.))/beta) > MassRange_pion_min) &&
       ((momentum*sqrt(1-pow(beta, 2.))/beta) < MassRange_pion_max)){
      //  if((momentum*sqrt(1-pow(beta, 2.))/beta) > MassRange_pion[1]){
      return true;
    }
  }
  return false;

}
//____________________________________________________________________________________________
void RKVertex_mkparam::FitGaus(TH1F *h, double &gamin, double &gamax, double range,int itr){
	 gamin = h  ->GetBinCenter(h ->GetMaximumBin())-2.0;
	 gamax = h  ->GetBinCenter(h ->GetMaximumBin())+2.0;
	for(Int_t l=0; l<itr; l++){
	TF1 *ga = new TF1("ga","gaus");
      ga->SetParameter(2,(gamin+gamax)/2.);
      h  ->Fit(ga,"0QR","",gamin,gamax);
      gamin = ga->GetParameter(1) - ga->GetParameter(2)*range;
      gamax = ga->GetParameter(1) + ga->GetParameter(2)*range;
      ga->Clear();
	}
 return;
 }
//////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "output0000.dat";
  string ofname = "root/hoge.root";
  string pdfname = "pdf/track/hoge.pdf";
  string input_param  = "param/default.param";
  string output_param = "param/0000.param";
  int ch;
  int MaxNum = 0;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  bool skip_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:w:n:bco:p:i:s"))!=-1){
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
    case 'i':
      input_param = optarg;
      break;
    case 'o':
      output_param = optarg;
      break;
    case 's':
      skip_flag = true;
      break;
    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output root filename"<<endl;
      cout<<"-n : maximum number of analysed events"<<endl;
      cout<<"-i : input param file"<<endl;
      cout<<"-o : output param file"<<endl;
      cout<<"-s : skip mode(fitting part only)"<<endl;
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
  RKVertex_mkparam *mkparam = new RKVertex_mkparam();

  mkparam->SetPIDfunc();
  mkparam->SetEgamma();
  mkparam->SetInputParam(input_param);//Set parameter file used to make tree
  mkparam->SetOutputParam(output_param);
  mkparam->SetMaxEvent(MaxNum);
  if(skip_flag)mkparam->GetHist(ifname);
  else{
    mkparam->SetRoot(ifname);
    mkparam->makehist(ofname);
    mkparam->loop();
  }
  //mkparam->fit_TagB();
  //mkparam->fit_TDL();
  //mkparam->fit_TDL2_6();
  mkparam->fit_TDL7_9();
  //mkparam->fit_OHVFW();
  //mkparam->fit_OHHmid();
  //mkparam->fit_OHVBW();
  //mkparam->fit_OHHtil();
  mkparam->draw();
  mkparam->SetPdfFilename(pdfname);
  mkparam->savecanvas();
  mkparam->WriteParam();
  delete mkparam;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}

