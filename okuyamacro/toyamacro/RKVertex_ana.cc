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
#include "ToF_ParamMan.h"

#define OutTree
//#undef OutTree

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

const double MaxCoinTime =  1.3;//IH-Tagger Coincidence Timing[ns]
const double MinCoinTime = -1.3;//IH-Tagger Coincidence Timing[ns]
const double MinTagB_w   = 11.;//TagB min Width[ns]
const double TagOffs     =  0.;//for accidental B.G. estimate
const double MaxChi2 =  1000.;//chi square of tracking

double vx_min = -3.2;
double vx_max =  24.5;
//double vx_max =  4.5;
double vy_min = -4.0;
double vy_max =  4.0;

double x1[2] = {0 , 1.1};
double yy[2] = {6.,  0.};
double x2[2] = {1. , 2.};
double y2[2] = {0. , 8.};

double tb_offs[24]={ 0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.8       ,
                     0.787073  ,
                     0.771146  ,
                     0.948757  ,
                     0.817174  };
/////////
struct TreeBranch{
  double Inv_mass, Mom_pi, Miss_mass;
  double Flength, Ftime, Dtime; 
};
static TreeBranch tr;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class RKVertex_ana : public RKVertex_tr
{
 public:
         RKVertex_ana();
        ~RKVertex_ana();
  void makehist(string ofname);
  void SetBranch();
  void SetPIDfunc();
  void loop();
  void loop_def();
  void fit();
  void draw(); 
  void savecanvas(); 
  void SetRoot(string ifname); 
  void SetRootSingle(string ifname); 
  void GetHist(string ifname); 
  void SetMaxEvent( int N )  { ENumMax = N; }
  void SetPdfFilename(string ifname); 
  double MissingMass(double Eg, TLorentzVector P1, TLorentzVector P2);
  int TagBHit(double *TagBTime, double *TagB_w, double IHTime, int *seg);
  bool SetEgamma(void);

  bool Proton(int charge, double mom, double beta);
  bool PiMinus(int charge, double mom, double beta);
  bool PiPulse(int charge, double mom, double beta);
  Settings *set;
  ToF_ParamMan *tofParamMan;

  private:
    TFile *ofp ;
    TTree *tree_out;
    string pdf_name;
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
    double MassRange_Lam_min;
    double MassRange_Lam_max;
    double MaxValue_of_chi2;
    double MaxValue_of_dca;
    double TPE[40]; //Tagged Photon Energy
    double oa_min;
    double oa_max;
    double coin_min[6];
    double m_sqr[3];//mass square
    double a1,a2,b1,b2;//constant value for dedx flag
    double a3,b3;//for pion beta flag
    double RF_offset;//[ns]
    double RFtime[5];
    int bunch[40];
    double val_tb[40];

    //cut condition study//
    double chi_max[12];
    double dca_max[12];
    double MassRange_Kth[12];//Kaon mass
    double p_min[12];//proton mass
    ////////////////////////
    
    bool tb_flag[40], tb_rf_flag[40];
    bool treeout_flag;

    TLorentzVector P_inv;    //lorentz vector of invariant mass
    TLorentzVector P_pi_minus, P_pi_plus, P_p;//lorentz vector of pi^{-] and proton
    TVector3 vec[3];//mom(vector) of pi^{-] and proton

    //vertex//
    TH2F *h2_dvol, *h2_dvol_oa, *h2_dvol_ppi, *h2_dvol_cut, *h2_dvol_Lam, *h2_dvol_bg;
    TH2F *h2_im_vs_mm;
    TH2F *h2_im_pimom;
    TH1F *h_oa_ee, *h_oa_ppi, *h_dca, *h_dca_cut, *h_dca_ppi;
    TH1F *h_inv_mass, *h_invm_lpi, *h_invm_Lam, *h_invm_bg, *h_mis_mass, *h_mm_pipi;
    TH1F *h_invm_mm[12];//missing mass cut dependence
    TH1F *h_invm_eg1, *h_invm_eg2, *h_invm_eg3;
    TH1F *h_invm_chi[12];//chi square cut dependence
    TH1F *h_invm_p_min[12];//missing mass cut dependence
    TH1F *h_invm_dca[12];
    TH1F *h_inv_coin[12];//cointime cut dependence
    TH1F *h_mom_Lam;//mom of Lambda

    //track//
    TH1F *h_chi2; 
    TH1F *h_chi2_pi[4], *h_chi2_p[4]; 
    TH1F *h_msqr, *h_msqr_a;
    TH1F *h_ft_pi;

    TH2F *h2_pid_all, *h2_pid_add;
    TH2F *h2_pid_ppi;
    TH2F *h2_pid_decut;
    TH2F *h2_pid_tag;
    TH2F *h2_fl_ftpi,*h2_beta_ftpi;
    TH2F *h2_beta_beta_ppi, *h2_beta_beta_pipi;
    TH2F *h2_beta_beta_all, *h2_beta_beta_cut;
    TH2F *h2_ft_chi2_pLam, *h2_ft_chi2_piLam;

    //cointime//
    TH1F *h_ctime_allpi, *h_ctime_allpi_wcut, *h_ctime_allpi_rf, *h_ctime_fastpi_rf, *h_ctime_ppipi,*h_ctime_ppipi_cor;
    TH1F *h_ctime_allpi_zoom, *h_ctime_allpi_zoom_rf, *h_ctime_fastpi_zoom_rf;//for drawing
    TH1F *h_ctime_Lam, *h_ctime_Lam_cor, *h_ctime_BG, *h_ctime_BG_cor, *h_ctime_pipi,*h_ctime_pipi_cor;
    TH2F *h2_ctime_xpos_Lam;
    TH2F *h2_ctime_ftih, *h2_ctime_ftih_bw, *h2_ctime_flih;
    TH2F *h2_ctime_ftpi, *h2_ctime_flpi, *h2_flpi_flih;
    TH1F *h_ctime_seg[40];
    TH1F *h_ctime_tb_RF[40], *h_ctime_tb_RF_cut[40], *h_ctime_tb_RF_fcut[40], *h_time_tb_RF[40];
    TH1F *h_ctime_tb_RF_all;
    TH2F *h2_ctime_tb_RF[40], *h2_ctime_tb_RF_cut[40], *h2_ctime_tb_RF_fcut[40], *h2_time_tb_RF[40];
    TH1F *h_cointime_proton;
    TH1F *h_ftpi_tdll[10], *h_ftpi_tdlr[10];
    TH2F *h2_ctime_ftih_tdll[10], *h2_ctime_ftih_tdllr[10];
    TH2F *h2_cointime_tdlseg, *h2_ctime_tdlseg_Lam;
    TH2F *h2_cointime_tdluw, *h2_cointime_tdldw, *h2_cointime_tagw, *h2_coincoin, *h2_cointime_oa, *h2_cointime_ohde,*h2_cointime_beta,*h2_cointime_beta_pi,*h2_cointime_ohde_a;
    TH2F *h2_cointime_tdluw_aa, *h2_cointime_tdldw_aa, *h2_cointime_tagw_aa;
    TH2F *h2_ctime_proton_mom;
    TH2F *h2_cointime_beta_tdll[10], *h2_cointime_beta_tdlr[10], *h2_cointime_beta_tagb[40];
    TH2F *h2_cointime_width_tdllu[10], *h2_cointime_width_tdlru[10], *h2_cointime_width_tagb[40];
    TH2F *h2_cointime_width_tdlld[10], *h2_cointime_width_tdlrd[10];
    TH2F *h2_cointime_ftpi_tdll[10], *h2_cointime_ftpi_tdlr[10], *h2_cointime_ftpi_tagb[40];
    TH2F *h2_cointime_runnum;
    TH1F *h_cointime_tb[40], *h_cointime_tdlr[10], *h_cointime_tdll[10];
    TH1F *h_time_RF, *h_time_RFdiff[4];

    //time at vertex//
    TH1F *h_time_at_ver_tdlr[10], *h_time_at_ver_tdll[10];
    TH2F *h2_time_at_ver_mom, *h2_time_at_ver_mom_ppi, *h2_time_at_ver_mom_pipi;

    //pipi//
    TH1F *h_cointime_tb_pipi[40], *h_cointime_tdlr_pipi[10], *h_cointime_tdll_pipi[10];
    TH1F *h_ftpi_tdll_pipi[10], *h_ftpi_tdlr_pipi[10];
    TH2F *h2_cointime_tdluw_pipi, *h2_cointime_tdldw_pipi, *h2_cointime_tagw_pipi, *h2_coincoin_pipi, *h2_cointime_ohde_pipi,*h2_cointime_beta_pipi,*h2_cointime_ohde_a_pipi;
    TH2F *h2_cointime_beta_tdll_pipi[10], *h2_cointime_beta_tdlr_pipi[10], *h2_cointime_beta_tagb_pipi[40];
    TH2F *h2_cointime_width_tdllu_pipi[10], *h2_cointime_width_tdlru_pipi[10], *h2_cointime_width_tagb_pipi[40];
    TH2F *h2_cointime_width_tdlld_pipi[10], *h2_cointime_width_tdlrd_pipi[10];
    TH2F *h2_cointime_ftpi_tdll_pipi[10], *h2_cointime_ftpi_tdlr_pipi[10], *h2_cointime_ftpi_tagb_pipi[40];

    //ih//
    TH1F *h_tdlseg,*h_tdlseg_pipi, *h_tdlseg_ppi, *h_tdlseg_Lam;
    TH2F *h2_pid_ihl[10], *h2_pid_ihr[10];
    TH2F *h2_ihdedxbeta_all;
    TH2F *h2_dedxbeta_ihl[10], *h2_dedxbeta_ihr[10];

    //oh//
    TH2F *h2_pid_ohv;
    TH2F *h2_pid_ohh;
    TH2F *h2_ohdedxbeta_all;
    TH2F *h2_ohdedxbeta_decut;
    TH2F *h2_ohdedxbeta_ohv;
    TH2F *h2_ohdedxbeta_ohh;
    TH2F *h2_dedxbeta_ohvl[12], *h2_dedxbeta_ohvr[12], *h2_dedxbeta_ohhl[9], *h2_dedxbeta_ohhr[9];
    TH2F *h2_pid_ohvl[12], *h2_pid_ohvr[12], *h2_pid_ohhl[9], *h2_pid_ohhr[9];
    TH1F *h_fl_ohvl[12], *h_fl_ohvr[12], *h_fl_ohhl[9], *h_fl_ohhr[9];


    //tag//
    TH2F *h2_ctime_uwid;
    TH1F *h_tag_nhit;
    TH1F *h_tag_seg; 
    TH1F *h_tdlseg_aa, *h_tagseg_aa, *h_runnum_aa;
    TH2F *h2_tagtdlseg_aa;

    TH2F *h_frame;

    TLatex *tex;
    TLine *line, *lineLam;
    int run_num;
    double beta[3];
    TCanvas *c1,*c2,*c3,*c4,*c5;
    TCanvas *cc[30];
    TLatex *tex_fl[4][12];
    TF1 *f_pi, *f_k, *f_p;
    TF1 *f_pi_min, *f_p_min;
    TF1 *f_pi_max, *f_p_max;

    TF1 *f_l[10], *f_r[10], *f[5];
    double param_l[10],param_r[10];
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
RKVertex_ana::RKVertex_ana()
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

  ////////////////////////////
  //finalized cut parameters//
  ////////////////////////////
  //MassRange_proton_min = 0.55;
  MassRange_proton_min = 0.75;
  MassRange_proton_max = 1.5;
  MassRange_pion_min   = 0.12;
  MassRange_pion_max   = 0.7;       
  MassRange_Lam_min   =  1.105;
  MassRange_Lam_max   =  1.122;       
  MassRange_K_min = 0.42;
  MassRange_K_max = 0.75;
  oa_min = -1.00;
  oa_max =  0.92;
  RF_offset = 0.;
  //RF_offset = -165.950 * 2.;

  for(int i=0;i<30;i++){
  cc[i]= new TCanvas(Form("c_%d",i+1),Form("c_%d",i+1),1400,800 );
  }


  tex  = new TLatex(0.5,0.7,"cut condition");
  tex  -> SetTextSize(0.050);
  tex  -> SetTextAlign(22);

  line = new TLine(0,0,1,1);
  line ->SetLineColor(6);
  line ->SetLineWidth(1);
  line ->SetLineStyle(1);

  lineLam = new TLine(0,0,1,1);
  lineLam ->SetLineColor(1);
  lineLam ->SetLineWidth(1);
  lineLam ->SetLineStyle(1);

  
  set = new Settings();
  tofParamMan = new ToF_ParamMan();
  tofParamMan ->read_param("param/tof/tof.param");
    chi_max[0] =200.; dca_max[0] = 50.; MassRange_Kth[0] =0.30; p_min[0] = 0.600; 
    chi_max[1] =180.; dca_max[1] = 45.; MassRange_Kth[1] =0.32; p_min[1] = 0.625; 
    chi_max[2] =160.; dca_max[2] = 40.; MassRange_Kth[2] =0.34; p_min[2] = 0.650; 
    chi_max[3] =140.; dca_max[3] = 35.; MassRange_Kth[3] =0.36; p_min[3] = 0.675; 
    chi_max[4] =120.; dca_max[4] = 30.; MassRange_Kth[4] =0.38; p_min[4] = 0.700; 
    chi_max[5] =100.; dca_max[5] = 25.; MassRange_Kth[5] =0.40; p_min[5] = 0.725; 
    chi_max[6] = 80.; dca_max[6] = 20.; MassRange_Kth[6] =0.42; p_min[6] = 0.750; 
    chi_max[7] = 60.; dca_max[7] = 15.; MassRange_Kth[7] =0.44; p_min[7] = 0.775; 
    chi_max[8] = 40.; dca_max[8] = 10.; MassRange_Kth[8] =0.46; p_min[8] = 0.800; 
    chi_max[9] = 20.; dca_max[9] =  5.; MassRange_Kth[9] =0.48; p_min[9] = 0.825; 
    chi_max[10]= 10.; dca_max[10]=  3.; MassRange_Kth[10]=0.50; p_min[10]= 0.850; 
    chi_max[11]=  5.; dca_max[11]=  1.; MassRange_Kth[11]=0.52; p_min[11]= 0.875; 
    
    coin_min[0]=-0.10;
    coin_min[1]= 0.00;
    coin_min[2]= 0.30;
    coin_min[3]= 0.60;
    coin_min[4]= 0.90;
    coin_min[5]= 1.20;


  a1=(yy[0]-yy[1])/(x1[0]-x1[1]);
  b1=(yy[0]+yy[1]-(yy[0]-yy[1])/(x1[0]-x1[1])*(x1[0]+x1[1]))/2.;
  a2=(y2[0]-y2[1])/(x2[0]-x2[1]);
  b2=(y2[0]+y2[1]-(y2[0]-y2[1])/(x2[0]-x2[1])*(x2[0]+x2[1]))/2.;

  treeout_flag = false;
  
}
////////////////////////////////////////////////////////////////////////////
RKVertex_ana::~RKVertex_ana(){
  ofp->Close();
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_ana::SetPIDfunc(){
  f_pi = new TF1("f_pi","[0]/sqrt(x*x-1)", 1,7);
  f_k  = new TF1("f_k" ,"[0]/sqrt(x*x-1)", 1,7);
  f_p  = new TF1("f_p" ,"[0]/sqrt(x*x-1)", 1,7);
  f_pi->SetParameter(0, 0.001*Mpi);
  f_k ->SetParameter(0, 0.001*MK);
  f_p ->SetParameter(0, 0.001*Mp);
  set->SetTF1(f_pi ,6,1,1.5);
  set->SetTF1(f_k  ,3,1,1.5);
  set->SetTF1(f_p  ,2,1,1.5);

  f_pi_min = new TF1("f_pi_min","[0]/sqrt(x*x-0.4)", 0.63256,7);
  f_p_min  = new TF1("f_p_min" ,"[0]/sqrt(x*x-1)", 1,7);
  f_pi_min->SetParameter(0, -1.*MassRange_pion_min  );
  f_p_min ->SetParameter(0, MassRange_proton_min);
  set->SetTF1(f_pi_min ,6,2,0.5);
  set->SetTF1(f_p_min  ,2,2,0.5);

  f_pi_max = new TF1("f_pi_max","[0]/sqrt(x*x-1)", 1,7);
  f_p_max  = new TF1("f_p_max" ,"[0]/sqrt(x*x-1)", 1,7);
  f_pi_max->SetParameter(0, -1.*MassRange_pion_max  );
  f_p_max ->SetParameter(0, MassRange_proton_max);
  set->SetTF1(f_pi_max ,6,2,0.5);
  set->SetTF1(f_p_max  ,2,2,0.5);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void RKVertex_ana::makehist(string ofname){
cout<<"makehist"<<endl;
  ofp = new TFile(Form("%s",ofname.c_str()),"recreate");

  h_frame = new TH2F("h_frame","h_frame",10,0,1,10,0,1);
  set->SetTH2(h_frame,"","","");
  h_frame->GetXaxis()->SetNdivisions(000);
  h_frame->GetYaxis()->SetNdivisions(000);
  h_frame->SetStats(0);

  ofp->mkdir("vertex"); ofp->cd("vertex");
  h_inv_mass  = new TH1F("h_inv_mass" ,"h_inv_mass", 250,  1., 1.5);
  h_invm_lpi  = new TH1F("h_invm_lpi" ,"h_invm_lpi", 250,  1., 1.5);
  h_invm_Lam  = new TH1F("h_invm_Lam" ,"h_invm_Lam", 250,  1., 1.5);
  h_invm_bg   = new TH1F("h_invm_bg"  ,"h_invm_bg" , 250,  1., 1.5);
  h_invm_eg1  = new TH1F("h_invm_eg1" ,"h_invm_eg1", 250,  1., 1.5);
  h_invm_eg2  = new TH1F("h_invm_eg2" ,"h_invm_eg2", 250,  1., 1.5);
  h_invm_eg3  = new TH1F("h_invm_eg3" ,"h_invm_eg3", 250,  1., 1.5);
  h_mis_mass  = new TH1F("h_mis_mass" ,"h_mis_mass", 200,  0., 1.0);
  h_mm_pipi   = new TH1F("h_mm_pipi"  ,"h_mm_pipi" , 200,  0., 2.0);
  h_oa_ee     = new TH1F("h_oa_ee"    ,"h_oa_ee"   , 200, -1., 1.0);
  h_oa_ppi    = new TH1F("h_oa_ppi"   ,"h_oa_ppi"  , 200, -1., 1.0);
  h_dca       = new TH1F("h_dca"      ,"h_dca"     ,1000,  0., 100.0);
  h_dca_ppi   = new TH1F("h_dca_ppi"  ,"h_dca_ppi" ,1000,  0., 100.0);
  h_dca_cut   = new TH1F("h_dca_cut"  ,"h_dca_cut" ,1000,  0., 100.0);
  h_mom_Lam   = new TH1F("h_mom_Lam"  ,"h_mom_Lam" ,  50,  0.,   1.5);

  set->SetTH1(h_inv_mass,"Invariant Mass"                       ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.0MeV/#it{c}^{2}"       , 1, 3000, 0);
  set->SetTH1(h_invm_lpi,"Invariant Mass(#pi^{-} mom.<250MeV/c)","Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.0MeV/#it{c}^{2}"       , 4, 3000, 0);
  set->SetTH1(h_invm_Lam,"Invariant Mass(#Lambda)"              ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.0MeV/#it{c}^{2}"       , 4, 3001, 7);
  set->SetTH1(h_invm_bg ,"Invariant Mass(not #Lambda event)"    ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.0MeV/#it{c}^{2}"       , 1, 3001, 9);
  set->SetTH1(h_invm_eg1,"Invariant Mass(E_{#gamma}<1.0)"       ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.0MeV/#it{c}^{2}"       , 1, 3000, 0);
  set->SetTH1(h_invm_eg2,"Invariant Mass(E_{#gamma}=1.0-1.1)"   ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.0MeV/#it{c}^{2}"       , 1, 3000, 0);
  set->SetTH1(h_invm_eg3,"Invariant Mass(E_{#gamma}>1.1)"       ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.0MeV/#it{c}^{2}"       , 1, 3000, 0);
  set->SetTH1(h_mis_mass,"Missing Mass"                         ,"Missing Mass [GeV/#it{c}^{2}]"        ,"Counts"                         , 1, 3000, 0);
  set->SetTH1(h_mm_pipi ,"Missing Mass(#pi^{+} #pi^{-})"        ,"Missing Mass [GeV/#it{c}^{2}]"        ,"Counts"                         , 2, 3000, 0);
  set->SetTH1(h_oa_ee   ,"Opening Angle(e^{+} e{-})"            ,"cos#theta"                            ,"Counts"                         , 1, 3000, 0);
  set->SetTH1(h_oa_ppi  ,"Opening Angle(p #pi^{-})"             ,"cos#theta"                            ,"Counts"                         , 1, 3000, 0);
  set->SetTH1(h_dca     ,"DCA"                                  ,"dca[cm]"                              ,"Counts"                         , 1, 3000, 0);
  set->SetTH1(h_dca_ppi ,"DCA(p#pi)"                            ,"dca[cm]"                              ,"Counts"                         , 1, 3000, 0);
  set->SetTH1(h_dca_cut ,"DCA(p#pi && vertex region cut)"       ,"dca[cm]"                              ,"Counts"                         , 1, 3000, 0);
  set->SetTH1(h_mom_Lam ,"momentum of #Lambda"       ,"momentum[GeV/#it{c}]"                              ,"Counts"                         , 1, 3000, 0);

  h2_dvol     = new TH2F("h2_dvol"    ,"h2_dvol"     , 200, -10., 10., 106, -10.6, 10.6);
  h2_dvol_oa  = new TH2F("h2_dvol_oa" ,"h2_dvol_oa"  , 200, -10., 10., 106, -10.6, 10.6);
  h2_dvol_ppi = new TH2F("h2_dvol_ppi","h2_dvol_ppi" , 200, -10., 10., 106, -10.6, 10.6);
  h2_dvol_cut = new TH2F("h2_dvol_cut","h2_dvol_cut" , 200, -10., 10., 106, -10.6, 10.6);
  h2_dvol_Lam = new TH2F("h2_dvol_Lam","h2_dvol_Lam" , 200, -10., 10., 106, -10.6, 10.6);
  h2_dvol_bg  = new TH2F("h2_dvol_bg" ,"h2_dvol_bg"  , 200, -10., 10., 106, -10.6, 10.6);
  h2_im_vs_mm = new TH2F("h2_im_vs_mm","h2_im_vs_mm" ,  50,   1., 1.5,  35,    0.,  0.7);
  h2_im_pimom = new TH2F("h2_im_pimom","h2_im_pimom" ,  50,   1., 1.5,  35,    0.,  1.0);
  set->SetTH2(h2_dvol    ,"Vertex Position"                          ,"x [cm]","y [cm]");
  set->SetTH2(h2_dvol_oa ,"Vertex Position(p-#pi^{-} && OA cut)"     ,"x [cm]","y [cm]");
  set->SetTH2(h2_dvol_ppi,"Vertex Position(p-#pi^{-})"               ,"x [cm]","y [cm]");
  set->SetTH2(h2_dvol_cut,"Vertex Position(p-#pi^{-}&&region select)","x [cm]","y [cm]");
  set->SetTH2(h2_dvol_Lam,"Vertex Position(#Lambda)"                 ,"x [cm]","y [cm]");
  set->SetTH2(h2_dvol_bg ,"Vertex Position(b.g.)"                    ,"x [cm]","y [cm]");
  set->SetTH2(h2_im_vs_mm,"M_{inv} vs M_{x} p#pi^{-}"                ,"Invariant Mass [GeV/#it{c}^{2}]","Missing Mass [GeV/#it{c}^{2}]");
  set->SetTH2(h2_im_pimom,"M_{inv} vs #pi^{-} mom"                   ,"Invariant Mass [GeV/#it{c}^{2}]","#pi^{-} momentum[GeV/#it{c}]");

  for(int i=0;i<12;i++){
    h_invm_chi[i]    = new TH1F(Form("h_invm_chi%d",i+1)    ,Form("h_invm_chi%d",i+1)   , 250, 1., 1.5);
    h_inv_coin[i]    = new TH1F(Form("h_inv_coin%d",i+1)    ,Form("h_inv_coin%d",i+1)   , 250, 1., 1.5);
    h_invm_dca[i]    = new TH1F(Form("h_invm_dca%d",i+1)    ,Form("h_invm_dca%d",i+1)   , 250, 1., 1.5);
    h_invm_mm[i]     = new TH1F(Form("h_invm_mm%d",i+1)     ,Form("h_invm_mm%d",i+1)    , 250, 1., 1.5);
    h_invm_p_min[i]  = new TH1F(Form("h_invm_p_min%d",i+1)  ,Form("h_invm_p_min%d",i+1) , 250, 1., 1.5);
    set->SetTH1(h_invm_chi[i]   ,Form("Invariant Mass(#chi^{2}<%.0lf)",chi_max[i])               ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.0MeV/#it{c}^{2}"       , 2, 3000, 2);
    set->SetTH1(h_inv_coin[i]   ,Form("Invariant Mass(coin. time>%.02lfns)",coin_min[i])         ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.0MeV/#it{c}^{2}"       , 4, 3000, 2);
    set->SetTH1(h_invm_dca[i]   ,Form("Invariant Mass(#dca<%.0lf)",dca_max[i])                   ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.0MeV/#it{c}^{2}"       , 4, 3000, 0);
    set->SetTH1(h_invm_mm[i]    ,Form("Invariant Mass(MM>%.2lfGeV/#it{c}^{2})",MassRange_Kth[i]) ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.0MeV/#it{c}^{2}"       , 1, 3000, 2);
    set->SetTH1(h_invm_p_min[i] ,Form("Invariant Mass(M_{p}%.2lfGeV/#it{c}^{2})",p_min[i])       ,"Invariant Mass [GeV/#it{c}^{2}]"      ,"Counts/2.0MeV/#it{c}^{2}"       , 1, 3000, 0);
  }
  ofp->cd();

  ofp->mkdir("track"); ofp->cd("track");

  h_chi2   = new TH1F("h_chi2"  ,"h_chi2"  ,1000, 0,200);
  h_msqr   = new TH1F("h_msqr"  ,"h_msqr"  , 400,-2,  2);
  h_msqr_a = new TH1F("h_msqr_a","h_msqr_a", 400,-2,  2);
  h_ft_pi  = new TH1F("h_ft_pi" ,"h_ft_pi" , 100,-1,  1);
  set->SetTH1(h_chi2    ,"#chi^{2}"                ,"#chi^{2}"                  ,"counts / 0.1"            , 1, 3000, 0);
  set->SetTH1(h_msqr    ,"mass square"             ,"mass square[GeV^{2}/c^{4}]","counts / 10MeV^{2}/c^{4}", 4, 3000, 0);
  set->SetTH1(h_msqr_a  ,"mass square(mom<0.15)"   ,"mass square[GeV^{2}/c^{4}]","counts / 10MeV^{2}/c^{4}", 4, 3000, 0);
  set->SetTH1(h_ft_pi   ,"Flight time of #pi^{-}"  ,"flight time[ns]"           ,"counts"                  , 1, 3000, 0);

  h2_pid_all         = new TH2F("h2_pid_all"           ,"h2_pid_all"         , 500, 0,   7, 1000, -1,   1);
  h2_pid_add         = new TH2F("h2_pid_add"           ,"h2_pid_add"         , 500, 0,   7, 1000, -1,   1);
  h2_pid_ppi         = new TH2F("h2_pid_ppi"           ,"h2_pid_ppi"         , 500, 0,   7, 1000, -1,   1);
  h2_pid_decut       = new TH2F("h2_pid_decut"         ,"h2_pid_decut"       , 500, 0,   7, 1000, -1,   1);
  h2_pid_tag         = new TH2F("h2_pid_tag"           ,"h2_pid_tag"         , 500, 0,   7, 1000, -1,   1);
  h2_fl_ftpi         = new TH2F("h2_fl_ftpi"           ,"h2_fl_ftpi"         , 160, 0,  10,  200, -1,   1);
  h2_beta_ftpi       = new TH2F("h2_beta_ftpi"         ,"h2_beta_ftpi"       , 100, 0,   1,  200,  0, 0.5);
  h2_beta_beta_pipi  = new TH2F("h2_beta_beta_pipi"    ,"h2_beta_beta_pipi"  , 500, 0, 1.5,  500,  0, 1.1);
  h2_beta_beta_ppi   = new TH2F("h2_beta_beta_ppi"     ,"h2_beta_beta_ppi"   , 500, 0, 1.5,  500,  0, 1.1);
  h2_beta_beta_all   = new TH2F("h2_beta_beta_all"     ,"h2_beta_beta_all"   , 500, 0, 1.5,  500,  0, 1.1);
  h2_beta_beta_cut   = new TH2F("h2_beta_beta_cut"     ,"h2_beta_beta_cut"   , 500, 0, 1.5,  500,  0, 1.1);
  h2_ft_chi2_pLam    = new TH2F("h2_ft_chi2_pLam "     ,"h2_ft_chi2_pLam "   ,  20, 0,  2.,  100,  0, 300);
  h2_ft_chi2_piLam   = new TH2F("h2_ft_chi2_piLam"     ,"h2_ft_chi2_piLam"   ,  20, 0,  2.,  100,  0, 300);
  set->SetTH2(h2_pid_all         ,"PID plot all"                                       ,"1/#beta"                ,"momentum[GeV/#it{c}]");
  set->SetTH2(h2_pid_add         ,"PID plot cointime>0.5"                              ,"1/#beta"                ,"momentum[GeV/#it{c}]");
  set->SetTH2(h2_pid_ppi         ,"PID plot ppi"                                       ,"1/#beta"                ,"momentum[GeV/#it{c}]");
  set->SetTH2(h2_pid_decut       ,"PID plot all w/ dEdx cut"                           ,"1/#beta"                ,"momentum[GeV/#it{c}]");
  set->SetTH2(h2_pid_tag         ,"PID plot all (OH-TagB)"                             ,"1/#beta"                ,"momentum[GeV/#it{c}]");
  set->SetTH2(h2_fl_ftpi         ,"flight time v.s. flight length"                     ,"flight length[cm]"      ,"flight time[ns]"     );
  set->SetTH2(h2_beta_ftpi       ,"#beta(from mom.)v.s. flight time"                   ,"#beta"                  ,"flight time[ns]"     );
  set->SetTH2(h2_beta_beta_pipi  ,"#beta of #pi^{-} (#pi^{+}#pi^{-})"                  ,"#beta_{#pi} from ToF"   ,"#beta_{#pi} from mom");
  set->SetTH2(h2_beta_beta_ppi   ,"#beta of #pi^{-} (p#pi^{-})"                        ,"#beta_{#pi} from ToF"   ,"#beta_{#pi} from mom");
  set->SetTH2(h2_beta_beta_all   ,"#beta of #pi^{-} all"                               ,"#beta_{#pi} from ToF"   ,"#beta_{#pi} from mom");
  set->SetTH2(h2_beta_beta_cut   ,"#beta of #pi^{-} after cut"                         ,"#beta_{#pi} from ToF"   ,"#beta_{#pi} from mom");
  set->SetTH2(h2_ft_chi2_pLam    ,"flight time vs #chi^{2} (p,#Lambda ev.)"            ,"#Flight time_{#pi}[ns]" ,"#chi^{2}_{p}"        );
  set->SetTH2(h2_ft_chi2_piLam   ,"flight time vs #chi^{2} (#pi^{-},#Lambda ev.)"      ,"#Flight time_{#pi}[ns]" ,"#chi^{2}_{#pi}"      );

  for(int i =0;i<4;i++){
    h_chi2_pi[i] = new TH1F(Form("h_chi2_pi%d",i+1),Form("h_chi2_pi%d",i+1), 1000,0,200);
    h_chi2_p[i]  = new TH1F(Form("h_chi2_p%d",i+1) ,Form("h_chi2_p%d",i+1) , 1000,0,200);
    set->SetTH1(h_chi2_pi[i]    ,Form("#chi^{2}(#pi^{-} region %d)",i+1)      ,"#chi^{2}","counts / 0.1", 1, 3000, 0);
    set->SetTH1(h_chi2_p[i]     ,Form("#chi^{2}(p region %d)",i+1)            ,"#chi^{2}","counts / 0.1", 1, 3000, 0);
  }
  ofp->cd();

  ofp->mkdir("cointime"); ofp->cd("cointime");
  for(int i=0;i<10;i++){//TDL seg
    h_cointime_tdll[i] = new TH1F(Form("h_cointime_tdll%d",i+1),Form("h_cointime_tdll%d",i+1) , 200,-10,  10);
    h_cointime_tdlr[i] = new TH1F(Form("h_cointime_tdlr%d",i+1),Form("h_cointime_tdlr%d",i+1) , 200,-10,  10);
    h_ftpi_tdll[i]     = new TH1F(Form("h_ftpi_tdll%d",i+1)    ,Form("h_ftpi_tdll%d",i+1)     , 100,  0,1.00);
    h_ftpi_tdlr[i]     = new TH1F(Form("h_ftpi_tdlr%d",i+1)    ,Form("h_ftpi_tdlr%d",i+1)     , 100,  0,1.00);
    set->SetTH1(h_cointime_tdll[i],Form("cointime TDLL%d",i+1)  ,"cointime[ns]"         ,"Counts"       , 1, 3001, 2);
    set->SetTH1(h_cointime_tdlr[i],Form("cointime TDLR%d",i+1)  ,"cointime[ns]"         ,"Counts"       , 1, 3001, 2);
    set->SetTH1(h_ftpi_tdll[i]    ,Form("ftpi TDLL%d",i+1)      ,"flight time #pi[ns]"  ,"Counts"       , 1, 3001, 2);
    set->SetTH1(h_ftpi_tdlr[i]    ,Form("ftpi TDLR%d",i+1)      ,"flight time #pi[ns]"  ,"Counts"       , 1, 3001, 2);

    h2_cointime_beta_tdll[i]   = new TH2F(Form("h2_cointime_beta_tdll%d",i+1)   ,Form("h2_cointime_beta_tdll%d",i+1)   , 200, -3, 3, 100, 0.3,  1.4);
    h2_cointime_beta_tdlr[i]   = new TH2F(Form("h2_cointime_beta_tdlr%d",i+1)   ,Form("h2_cointime_beta_tdlr%d",i+1)   , 200, -3, 3, 100, 0.3,  1.4);
    h2_cointime_width_tdllu[i] = new TH2F(Form("h2_cointime_width_tdllu%d",i+1) ,Form("h2_cointime_width_tdllu%d",i+1) , 200, -3, 3, 100,  0.,  54.);
    h2_cointime_width_tdlru[i] = new TH2F(Form("h2_cointime_width_tdlru%d",i+1) ,Form("h2_cointime_width_tdlru%d",i+1) , 200, -3, 3, 100,  0.,  54.);
    h2_cointime_width_tdlld[i] = new TH2F(Form("h2_cointime_width_tdlld%d",i+1) ,Form("h2_cointime_width_tdlld%d",i+1) , 200, -3, 3, 100,  0.,  54.);
    h2_cointime_width_tdlrd[i] = new TH2F(Form("h2_cointime_width_tdlrd%d",i+1) ,Form("h2_cointime_width_tdlrd%d",i+1) , 200, -3, 3, 100,  0.,  54.);
    h2_cointime_ftpi_tdll[i]   = new TH2F(Form("h2_cointime_ftpi_tdll%d",i+1)   ,Form("h2_cointime_ftpi_tdll%d",i+1)   , 200, -2, 2, 100,  0., 0.45);
    h2_cointime_ftpi_tdlr[i]   = new TH2F(Form("h2_cointime_ftpi_tdlr%d",i+1)   ,Form("h2_cointime_ftpi_tdlr%d",i+1)   , 200, -2, 2, 100,  0., 0.45);
    set->SetTH2(h2_cointime_beta_tdll[i]  ,Form("cointime TDLL%d vs #beta",i+1)             ,"cointime[ns]"      ,"#beta"          );
    set->SetTH2(h2_cointime_beta_tdlr[i]  ,Form("cointime TDLR%d vs #beta",i+1)             ,"cointime[ns]"      ,"#beta"          );
    set->SetTH2(h2_cointime_width_tdllu[i],Form("cointime TDLL%d vs WidthU",i+1)            ,"cointime[ns]"      ,"Width[ns]"      );
    set->SetTH2(h2_cointime_width_tdlru[i],Form("cointime TDLR%d vs WidthU",i+1)            ,"cointime[ns]"      ,"Width[ns]"      );
    set->SetTH2(h2_cointime_width_tdlld[i],Form("cointime TDLL%d vs WidthD",i+1)            ,"cointime[ns]"      ,"Width[ns]"      );
    set->SetTH2(h2_cointime_width_tdlrd[i],Form("cointime TDLR%d vs WidthD",i+1)            ,"cointime[ns]"      ,"Width[ns]"      );
    set->SetTH2(h2_cointime_ftpi_tdll[i]  ,Form("cointime TDLL%d vs flight time(#pi)",i+1)  ,"cointime[ns]"      ,"flight time[ns]");
    set->SetTH2(h2_cointime_ftpi_tdlr[i]  ,Form("cointime TDLR%d vs flight time(#pi)",i+1)  ,"cointime[ns]"      ,"flight time[ns]");
  }
  for(int i=0;i<40;i++){//TagB seg
    h_cointime_tb[i]          = new TH1F(Form("h_cointime_tb%d",i+1)         ,Form("h_cointime_tb%d",i+1)          , 200,-10,10);
    set->SetTH1(h_cointime_tb[i] , Form("cointime TagB%d",i+1)  ,"cointime[ns]"      ,"Counts"       , 1, 3001, 4);

    h2_cointime_beta_tagb[i]  = new TH2F(Form("h2_cointime_beta_tagb%d",i+1)  ,Form("h2_cointime_beta_tagb%d",i+1)  , 200, -10, 10, 100, 0.3,  1.4);
    h2_cointime_width_tagb[i] = new TH2F(Form("h2_cointime_width_tagb%d",i+1) ,Form("h2_cointime_width_tagb%d",i+1) , 200, -10, 10, 100,  0.,  50.);
    h2_cointime_ftpi_tagb[i]  = new TH2F(Form("h2_cointime_ftpi_tagb%d",i+1)  ,Form("h2_cointime_ftpi_tagb%d",i+1)  , 200,  -2,  2, 100,  0., 0.45);
    set->SetTH2(h2_cointime_beta_tagb[i]  ,Form("cointime TagB%d vs #beta",i+1)             ,"cointime[ns]"      ,"#beta"          );
    set->SetTH2(h2_cointime_width_tagb[i] ,Form("cointime TagB%d vs width",i+1)             ,"cointime[ns]"      ,"Width[ns]"      );
    set->SetTH2(h2_cointime_ftpi_tagb[i]  ,Form("cointime TagB%d vs flight time(#pi)",i+1)  ,"cointime[ns]"      ,"flight time[ns]");
  }

  h_ctime_allpi        = new TH1F("h_ctime_allpi"     ,"h_ctime_allpi"     ,2400,-30,30);
  h_ctime_allpi_wcut   = new TH1F("h_ctime_allpi_wcut","h_ctime_allpi_wcut",2400,-30,30);
  h_ctime_allpi_rf     = new TH1F("h_ctime_allpi_rf"  ,"h_ctime_allpi_rf"  ,2400,-30,30);
  h_ctime_fastpi_rf    = new TH1F("h_ctime_fastpi_rf" ,"h_ctime_fastpi_rf" ,2400,-30,30);
  h_ctime_ppipi        = new TH1F("h_ctime_ppipi"     ,"h_ctime_ppipi"     , 160, -2, 2);
  h_ctime_ppipi_cor    = new TH1F("h_ctime_ppipi_cor" ,"h_ctime_ppipi_cor" , 160, -2, 2);
  h_ctime_Lam          = new TH1F("h_ctime_Lam"       ,"h_ctime_Lam"       , 160, -2, 2);
  h_ctime_Lam_cor      = new TH1F("h_ctime_Lam_cor"   ,"h_ctime_Lam_cor"   , 160, -2, 2);
  h_ctime_BG           = new TH1F("h_ctime_BG"        ,"h_ctime_BG"        , 160, -2, 2);
  h_ctime_BG_cor       = new TH1F("h_ctime_BG_cor"    ,"h_ctime_BG_cor"    , 160, -2, 2);
  h_time_RF            = new TH1F("h_time_RF"         ,"h_time_RF"         ,10000, 0,3200);
  h_cointime_proton    = new TH1F("h_cointime_proton" ,"h_cointime_proton" , 160, -2,2);
  set->SetTH1(h_ctime_allpi     ,"TDL#otimesTag coin. time(all: #pi^{-})"                        ,"cointime (ns)" ,"counts / 25ps", 1, 3000, 0);
  set->SetTH1(h_ctime_allpi_wcut,"TDL#otimesTag coin. time(all: #pi^{-} TDL width cut)"          ,"cointime (ns)" ,"counts / 25ps", 2, 3000, 0);
  set->SetTH1(h_ctime_allpi_rf  ,"TDL#otimesTag coin. time(all: #pi^{-} RF cut)"                 ,"cointime (ns)" ,"counts / 25ps", 4, 3000, 0);
  set->SetTH1(h_ctime_fastpi_rf ,"TDL#otimesTag coin. time(#beta_{#pi}>0.85 RF cut)"             ,"cointime (ns)" ,"counts / 25ps", 2, 3000, 0);
  set->SetTH1(h_ctime_ppipi     ,"TDL#otimesTag coin. time(p&#pi^{-}:#pi^{-})"                   ,"cointime (ns)" ,"counts / 25ps", 1, 3001, 2);
  set->SetTH1(h_ctime_ppipi_cor ,"TDL#otimesTag coin. time w/ ToF correction(p&#pi^{-}:#pi^{-})" ,"cointime (ns)" ,"counts / 25ps", 1, 3001, 5);
  set->SetTH1(h_ctime_Lam       ,"TDL#otimesTag coin. time(#Lambda)"                             ,"cointime (ns)" ,"counts / 25ps", 1, 3000, 0);
  set->SetTH1(h_ctime_Lam_cor   ,"TDL#otimesTag coin. time(#Lambda w/ ToF cor)"                  ,"cointime (ns)" ,"counts / 25ps", 2, 3000, 0);
  set->SetTH1(h_ctime_BG        ,"TDL#otimesTag coin. time(B.G.)"                                ,"cointime (ns)" ,"counts / 25ps", 1, 3000, 0);
  set->SetTH1(h_ctime_BG_cor    ,"TDL#otimesTag coin. time(B.G. w/ ToF cor)"                     ,"cointime (ns)" ,"counts / 25ps", 1, 3000, 0);
  set->SetTH1(h_cointime_proton ,"TDL#otimesTag coin. time(proton)"                              ,"cointime (ns)" ,"counts / 25ps", 1, 3001, 2);

  h2_cointime_tdlseg   = new TH2F("h2_cointime_tdlseg"    ,"h2_cointime_tdlseg"    ,  300, -10, 10,  21,   -10.5,    10.5);
  h2_cointime_tdluw    = new TH2F("h2_cointime_tdluw"     ,"h2_cointime_tdluw"     , 1000, -10, 10, 500,       0,      50);
  h2_cointime_tdldw    = new TH2F("h2_cointime_tdldw"     ,"h2_cointime_tdldw"     , 1000, -10, 10, 500,       0,      50);
  h2_cointime_tagw     = new TH2F("h2_cointime_tagw"      ,"h2_cointime_tagw"      , 1000, -10, 10, 500,       0,      50);
  h2_cointime_runnum   = new TH2F("h2_cointime_runnum"    ,"h2_cointime_runnum"    , 1000, -10, 10, 125, 10225.5, 10350.5);
  h2_cointime_tdluw_aa = new TH2F("h2_cointime_tdluw_aa"  ,"h2_cointime_tdluw_aa"  , 1000, -10, 10, 500,       0,      50);
  h2_cointime_tdldw_aa = new TH2F("h2_cointime_tdldw_aa"  ,"h2_cointime_tdldw_aa"  , 1000, -10, 10, 500,       0,      50);
  h2_cointime_tagw_aa  = new TH2F("h2_cointime_tagw_aa"   ,"h2_cointime_tagw_aa"   , 1000, -10, 10, 500,       0,      50);
  h2_ctime_proton_mom  = new TH2F("h2_ctime_proton_mom"   ,"h2_ctime_proton_mom"   ,  160,  -2,  2, 200,       0,       1);
  h2_coincoin          = new TH2F("h2_coincoin"           ,"h2_coincoin"           , 1000, -10, 10,1000,     -10,      10);
  h2_cointime_ohde     = new TH2F("h2_cointime_ohde"      ,"h2_cointime_ohde"      , 1000, -10, 10, 100,       0,      10);
  h2_cointime_ohde_a   = new TH2F("h2_cointime_ohde_a"    ,"h2_cointime_ohde_a"    , 1000, -10, 10, 100,       0,      10);
  h2_cointime_oa       = new TH2F("h2_cointime_oa"        ,"h2_cointime_oa"        , 1000, -10, 10, 100,      -1,       1);
  h2_cointime_beta     = new TH2F("h2_cointime_beta"      ,"h2_cointime_beta"      , 1000, -10, 10, 100,     0.3,     1.4);
  h2_cointime_beta_pi  = new TH2F("h2_cointime_beta_pi"   ,"h2_cointime_beta_pi"   , 1000, -10, 10, 100,     0.3,     1.4);
  h2_ctime_ftih        = new TH2F("h2_ctime_ftih"         ,"h2_ctime_ftih"         ,  500,  -2,  2, 200,      -1,    1.50);
  h2_ctime_xpos_Lam    = new TH2F("h2_ctime_xpos_Lam"     ,"h2_ctime_xpos_Lam"     ,  500,  -1,  2, 200,     -10,    15.0);
  h2_ctime_ftih_bw     = new TH2F("h2_ctime_ftih_bw"      ,"h2_ctime_ftih_bw"      ,  500,  -2,  2, 200,      -1,    1.50);
  h2_ctime_flih        = new TH2F("h2_ctime_flih"         ,"h2_ctime_flih"         ,  200,  -2,  2, 100,       7,      14);
  h2_ctime_flpi        = new TH2F("h2_ctime_flpi"         ,"h2_ctime_flpi"         ,  200,  -2,  2, 200,       0,      20);
  h2_flpi_flih         = new TH2F("h2_flpi_flih"          ,"h2_flpi_flih"          ,  200,   0, 20, 200,       0,      20);
  set->SetTH2(h2_cointime_tdlseg    ,"coin. time vs TDL seg"                                 ,"cointime[ns]"     ,"TDL seg"             );
  set->SetTH2(h2_cointime_tdluw     ,"coin. time vs TDLU width"                              ,"cointime[ns]"     ,"TDL pulse width[ns]" );
  set->SetTH2(h2_cointime_tdldw     ,"coin. time vs TDLD width"                              ,"cointime[ns]"     ,"TDL pulse width[ns]" );
  set->SetTH2(h2_cointime_tagw      ,"coin. time vs TagB width"                              ,"cointime[ns]"     ,"TagB pulse width[ns]");
  set->SetTH2(h2_cointime_runnum    ,"coin. time vs RunNum"                                  ,"cointime[ns]"     ,"runnum[ns]"          );
  set->SetTH2(h2_cointime_tdluw_aa  ,"coin. time vs TDLU width (cluster selected)"           ,"cointime[ns]"     ,"TDL pulse width[ns]" );
  set->SetTH2(h2_cointime_tdldw_aa  ,"coin. time vs TDLD width (cluster selected)"           ,"cointime[ns]"     ,"TDL pulse width[ns]" );
  set->SetTH2(h2_cointime_tagw_aa   ,"coin. time vs TagB width (cluster selected)"           ,"cointime[ns]"     ,"TagB pulse width[ns]");
  set->SetTH2(h2_ctime_proton_mom   ,"TDL#otimesTag coin. time(proton) vs mom"               ,"cointime [ns]"    ,"mom[GeV/#it{c}]"     );
  set->SetTH2(h2_coincoin           ,"coin. time(TDL#otimesTag) vs coin. time(OH#otimesTag)" ,"TDL#otimesTag[ns]","OH#otimesTag[ns]"    );
  set->SetTH2(h2_cointime_ohde      ,"coin. time vs OH dE"                                   ,"cointime[ns]"     ,"OH dE[MeV]"          );
  set->SetTH2(h2_cointime_ohde_a    ,"coin. time vs OH dE(cluster selected)"                 ,"cointime[ns]"     ,"OH dE[MeV]"          );
  set->SetTH2(h2_cointime_oa        ,"coin. time vs Opening Angle"                           ,"cointime[ns]"     ,"cos#theta"           );
  set->SetTH2(h2_cointime_beta      ,"coin. time vs #beta"                                   ,"cointime[ns]"     ,"#beta"               );
  set->SetTH2(h2_cointime_beta_pi   ,"coin. time vs #beta (#pi^{-} strict cut)"              ,"cointime[ns]"     ,"#beta"               );
  set->SetTH2(h2_ctime_ftih         ,"Cointime vs flight time to TDL(#beta from mom.)"       ,"cointime [ns]"    ,"flight time[ns]"     );
  set->SetTH2(h2_ctime_xpos_Lam     ,"Cointime vs Vertex x pos.(#Lambda)"                    ,"time at ver[ns]"  ,"x pos[cm]"     );
  set->SetTH2(h2_ctime_ftih_bw      ,"Cointime vs flight time to TDL(BW)"                    ,"cointime [ns]"    ,"flight time[ns]"     );
  set->SetTH2(h2_ctime_flih         ,"Cointime vs flight length to TDL"                      ,"cointime [ns]"    ,"flight length[cm]"   );
  set->SetTH2(h2_ctime_flpi         ,"Cointime vs flight length to TDL(#pi)"                 ,"cointime [ns]"    ,"flight length[cm]"   );
  set->SetTH2(h2_flpi_flih          ,"Distance vs Flight length to TDL(#pi)"                 ,"distance[cm]"     ,"flight length[cm]"   );
  ofp->cd();

  ofp->mkdir("pipi"); ofp->cd("pipi");
  for(int i=0;i<10;i++){//TDL seg
    h_cointime_tdll_pipi[i] = new TH1F(Form("h_cointime_tdll_pipi%d",i+1),Form("h_cointime_tdll_pipi%d",i+1) , 200,-10,  10);
    h_cointime_tdlr_pipi[i] = new TH1F(Form("h_cointime_tdlr_pipi%d",i+1),Form("h_cointime_tdlr_pipi%d",i+1) , 200,-10,  10);
    h_ftpi_tdll_pipi[i]     = new TH1F(Form("h_ftpi_tdll_pipi%d",i+1)    ,Form("h_ftpi_tdll_pipi%d",i+1)     , 100,  0,1.00);
    h_ftpi_tdlr_pipi[i]     = new TH1F(Form("h_ftpi_tdlr_pipi%d",i+1)    ,Form("h_ftpi_tdlr_pipi%d",i+1)     , 100,  0,1.00);
    set->SetTH1(h_cointime_tdll_pipi[i],Form("cointime TDLL%d",i+1)  ,"cointime[ns]"             ,"Counts"       , 1, 3001, 2);
    set->SetTH1(h_cointime_tdlr_pipi[i],Form("cointime TDLR%d",i+1)  ,"cointime[ns]"             ,"Counts"       , 1, 3001, 2);
    set->SetTH1(h_ftpi_tdll_pipi[i]    ,Form("ftpi TDLL%d",i+1)      ,"flight time #pi[ns]"      ,"Counts"       , 1, 3001, 2);
    set->SetTH1(h_ftpi_tdlr_pipi[i]    ,Form("ftpi TDLR%d",i+1)      ,"flight time #pi[ns]"      ,"Counts"       , 1, 3001, 2);

    h2_cointime_beta_tdll_pipi[i]   = new TH2F(Form("h2_cointime_beta_tdll_pipi%d",i+1)  ,Form("h2_cointime_beta_tdll_pipi%d",i+1)   , 200,-3, 3, 100, 0.3, 1.4);
    h2_cointime_beta_tdlr_pipi[i]   = new TH2F(Form("h2_cointime_beta_tdlr_pipi%d",i+1)  ,Form("h2_cointime_beta_tdlr_pipi%d",i+1)   , 200,-3, 3, 100, 0.3, 1.4);
    h2_cointime_width_tdllu_pipi[i] = new TH2F(Form("h2_cointime_width_tdllu_pipi%d",i+1),Form("h2_cointime_width_tdllu_pipi%d",i+1) , 200,-3, 3, 100,  0., 54.);
    h2_cointime_width_tdlru_pipi[i] = new TH2F(Form("h2_cointime_width_tdlru_pipi%d",i+1),Form("h2_cointime_width_tdlru_pipi%d",i+1) , 200,-3, 3, 100,  0., 54.);
    h2_cointime_width_tdlld_pipi[i] = new TH2F(Form("h2_cointime_width_tdlld_pipi%d",i+1),Form("h2_cointime_width_tdlld_pipi%d",i+1) , 200,-3, 3, 100,  0., 54.);
    h2_cointime_width_tdlrd_pipi[i] = new TH2F(Form("h2_cointime_width_tdlrd_pipi%d",i+1),Form("h2_cointime_width_tdlrd_pipi%d",i+1) , 200,-3, 3, 100,  0., 54.);
    h2_cointime_ftpi_tdll_pipi[i]   = new TH2F(Form("h2_cointime_ftpi_tdll_pipi%d",i+1)  ,Form("h2_cointime_ftpi_tdll_pipi%d",i+1)   , 200,-2, 2, 100,  0.,0.45);
    h2_cointime_ftpi_tdlr_pipi[i]   = new TH2F(Form("h2_cointime_ftpi_tdlr_pipi%d",i+1)  ,Form("h2_cointime_ftpi_tdlr_pipi%d",i+1)   , 200,-2, 2, 100,  0.,0.45);
    set->SetTH2(h2_cointime_beta_tdll_pipi[i]  ,Form("cointime TDLL%d vs #beta(#pi#pi)",i+1)     ,"cointime[ns]"      ,"#beta"          );
    set->SetTH2(h2_cointime_beta_tdlr_pipi[i]  ,Form("cointime TDLR%d vs #beta(#pi#pi)",i+1)     ,"cointime[ns]"      ,"#beta"          );
    set->SetTH2(h2_cointime_width_tdllu_pipi[i],Form("cointime TDLL%d vs WidthU(#pi#pi)",i+1)    ,"cointime[ns]"      ,"Width[ns]"      );
    set->SetTH2(h2_cointime_width_tdlru_pipi[i],Form("cointime TDLR%d vs WidthD(#pi#pi)",i+1)    ,"cointime[ns]"      ,"Width[ns]"      );
    set->SetTH2(h2_cointime_width_tdlld_pipi[i],Form("cointime TDLL%d vs WidthU(#pi#pi)",i+1)    ,"cointime[ns]"      ,"Width[ns]"      );
    set->SetTH2(h2_cointime_width_tdlrd_pipi[i],Form("cointime TDLR%d vs WidthD(#pi#pi)",i+1)    ,"cointime[ns]"      ,"Width[ns]"      );
    set->SetTH2(h2_cointime_ftpi_tdll_pipi[i]  ,Form("cointime TDLL%d vs flight time(#pi)",i+1)  ,"cointime[ns]"      ,"flight time[ns]");
    set->SetTH2(h2_cointime_ftpi_tdlr_pipi[i]  ,Form("cointime TDLR%d vs flight time(#pi)",i+1)  ,"cointime[ns]"      ,"flight time[ns]");
  }

  for(int i=0;i<40;i++){//TagB seg
    h_cointime_tb_pipi[i] = new TH1F(Form("h_cointime_tb_pipi%d",i+1),Form("h_cointime_tb_pipi%d",i+1) , 200,-10,10);
    set->SetTH1(h_cointime_tb_pipi[i],Form("cointime TagB%d(#pi#pi)",i+1)  ,"cointime[ns]"      ,"Counts"       , 1, 3001, 4);

    h2_cointime_beta_tagb_pipi[i]  = new TH2F(Form("h2_cointime_beta_tagb_pipi%d",i+1)  ,Form("h2_cointime_beta_tagb_pipi%d",i+1)  , 200,-10,10, 100,0.3, 1.4);
    h2_cointime_width_tagb_pipi[i] = new TH2F(Form("h2_cointime_width_tagb_pipi%d",i+1) ,Form("h2_cointime_width_tagb_pipi%d",i+1) , 200,-10,10, 100, 0., 50.);
    h2_cointime_ftpi_tagb_pipi[i]  = new TH2F(Form("h2_cointime_ftpi_tagb_pipi%d",i+1)  ,Form("h2_cointime_ftpi_tagb_pipi%d",i+1)  , 200, -2, 2, 100, 0.,0.45);
    set->SetTH2(h2_cointime_beta_tagb_pipi[i]  ,Form("cointime TagB%d vs #beta(#pi#pi)",i+1)     ,"cointime[ns]"      ,"#beta"          );
    set->SetTH2(h2_cointime_width_tagb_pipi[i] ,Form("cointime TagB%d vs width(#pi#pi)",i+1)     ,"cointime[ns]"      ,"Width[ns]"      );
    set->SetTH2(h2_cointime_ftpi_tagb_pipi[i]  ,Form("cointime TagB%d vs flight time(#pi)",i+1)  ,"cointime[ns]"      ,"flight time[ns]");

  }
  h_ctime_pipi           = new TH1F("h_ctime_pipi"            ,"h_ctime_pipi"     , 160, -2, 2);
  h_ctime_pipi_cor       = new TH1F("h_ctime_pipi_cor"        ,"h_ctime_pipi_cor" , 160, -2, 2);
  set->SetTH1(h_ctime_pipi           ,"TDL#otimesTag coin. time(#pi^{+}&#pi^{-}:#pi^{-})"                   ,"cointime (ns)","counts / 25ps", 1, 3001, 6);
  set->SetTH1(h_ctime_pipi_cor       ,"TDL#otimesTag coin. time w/ ToF correction(#pi^{+}&#pi^{-}:#pi^{-})" ,"cointime (ns)","counts / 25ps", 1, 3001, 8);

  h2_cointime_beta_pipi  = new TH2F("h2_cointime_beta_pipi"   ,"h2_cointime_beta_pipi"   , 1000, -10, 10, 100, 0.3,1.4);
  set->SetTH2(h2_cointime_beta_pipi  ,"coin. time vs #beta(#pi#pi)"                                    ,"cointime[ns]"     ,"#beta"         );
  ofp->cd();

  ofp->mkdir("time_at_vertex"); ofp->cd("time_at_vertex");
  for(int i =0;i<40;i++){
    h_time_at_ver_tdlr[i] = new TH1F(Form("h_time_at_ver_tdlr%d",i+1),Form("h_time_at_ver_tdlr%d",i+1) , 200,-10,  10);
    h_time_at_ver_tdll[i] = new TH1F(Form("h_time_at_ver_tdll%d",i+1),Form("h_time_at_ver_tdll%d",i+1) , 200,-10,  10);
  }
  h2_time_at_ver_mom      = new TH2F("h2_time_at_ver_mom"      , "h2_time_at_ver_mom"      ,200,-3,3,100,0,1.2);
  h2_time_at_ver_mom_ppi  = new TH2F("h2_time_at_ver_mom_ppi"  , "h2_time_at_ver_mom_ppi"  ,200,-3,3,100,0,1.2);
  h2_time_at_ver_mom_pipi = new TH2F("h2_time_at_ver_mom_pipi" , "h2_time_at_ver_mom_pipi" ,200,-3,3,100,0,1.2);
  set->SetTH2(h2_time_at_ver_mom         ,"time at vertex vs mom.()"                  ,"time at vertex[ns]"           ,"mom[GeV/#it{c}]"  );
  set->SetTH2(h2_time_at_ver_mom_ppi     ,"time at vertex vs mom.(p&#pi^{-})"         ,"time at vertex[ns]"           ,"mom[GeV/#it{c}]"  );
  set->SetTH2(h2_time_at_ver_mom_pipi    ,"time at vertex vs mom.(#pi^{+}&#pi^{-})"   ,"time at vertex[ns]"           ,"mom[GeV/#it{c}]"  );
  ofp->cd();


  ofp->mkdir("tag"); ofp->cd("tag");
  h_tag_nhit         = new TH1F("h_tag_nhit"       ,  "h_tag_nhit"     ,  15,      0.,     15.);
  h_tag_seg          = new TH1F("h_tag_seg"        ,  "h_tag_seg"      ,  45,      0.,     45.);
  h_tdlseg_aa        = new TH1F("h_tdlseg_aa"      ,  "h_tdlseg_aa"    ,  21,   -10.5,    10.5);
  h_tagseg_aa        = new TH1F("h_tagseg_aa"      ,  "h_tagseg_aa"    ,  40,     0.5,    40.5);
  h_runnum_aa        = new TH1F("h_runnum_aa"      ,  "h_runnum_aa"    , 125, 10225.5, 10350.5);
  h_ctime_tb_RF_all  = new TH1F("h_ctime_tb_RF_all","h_ctime_tb_RF_all",1000,      -3,       3);
  set->SetTH1(h_tag_nhit         ,"TagB Multiplicity"               ,"Num. of Hit"  ,"Counts"       , 1, 3000, 0);
  set->SetTH1(h_tag_seg          ,"TagB Segment"                    ,"Segment"      ,"Counts"       , 1, 3000, 0);
  set->SetTH1(h_ctime_tb_RF_all  ,"TagB - RF"                       ,"TagB - RF[ns]","Counts"       , 1, 3001, 9);

  h2_ctime_uwid   = new TH2F("h2_ctime_uwid"   ,"h2_ctime_uwid"   ,200,  -2,    2, 200,     0,   50);
  h2_tagtdlseg_aa = new TH2F("h2_tagtdlseg_aa" ,"h2_tagtdlseg_aa" , 40, 0.5, 40.5,  21, -10.5, 10.5);
  set->SetTH2(h2_ctime_uwid,"Cointime vs Width"               ,"cointime (ns)","TDLU width (ns)");

  for(int i =0;i<40;i++){
    h_ctime_seg[i]    = new TH1F(Form("h_ctime_seg%d",i+1)      , Form("h_ctime_seg%d",i+1)    ,2400,-30,30);
    h_ctime_tb_RF[i]  = new TH1F(Form("h_ctime_tb_RF%d",i+1)    , Form("h_ctime_tb_RF%d",i+1)  ,1000, -3, 3);
    set->SetTH1(h_ctime_seg[i]     ,Form("IH#otimesTagB%d coin. time(all:#pi^{-})",i+1)      ,"cointime (ns)","counts / 25ps", 1, 3000, 0);
    set->SetTH1(h_ctime_tb_RF[i]   ,Form("TagB%d - RF(w/ TWC)"             ,i+1)             ,"TagB - RF[ns]","Counts"       , 1, 3000, 0);

    h2_ctime_tb_RF[i]      =new TH2F(Form("h2_ctime_tb_RF%d",i+1)     ,Form("h2_ctime_tb_RF%d",i+1)     ,1000, -3, 3,1000,0,50);
    set->SetTH2(h2_ctime_tb_RF[i]     ,Form("TagB%d - RF(w/ TWC)"             ,i+1),"TagB - RF[ns]","Width[ns]");
  }
  ofp->cd();

  ofp->mkdir("ih"); ofp->cd("ih");
  h_tdlseg        = new TH1F("h_tdlseg"          ,"h_tdlseg",20,-10.5,10.5);
  h_tdlseg_pipi   = new TH1F("h_tdlseg_pipi"     ,"h_tdlseg_pipi",21,-10.5,10.5);
  h_tdlseg_ppi    = new TH1F("h_tdlseg_ppi"      ,"h_tdlseg_ppi" ,21,-10.5,10.5);
  h_tdlseg_Lam    = new TH1F("h_tdlseg_Lam"      ,"h_tdlseg_Lam" ,21,-10.5,10.5);
  h2_ctime_tdlseg_Lam    = new TH2F("h2_ctime_tdlseg_Lam"      ,"h2_ctime_tdlseg_Lam" ,20,-2,2,21,-10.5,10.5);
  set->SetTH1(h_tdlseg               ,"TDL seg(all #pi^{^})"                   ,"TDL seg","Counts"       , 1, 3000, 0);
  set->SetTH1(h_tdlseg_pipi          ,"TDL seg(#pi^{+}#pi^{-})"                ,"TDL seg","Counts"       , 2, 3000, 0);
  set->SetTH1(h_tdlseg_ppi           ,"TDL seg(p#pi^{-})"                      ,"TDL seg","Counts"       , 4, 3000, 0);
  set->SetTH1(h_tdlseg_Lam           ,"TDL seg(#Lambda)"                       ,"TDL seg","Counts"       , 6, 3000, 0);
  set->SetTH2(h2_ctime_tdlseg_Lam    ,"coin. time(ToF cor.) vs TDL seg(#Lambda)"         ,"cointime[ns]"           ,"TDL seg"  );

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

    h2_pid_ohv         = new TH2F("h2_pid_ohv"          ,"h2_pid_ohv"           , 500, 0, 7, 1000, -1, 1);
    h2_pid_ohh         = new TH2F("h2_pid_ohh"          ,"h2_pid_ohh"           , 500, 0, 7, 1000, -1, 1);
    h2_ohdedxbeta_all  = new TH2F("h2_ohdedxbeta_all"   ,"h2_ohdedxbeta_all"    , 500, 0, 2, 1000,  0, 9);
    h2_ohdedxbeta_decut= new TH2F("h2_ohdedxbeta_decut" ,"h2_ohdedxbeta_decut"  , 500, 0, 2, 1000,  0, 9);
    h2_ohdedxbeta_ohv  = new TH2F("h2_ohdedxbeta_ohv"   ,"h2_ohdedxbeta_ohv"    , 500, 0, 2, 1000,  0, 9);
    h2_ohdedxbeta_ohh  = new TH2F("h2_ohdedxbeta_ohh"   ,"h2_ohdedxbeta_ohh"    , 500, 0, 2, 1000,  0, 9);
    set->SetTH2(h2_pid_ohv         ,"PID plot OHV"                   ,"1/#beta","momentum[GeV/#it{c}]");
    set->SetTH2(h2_pid_ohh         ,"PID plot OHH"                   ,"1/#beta","momentum[GeV/#it{c}]");
    set->SetTH2(h2_ohdedxbeta_all  ,"de/dx(OH) vs #beta all"         ,"#beta"  ,"momentum[GeV/#it{c}]");
    set->SetTH2(h2_ohdedxbeta_decut,"de/dx(OH) vs #beta all(dE cut)" ,"#beta"  ,"momentum[GeV/#it{c}]");
    set->SetTH2(h2_ohdedxbeta_ohv  ,"de/dx(OH) vs #beta OHV1-8"      ,"#beta"  ,"momentum[GeV/#it{c}]");
    set->SetTH2(h2_ohdedxbeta_ohh  ,"de/dx(OH) vs #beta OHH4-6"      ,"#beta"  ,"momentum[GeV/#it{c}]");


  for(int i=0;i<12;i++){
    h_fl_ohvl[i]  = new TH1F(Form("h_fl_ohvl%d",i+1) ,Form("h_fl_ohvl%d",i+1) , 1000, 0, 1200);
    h_fl_ohvr[i]  = new TH1F(Form("h_fl_ohvr%d",i+1) ,Form("h_fl_ohvr%d",i+1) , 1000, 0, 1200);
    set->SetTH1(h_fl_ohvl[i] ,Form("Flight length OHVL%d",i+1) ,"flight length[cm]","counts", 1, 3001, 3);
    set->SetTH1(h_fl_ohvr[i] ,Form("Flight length OHVR%d",i+1) ,"flight length[cm]","counts", 1, 3001, 4);

    h2_pid_ohvl[i] = new TH2F(Form("h2_pid_ohvl%d",i+1) ,Form("h2_pid_ohvl%d",i+1) , 500, 0, 9, 1000, -1, 1);
    h2_pid_ohvr[i] = new TH2F(Form("h2_pid_ohvr%d",i+1) ,Form("h2_pid_ohvr%d",i+1) , 500, 0, 9, 1000, -1, 1);
    set->SetTH2(h2_pid_ohvl[i],Form("PID plot OHVL%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
    set->SetTH2(h2_pid_ohvr[i],Form("PID plot OHVR%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
  }
  for(int i=0;i<9;i++){
    h_fl_ohhl[i]  = new TH1F(Form("h_fl_ohhl%d",i+1) ,Form("h_fl_ohhl%d",i+1) , 1000, 0, 1200);
    h_fl_ohhr[i]  = new TH1F(Form("h_fl_ohhr%d",i+1) ,Form("h_fl_ohhr%d",i+1) , 1000, 0, 1200);
    set->SetTH1(h_fl_ohhl[i] ,Form("Flight length OHHL%d",i+1) ,"flight length[cm]","counts", 1, 3001, 8);
    set->SetTH1(h_fl_ohhr[i] ,Form("Flight length OHHR%d",i+1) ,"flight length[cm]","counts", 1, 3001, 7);

    h2_pid_ohhl[i] = new TH2F(Form("h2_pid_ohhl%d",i+1) ,Form("h2_pid_ohhl%d",i+1) , 500, 0, 9, 1000, -1, 1);
    h2_pid_ohhr[i] = new TH2F(Form("h2_pid_ohhr%d",i+1) ,Form("h2_pid_ohhr%d",i+1) , 500, 0, 9, 1000, -1, 1);;
    set->SetTH2(h2_pid_ohhl[i],Form("PID plot OHHL%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
    set->SetTH2(h2_pid_ohhr[i],Form("PID plot OHHR%d",i+1) ,"1/#beta","momentum[GeV/#it{c}]");
  }

}
////////////////////////////////////////////////////////////////////////////
void RKVertex_ana::SetBranch(){
  cout<<"make branch"<<endl;
  //ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  ofp->cd();
  tree_out = new TTree("tree","tree");
  tree_out->Branch("inv_mass"   , &tr.Inv_mass     , "inv_mass/D"        );
  tree_out->Branch("mom_pi"     , &tr.Mom_pi       , "mom_pi/D"          );
  tree_out->Branch("miss_mass"  , &tr.Miss_mass    , "miss_mass/D"       );
  tree_out->Branch("flength"    , &tr.Flength      , "flength/D"         );
  tree_out->Branch("ftime"      , &tr.Ftime        , "ftime/D"           );
  tree_out->Branch("dtime"      , &tr.Dtime        , "dtime/D"           );
  treeout_flag = true;
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_ana::loop(){
cout<<"start loop"<<endl;
runmin=99999;runmax=-99;

  ofstream fout;

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){
    tree->GetEntry(n);

    if(charge[0]*charge[1]<0){
      if( fout.is_open() ) fout.close();
      fout.open(Form("conf/run%05dplus_minus_event.dat",runnum), ios::out|ios::app);
      fout.setf(ios_base::fixed);
      fout<<"SAVE: "<<evnum<<endl;
      //cout<<"SAVE: "<<evnum<<endl;
    }
  }
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_ana::loop_def(){
cout<<"start loop"<<endl;
runmin=99999;runmax=-99;

int ppicounter=0;
int pcounter=0;
int picounter=0;
bool ppi_flag;//proton & pi^{-} vertex
bool pipi_flag;//pi^{+} & pi^{-} vertex
bool missing_proton_flag;//pi^{+} & pi^{-} vertex
bool pibeta_flag;//pi^{-} beta is consistent with momentum
bool pi_plus_beta_flag;//pi^{+} beta is consistent with momentum
bool oa_flag, tag_flag, dedx_flag, ohseg_flag, vpos_flag, width_flag,p_mass_flag;
bool chi2_flag, dca_flag, mm_flag;
bool de_flag[3];
int pi_m_id, pi_p_id, p_id;//track ID number for pi- and proton
double Egamma, mm, inv_mass;
double cos_oa, cos_pipi;
double fl_pi;//flight length of pion from vertex to TDL
double beta_pi;//pion velocity calcurated from momentum
double ft_pi;//flight time of pion from vertex to TDL
double tmp_tb;

  ofstream fout;
  if( fout.is_open() ) fout.close();
  fout.open("list/plus_minus_event.dat", ios::out|ios::trunc);
  //fout.open("list/ppi_event.dat", ios::out|ios::trunc);
  fout.setf(ios_base::fixed);

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){
    tree->GetEntry(n);
    if(runnum>runmax)runmax=runnum;
    if(runnum<runmin)runmin=runnum;

    //////////////
    //Initialize//
    //////////////
    ppi_flag=pipi_flag=pibeta_flag=pi_plus_beta_flag=false;
    missing_proton_flag = false;
    pi_m_id = pi_p_id = p_id =-1;
    oa_flag = tag_flag = dedx_flag=ohseg_flag=vpos_flag=false;
    chi2_flag = dca_flag = mm_flag = false;
    width_flag=false;
    for(int i=0;i<3;i++){
      de_flag[i]=false;
    }

    for(int i=0;i<40;i++){
      tb_flag[i] = tb_rf_flag[i]=false;
    }

    cos_oa = cos_pipi = mm = -2;
    if(n%100000==0) cout<<n<<" / "<<ENum<<endl;

    h2_dvol    ->Fill(vertex[0],vertex[1]);//x,y
    h_dca      ->Fill(dca);//x,y

    if(vertex[0]>vx_min && vertex[0]<vx_max && vertex[1]>vy_min && vertex[1]<vy_max)vpos_flag=true;
    //h2_dvol_oa ->Fill(vertex[0],vertex[1]);

    //////
    //RF//
    //////
    for(int i=0;i<5;i++){
      RFtime[i] = tagbtime_l[39][i]+RF_offset;
    }
      h_time_RF -> Fill(RFtime[0]);
    for(int i=0;i<24;i++){
      bunch[i]  =0;
      val_tb[i] =1000.;
      if(tagbctime[i]>-100){
        for(int k=0;k<166;k++){
          tmp_tb = RFtime[3]+k*(RFtime[2]-RFtime[3])/83. - tagbctime[i]-tb_offs[i];
          //tmp_tb = RFtime[1]+k*(RFtime[0]-RFtime[1])/83. - tagbctime[i]-tb_offs[i];
          if(runnum<10229)tmp_tb += tb_offs[i];
          if(abs(tmp_tb)<abs(val_tb[i])){val_tb[i]=tmp_tb;bunch[i]=k;}
        }
       }
       h_ctime_tb_RF[i]  ->Fill(val_tb[i]);
       h_ctime_tb_RF_all ->Fill(val_tb[i]);
       h2_ctime_tb_RF[i] ->Fill(val_tb[i], tagbtime_w[i]);
       if(abs(val_tb[i])<0.5)tb_rf_flag[i]=true;
     }
    //h_time_RFdiff -> Fill(RFtime[0]-RFtime[1]);

      ///////////////
      //search pion//
      ///////////////
    for(int i=0;i<ntr;i++){
      if(ass_noh[i]<1||ass_nih[i]<1)continue;//if no hit in hodoscopes
      beta[i] = (fl_oh[i] - fl_ih[i])/(ohct[i] - ihct[i])/(c*100.);
      h2_pid_all->Fill(1./beta[i], mom[i]*charge[i]);
      if(PiMinus(charge[i],mom[i],beta[i])) {
        pi_m_id=i;
      }
      if(PiPulse(charge[i],mom[i],beta[i])) {
        pi_p_id=i;
      }
      vec[i].SetXYZ(momx[i], momy[i], momz[i]);

        if(ihlr[i]==-1){//IHL
          h2_pid_ihl[ihseg[i]-1]  ->Fill(1./beta[i], mom[i]*charge[i]);
        }     
        if(ihlr[i]== 1){//IHR
          h2_pid_ihr[ihseg[i]-1]  ->Fill(1./beta[i], mom[i]*charge[i]);
        }
    }//for each track


    //////////////
    //TagB info.//
    //////////////
    int seg[3];int NHit;
    double tdlmctime,tdltime_uw,tdltime_dw,taggertime;
    double tdl_uw_min,tdl_dw_min, tdl_uw_max,tdl_dw_max;
    if(pi_m_id>-1){
      h_tdlseg ->Fill(ihseg[pi_m_id]*ihlr[pi_m_id]);
      tdl_uw_min = tofParamMan->GetWidthLimit(ihlr[pi_m_id],15,1,0,ihseg[pi_m_id]);
      tdl_uw_max = tofParamMan->GetWidthLimit(ihlr[pi_m_id],15,2,0,ihseg[pi_m_id]);
      tdl_dw_min = tofParamMan->GetWidthLimit(ihlr[pi_m_id],15,1,1,ihseg[pi_m_id]);
      tdl_dw_max = tofParamMan->GetWidthLimit(ihlr[pi_m_id],15,2,1,ihseg[pi_m_id]);

      if(ihlr[pi_m_id]==-1){
        tdlmctime  = tdllmctime[ihseg[pi_m_id]-1];
        tdltime_uw = tdllutime_w[ihseg[pi_m_id]-1];
        tdltime_dw = tdlldtime_w[ihseg[pi_m_id]-1];
      }
      else if(ihlr[pi_m_id]==1){
        tdlmctime  = tdlrmctime[ihseg[pi_m_id]-1];
        tdltime_uw = tdlrutime_w[ihseg[pi_m_id]-1];
        tdltime_dw = tdlrdtime_w[ihseg[pi_m_id]-1];
      }
      if(tdltime_uw<tdl_uw_max&&tdltime_dw<tdl_dw_max&&tdltime_uw>tdl_uw_min&&tdltime_dw>tdl_dw_min)width_flag=true;
      NHit = TagBHit(tagbctime, tagbtime_w, tdlmctime, seg);
      h_tag_nhit ->Fill(NHit);
      taggertime=tagbctime[seg[0]-1];
      Egamma = -1.; 
      for(int i=0;i<24;i++){
        if(tagbctime[i]>-100&&tagbtime_w[i]>11&&vpos_flag){
          h_ctime_allpi ->Fill(tagbctime[i]-tdlmctime);

          if(tb_rf_flag[i]){
            h_ctime_allpi_rf ->Fill(tagbctime[i]-tdlmctime);
            if(width_flag)h_ctime_allpi_wcut ->Fill(tagbctime[i]-tdlmctime);
            if(beta_pi>0.85)h_ctime_fastpi_rf ->Fill(tagbctime[i]-tdlmctime);
            //h2_ctime_ftih ->Fill(tdlmctime-tagbctime[i],fl_ih[pi_m_id]/(beta[pi_m_id]*c*100.) );
            h2_ctime_flih          ->Fill(tagbctime[i]-tdlmctime,fl_ih[pi_m_id]);
            h2_cointime_tdlseg     ->Fill(tagbctime[i]-tdlmctime,ihlr[pi_m_id]*ihseg[pi_m_id]);
            h2_cointime_runnum     ->Fill(tagbctime[i]-tdlmctime,runnum);
            h_runnum_aa     ->Fill(runnum);
            if(ihlr[pi_m_id]==-1){
              h_cointime_tdll[ihseg[pi_m_id]-1]        ->Fill(tagbctime[i]-tdlmctime);
              h2_cointime_width_tdllu[ihseg[pi_m_id]-1]->Fill(tagbctime[i]-tdllmctime[ihseg[pi_m_id]-1],tdllutime_w[ihseg[pi_m_id]-1]);
              h2_cointime_width_tdlld[ihseg[pi_m_id]-1]->Fill(tagbctime[i]-tdllmctime[ihseg[pi_m_id]-1],tdlldtime_w[ihseg[pi_m_id]-1]);
            }
            else if(ihlr[pi_m_id]==1){
              h_cointime_tdlr[ihseg[pi_m_id]-1]        ->Fill(tagbctime[i]-tdlmctime);
              h2_cointime_width_tdlru[ihseg[pi_m_id]-1]->Fill(tagbctime[i]-tdlrmctime[ihseg[pi_m_id]-1],tdlrutime_w[ihseg[pi_m_id]-1]);
              h2_cointime_width_tdlrd[ihseg[pi_m_id]-1]->Fill(tagbctime[i]-tdlrmctime[ihseg[pi_m_id]-1],tdlrdtime_w[ihseg[pi_m_id]-1]);
            }
          }
        }

        h2_cointime_tdluw ->Fill(tagbctime[i]-tdlmctime,tdltime_uw);
        h2_cointime_tdldw ->Fill(tagbctime[i]-tdlmctime,tdltime_dw);
        h2_cointime_tagw  ->Fill(tagbctime[i]-tdlmctime,tagbtime_w[i]);

        for(int j=0;j<ntr;j++){
          if(ass_noh[j]<1)continue;//if no hit in OH
          h2_coincoin      ->Fill(tagbctime[i]-tdlmctime,ohct[j]-tagbctime[i]);

          if(tagbctime[i]>-100&&tagbtime_w[i]>11&&tdltime_uw>10&&tdltime_dw>10&&vpos_flag){
            h2_cointime_ohde ->Fill(tagbctime[i]-tdlmctime,ohde[j]);
          }

          if(tagbctime[i]-tdlmctime>0.5&&tagbctime[i]-tdlmctime<1.5&&ohde[j]>3.&&ohde[j]<6.){
            h2_cointime_ohde_a ->Fill(tagbctime[i]-tdlmctime,ohde[j]);
            h2_pid_add      ->Fill(1./beta[j], mom[j]*charge[j]);
            h_tdlseg_aa     ->Fill(ihlr[j]*ihseg[j]);
            h_tagseg_aa     ->Fill(i+1);
            h2_tagtdlseg_aa ->Fill(i+1,ihlr[j]*ihseg[j]);
            h2_cointime_tdluw_aa ->Fill(tagbctime[i]-tdlmctime,tdltime_uw);
            h2_cointime_tdldw_aa ->Fill(tagbctime[i]-tdlmctime,tdltime_dw);
            h2_cointime_tagw_aa  ->Fill(tagbctime[i]-tdlmctime,tagbtime_w[i]);
          }
        }
      }


      if(NHit == 1){
        Egamma = TPE[seg[0]-1];
        //h_ctime_allpi ->Fill(tagbctime[seg[0]]-ihct[pi_m_id]);
        //h_ctime_seg[seg[0]-1]->Fill(tagbctime[seg[0]]-ihct[pi_m_id]);
        h_cointime_tb[seg[0]-1] ->Fill(taggertime-tdlmctime);
        //if(ihlr[pi_m_id]==-1)    h_cointime_tdll[ihseg[pi_m_id]-1] ->Fill(taggertime-tdlmctime);
        //else if(ihlr[pi_m_id]==1)h_cointime_tdlr[ihseg[pi_m_id]-1] ->Fill(taggertime-tdlmctime);
        //if(ppi_flag)h_ctime_ppipi ->Fill(tdlmctime-taggertime);
        tag_flag = true;
      }
      else if(NHit == 2 && abs(seg[0]-seg[1])==1){
        Egamma = (TPE[seg[0]-1]+TPE[seg[1]-1])/2; 
        //h_ctime_allpi ->Fill(tagbctime[seg[0]]-ihct[pi_m_id]);
        //h_ctime_seg[seg[0]-1]->Fill(tagbctime[seg[0]]-ihct[pi_m_id]);
        //h_ctime_seg[seg[1]-1]->Fill(tagbctime[seg[1]]-ihct[pi_m_id]);
        h_cointime_tb[seg[0]-1] ->Fill(taggertime-tdlmctime);
        h_cointime_tb[seg[1]-1] ->Fill(tagbctime[seg[1]-1]-tdlmctime);
        //if(ihlr[pi_m_id]==-1)    h_cointime_tdll[ihseg[pi_m_id]-1] ->Fill(taggertime-tdlmctime);
        //else if(ihlr[pi_m_id]==1)h_cointime_tdlr[ihseg[pi_m_id]-1] ->Fill(taggertime-tdlmctime);

        //if(ppi_flag)h_ctime_ppipi ->Fill(tagbctime[seg[0]]-ihct[pi_m_id]);
        tag_flag = true;
      }
      if(NHit == 1){h_tag_seg -> Fill(seg[0]);}
      else if((NHit == 2) && (abs(seg[0]-seg[1])==1)){h_tag_seg -> Fill((seg[0]+seg[1])/2.);}
    }


      //////////////
      //each track//
      //////////////
    for(int i=0;i<ntr;i++){
      if(ass_noh[i]<1||!tag_flag)continue;//if no hit in hodoscopes
       h_chi2 ->Fill(chi2[i]);
      if(chi2[i]>MaxChi2)continue;
      beta[i] = (fl_oh[i])/(ohct[i] + taggertime)/(c*100.);
      h2_pid_tag  ->Fill(1./beta[i], mom[i]*charge[i]);
      m_sqr[i] = (1./(beta[i]*beta[i])-1.)*mom[i]*mom[i];
      h_msqr      ->Fill(m_sqr[i]);
      if(fabs(mom[i])<0.15)h_msqr_a    ->Fill(m_sqr[i]);
      h2_ohdedxbeta_all ->Fill(beta[i],ohdedx[i]);
      if(ohde[i]/2.> a1*beta[i]+b1 &&ohde[i]/2.>a2*beta[i]+b2){
        h2_pid_decut        ->Fill(1./beta[i],mom[i]*charge[i]);
        h2_ohdedxbeta_decut ->Fill(beta[i],ohdedx[i]);
        de_flag[i]=true;
      }
      if(Proton(charge[i],mom[i],beta[i]) ) {
        p_id =i; pcounter++; 
        h_chi2_p[0]  ->Fill(chi2[i]);
        if(ass_nih[i]>0){
          if(ihlr[i]==-1){
            h_cointime_proton -> Fill(taggertime - tdllmctime[ihseg[i]-1]);
            h2_ctime_proton_mom -> Fill(taggertime - tdllmctime[ihseg[i]-1], mom[p_id]);
          }
          else if(ihlr[i]==1){
            h_cointime_proton -> Fill(taggertime - tdlrmctime[ihseg[i]-1]);
            h2_ctime_proton_mom -> Fill(taggertime - tdlrmctime[ihseg[i]-1], mom[p_id]);
          }
        }
      }

      if(PiMinus(charge[i],mom[i],beta[i])&&ass_nih[i]>0) {
        pi_m_id=i;
        picounter++;
        h_chi2_pi[0] ->Fill(chi2[i] ); 
        fl_pi=sqrt( (vertex[0]-ih_x[i])*(vertex[0]-ih_x[i]) + (vertex[1]-ih_y[i])*(vertex[1]-ih_y[i])  );
        beta_pi = mom[i]/sqrt(mom[i]*mom[i]+Mpi*Mpi*1e-6);
        //ft_pi = fl_pi/(c*100.);
        //ft_pi = fl_pi/(beta[i]*c*100.);
        ft_pi = fl_pi/(beta_pi*c*100.) + tofParamMan->GetToF(ihlr[pi_m_id],15,ihseg[pi_m_id]);
      //cout<<tofParamMan->GetToF(ihlr[pi_m_id],15,ihseg[pi_m_id])<<endl;

         //if(beta_pi>0.9)cout<<beta[i]<<endl;
        if((beta[i]>(beta_pi-0.2))&&(beta[i]<(beta_pi+0.2))){pibeta_flag = true;
          if(vpos_flag){
            h2_ctime_ftih ->Fill(taggertime-tdlmctime,ft_pi );
            if(ihseg[i]>5)h2_ctime_ftih_bw ->Fill(taggertime-tdlmctime,ft_pi );
          }
        }

        if(vpos_flag){
          h_ft_pi       ->Fill(ft_pi);
          h2_fl_ftpi    ->Fill(fl_pi,ft_pi);
          h2_beta_ftpi  ->Fill(beta_pi,ft_pi);
        }
      }
      if(PiPulse(charge[i],mom[i],beta[i])) {
        pi_p_id=i;
        if((beta[i]>(beta_pi-0.2))&&(beta[i]<(beta_pi+0.2)))pi_plus_beta_flag = true;
      }

        if(ohlr[i]==-1){//OHVL
          h2_pid_ohv              ->Fill(1./beta[i], mom[i]*charge[i]);
          h2_pid_ohvl[ohseg[i]-1] ->Fill(1./beta[i], mom[i]*charge[i]);
        }     
        if(ohlr[i]== 1){//OHVR
          h2_pid_ohv              ->Fill(1./beta[i], mom[i]*charge[i]);
          h2_pid_ohvr[ohseg[i]-1] ->Fill(1./beta[i], mom[i]*charge[i]);
        }
        if(ohlr[i]==-2){//OHHL
          h2_pid_ohh              ->Fill(1./beta[i], mom[i]*charge[i]);
          h2_pid_ohhl[ohseg[i]-1] ->Fill(1./beta[i], mom[i]*charge[i]);
        }
        if(ohlr[i]== 2){//OHHR
          h2_pid_ohh              ->Fill(1./beta[i], mom[i]*charge[i]);
          h2_pid_ohhr[ohseg[i]-1] ->Fill(1./beta[i], mom[i]*charge[i]);
        }
      vec[i].SetXYZ(momx[i], momy[i], momz[i]);
    }//for each track
    cos_oa = vec[0]*vec[1]/(vec[0].Mag()*vec[1].Mag());


    //////////////
    //make flags//
    //////////////
    if(pi_m_id>-1&&pi_p_id>-1){
      cos_pipi = vec[pi_m_id]*vec[pi_p_id]/(vec[pi_m_id].Mag()*vec[pi_p_id].Mag());
      h2_cointime_oa                             ->Fill(taggertime-tdlmctime,cos_pipi);
      if(cos_pipi>-1.0&&cos_pipi<0.8&&vpos_flag&&pibeta_flag&&pi_plus_beta_flag)pipi_flag=true;
      P_pi_minus.SetPxPyPzE(momx[pi_m_id], momy[pi_m_id], momz[pi_m_id], sqrt(0.000001*Mpi*Mpi +mom[pi_m_id]*mom[pi_m_id] ));
      P_pi_plus.SetPxPyPzE( momx[pi_p_id], momy[pi_p_id], momz[pi_p_id], sqrt(0.000001*Mpi*Mpi +mom[pi_p_id]*mom[pi_p_id] ));
      P_inv = P_pi_minus + P_pi_plus;
      inv_mass = P_inv.M();
      mm = MissingMass(Egamma, P_pi_minus, P_pi_plus);
      if(mm>0.88&&mm<1.05)missing_proton_flag=true;
    }

    if(pi_m_id>-1&&p_id>-1&&ass_nih[pi_m_id]>0){
      if(de_flag[pi_m_id]&&de_flag[p_id])dedx_flag=true;
      ppi_flag=true;
      P_pi_minus.SetPxPyPzE(momx[pi_m_id], momy[pi_m_id], momz[pi_m_id], sqrt(0.000001*Mpi*Mpi +mom[pi_m_id]*mom[pi_m_id] ));
      P_p.SetPxPyPzE( momx[p_id] , momy[p_id] , momz[p_id] , sqrt(0.000001*Mp*Mp   +mom[p_id] *mom[p_id]    ));
      P_inv = P_p + P_pi_minus;
      inv_mass = P_inv.M();
      mm = MissingMass(Egamma, P_pi_minus, P_p);
    }

    if(mm<MassRange_K_max&&mm>MassRange_K_min)mm_flag=true;
    if(dca<20.)dca_flag=true;
    if(chi2[0]<140.&&chi2[1]<140.)chi2_flag=true;
    if(cos_oa>oa_min&&cos_oa<oa_max)oa_flag=true;
    if(p_id>0){if(m_sqr[p_id]>0.75*0.75)p_mass_flag=true;}
    
    /////////////////////////////////////
    //analysis of proton & pi^{-} event//
    /////////////////////////////////////
    if(ppi_flag&&dedx_flag){
      ppicounter++;
      h_oa_ppi   ->Fill(cos_oa);
      h2_dvol_ppi->Fill(vertex[0],vertex[1]);

      if(oa_flag){
        h2_dvol_oa->Fill(vertex[0],vertex[1]);
        h2_pid_ppi ->Fill(1./beta[pi_m_id], mom[pi_m_id]*charge[pi_m_id]);
        h2_pid_ppi ->Fill(1./beta[p_id] , mom[p_id]*charge[p_id]  );
        h_dca_ppi  ->Fill(dca);//

        if(vertex[0]>1.&&vertex[0]<6.&&vertex[1]>-2.&&vertex[1]<-1){
          h_dca_cut  ->Fill(dca);//
          h2_dvol_cut->Fill(vertex[0],vertex[1]);
        }
      }
    }//ppi_flag

    if(vpos_flag&&oa_flag&&pi_m_id>-1){
      h2_beta_beta_ppi -> Fill(beta[pi_m_id],beta_pi);
      h2_beta_beta_all -> Fill(beta[pi_m_id],beta_pi);
      if(pibeta_flag)h2_beta_beta_cut -> Fill(beta[pi_m_id],beta_pi);

      if(ihlr[pi_m_id]==-1){
        h_ftpi_tdll[ihseg[pi_m_id]-1]            ->Fill(ft_pi);
        h2_cointime_beta_tdll[ihseg[pi_m_id]-1]  ->Fill(taggertime-tdlmctime-ft_pi,beta[pi_m_id]);
        h2_cointime_ftpi_tdll[ihseg[pi_m_id]-1]  ->Fill(taggertime-tdlmctime,ft_pi);
      }
      if(ihlr[pi_m_id]== 1){
        h_ftpi_tdlr[ihseg[pi_m_id]-1]            ->Fill(ft_pi);
        h2_cointime_beta_tdlr[ihseg[pi_m_id]-1]  ->Fill(taggertime-tdlmctime-ft_pi,beta[pi_m_id]);
        h2_cointime_ftpi_tdlr[ihseg[pi_m_id]-1]  ->Fill(taggertime-tdlmctime,ft_pi);
      }

      h2_cointime_beta ->Fill(taggertime-tdlmctime,beta_pi);
      //h2_cointime_beta ->Fill(taggertime-tdlmctime-ft_pi,beta[pi_m_id]);
      if(pibeta_flag){
        h2_cointime_beta_pi ->Fill(taggertime-tdlmctime,beta[pi_m_id]);
        h_tdlseg_ppi ->Fill(ihseg[pi_m_id]*ihlr[pi_m_id]);
      }
    }

   if(tag_flag&&ppi_flag&&oa_flag&&dedx_flag&&vpos_flag&&width_flag){
     for(int i=0;i<12;i++){
       //all cut : chi2_flag, dca_flag, mm_flag, p_mass_flag
       //if(chi2_flag && dca_flag && mm_flag && p_mass_flag)
       if(mm<MassRange_K_max&&mm>MassRange_Kth[i]&&mom[pi_m_id]<0.25&&dca_flag&&chi2_flag&&p_mass_flag)         h_invm_mm[i]  -> Fill(inv_mass);
       if(chi2[pi_m_id]<chi_max[i]&&chi2[p_id]<chi_max[i]&&mom[pi_m_id]<0.25&&mm_flag&&dca_flag&&p_mass_flag)   h_invm_chi[i] -> Fill(inv_mass);
       if(chi2_flag&&dca<dca_max[i]&&mom[pi_m_id]<0.25&&p_mass_flag && mm_flag)                                 h_invm_dca[i] -> Fill(inv_mass);
       if(chi2_flag && dca_flag && mm_flag && m_sqr[p_id]>p_min[i]*p_min[i]&&mom[pi_m_id]<0.25)                 h_invm_p_min[i] ->Fill(inv_mass);
       if((taggertime-tdlmctime)>coin_min[i]&&mom[pi_m_id]<0.25&&p_mass_flag)h_inv_coin[i] -> Fill(inv_mass);
     }

     if(mm_flag && p_mass_flag && dca_flag && chi2_flag &&mom[pi_m_id]<0.25){
       h_inv_mass  -> Fill(inv_mass);
       h_invm_lpi  -> Fill(inv_mass);

       if(inv_mass>MassRange_Lam_min &&inv_mass<MassRange_Lam_max){
         h_invm_Lam        ->Fill(inv_mass);
         h2_dvol_Lam       ->Fill(vertex[0],vertex[1]);
         h2_ctime_xpos_Lam ->Fill(taggertime-tdlmctime-ft_pi,vertex[0]);
         h_mom_Lam         ->Fill(P_inv.Vect().Mag());
          //if(taggertime-tdlmctime-ft_pi<-0.7)cout<<Form("TDL lr%d seg%d ft=%.2lf beta=%.2lf verx %.2lf very %.2lf",ihlr[pi_m_id],ihseg[pi_m_id],ft_pi,beta_pi, vertex[0] ,vertex[1] )<<endl;
         h_ctime_Lam             ->Fill(taggertime-tdlmctime);
         h_ctime_Lam_cor         ->Fill(taggertime-tdlmctime-ft_pi);
         h2_ft_chi2_pLam         ->Fill(ft_pi,chi2[p_id]);
         h2_ft_chi2_piLam        ->Fill(ft_pi,chi2[pi_m_id]);
         h_tdlseg_Lam            ->Fill(ihseg[pi_m_id]*ihlr[pi_m_id]);
         h2_ctime_tdlseg_Lam     ->Fill(taggertime-tdlmctime-ft_pi, ihseg[pi_m_id]*ihlr[pi_m_id]);
       }

       if( (inv_mass<MassRange_Lam_min || inv_mass>MassRange_Lam_max) ){
         h_invm_bg  ->Fill(inv_mass);
         h2_dvol_bg ->Fill(vertex[0],vertex[1]);
         h_ctime_BG->Fill(taggertime-tdlmctime);
         h_ctime_BG_cor  ->Fill(taggertime-tdlmctime-ft_pi);
       }
       if(Egamma<1.0)            h_invm_eg1  -> Fill(inv_mass);
       if(Egamma<1.1&&Egamma>1.0)h_invm_eg2  -> Fill(inv_mass);
       if(Egamma>1.1)            h_invm_eg3  -> Fill(inv_mass);
       if(inv_mass<1.1){
        h_chi2_pi[1] ->Fill(chi2[pi_m_id]);
        h_chi2_p[1]  ->Fill(chi2[p_id] );
       }
       if(inv_mass>1.1 && inv_mass<1.2){
        h_chi2_pi[2] ->Fill(chi2[pi_m_id]);
        h_chi2_p[2]  ->Fill(chi2[p_id] );
       }
       if(inv_mass>1.2){
        h_chi2_pi[3] ->Fill(chi2[pi_m_id]);
        h_chi2_p[3]  ->Fill(chi2[p_id] );
       }
       if(treeout_flag){
         tr.Inv_mass = inv_mass; tr.Miss_mass = mm; tr.Mom_pi = mom[pi_m_id];
         tr.Flength  = fl_pi   ; tr.Ftime     = ft_pi; tr.Dtime = taggertime-tdlmctime-ft_pi;
         tree_out->Fill();
       }
     }
     h_mis_mass  -> Fill(mm);
     h2_im_vs_mm -> Fill(inv_mass, mm);
     h2_im_pimom -> Fill(inv_mass, mom[pi_m_id]);
     if(p_mass_flag && dca_flag && chi2_flag && width_flag){
       h_ctime_ppipi    ->Fill(taggertime-tdlmctime);
       h_ctime_ppipi_cor->Fill(taggertime-tdlmctime-ft_pi);
     }
     h2_cointime_beta_tagb[seg[0]-1]       ->Fill(taggertime-tdlmctime,beta[pi_m_id]);
     h2_cointime_width_tagb[seg[0]-1]      ->Fill(taggertime-tdlmctime,tagbtime_w[seg[0]-1]);
     h2_cointime_ftpi_tagb[seg[0]-1]       ->Fill(taggertime-tdlmctime,ft_pi);
     h2_ctime_flpi                         ->Fill(taggertime-tdlmctime,fl_pi);
     h2_flpi_flih                          ->Fill(fl_pi,fl_ih[pi_m_id]);

     if(chi2_flag && dca_flag && mm_flag && p_mass_flag){
       h2_time_at_ver_mom_ppi ->Fill(taggertime-tdlmctime-ft_pi,mom[pi_m_id]);
     }
   }//tag_flag&&ppi_flag&&oa_flag&&dedx_flag&&vpos_flag&&width_flag

    //analysis of pi^{+} & pi^{-} event
    if(pipi_flag&&de_flag[pi_m_id]&&de_flag[pi_p_id]&&tag_flag&&width_flag&&chi2_flag && dca_flag){
      h_mm_pipi  -> Fill(mm);

      h2_beta_beta_pipi -> Fill(beta[pi_m_id],beta_pi);
      h2_beta_beta_all  -> Fill(beta[pi_m_id],beta_pi);
      if(pibeta_flag)h2_beta_beta_cut -> Fill(beta[pi_m_id],beta_pi);

      h2_cointime_beta ->Fill(taggertime-tdlmctime,beta_pi);
      if(pibeta_flag){
        h2_cointime_beta_pi ->Fill(taggertime-tdlmctime,beta[pi_m_id]);
        h_tdlseg_pipi ->Fill(ihseg[pi_m_id]*ihlr[pi_m_id]);
      }

      if(ihlr[pi_m_id]==-1){
        h_ftpi_tdll_pipi[ihseg[pi_m_id]-1]            ->Fill(ft_pi);
        h2_cointime_beta_tdll_pipi[ihseg[pi_m_id]-1]  ->Fill(taggertime-tdlmctime,beta_pi);
        h2_cointime_ftpi_tdll_pipi[ihseg[pi_m_id]-1]  ->Fill(taggertime-tdlmctime,ft_pi);
        h2_cointime_width_tdllu_pipi[ihseg[pi_m_id]-1]->Fill(taggertime-tdlmctime,tdllutime_w[ihseg[pi_m_id]-1]);
        h2_cointime_width_tdlld_pipi[ihseg[pi_m_id]-1]->Fill(taggertime-tdlmctime,tdlldtime_w[ihseg[pi_m_id]-1]);
      }
      if(ihlr[pi_m_id]== 1){
        h_ftpi_tdlr_pipi[ihseg[pi_m_id]-1]            ->Fill(ft_pi);
        h2_cointime_beta_tdlr_pipi[ihseg[pi_m_id]-1]  ->Fill(taggertime-tdlmctime,beta_pi);
        h2_cointime_ftpi_tdlr_pipi[ihseg[pi_m_id]-1]  ->Fill(taggertime-tdlmctime,ft_pi);
        h2_cointime_width_tdlru_pipi[ihseg[pi_m_id]-1]->Fill(taggertime-tdlmctime,tdlrutime_w[ihseg[pi_m_id]-1]);
        h2_cointime_width_tdlrd_pipi[ihseg[pi_m_id]-1]->Fill(taggertime-tdlmctime,tdlrdtime_w[ihseg[pi_m_id]-1]);
      }

      if(missing_proton_flag){
        h2_cointime_beta_pipi                      ->Fill(taggertime-tdlmctime,beta_pi);
        h_ctime_pipi                               ->Fill(taggertime-tdlmctime);
        h_ctime_pipi_cor                           ->Fill(taggertime-tdlmctime-ft_pi);
        h2_time_at_ver_mom_pipi                    ->Fill(taggertime-tdlmctime-ft_pi,mom[pi_m_id]);
        h2_cointime_beta_tagb_pipi[seg[0]-1]       ->Fill(taggertime-tdlmctime,beta_pi);
        h2_cointime_width_tagb_pipi[seg[0]-1]      ->Fill(taggertime-tdlmctime,tagbtime_w[seg[0]-1]);
        h2_cointime_ftpi_tagb_pipi[seg[0]-1]       ->Fill(taggertime-tdlmctime,ft_pi);
        //h2_ctime_ftih_pipi                         ->Fill(taggertime-tdlmctime,ft_pi );
        //h2_ctime_flpi_pipi                         ->Fill(taggertime-tdlmctime,fl_pi);
      }
    }

  }//event loop

cout<<"pi- :"<<picounter<<"  p :"<<pcounter<<"  p&pi- :"<<ppicounter<<endl;
  ofp->Write();
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_ana::fit(){
  h_ctime_allpi_zoom = (TH1F*)h_ctime_allpi->Clone();
  h_ctime_allpi_zoom -> GetXaxis()->SetRangeUser(-3,3);

  h_ctime_allpi_zoom_rf = (TH1F*)h_ctime_allpi_rf->Clone();
  h_ctime_allpi_zoom_rf -> GetXaxis()->SetRangeUser(-3,3);

  h_ctime_fastpi_zoom_rf = (TH1F*)h_ctime_fastpi_rf->Clone();
  h_ctime_fastpi_zoom_rf -> GetXaxis()->SetRangeUser(-3,3);

  //TF1 *f1 = new TF1("f1","gaus+pol0[3]",-5,5);
  //set->SetTF1(f1,2,1,1.5);
  //f1->SetParameter(0,100);
  //f1->SetParameter(1,0);
  //f1->SetParameter(2,0.2);
  //f1->SetParameter(3,10);

  //h_ctime_allpi_zoom -> Fit(f1,"QR","",-1.2,1.2);
  for(int i=0;i<5;i++){
    f[i] = new TF1(Form("f%d",i+1),"gaus",-5,5);
    set->SetTF1(f[i],2,1,1.5);
  }

  h_ctime_allpi_zoom_rf -> Fit(f[0],"QR","",-1.2,1.2);
  h_ctime_ppipi      -> Fit(f[1],"QR","",-1.2,1.2);
  h_ctime_ppipi_cor  -> Fit(f[2],"QR","",-1.2,1.2);

  h_ctime_pipi      -> Fit(f[3],"QR","",-1.2,1.2);
  h_ctime_pipi_cor  -> Fit(f[4],"QR","",-1.2,1.2);

  //for(int i=1;i<9;i++){
  //  f_l[i] = new TF1(Form("f_l%d",i+1),"gaus",0,1.0);
  //  set->SetTF1(f_l[i],2,1,1.5);
  //  h_ftpi_tdll[i]->Fit(f_l[i],"QR","",0,1);
  //  param_l[i] = f_l[i]->GetParameter(1);

  //  f_r[i] = new TF1(Form("f_r%d",i+1),"gaus",0,1.0);
  //  set->SetTF1(f_r[i],2,1,1.5);
  //  h_ftpi_tdlr[i]->Fit(f_r[i],"QR","",0,1);
  //  param_r[i] = f_r[i]->GetParameter(1);
  //}

  //for(int i=1;i<9;i++){
  // cout<<"TDLL"<<i+1<<"  "<<param_l[i]<<endl;
  //}
  //for(int i=1;i<9;i++){
  // cout<<"TDLR"<<i+1<<"  "<<param_r[i]<<endl;
  //}

}
////////////////////////////////////////////////////////////////////////////
void RKVertex_ana::draw(){

  for(int i=0;i<30;i++){
    cc[i]->Clear();
  }

  cc[0]->Divide(5,3);
  cc[0]->cd(1);h_frame->Draw("");
                tex->DrawLatex(0.5,0.9,Form("run info.: %d - %d"               ,runmin,runmax));
                //tex->DrawLatex(0.5,0.7,Form("vertex xpos:%.02lf - %.02lf cm"           , vxpos_min,         vxpos_max )          );
                tex->DrawLatex(0.5,0.5,Form("Missing mass:%.02lf - %.02lf GeV/#it{c}^{2}",MassRange_K_min,     MassRange_K_max)      );
                //tex->DrawLatex(0.5,0.3,Form("#Lambda mass:%.02lf - %.02lf GeV/#it{c}^2",MassRange_Lambda[0],MassRange_Lambda[1]) );
                tex->DrawLatex(0.5,0.1,Form("opening angle(cos#theta):%.02lf - %.02lf" ,oa_min , oa_max)               );
  cc[0]->cd(2);gPad->SetLogz(0);h2_dvol    ->Draw("colz");
  cc[0]->cd(3);gPad->SetLogz(0);h2_dvol_ppi->Draw("colz");
  cc[0]->cd(4);gPad->SetLogz(0);h2_dvol_Lam->Draw("colz");
  cc[0]->cd(5);gPad->SetLogz(1);h2_pid_tag ->Draw("colz");
                                f_pi->Draw("same");
                                f_k ->Draw("same");
                                f_p ->Draw("same");
  cc[0]->cd(6);gPad->SetLogz(1);h2_pid_all ->Draw("colz");
                                f_pi_min->Draw("same");
                                f_p_min ->Draw("same");
                                f_pi_max->Draw("same");
                                f_p_max ->Draw("same");
  cc[0]->cd(7);gPad->SetLogz(1);h2_pid_ppi ->Draw("colz");
  //cc[0]->cd(8);gPad->SetLogz(1);h2_pid_add ->Draw("colz");
                                //f_pi->Draw("same");
                                //f_k ->Draw("same");
                                //f_p ->Draw("same");
  cc[0]->cd(8); gPad->SetLogz(1);h2_pid_decut       ->Draw("colz");
                                f_pi_min->Draw("same");
                                f_p_min ->Draw("same");
                                f_pi_max->Draw("same");
                                f_p_max ->Draw("same");
  cc[0]->cd(9); gPad->SetLogz(1);h2_ohdedxbeta_all  ->Draw("colz");
  cc[0]->cd(10);gPad->SetLogz(1);h2_ohdedxbeta_decut->Draw("colz");
                                 line->DrawLine(x1[0],yy[0],x1[1],yy[1]);
                                 line->DrawLine(x2[0],y2[0],x2[1],y2[1]);
  cc[0]->cd(11);gPad->SetLogy(1);h_msqr            ->Draw("");
  cc[0]->cd(12);gPad->SetLogz(1);h2_pid_add        ->Draw("colz");
  cc[0]->cd(13);gPad->SetLogz(1);h2_ctime_xpos_Lam ->Draw("colz");

  cc[1]->Divide(5,3);
  cc[1]->cd(1);gPad->SetLogy(0);h_inv_mass      ->Draw();
                                //h_invm_lpi      ->Draw("same");
                                h_invm_Lam      ->Draw("same");
                                line->DrawLine(MassRange_Lam_min,0,MassRange_Lam_min,20);
                                line->DrawLine(MassRange_Lam_max,0,MassRange_Lam_max,20);
  cc[1]->cd(2); gPad->SetLogy(0);h_mis_mass      ->Draw();
  cc[1]->cd(3); gPad->SetLogz(0);h2_im_vs_mm     ->Draw("colz");line->DrawLine(1.0,MassRange_K_min,1.5,MassRange_K_min);line->DrawLine(1.0,MassRange_K_max,1.5,MassRange_K_max);
  cc[1]->cd(4); gPad->SetLogz(0);h2_im_pimom     ->Draw("colz");//h_invm_eg1      ->Draw();
  cc[1]->cd(5); gPad->SetLogy(0);h_invm_eg2      ->Draw();
  cc[1]->cd(6); gPad->SetLogy(0);h_invm_eg3      ->Draw();
  cc[1]->cd(7); gPad->SetLogy(0);h_ctime_Lam     ->Draw();h_ctime_Lam_cor     ->Draw("same");
  cc[1]->cd(8); gPad->SetLogz(0);h2_ft_chi2_pLam ->Draw("colz");
  cc[1]->cd(9); gPad->SetLogz(0);h2_ft_chi2_piLam->Draw("colz");
  //cc[1]->cd(10);gPad->SetLogy(0);h_inv_coin[0]   ->Draw();
  cc[1]->cd(10);gPad->SetLogy(0);h_mm_pipi       ->Draw();
  cc[1]->cd(11);gPad->SetLogy(0);h_inv_coin[1]   ->Draw();
  cc[1]->cd(12);gPad->SetLogy(0);h_inv_coin[2]   ->Draw();//h_tdlseg_aa     ->Draw("");
  cc[1]->cd(13);gPad->SetLogy(0);h_inv_coin[3]   ->Draw();//h_tagseg_aa     ->Draw("");
  cc[1]->cd(14);gPad->SetLogy(0);h_inv_coin[4]   ->Draw();//h2_tagtdlseg_aa ->Draw("colz");
  cc[1]->cd(15);gPad->SetLogy(0);h_inv_coin[5]   ->Draw();//h2_tagtdlseg_aa ->Draw("colz");


  cc[2]->Divide(5,3);
  double maxy=0.9*(h_ctime_allpi ->GetBinContent(h_ctime_allpi ->GetMaximumBin()) );
  cc[2]->cd(1); gPad->SetLogy(0);h_ctime_allpi   ->Draw("");
                                 line->DrawLine(MaxCoinTime,0,MaxCoinTime,maxy);
                                 line->DrawLine(MinCoinTime,0,MinCoinTime,maxy);
  cc[2]->cd(2); gPad->SetLogy(1);h_ctime_allpi_zoom      ->Draw("");
                                 h_ctime_allpi_zoom_rf   ->Draw("same");
                                 h_ctime_allpi_wcut      ->Draw("same");
                                 //h_ctime_fastpi_zoom_rf  ->Draw("same");
                                 line->DrawLine(MaxCoinTime,0,MaxCoinTime,maxy);
                                 line->DrawLine(MinCoinTime,0,MinCoinTime,maxy);
  cc[2]->cd(3); gPad->SetLogy(1);h_ctime_ppipi    ->Draw("");
  cc[2]->cd(4); gPad->SetLogy(1);h_ctime_ppipi_cor->Draw("");
  //cc[2]->cd(2); gPad->SetLogy(0);h_tag_nhit      ->Draw("");
  //cc[2]->cd(3); gPad->SetLogy(0);h_tag_seg       ->Draw("");
  cc[2]->cd(5); gPad->SetLogz(1);h2_cointime_tagw ->Draw("colz");
  cc[2]->cd(6); gPad->SetLogy(1);h_oa_ppi         ->Draw("");
  cc[2]->cd(7); gPad->SetLogz(1);h2_cointime_tdluw->Draw("colz");
  cc[2]->cd(8); gPad->SetLogz(1);h2_cointime_tdldw->Draw("colz");
  //cc[2]->cd(9); gPad->SetLogz(1);h2_coincoin      ->Draw("colz");
  cc[2]->cd(9); gPad->SetLogz(1);h2_cointime_ohde   ->Draw("colz");
  //cc[2]->cd(11);gPad->SetLogz(1);h2_ctime_flpi    ->Draw("colz");
  cc[2]->cd(10);gPad->SetLogz(1);h2_ctime_ftih      ->Draw("colz");
  //cc[2]->cd(11);gPad->SetLogz(1);h2_ctime_flih    ->Draw("colz");
  cc[2]->cd(11);gPad->SetLogz(1);h2_cointime_beta    ->Draw("colz");
  cc[2]->cd(12);gPad->SetLogz(1);h2_cointime_beta_pi ->Draw("colz");
  //cc[2]->cd(5); gPad->SetLogz(1);h2_cointime_tagw_aa ->Draw("colz");
  //cc[2]->cd(7); gPad->SetLogz(1);h2_cointime_tdluw_aa->Draw("colz");
  //cc[2]->cd(8); gPad->SetLogz(1);h2_cointime_tdldw_aa->Draw("colz");

  //cc[2]->cd(13);gPad->SetLogy(0);h_tdlseg_aa     ->Draw("");
  //cc[2]->cd(14);gPad->SetLogy(0);h_tagseg_aa     ->Draw("");
  cc[2]->cd(13); gPad->SetLogy(1);h_cointime_proton    ->Draw("");
  cc[2]->cd(14); gPad->SetLogz(1);h2_ctime_proton_mom  ->Draw("colz");
  cc[2]->cd(15);gPad->SetLogz(1);h2_cointime_runnum    ->Draw("colz");//h_runnum_aa     ->Draw("");

  cc[3]->Divide(5,4);
  cc[3]->cd(1) ;gPad->SetLogy(0);h_invm_chi[0]      ->Draw("");
  cc[3]->cd(2) ;gPad->SetLogy(0);h_invm_chi[1]      ->Draw("");
  cc[3]->cd(3) ;gPad->SetLogy(0);h_invm_chi[2]      ->Draw("");
  cc[3]->cd(4) ;gPad->SetLogy(0);h_invm_chi[3]      ->Draw("");
  cc[3]->cd(5) ;gPad->SetLogy(0);h_invm_chi[4]      ->Draw("");
  cc[3]->cd(6) ;gPad->SetLogy(0);h_invm_chi[6]      ->Draw("");
  cc[3]->cd(7) ;gPad->SetLogy(0);h_invm_chi[7]      ->Draw("");
  cc[3]->cd(8) ;gPad->SetLogy(0);h_invm_chi[8]      ->Draw("");
  cc[3]->cd(9) ;gPad->SetLogy(0);h_invm_chi[9]      ->Draw("");
  cc[3]->cd(10);gPad->SetLogy(0);h_invm_chi[10]     ->Draw("");
  cc[3]->cd(11);gPad->SetLogy(0);h_invm_mm[0]       ->Draw("");
  cc[3]->cd(12);gPad->SetLogy(0);h_invm_mm[1]       ->Draw("");
  cc[3]->cd(13);gPad->SetLogy(0);h_invm_mm[2]       ->Draw("");
  cc[3]->cd(14);gPad->SetLogy(0);h_invm_mm[3]       ->Draw("");
  cc[3]->cd(15);gPad->SetLogy(0);h_invm_mm[4]       ->Draw("");
  cc[3]->cd(16);gPad->SetLogy(1);h_dca              ->Draw("");
  cc[3]->cd(17);gPad->SetLogy(1);h_chi2             ->Draw("");
  cc[3]->cd(18);gPad->SetLogz(1);h2_time_at_ver_mom_pipi ->Draw("colz");
  cc[3]->cd(19);gPad->SetLogz(1);h2_time_at_ver_mom_ppi  ->Draw("colz");
  cc[3]->cd(20);gPad->SetLogy(0);h_mom_Lam          ->Draw();

  //cc[4]->Divide(5,2);
  //cc[4]->cd(2) ;gPad->SetLogy(1);h_chi2_pi[0]->Draw("");
  //cc[4]->cd(3) ;gPad->SetLogy(1);h_chi2_pi[1]->Draw("");
  //cc[4]->cd(4) ;gPad->SetLogy(1);h_chi2_pi[2]->Draw("");
  //cc[4]->cd(5) ;gPad->SetLogy(1);h_chi2_pi[3]->Draw("");
  //cc[4]->cd(7) ;gPad->SetLogy(1);h_chi2_p[0] ->Draw("");
  //cc[4]->cd(8) ;gPad->SetLogy(1);h_chi2_p[1] ->Draw("");
  //cc[4]->cd(9) ;gPad->SetLogy(1);h_chi2_p[2] ->Draw("");
  //cc[4]->cd(10);gPad->SetLogy(1);h_chi2_p[3] ->Draw("");

  cc[4]->Divide(5,4);
  cc[4]->cd(1) ;gPad->SetLogy(1);h_ctime_pipi         ->Draw("");
  cc[4]->cd(2) ;gPad->SetLogy(1);h_ctime_pipi_cor     ->Draw("");
  cc[4]->cd(3) ;gPad->SetLogz(1);h2_cointime_beta_pipi->Draw("colz");
  cc[4]->cd(4) ;gPad->SetLogz(1);h2_cointime_oa       ->Draw("colz");
  cc[4]->cd(5) ;gPad->SetLogz(1);h2_flpi_flih         ->Draw("colz");
  cc[4]->cd(6) ;gPad->SetLogz(1);h2_ctime_ftih        ->Draw("colz");
  cc[4]->cd(7) ;gPad->SetLogz(1);h2_ctime_flih        ->Draw("colz");
  cc[4]->cd(8) ;gPad->SetLogz(1);h2_fl_ftpi           ->Draw("colz");
  cc[4]->cd(9) ;gPad->SetLogz(1);h2_beta_ftpi         ->Draw("colz");
  cc[4]->cd(10);gPad->SetLogz(1);h2_cointime_tdlseg   ->Draw("colz");
  cc[4]->cd(11);gPad->SetLogy(1);h_ctime_tb_RF_all    ->Draw("");
  cc[4]->cd(12);gPad->SetLogz(1);h2_ctime_ftih_bw     ->Draw("colz");
  cc[4]->cd(13);gPad->SetLogz(1);h2_beta_beta_pipi    ->Draw("colz");
  cc[4]->cd(14);gPad->SetLogz(1);h2_beta_beta_ppi     ->Draw("colz");
  cc[4]->cd(15);gPad->SetLogz(1);h2_beta_beta_cut     ->Draw("colz");
  cc[4]->cd(16);gPad->SetLogy(0);h_tdlseg             ->Draw("");
  cc[4]->cd(17);gPad->SetLogy(0);h_tdlseg_pipi        ->Draw("");
  cc[4]->cd(18);gPad->SetLogy(0);h_tdlseg_ppi         ->Draw("");
  cc[4]->cd(19);gPad->SetLogy(0);h_tdlseg_Lam         ->Draw("");


//PID                                  
  cc[5]->Divide(4,3);
  for(int i=0;i<12;i++){
  cc[5]->cd(i+1);gPad->SetLogz(1);h2_pid_ohvl[i]->Draw("colz");
  }
  cc[6]->Divide(4,3);
  for(int i=0;i<12;i++){
  cc[6]->cd(i+1);gPad->SetLogz(1);h2_pid_ohvr[i]->Draw("colz");
  }
  
  cc[7]->Divide(4,3);
  for(int i=0;i<8;i++){
  cc[7]->cd(i+1);gPad->SetLogz(1);h2_pid_ohhl[i]->Draw("colz");
  }
  
  cc[8]->Divide(4,3);
  for(int i=0;i<9;i++){
  cc[8]->cd(i+1);gPad->SetLogz(1);h2_pid_ohhr[i]->Draw("colz");
  }

  cc[9]->Divide(4,3);
  for(int i=0;i<10;i++){
  cc[9]->cd(i+1);gPad->SetLogz(1);h2_pid_ihl[i]->Draw("colz");
  }

  cc[10]->Divide(4,3);
  for(int i=0;i<10;i++){
  cc[10]->cd(i+1);gPad->SetLogz(1);h2_pid_ihr[i]->Draw("colz");
  }

//coin time
  cc[11]->Divide(6,4);
  for(int i=0;i<24;i++){
  cc[11]->cd(i+1);gPad->SetLogy(1);h_cointime_tb[i]->Draw("");
  }

  cc[12]->Divide(6,4);
  for(int i=0;i<24;i++){
  cc[12]->cd(i+1);gPad->SetLogz(1);h2_cointime_width_tagb[i]->Draw("colz");//h2_cointime_beta_tagb[i]->Draw("colz");
  }

  cc[13]->Divide(6,4);
  for(int i=0;i<24;i++){
  cc[13]->cd(i+1);gPad->SetLogz(1);h2_cointime_ftpi_tagb[i]->Draw("colz");
  }

  cc[14]->Divide(5,4);
  for(int i=0;i<10;i++){
  cc[14]->cd(i+1); gPad->SetLogy(1);h_cointime_tdll[i]->Draw("");//h_ftpi_tdll[i]->Draw("");
  cc[14]->cd(i+11);gPad->SetLogy(1);h_cointime_tdlr[i]->Draw("");//h_ftpi_tdlr[i]->Draw("");
  }

  cc[15]->Divide(5,4);
  for(int i=1;i<9;i++){
  cc[15]->cd(i+1); gPad->SetLogz(1);h2_cointime_beta_tdll[i]->Draw("colz");
  cc[15]->cd(i+11);gPad->SetLogz(1);h2_cointime_beta_tdlr[i]->Draw("colz");
  }

  cc[16]->Divide(5,4);
  for(int i=1;i<9;i++){
  cc[16]->cd(i+1); gPad->SetLogz(1);h2_cointime_width_tdllu[i]->Draw("colz");//h2_cointime_ftpi_tdll_pipi[i]->Draw("colz");
  cc[16]->cd(i+11);gPad->SetLogz(1);h2_cointime_width_tdlld[i]->Draw("colz");//h2_cointime_ftpi_tdlr_pipi[i]->Draw("colz");
  }

  cc[17]->Divide(5,4);
  for(int i=1;i<9;i++){
  cc[17]->cd(i+1); gPad->SetLogz(1);h2_cointime_width_tdlru[i]->Draw("colz");
  cc[17]->cd(i+11);gPad->SetLogz(1);h2_cointime_width_tdlrd[i]->Draw("colz");
  }

  cc[18]->Divide(5,4);
  for(int i=1;i<9;i++){
  cc[18]->cd(i+1); gPad->SetLogz(1);h2_cointime_ftpi_tdll[i]->Draw("colz");
  cc[18]->cd(i+11);gPad->SetLogz(1);h2_cointime_ftpi_tdlr[i]->Draw("colz");
  }

  cc[19]->Divide(6,4);
  for(int i=0;i<24;i++){
    cc[19]->cd(i+1);gPad->SetLogy(1);h_ctime_tb_RF[i]->Draw("");
  }
  cc[20]->Divide(6,4);
  for(int i=0;i<24;i++){
    cc[20]->cd(i+1);gPad->SetLogz(1);h2_ctime_tb_RF[i]->Draw("colz");
  }

  //cc[19]->Divide(5,4);
  //for(int i=0;i<10;i++){
  //cc[19]->cd(i+1); gPad->SetLogy(1);h_ftpi_tdll[i]->Draw("");
  //cc[19]->cd(i+11);gPad->SetLogy(1);h_ftpi_tdlr[i]->Draw("");
  //}

}
////////////////////////////////////////////////////////////////////////////
void RKVertex_ana::SetPdfFilename(string ifname){
pdf_name = ifname;
} 
////////////////////////////////////////////////////////////////////////////
void RKVertex_ana::savecanvas(){
  cc[0] ->Print(Form("%s[",pdf_name.c_str()) );
  cc[0] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[1] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[2] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[3] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[4] ->Print(Form("%s" ,pdf_name.c_str()) );
  cc[4] ->Print(Form("%s]",pdf_name.c_str()) );

 pdf_name.erase(pdf_name.size()-4);
  cc[5] ->Print(Form("%s_pid.pdf[",pdf_name.c_str()) );
  cc[5] ->Print(Form("%s_pid.pdf" ,pdf_name.c_str()) );
  cc[6] ->Print(Form("%s_pid.pdf" ,pdf_name.c_str()) );
  cc[7] ->Print(Form("%s_pid.pdf" ,pdf_name.c_str()) );
  cc[8] ->Print(Form("%s_pid.pdf" ,pdf_name.c_str()) );
  cc[9] ->Print(Form("%s_pid.pdf" ,pdf_name.c_str()) );
  cc[10]->Print(Form("%s_pid.pdf" ,pdf_name.c_str()) );
  cc[10]->Print(Form("%s_pid.pdf]",pdf_name.c_str()) );

  cc[11]->Print(Form("%s_coin.pdf[",pdf_name.c_str()) );
  cc[11]->Print(Form("%s_coin.pdf" ,pdf_name.c_str()) );
  cc[12]->Print(Form("%s_coin.pdf" ,pdf_name.c_str()) );
  cc[13]->Print(Form("%s_coin.pdf" ,pdf_name.c_str()) );
  cc[14]->Print(Form("%s_coin.pdf" ,pdf_name.c_str()) );
  cc[15]->Print(Form("%s_coin.pdf" ,pdf_name.c_str()) );
  cc[16]->Print(Form("%s_coin.pdf" ,pdf_name.c_str()) );
  cc[17]->Print(Form("%s_coin.pdf" ,pdf_name.c_str()) );
  cc[18]->Print(Form("%s_coin.pdf" ,pdf_name.c_str()) );
  cc[19]->Print(Form("%s_coin.pdf" ,pdf_name.c_str()) );
  cc[20]->Print(Form("%s_coin.pdf" ,pdf_name.c_str()) );
  cc[20]->Print(Form("%s_coin.pdf]",pdf_name.c_str()) );
  cout<<pdf_name<<".pdf saved!"<<endl;
}
////////////////////////////////////////////////
void RKVertex_ana::SetRoot(string ifname)
{
  cout<<"SetRoot"<<endl;
  std::ifstream fp(ifname.c_str());
  if(fp.fail()){ cout<<"file open fail!"<<endl; exit(1); }
  vector<string> runname;
  string line, aa;
  while(1){
    getline(fp,line);
    if(line[0]=='#') continue;
    if( fp.eof() ) break;
    istringstream sline(line);
    sline >> aa >> ws;
    runname.push_back(line);
  }

  cout<<"Run list"<<endl;
  for(int i=0;i<(int)runname.size();i++){
    cout<<runname[i]<<endl;
    add(Form("%s",runname[i].c_str()));
  }

   readtree();
  ENum = GetEntries();
}
////////////////////////////////////////////////
void RKVertex_ana::SetRootSingle(string ifname)
{
  cout<<"SetRootSingle"<<endl;
  add(Form("%s",ifname.c_str()));

  readtree();
  ENum = GetEntries();
}
////////////////////////////////////////////////
void RKVertex_ana::GetHist(string ifname)
{
  cout<<"Get Histogram from root file"<<endl;
  TFile *file = new TFile(ifname.c_str() );
}

////////////////////////////////////////////////
int RKVertex_ana::TagBHit(double *TagBTime, double *TagB_w, double IHTime, int *seg){

  seg[0] = seg[1] =seg[2] = -1;
  int NHit = 0;

  for(int i=0; i<27; i++){
    TagBTime[i]+=TagOffs;
    if( (IHTime-TagBTime[i])>MinCoinTime && (IHTime-TagBTime[i])<MaxCoinTime && TagB_w[i]>MinTagB_w&&tb_rf_flag[i]){
      if(NHit<2)seg[NHit] = i+1;
      NHit++;
    }
  }

  return NHit;
}
////////////////////////////////////////////////

bool RKVertex_ana::SetEgamma(void){

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

double RKVertex_ana::MissingMass(double Eg, TLorentzVector P1, TLorentzVector P2){
  TLorentzVector P_gamma(Eg, 0., 0., Eg);
  TLorentzVector P_target(0., 0., 0., 0.001*Mp);

  TLorentzVector P_X = P_gamma + P_target - P1 - P2;
  return P_X.M();

}

////////////////////////////////////////////////

bool RKVertex_ana::Proton(int charge, double mom, double beta){

  if(charge < 0){return false;}
  if(((mom*sqrt(1-pow(beta, 2.))/beta) > MassRange_proton_min) &&
     ((mom*sqrt(1-pow(beta, 2.))/beta) < MassRange_proton_max)){
    return true;
  }
  else return false;

}

////////////////////////////////////////////////

bool RKVertex_ana::PiMinus(int charge, double mom, double beta){

  //mom>0
  if(charge > 0){return false;}
  if(1./(beta*beta)<0.4){return false;}
  else{
    if(((mom*sqrt(1/pow(beta, 2.)-0.4)) > MassRange_pion_min) &&
       ((mom*sqrt(1/pow(beta, 2.)-1.)) < MassRange_pion_max)){
      //  if((mom*sqrt(1-pow(beta, 2.))/beta) > MassRange_pion[1]){
      return true;
    }
    else if(beta>1.&&(mom*sqrt(1/pow(beta, 2.)-0.4)) > MassRange_pion_min){
      return true;
    }
  }
  return false;

}

////////////////////////////////////////////////

bool RKVertex_ana::PiPulse(int charge, double mom, double beta){

  //mom>0
  if(charge < 0){return false;}
  if(1./(beta*beta)<0.4){return false;}
  else{
    if(((mom*sqrt(1/pow(beta, 2.)-0.4)) > MassRange_pion_min) &&
       ((mom*sqrt(1/pow(beta, 2.)-1.)) < MassRange_pion_max)){
      //  if((mom*sqrt(1-pow(beta, 2.))/beta) > MassRange_pion[1]){
      return true;
    }
    else if(beta>1.&&(mom*sqrt(1/pow(beta, 2.)-0.4)) > MassRange_pion_min){
      return true;
    }
  }
  return false;

}

////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "output0000.dat";
  string ifname_root = "rkvertex.root";
  string ofname = "root/hoge.root";
  string pdfname = "pdf/track/hoge.pdf";
  int ch;
  int MaxNum = 0;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  bool single_rootfile_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:s:w:n:bcop:"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename list: "<<ifname<<endl;
      break;
    case 's':
      ifname_root = optarg;
      single_rootfile_flag = true;
      cout<<"input root filename: "<<ifname_root<<endl;
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
      cout<<"-f : input runlist"<<endl;
      cout<<"-s : input root filename"<<endl;
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
  RKVertex_ana *ana = new RKVertex_ana();

  ana->SetPIDfunc();
  ana->SetEgamma();
  ana->SetMaxEvent(MaxNum);
  if(!single_rootfile_flag)ana->SetRoot(ifname);
  if( single_rootfile_flag)ana->SetRootSingle(ifname_root);
  ana->makehist(ofname);
  #ifdef OutTree
    ana->SetBranch();
  #endif
  ana->loop_def();
  ana->fit();
  ana->draw();
  ana->SetPdfFilename(pdfname);
  ana->savecanvas();
  delete ana;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}

