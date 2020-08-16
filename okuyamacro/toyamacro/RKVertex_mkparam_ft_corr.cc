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

static const double PI = 4.0*atan(1.);
static const double mrad_to_deg = 1./1000*180./PI;
const double Mp = 938.272046;          // proton       mass (MeV/c2)
const double Mpi = 139.57018;          // charged pion mass (MeV/c2)
const double MK = 493.677;             // charged Kaon mass (MeV/c2)
const double c = 0.299792458;          // speed of light in vacuum (m/ns)

double vx_min = -3.2;
double vx_max =  4.5;
double vy_min = -4.0;
double vy_max =  4.0;

const double MaxCoinTime =  2.0;//IH-Tagger Coincidence Timing[ns]
const double MinCoinTime = -2.0;//IH-Tagger Coincidence Timing[ns]
const int NCanvas=8; //num of TCanvas
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class RKVertex_mkparam_ft_corr : public RKVertex_tr
{
 public:
         RKVertex_mkparam_ft_corr();
        ~RKVertex_mkparam_ft_corr();
  void makehist(string ofname);
  void SetPIDfunc();
  void loop();
  void fit_TDL();
  void draw(); 
  void savecanvas(); 
  void SetRoot(string ifname); 
  void SetRootSingle(string ifname); 
  void GetHist(string ifname); 
  void SetMaxEvent( int N )  { ENumMax = N; }
  void SetPdfFilename(string ifname); 
  void SetInputParam(string ifname); 
  void SetOutputParam(string ifname); 
  void WriteParam();
  void FitGaus(TH1F *h, double &gamin, double &gamax, double range=2.5,int itr = 10);

  bool PiMinus(int charge, double momentum, double beta);
  bool PiPlus(int charge, double momentum, double beta);
  int TagBHit(double *TagBTime, double IHTime, int *seg);
  bool SetEgamma(void);
  Settings *set;
  ToF_ParamMan *tofParamMan;

  private:
    TFile *ofp;
    string pdf_name;
    string input_param;
    string output_param;
    int GetMaxEvent() { return ENumMax; }
    int ENumMax;
    int ENum;
    int runmin,runmax;
    //PID Mass Range
    double MassRange_pion_min;
    double MassRange_pion_max;
    double oa_min;
    double oa_max;
    double TPE[40];

    double gamin, gamax;
    double mean,e_mean,p0;//p0->default param
    

    TLorentzVector P_inv;    //lorentz vector of invariant mass
    TLorentzVector P_pi, P_p;//lorentz vector of pi^{-] and proton
    TVector3 vec[3];//mom(vector) of pi^{-] and proton

    TH1F *h_cointime_tb[40], *h_cointime_tdlr[10], *h_cointime_tdll[10];
    TH1F *h_time_at_ver_tdlr[10], *h_time_at_ver_tdll[10];
    TH1F *h_flight_time_tdlr[10], *h_flight_time_tdll[10];

    TH2F *h2_cointime_mom_tdlr[10]   , *h2_cointime_mom_tdll[10];
    TH2F *h2_time_at_ver_mom_tdlr[10], *h2_time_at_ver_mom_tdll[10];
    TH2F *h2_ftime_mom_tdlr[10], *h2_ftime_mom_tdll[10];
    TH2F *h2_ftime_beta_tdlr[10], *h2_ftime_beta_tdll[10];

    TH2F *h2_dvol, *h2_dvol_oa, *h2_dvol_ppi;
    TH2F *h2_pid_all, *h2_pid_add;
    TH1F *h_beta, *h_oa, *h_oa_ee, *h_oa_ppi;

    TLatex *tex;
    TLine *line;
    int run_num;
    double beta[3];
    double beta_pi, fl_pi, ft_pi;
    TCanvas *cc[NCanvas];
    TLatex *tex_fl[4][12];
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
RKVertex_mkparam_ft_corr::RKVertex_mkparam_ft_corr()
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

  MassRange_pion_min   = 0.0;
  MassRange_pion_max   = 0.8;       

  for(int i=0;i<NCanvas;i++){
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
  tofParamMan = new ToF_ParamMan();
}
////////////////////////////////////////////////////////////////////////////
RKVertex_mkparam_ft_corr::~RKVertex_mkparam_ft_corr(){
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam_ft_corr::makehist(string ofname){
cout<<"makehist"<<endl;
  ofp = new TFile(Form("%s",ofname.c_str()),"recreate");

  ofp->mkdir("cointime"); ofp->cd("cointime");
  for(int i=0;i<10;i++){//TDL seg
    h_cointime_tdll[i]    = new TH1F(Form("h_cointime_tdll%d",i+1)   ,Form("h_cointime_tdll%d",i+1)    , 100, -2, 2);
    h_cointime_tdlr[i]    = new TH1F(Form("h_cointime_tdlr%d",i+1)   ,Form("h_cointime_tdlr%d",i+1)    , 100, -2, 2);
    h_time_at_ver_tdll[i] = new TH1F(Form("h_time_at_ver_tdll%d",i+1),Form("h_time_at_ver_tdll%d",i+1) , 100, -2, 2);
    h_time_at_ver_tdlr[i] = new TH1F(Form("h_time_at_ver_tdlr%d",i+1),Form("h_time_at_ver_tdlr%d",i+1) , 100, -2, 2);
    h_flight_time_tdll[i] = new TH1F(Form("h_flight_time_tdll%d",i+1),Form("h_flight_time_tdll%d",i+1) , 200,  0, 3);
    h_flight_time_tdlr[i] = new TH1F(Form("h_flight_time_tdlr%d",i+1),Form("h_flight_time_tdlr%d",i+1) , 200,  0, 3);

    h2_cointime_mom_tdll[i]   = new TH2F(Form("h2_cointime_mom_tdll%d",i+1)   ,Form("h2_cointime_mom_tdll%d",i+1)    , 100, -2, 2,100,0,1);
    h2_cointime_mom_tdlr[i]   = new TH2F(Form("h2_cointime_mom_tdlr%d",i+1)   ,Form("h2_cointime_mom_tdlr%d",i+1)    , 100, -2, 2,100,0,1);
    h2_time_at_ver_mom_tdll[i]= new TH2F(Form("h2_time_at_ver_mom_tdll%d",i+1),Form("h2_time_at_ver_mom_tdll%d",i+1) , 100, -2, 2,100,0,1);
    h2_time_at_ver_mom_tdlr[i]= new TH2F(Form("h2_time_at_ver_mom_tdlr%d",i+1),Form("h2_time_at_ver_mom_tdlr%d",i+1) , 100, -2, 2,100,0,1);
    h2_ftime_mom_tdll[i]  = new TH2F(Form("h2_ftime_mom_tdll%d",i+1) ,Form("h2_ftime_mom_tdll%d",i+1)  , 200,0,3,100,0,1);
    h2_ftime_mom_tdlr[i]  = new TH2F(Form("h2_ftime_mom_tdlr%d",i+1) ,Form("h2_ftime_mom_tdlr%d",i+1)  , 200,0,3,100,0,1);
    h2_ftime_beta_tdll[i] = new TH2F(Form("h2_ftime_beta_tdll%d",i+1),Form("h2_ftime_beta_tdll%d",i+1) , 200,0,3,100,0,1);
    h2_ftime_beta_tdlr[i] = new TH2F(Form("h2_ftime_beta_tdlr%d",i+1),Form("h2_ftime_beta_tdlr%d",i+1) , 200,0,3,100,0,1);

    set->SetTH1(h_cointime_tdll[i],Form("cointime TDLL%d",i+1)  ,"cointime[ns]"      ,"Counts"       , 1, 3001, 2);
    set->SetTH1(h_cointime_tdlr[i],Form("cointime TDLR%d",i+1)  ,"cointime[ns]"      ,"Counts"       , 1, 3001, 2);
    set->SetTH1(h_time_at_ver_tdll[i],Form("Time at vertex TDLL%d",i+1)  ,"cointime[ns]"      ,"Counts"       , 1, 3001, 2);
    set->SetTH1(h_time_at_ver_tdlr[i],Form("Time at vertex TDLR%d",i+1)  ,"cointime[ns]"      ,"Counts"       , 1, 3001, 2);
    set->SetTH1(h_flight_time_tdll[i],Form("Flight time TDLL%d",i+1)  ,"flight_time[ns]"      ,"Counts"       , 1, 3001, 2);
    set->SetTH1(h_flight_time_tdlr[i],Form("Flight time TDLR%d",i+1)  ,"flight_time[ns]"      ,"Counts"       , 1, 3001, 2);
  }
  for(int i=0;i<40;i++){//TagB seg
    h_cointime_tb[i] = new TH1F(Form("h_cointime_tb%d",i+1),Form("h_cointime_tb%d",i+1) , 200,-10,10);
    set->SetTH1(h_cointime_tb[i],Form("cointime TagB%d",i+1)  ,"cointime[ns]"      ,"Counts"       , 1, 3001, 4);
  }

  ofp->cd();

  ofp->mkdir("vertex"); ofp->cd("vertex");
  h2_pid_all  = new TH2F("h2_pid_all" ,"h2_pid_all"  , 500,    0,   7,1000,    -1,    1);
  h2_dvol     = new TH2F("h2_dvol"    ,"h2_dvol"     , 200, -10., 10., 106, -10.6, 10.6);
  h2_dvol_oa  = new TH2F("h2_dvol_oa" ,"h2_dvol_oa"  , 200, -10., 10., 106, -10.6, 10.6);
  h2_dvol_ppi = new TH2F("h2_dvol_ppi","h2_dvol_ppi" , 200, -10., 10., 106, -10.6, 10.6);
  h_oa        = new TH1F("h_oa"     ,"h_oa"    , 200, -1., 1.);
  h_oa_ee     = new TH1F("h_oa_ee"  ,"h_oa_ee" , 200, -1., 1.);
  h_oa_ppi    = new TH1F("h_oa_ppi" ,"h_oa_ppi", 200, -1., 1.);
  h_beta    = new TH1F("h_beta" ,"h_beta", 100, 0., 1.5);
  set->SetTH2(h2_dvol    ,"Vertex Position"           ,"x [cm]","y [cm]");
  set->SetTH2(h2_dvol_oa ,"Vertex Position(OA cut)"   ,"x [cm]","y [cm]");
  set->SetTH2(h2_dvol_ppi,"Vertex Position(p-#pi^{-})","x [cm]","y [cm]");
  set->SetTH1(h_oa      ,"Opening Angle"                 ,"cos#theta"        ,"Counts"       , 1, 3000, 1);
  set->SetTH1(h_oa_ee   ,"Opening Angle(e^{+} e{-})"     ,"cos#theta"        ,"Counts"       , 1, 3000, 1);
  set->SetTH1(h_oa_ppi  ,"Opening Angle(p #pi^{-})"      ,"cos#theta"        ,"Counts"       , 1, 3000, 1);
  set->SetTH1(h_beta    ,"#beta(#pi)"                    ,"#beta"            ,"Counts"       , 1, 3000, 1);
  ofp->cd();


}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam_ft_corr::loop(){
cout<<"start loop"<<endl;
runmin=99999;runmax=-99;

int ppicounter=0;
int pcounter=0;
int picounter=0;
bool ppi_flag;//proton & pi^{-} vertex
bool pibeta_flag;
bool oa_flag, tag_flag, dedx_flag, ohseg_flag, vpos_flag, width_flag;
bool chi2_flag, dca_flag, mm_flag;
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
    oa_flag = tag_flag = dedx_flag=ohseg_flag=vpos_flag=false;
    chi2_flag = dca_flag = mm_flag = false;

    if(n%100000==0) cout<<n<<" / "<<ENum<<endl;

    for(int i=0;i<ntr;i++){
      vec[i].SetXYZ(momx[i], momy[i], momz[i]);
    }//for each track
    cos_oa = vec[0]*vec[1]/(vec[0].Mag()*vec[1].Mag());

    if(vertex[0]>vx_min && vertex[0]<vx_max && vertex[1]>vy_min && vertex[1]<vy_max)vpos_flag=true;
    h2_dvol    ->Fill(vertex[0],vertex[1]);//x,y
    if(dca<10.)dca_flag=true;
    if(chi2[0]<100.&&chi2[1]<100.)chi2_flag=true;
    //if(cos_oa>oa_min&&cos_oa<oa_max)oa_flag=true;
    //h2_dvol_oa ->Fill(vertex[0],vertex[1]);

    

      ///////////////
      //pion       //
      ///////////////
    bool really_pion=false;
    for(int i=0;i<ntr;i++){
      pibeta_flag = false;
      if(ass_noh[i]<1||ass_nih[i]<1)continue;//if no hit in hodoscopes
      beta[i] = (fl_oh[i] - fl_ih[i])/(ohct[i] - ihct[i])/(c*100.);
      h_beta->Fill(beta[i]);


      if(cos_oa<0.8){
        if(PiMinus(charge[i],mom[i],beta[i])||PiPlus(charge[i],mom[i],beta[i])) {
          pi_id=i;
          beta_pi = mom[i]/sqrt(mom[i]*mom[i]+Mpi*Mpi*1e-6);
          if((beta[i]>(beta_pi-0.2))&&(beta[i]<(beta_pi+0.2))){pibeta_flag = true;}
          fl_pi=sqrt( (vertex[0]-ih_x[i])*(vertex[0]-ih_x[i]) + (vertex[1]-ih_y[i])*(vertex[1]-ih_y[i])  );
          //ft_pi = fl_pi/(beta[i]*c*100.);
          ft_pi = fl_pi/(beta_pi*c*100.);
          if(pibeta_flag&&dca_flag&&chi2_flag&&vpos_flag){
          //////////////
          //TagB info.//
          //////////////
            double tdlmctime,taggertime=-999.;
            if(ihlr[pi_id]==-1)    tdlmctime=tdllmctime[ihseg[pi_id]-1];
            else if(ihlr[pi_id]==1)tdlmctime=tdlrmctime[ihseg[pi_id]-1];
            int seg[3];
            int NHit = TagBHit(tagbctime, tdlmctime, seg);
            taggertime=tagbctime[seg[0]-1];
            Egamma = -1.; 
            for(int i=0;i<24;i++){
              if(tagbctime[i]>-100){
                if(ihlr[pi_id]==-1)    {
                  h_cointime_tdll[ihseg[pi_id]-1]         ->Fill(tagbctime[i]-tdlmctime);
                  h_time_at_ver_tdll[ihseg[pi_id]-1]      ->Fill(tagbctime[i]-tdlmctime -ft_pi);
                  h2_cointime_mom_tdll[ihseg[pi_id]-1]    ->Fill(tagbctime[i]-tdlmctime       , mom[pi_id]);
                  h2_time_at_ver_mom_tdll[ihseg[pi_id]-1] ->Fill(tagbctime[i]-tdlmctime -ft_pi, mom[pi_id]);
                }
                else if(ihlr[pi_id]==1){
                  h_cointime_tdlr[ihseg[pi_id]-1]         ->Fill(tagbctime[i]-tdlmctime);
                  h_time_at_ver_tdlr[ihseg[pi_id]-1]      ->Fill(tagbctime[i]-tdlmctime -ft_pi);
                  h2_cointime_mom_tdlr[ihseg[pi_id]-1]    ->Fill(tagbctime[i]-tdlmctime       , mom[pi_id]);
                  h2_time_at_ver_mom_tdlr[ihseg[pi_id]-1] ->Fill(tagbctime[i]-tdlmctime -ft_pi, mom[pi_id]);
                }
              }
            }
            if(NHit == 1){
              Egamma = TPE[seg[0]-1];
              h_cointime_tb[seg[0]-1] ->Fill(tagbctime[seg[0]-1]-tdlmctime);
              tag_flag = true;
            }
            else if(NHit == 2 && abs(seg[0]-seg[1])==1){
              Egamma = (TPE[seg[0]-1]+TPE[seg[1]-1])/2; 
              tag_flag = true;
            }
          //////////////

            really_pion=true;
            h_oa ->Fill(cos_oa);
            h2_pid_all->Fill(1./beta[i], mom[i]*charge[i]);
            if(ihlr[i]==-1){
              h_flight_time_tdll[ihseg[i]-1] ->Fill(ft_pi);
              h2_ftime_mom_tdll[ihseg[i]-1]  ->Fill(ft_pi,mom[i]);
              h2_ftime_beta_tdll[ihseg[i]-1] ->Fill(ft_pi,beta_pi);
              //h2_ftime_beta_tdll[ihseg[i]-1] ->Fill(ft_pi,beta[i]);
            }
            else if(ihlr[i]==1){
              h_flight_time_tdlr[ihseg[i]-1] ->Fill(ft_pi);
              h2_ftime_mom_tdlr[ihseg[i]-1]  ->Fill(ft_pi,mom[i]);
              h2_ftime_beta_tdlr[ihseg[i]-1] ->Fill(ft_pi,beta_pi);
              //h2_ftime_beta_tdlr[ihseg[i]-1] ->Fill(ft_pi,beta[i]);
            }
          }
        }
      }
    }//for each track


    if(pi_id>-1){
    }

    if(pi_id>-1&&p_id>-1)ppi_flag=true;
    

  }//event loop
  ofp->Write();
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam_ft_corr::fit_TDL(){
int lr,cid;
cid=15;
  for(int i=1; i<9; i++){//Fitting TDLL
    lr=-1;
    FitGaus(h_time_at_ver_tdll[i],gamin,gamax,2.0);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_time_at_ver_tdll[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    if(e_mean<1.)tofParamMan->SetToF( lr,cid,i+1,mean);
  }

  for(int i=1; i<9; i++){//Fitting TDLR
    lr= 1;
    FitGaus(h_time_at_ver_tdlr[i],gamin,gamax,2.0);
    TF1 * fit = new TF1("fit", "gaus", gamin, gamax);
    set->SetTF1(fit, 2, 1,1);
    h_time_at_ver_tdlr[i]->Fit(fit,"RQ");
    mean  =fit->GetParameter(1);
    e_mean=fit->GetParError(1);
    if(e_mean<1.)tofParamMan->SetToF( lr,cid,i+1,mean);
  }

}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam_ft_corr::draw(){
  cc[0]->Clear();
  cc[0]->Divide(5,4);
  cc[0]->cd(1); gPad->SetLogz(1);h2_pid_all->Draw("colz");
  cc[0]->cd(2); gPad->SetLogy(1);h_oa ->Draw();

  cc[1]->Clear();
  cc[1]->Divide(5,4);
  for(int i=0;i<10;i++){
    cc[1]->cd(i+1); h_flight_time_tdll[i]->Draw("");
    cc[1]->cd(i+11);h_flight_time_tdlr[i]->Draw("");
  }

  cc[2]->Clear();
  cc[2]->Divide(5,4);
  for(int i=0;i<10;i++){
    cc[2]->cd(i+1); h_cointime_tdll[i]->Draw("");
    cc[2]->cd(i+11);h_cointime_tdlr[i]->Draw("");
  }

  cc[3]->Clear();
  cc[3]->Divide(5,4);
  for(int i=0;i<10;i++){
    cc[3]->cd(i+1); h_time_at_ver_tdll[i]->Draw("");
    cc[3]->cd(i+11);h_time_at_ver_tdlr[i]->Draw("");
  }

  cc[4]->Clear();
  cc[4]->Divide(5,4);
  for(int i=0;i<10;i++){
    cc[4]->cd(i+1); gPad->SetLogz(1);h2_ftime_mom_tdll[i]->Draw("colz");
    cc[4]->cd(i+11);gPad->SetLogz(1);h2_ftime_mom_tdlr[i]->Draw("colz");
  }

  //cc[5]->Clear();
  //cc[5]->Divide(5,4);
  //for(int i=0;i<10;i++){
  //  cc[5]->cd(i+1); gPad->SetLogz(1);h2_ftime_beta_tdll[i]->Draw("colz");
  //  cc[5]->cd(i+11);gPad->SetLogz(1);h2_ftime_beta_tdlr[i]->Draw("colz");
  //}

  cc[5]->Clear();
  cc[5]->Divide(5,4);
  for(int i=0;i<10;i++){
    cc[5]->cd(i+1); gPad->SetLogz(1);h2_cointime_mom_tdll[i]->Draw("colz");
    cc[5]->cd(i+11);gPad->SetLogz(1);h2_cointime_mom_tdlr[i]->Draw("colz");
  }

  cc[6]->Clear();
  cc[6]->Divide(5,4);
  for(int i=0;i<10;i++){
    cc[6]->cd(i+1); gPad->SetLogz(1);h2_time_at_ver_mom_tdll[i]->Draw("colz");
    cc[6]->cd(i+11);gPad->SetLogz(1);h2_time_at_ver_mom_tdlr[i]->Draw("colz");
  }

}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam_ft_corr::SetPdfFilename(string ifname){
  pdf_name = ifname;
} 
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam_ft_corr::SetInputParam(string ifname){
  input_param = ifname;
} 
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam_ft_corr::SetOutputParam(string ifname){
  output_param = ifname;
} 
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam_ft_corr::WriteParam(){
  tofParamMan -> write_param(output_param.c_str());
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_mkparam_ft_corr::savecanvas(){
  cc[0] ->Print(Form("%s[",pdf_name.c_str()) );
  for(int i=0;i<NCanvas;i++){
  cc[i] ->Print(Form("%s" ,pdf_name.c_str()) );
  }
  cc[NCanvas-1] ->Print(Form("%s]",pdf_name.c_str()) );
  cout<<pdf_name<<" saved!"<<endl;
}
////////////////////////////////////////////////
void RKVertex_mkparam_ft_corr::SetRoot(string ifname)
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
void RKVertex_mkparam_ft_corr::SetRootSingle(string ifname)
{
  cout<<"SetRoot single"<<endl;
  add(ifname);
  readtree();
  ENum = GetEntries();
}
////////////////////////////////////////////////
void RKVertex_mkparam_ft_corr::GetHist(string ifname)
{
  cout<<"Get Histogram from root file"<<endl;
  TFile *file = new TFile(ifname.c_str() );
  for(int i=0;i<10;i++){//TDL seg
    h_cointime_tdll[i]    = (TH1F*)file->Get(Form("cointime/h_cointime_tdll%d",i+1));
    h_cointime_tdlr[i]    = (TH1F*)file->Get(Form("cointime/h_cointime_tdlr%d",i+1));
    h_time_at_ver_tdll[i] = (TH1F*)file->Get(Form("cointime/h_time_at_ver_tdll%d",i+1));
    h_time_at_ver_tdlr[i] = (TH1F*)file->Get(Form("cointime/h_time_at_ver_tdlr%d",i+1));
    h_flight_time_tdll[i] = (TH1F*)file->Get(Form("cointime/h_flight_time_tdll%d",i+1));
    h_flight_time_tdlr[i] = (TH1F*)file->Get(Form("cointime/h_flight_time_tdlr%d",i+1));

    h2_cointime_mom_tdll[i]   = (TH2F*)file->Get(Form("cointime/h2_cointime_mom_tdll%d",i+1)   );
    h2_cointime_mom_tdlr[i]   = (TH2F*)file->Get(Form("cointime/h2_cointime_mom_tdlr%d",i+1)   );
    h2_time_at_ver_mom_tdll[i]= (TH2F*)file->Get(Form("cointime/h2_time_at_ver_mom_tdll%d",i+1));
    h2_time_at_ver_mom_tdlr[i]= (TH2F*)file->Get(Form("cointime/h2_time_at_ver_mom_tdlr%d",i+1));
    h2_ftime_mom_tdll[i]      = (TH2F*)file->Get(Form("cointime/h2_ftime_mom_tdll%d",i+1) );
    h2_ftime_mom_tdlr[i]      = (TH2F*)file->Get(Form("cointime/h2_ftime_mom_tdlr%d",i+1) );
    h2_ftime_beta_tdll[i]     = (TH2F*)file->Get(Form("cointime/h2_ftime_beta_tdll%d",i+1));
    h2_ftime_beta_tdlr[i]     = (TH2F*)file->Get(Form("cointime/h2_ftime_beta_tdlr%d",i+1));
  }
  for(int i=0;i<40;i++){//TagB seg
    h_cointime_tb[i] = (TH1F*)file->Get(Form("cointime/h_cointime_tb%d",i+1));
  }
  h2_pid_all  = (TH2F*)file->Get("vertex/h2_pid_all" );
  h2_dvol     = (TH2F*)file->Get("vertex/h2_dvol"    );
  h2_dvol_oa  = (TH2F*)file->Get("vertex/h2_dvol_oa" );
  h2_dvol_ppi = (TH2F*)file->Get("vertex/h2_dvol_ppi");
  h_oa        = (TH1F*)file->Get("vertex/h_oa"       );
  h_oa_ee     = (TH1F*)file->Get("vertex/h_oa_ee"    );
  h_oa_ppi    = (TH1F*)file->Get("vertex/h_oa_ppi"   );
  h_beta      = (TH1F*)file->Get("vertex/h_beta"     );
cout<<"end of gethist"<<endl;
}

////////////////////////////////////////////////
int RKVertex_mkparam_ft_corr::TagBHit(double *TagBTime, double IHTime, int *seg){

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

bool RKVertex_mkparam_ft_corr::SetEgamma(void){

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

bool RKVertex_mkparam_ft_corr::PiMinus(int charge, double momentum, double beta){

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

bool RKVertex_mkparam_ft_corr::PiPlus(int charge, double momentum, double beta){

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
void RKVertex_mkparam_ft_corr::FitGaus(TH1F *h, double &gamin, double &gamax, double range,int itr){
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
  string input_param  = "param/tof/default.param";
  string output_param = "param/tof/default.param";
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
      cout<<"-f : input root filename list"<<endl;
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
  RKVertex_mkparam_ft_corr *mkparam = new RKVertex_mkparam_ft_corr();

  mkparam->SetEgamma();
  mkparam->SetInputParam(input_param);//Set parameter file used to make tree
  mkparam->SetOutputParam(output_param);
  mkparam->SetMaxEvent(MaxNum);
  if(skip_flag)mkparam->GetHist(ifname);
  else{
    //mkparam->SetRootSingle(ifname);
    mkparam->SetRoot(ifname);
    mkparam->makehist(ofname);
    mkparam->loop();
  }
  mkparam->fit_TDL();
  mkparam->draw();
  mkparam->SetPdfFilename(pdfname);
  mkparam->savecanvas();
  mkparam->WriteParam();
  delete mkparam;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}

