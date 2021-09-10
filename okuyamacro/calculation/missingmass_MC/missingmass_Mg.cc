#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <string>
#include <numeric>//for accumulate()
using namespace std;

#include "TApplication.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TArrow.h"
#include "TDatime.h"
#include "TMarker.h"
#include "TRandom.h"
#include "Math/Integrator.h"
double fmm_total(double *x, double *par){
  double val = par[0] * TMath::Gaus(x[0],par[1],par[2],kTRUE);//s-orbit
  val += par[3] * TMath::Gaus(x[0],par[4],par[5], kTRUE);//p-orbit
  val += par[6] * TMath::Gaus(x[0],par[7],par[8], kTRUE);//p-orbit
  val += par[9] * TMath::Gaus(x[0],par[10],par[11],kTRUE);//p-orbit
  val += par[12];//Accidentals
  if(x[0]>0)val += par[13] + par[14]*x[0] + par[15]*x[0];//Accidentals

  return val;
}

static const double PI = 4.*atan(1.);
static const double deg_to_rad = PI / 180.;//degree -> radian
static const double rad_to_deg = 180. / PI;//radian -> degree
static const double sigma_to_fwhm = 2.*sqrt(2.*log(2.));//sigma -> FWHM
static const double fwhm_to_sigma = 1./sigma_to_fwhm;//FWHM -> sigma
static const double cm_to_barn = 1e+24;// cm^2 --> b
static const double alpha = 1./137.035999074; // fine structure constant
static const double hbar = 6.58211928*1e-22;  // Planc constant (reduced) (MeV x s)
static const double hbarc = 197.3269718;      // conversion constant (MeV x fm)
static const double kb = 8.6173324*1e-5;      // Boltzmann constant
static const double e = 1.602176565*1e-19;    // electron charge magnitude (C)
static const double c = 0.299792458;          // speed of light in vacuum (m/ns)
static const double re = 2.817e-13;           // classical electron radius (cm)
static const double Na = 6.02214129*1e+23;    // Avogadro constant
static const double Me = 0.510998928;         // electron     mass (MeV/c2)
static const double Mmu = 105.6583715;        // muon         mass (MeV/c2)
static const double Mpi = 139.57018;          // charged pion mass (MeV/c2)
static const double Mpi0 = 134.9766;          // charged pion mass (MeV/c2)
static const double MK = 493.677;             // charged Kaon mass (MeV/c2)
static const double Mp = 938.272046;          // proton       mass (MeV/c2)
static const double Mn = 939.565379;          // proton       mass (MeV/c2)
static const double Mu = 931.494061;          // proton       mass (MeV/c2)
static const double ML = 1115.683;            // Lambda       mass (MeV/c2)
static const double MS0 = 1192.642;           // Sigma Zero   mass (MeV/c2)
static const double MSm = 1197.449;           // Sigma Minus  mass (MeV/c2)
static const double MSp = 1189.37;            // Sigma Plus   mass (MeV/c2)


double land(double *x, double *par){
  double val;
  double PI = 4.*atan(1.);

  val =  exp(-x[0]*log(x[0])-par[0]*x[0]) * sin(PI*x[0]);

  return val;
}

//Leo: Landau's theory
double dEdx(double Delta, double thick, double Minc, double mom){ // dE, target thick, mass, momentum
  gErrorIgnoreLevel = kError;

  double cons = 2.*PI*Na*re*re*Me;
  double ro = 1.0;              // density (g/cm3)
//  double thick = 0.1;  // thickness (cm)
  double Z = 1;                 // atomic number of absorbing material
  double A = 1;               // atomic weight of absorbing material
  double M = Minc;               // incident particle mass (MeV/c2)
  double z = 1.;                // charge of incident particle
//  double ss = Me/M;
  double delta;                          // density correction
  double C = 0;                              // shell correction
//  double W;                              // maximum energy transfer in a single collision
  double I;                              // mean excitation potential
  double C0, a, m, X0, X1;               // parameters of delta
  double Ne = Na * ro * Z / A;           // density of electrons
  double mup = sqrt(80.617 * 1e+6 * Ne); // plasma frequency (cm^3/2 x Hz)

  if(Z<13) I = (12.*Z + 7.) * 1e-6;
  else I = (9.76 + 58.8*pow(Z,-1.19)) * Z * 1e-6; // MeV
  C0 = - (2. * log(I/(2.*PI) / hbar / mup) + 1.); // MeV
  a = 0.1;
  m = 3.0;
  int tar = 4; // C:1, Ca40:2, Scin:3, Al:4
  switch(tar){
  case 1:  // C
    ro = 1.80;  Z = 6;  A = 12.011;
    I = 78 * 1e-6;  C0 = -2.99;
    a = 0.2024;  m = 3.00;
    X1 = 2.486;  X0 = -0.0351;
    break;
  case 2: // Ca40
    ro = 1.58;  Z = 20;  A = 40;
    Ne = Na * ro * Z / A;  mup = sqrt(80.617 * 1e+6 * Ne);
    I = 64.7 * 1e-6;  C0 = -3.20;
    a = 0.1610;  m = 3.24;
    X1 = 2.49;  X0 = 0.1464;
    break;
  case 3: // Scinti
    ro = 1.032;  Z = 6.5;  A = 13;
    Ne = Na * ro * Z / A;  mup = sqrt(80.617 * 1e+6 * Ne);
    I = 64.7 * 1e-6;  C0 = -3.20;
    a = 0.1610;  m = 3.24;
    X1 = 2.49;  X0 = 0.1464;
    break;
  case 4: // Al (Ref. Leo text)
    ro = 2.699;  Z = 13;  A = 27;
    Ne = Na * ro * Z / A;  mup = sqrt(80.617 * 1e+6 * Ne);
    I = 166.0 * 1e-6;  C0 = -4.24;
    a = 0.0802;  m = 3.63;
    X1 = 3.01;  X0 = 0.1708;
    break;
  default:
    std::cout<<"unknown target! break!"<<std::endl;
    break;
  }

  double Fai;

  double p = mom;//momentum
  double beta = p / sqrt(p*p + M*M);
  double gamma = sqrt(p*p + M*M) / M;
  double eta = p / M;
  double ss = Me / M;
  double X = log10(beta*gamma);

//  double W = 2.*Me*eta*eta / (1.+2.*ss*sqrt(1+eta*eta)+ss*ss);
  double nonrela = cons * ro * Z/A * z*z/beta/beta;
//  double rela = log(2.*Me*gamma*gamma*beta*beta*W/I/I)-2.*beta*beta;
  if(X<X0) delta = 0;
  else if(X>=X0 && X<X1) delta = 4.6052*X + C0 + a*pow(X1-X,m);
  else delta = 4.6052*X + C0;
  if(eta<0.1) C = 0;
  else C = (0.422377*pow(eta,-2) + 0.0304043*pow(eta,-4) - 0.00038106*pow(eta,-6)) * 1e-6 * pow(I*1e+6,2)
          +(3.850190*pow(eta,-2) - 0.1667989*pow(eta,-4) + 0.00159755*pow(eta,-6)) * 1e-9 * pow(I*1e+6,3);
  double Delta_b = nonrela * thick;//bar{Delta}

  double LnEps = log((1.-beta*beta)*I*I / (2.*Me*beta*beta)) + beta*beta;//ln(eps)
  double Xi = Delta_b * (log(Delta_b) - LnEps + 0.198 -delta);//Delta_mp (most probable)

  TF1 *f;

  double Lambda = 1./Delta_b * (Xi - Delta_b * (log(Delta_b)- LnEps + 1. - 0.577));//lambda
  if(Lambda<-4) Fai = 0;
  else{
    f = new TF1("f",land,0,1e+4,1);
    f->SetParameter(0,Lambda);
    Fai = f->Integral(0,1e+2,1e-6) / PI;//varphi(lambda)
    f->Clear();
  }
  double Func_max = Fai / Delta_b;//f(thick,dE)


  Lambda = 1./Delta_b * (Delta - Delta_b * (log(Delta_b)- LnEps + 1. - 0.577));
  if(Lambda<-4) Fai = 0;
  else{
    f = new TF1("f",land,0,1e+4,1);
    f->SetParameter(0,Lambda);
    Fai = f->Integral(0,1e+2,1e-6) / PI;
    f->Clear();
  }

  double Func = Fai / Delta_b;//
//  cout<<Delta<<"  "<<Xi<<"  "<<Func_max<<"  "<<Func<<endl;
//  Func /= Func_max;

  return Func;
}

void SetTH1(TH1 *h, TString hname, TString xname, TString yname, int LColor, int LStyle, int FStyle, int FColor){
  h->SetTitle(hname);
  h->SetLineColor(LColor);
  h->SetLineStyle(LStyle);
  //h->SetLineWidth(1);
  h->SetLineWidth(2);
  h->SetFillStyle(FStyle);
  h->SetFillColor(FColor);
//  h->SetMinimum(0.8);
  h->SetMarkerColor(FColor);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(3);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double p2E(double p, double M){
  return sqrt(p*p + M*M);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double E2p(double E, double M){
  return sqrt(E*E - M*M);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double GetQ2(double Ee, double Eep, double theta){  // Q2 (MeV^2)
// Q2 = - q^{\mu} * q_{\mu}
  TVector3 v_pe;  v_pe.SetMagThetaPhi(E2p(Ee,Me),0,0);
  TVector3 v_pep; v_pep.SetMagThetaPhi(E2p(Eep,Me),theta,0);
  TLorentzVector lv_e(v_pe,Ee);
  TLorentzVector lv_ep(v_pep,Eep);
//  return 2.*Ee*Eep*(1.-cos(theta)); // URA
//  return 2.*(Ee*Eep - Me*Me - E2p(Ee,Me)*E2p(Eep,Me)*cos(theta)); // Exact
  return -(lv_e - lv_ep)*(lv_e - lv_ep);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double GetVPEnergy(double Ee, double Eep){
  return Ee - Eep;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double GetVPMomentum(double Ee, double Eep, double theta){
  return sqrt( GetQ2(Ee,Eep,theta) + GetVPEnergy(Ee,Eep)*GetVPEnergy(Ee,Eep) );
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double GetVPInvariantEnergy(double Ee, double Mt, double Eep, double theta){
  TVector3 v_pt; v_pt.SetMagThetaPhi(0,0,0);
  TVector3 v_pb; v_pb.SetMagThetaPhi(GetVPMomentum(Ee,Eep,theta),PI,0);
  TLorentzVector lv_t(v_pt,Mt);
  TLorentzVector lv_b(v_pb,(Double_t)GetVPEnergy(Ee,Eep));
//  return sqrt( Mt*Mt - GetQ2(Ee,Eep,theta) + 2.*Mt*GetVPEnergy(Ee,Eep) );
  return sqrt((lv_t + lv_b)*(lv_t + lv_b));
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double GetVPEquivalentEnergy(double Ee, double Mt, double Eep, double theta){
  double s = GetVPInvariantEnergy(Ee,Mt,Eep,theta);
  return ( s*s - Mt*Mt ) / (2.*Mt);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double GetVPTheta(double E, double Ep, double theta){
  double Q2    = GetQ2(E, Ep, theta);
  double omega = E - Ep;
  double vptheta = acos((E-Ep*(1.-Q2/2./E/Ep))/sqrt(omega*omega+Q2));

  return vptheta;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main(int argc, char** argv){
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  TApplication theApp("App", &argc, argv);
  gRandom->SetSeed(0);
ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6); 

// Inputs
  double E_beam     = 4240;  // 4240 2344 4318
  double M_target   = 25133.144;// t: 2808.921 C12:11174.864  Ca40:37214.521  Ca48:44657.300 Al27:25133.144
// CS from M. Iodice et al., Phys. Rev. Lett. 99 (2007), 052501.
//  double Ex[11] = {0.00, 0.14, 2.67, 5.74, 5.85, 10.48, 10.52, 10.98, 11.05, 12.95, 13.05};
//  double CS[11] = {1.02, 3.66, 1.54, 0.58, 0.18,  0.24,  0.12,  1.43,  2.19,  0.91,  0.27};
//  for(int i=0;i<11;i++){ Ex[i] -= 11.5; }
// CS from L. Tang et al.,
//  double Ex[9] = {-11.529, -11.348, -8.425, -5.488, -2.499, -1.220, -0.524, 0.223, 1.047};
//  double CS[9] = {   22.2,    77.7,   33.5,   26.0,   20.5,   31.5,   87.7,  46.3,  28.5};
// CS 
//  double Ex[1] = {  0.0};
//  double CS[1] = {   1.};
// CS P.Bydzovsky NPA881(2012)119. T.Motoba PTPS185(2010)224.
//  double Ex[30] = {-18.1, -18.0, -10.9, -10.8, -10.7, -7.0, -3.0, -3.1,  -3.2, -2.5, 4.0,  4.4,  4.8, -14.8, -8.3, -8.2, -0.8, -12.8, -11.8, -5.8, -4.5,  2, 3};
//  double CS[30] = { 13.7,  41.0,   8.1,  15.3,  90.8,  5.0, 30.2, 27.2, 108.1, 47.1, 50.6, 36.9, 86.3, 53.0, 11.2, 33.3,  5.2,    10,    10,   30,   30, 60,60};
// CS from 28Si(e,e'K+)28LAl 2021/8/18
//					1/2+(s)  1/2-,3/2-(p) 125 nb/sr -> divided by three.
  double CSp_tot = 125.*5./6.;//Lambda in p-orbit
  double CSpval1 = CSp_tot/3.;
  double CSpval2 = CSp_tot/3.;
  double CSpval3 = CSp_tot/3.;
  //CSpval1 = CSp_tot;
  //CSpval2 = 0.;
  //CSpval3 = 0.;
  double Expval1 = -6.3;
  double Expval2 = -3.3;
  double Expval3 = -0.8;
  //double Ex[1]={-16.3};
  //double Cs[1]={ 79.0};
  double Ex[19] = { -16.3, -15.2, -13.5, -12.0, Expval1, Expval1*7.5/8.3, Expval1*5.5/8.3, Expval1*4.0/8.3, Expval1*2.8/8.3, Expval2, (Expval2-Expval1)+Expval1*7.5/8.3, (Expval2-Expval1)+Expval1*5.5/8.3, (Expval2-Expval1)+Expval1*4.0/8.3, (Expval2-Expval1)+Expval1*2.8/8.3, Expval3, (Expval3-Expval1)+Expval1*7.5/8.3, (Expval3-Expval1)+Expval1*5.5/8.3, (Expval3-Expval1)+Expval1*4.0/8.3, (Expval3-Expval1)+Expval1*2.8/8.3};
  double CS[19] = {  79.0*5./6.,  26.3,  13.2,  35.1, CSpval1, CSpval1*100/320, CSpval1*60./320, CSpval1*100/320, CSpval1*60./320, CSpval2, CSpval2*100/320, CSpval2*60./320, CSpval2*100/320, CSpval2*60./320, CSpval3, CSpval3*100/320, CSpval3*60./320, CSpval3*100/320, CSpval3*60./320};
  int    CS_num = sizeof(CS)/sizeof(CS[0]);
  double CS_sum = accumulate(CS, CS + CS_num, 0.0f);
cout<<"CS_sum = "<<CS_sum<<" [nb/sr]"<<endl;
  double P_scat     = 2740;  // 2740  844  2100
  double Theta_scat = 6.5 * PI/180.;
  double Phi_scat   = 0.0 * PI/180.;
  double Theta_K    = 11.5 * PI/180.;
  double Phi_K      = 180.0 * PI/180.;
  double M_hyper    = 25318.694; // nnL:2992.227 B12L:11356.861  K40L:37382.452  K48L:44832.515 Mg27L:25318.694
                                                       // E05  E12-15
  double Ereso_beam     = 1.0 * 1E-4 * fwhm_to_sigma;  // 1.0  1.0
  double Preso_scat     = 1.08 * 1E-4 * fwhm_to_sigma;  // 4.2  2.0
  double Thetareso_scat = 0.389 * 1E-3 * fwhm_to_sigma;                  // 0.5  0.4  // rad
  double Phireso_scat   = 0.8 * 1E-3;                  // 1.0  0.4  // rad
  double Preso_K        = 2.32 * 1E-4 * fwhm_to_sigma;  // 2.0  2.0
  double Thetareso_K    = 0.487 * 1E-3 * fwhm_to_sigma;                  // 0.4  0.4  // rad  ~dY'
  double Phireso_K      = 0.4 * 1E-3;                  // 1.0  0.4  // rad  ~dX'

  double thick = 0.037*2.;  // target thickness in cm (0.05cm = 500um)
                         // E05  E12-15
  //int roop = 100000;       // 1800  800
  //int roop = 314.5*14.5;       // 1800  800
  int roop = 2.*(CS_sum/79.)*14.5*24.*0.5;// 2 [counts/hour] * CS factor * 14.5 [day(w/ JLab eff.)] * 24 [hours/day] * 0.5 (Eff.)
  //int QF   = roop * 0.;  // (50) 2.0% sticking
  //int Acc  = roop * 0.; // 37    3.7
  //int QF   = roop * 50.;  // (50) 2.0% sticking
  int QF   = roop * 30.;  // (50) 2.0% sticking
  //int Acc  = roop * 3.7; // 37    3.7
  double rate_e = 16.6*1000.;//Hz
  double rate_k_e05 = 110.;//Hz (E05-115:Carbon, 7.78msr, 87.5mg/cm2, 19.3uA, Chiba M-thesis)
  double rate_k = rate_k_e05 * (4.98/8.50) * (100./87.5) * (20./19.3) * (pow(27.,0.8-1.)/pow(12.,0.8-1.));//dOmega,Nt,Ie,A^(0.8)
  double rate_acc = rate_e * rate_k * 2.0E-9;//Hz
cout<<"Kaon Rate = "<<rate_k<<" [Hz]"<<endl;
cout<<"Accidnetal Rate = "<<rate_acc<<" [Hz]"<<endl;
  //int Acc  = roop * (rate_acc/6.1E-3); //okuyama
  int Acc  = rate_acc * 14.5 * 24. * 3600. ;//Acc. rate [Hz] * 14.5 [day] * 24 [day/hour] * 3600 [sec/hour] * (real factor)
cout<<"Acc = "<<Acc<<endl;
  //int Acc  = roop * 1.0; // 37    3.7
  //int Acc  = 4170; //akiyama

  //TH1D *h1_mm_p = new TH1D("h1_mm_p","h1_mm_p",700,-40,100);
  //TH1D *h1_mm_q = new TH1D("h1_mm_q","h1_mm_q",700,-40,100);
  //TH1D *h1_mm_a = new TH1D("h1_mm_a","h1_mm_a",700,-40,100);
  //TH1D *h1_mm   = new TH1D("h1_mm"  ,"h1_mm"  ,700,-40,100);
  TH1D *h1_dEe = new TH1D("dEe","dEe",400,0.,10.);
  TH1D *h1_dEep = new TH1D("dEep","dEep",400,0.,10.);
  TH1D *h1_dEk = new TH1D("dEk","dEk",400,0.,10.);
  TH1D *h1_mm_p = new TH1D("h1_mm_p","h1_mm_p",300,-20,10);
  TH1D *h1_mm_q = new TH1D("h1_mm_q","h1_mm_q",300,-20,10);
  TH1D *h1_mm_a = new TH1D("h1_mm_a","h1_mm_a",300,-20,10);
TH1D *h1_mm_a2= new TH1D("h1_mm_a2","h1_mm_a2",300,-20,10);
  TH1D *h1_mm   = new TH1D("h1_mm"  ,"h1_mm"  ,300,-20,10);
  TH1D *h1_mm2  = new TH1D("h1_mm2" ,"h1_mm2" ,300,-20,10);
  TH1D *h1_mm3  = new TH1D("h1_mm3" ,"h1_mm3" ,300,-20,10);
  SetTH1(h1_mm_p  ,"" ,"-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.1 MeV/#it{c}^{2})",kBlack,1,3004,kBlack);
  SetTH1(h1_mm_q  ,"" ,"-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.1 MeV/#it{c}^{2})",kRed,1,3004,kRed);
  SetTH1(h1_mm_a  ,"" ,"-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.1 MeV/#it{c}^{2})",kGreen ,1,1001,kGreen);
  SetTH1(h1_mm_a2 ,"" ,"-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.1 MeV/#it{c}^{2})",kGreen ,1,1001,kGreen);
  SetTH1(h1_mm    ,"" ,"-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.1 MeV/#it{c}^{2})",kBlack ,1,3004,kBlack);
  SetTH1(h1_mm2   ,"" ,"-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.1 MeV/#it{c}^{2})",kBlack ,1,3004,kBlack);
  SetTH1(h1_mm3   ,"" ,"-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.1 MeV/#it{c}^{2})",kBlack ,1,3004,kBlack);
  //SetTH1(h1_mm_p  ,"" ,"-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.2 MeV/#it{c}^{2})",kBlack   ,1,3004,kPink  +7);
  //SetTH1(h1_mm_q  ,"" ,"-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.2 MeV/#it{c}^{2})",kBlack   ,1,3004,kViolet+7);
  //SetTH1(h1_mm_a  ,"" ,"-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.2 MeV/#it{c}^{2})",kBlack   ,1,1001,kBlack   );
  //SetTH1(h1_mm    ,"" ,"-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.2 MeV/#it{c}^{2})",kBlack   ,1,3004,kPink  +7);
  //SetTH1(h1_mm_p,"Missing Mass","-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.2MeV/#it{c}^{2})",1,3001,3);
  //SetTH1(h1_mm_q,"Missing Mass","-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.2MeV/#it{c}^{2})",1,3003,4);
  //SetTH1(h1_mm_a,"Missing Mass","-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.2MeV/#it{c}^{2})",1);
  //SetTH1(h1_mm  ,"Missing Mass","-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.2MeV/#it{c}^{2})");
  //SetTH1(h1_mm_p,"","-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.2MeV/#it{c}^{2})",1,3001,3);
  //SetTH1(h1_mm_q,"","-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.2MeV/#it{c}^{2})",1,3003,4);
  //SetTH1(h1_mm_a,"","-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.2MeV/#it{c}^{2})",1);
  //SetTH1(h1_mm  ,"","-B_{#Lambda} (MeV/#it{c}^{2})","Counts / (0.2MeV/#it{c}^{2})");
  //h1_mm->SetStats(kFALSE);
  h1_mm->SetStats(1111);
  h1_mm->SetMarkerStyle(20);
  h1_mm->SetMarkerSize(0.8);
  h1_mm2->SetStats(1111);
  h1_mm2->SetMarkerStyle(20);
  h1_mm2->SetMarkerSize(0.8);
  h1_mm3->SetStats(1111);
  h1_mm3->SetMarkerStyle(20);
  h1_mm3->SetMarkerSize(0.8);

  int gs_num = 0;
  for(int n=0;n<roop;n++){
if(n%500==0){cout<<n<<"/"<<roop<<endl;}
// electron beam
    double M1     = Me;
    double E1     = E_beam;
    double P1     = E2p(E1,M1);
    double Theta1 = 0.;
    double Phi1   = 0.;
    TVector3 p1;
    p1.SetMagThetaPhi(P1,Theta1,Phi1);
    TLorentzVector lp1(p1,E1);

// target

//-- State Selection --//
    double Y = gRandom->Uniform(0,CS_sum);
    int peak_num = 0;
    for(int i=0;i<CS_sum;i++){
      double yy = accumulate(CS, CS + (i+1), 0.0f);
      if(Y<yy){
        if(i==0){ gs_num++;} //cout<<n<<"  "<<gs_num<<endl; }
        //if(i==0 || i==1){ gs_num++; cout<<n<<"  "<<gs_num<<endl; }
        break;
      }
      peak_num += 1;
    }
//---------------------//

    double M2     = M_target;
    double P2     = 0.;
    double E2     = p2E(M2,P2);
    double Theta2 = 0.;
    double Phi2   = 0.;
    TVector3 p2;
    p2.SetMagThetaPhi(P2,Theta2,Phi2);
    TLorentzVector lp2(p2,E2);

// scattered electron
    double M3     = Me;
    double P3     = P_scat;
    double E3     = p2E(M3,P3);
    double Theta3 = Theta_scat;
    double Phi3   = Phi_scat;
    TVector3 p3;
    p3.SetMagThetaPhi(P3,Theta3,Phi3);
    TLorentzVector lp3(p3,E3);

// virtual photon
    double E4      = GetVPEnergy(E1,E3);
//    double P4      = GetVPMomentum(E1,E3,Theta3);
//    double Q2      = GetQ2(E1,E3,Theta3);
//    double k_gamma = GetVPEquivalentEnergy(E1,M2,E3,Theta3);
    TLorentzVector lp4;
    TVector3 p4;
    lp4           = lp1 - lp3;
    p4            = lp4.Vect();

// VP+target CM system transfer
    double M5     = MK;
    double M6     = M_hyper + Ex[peak_num];
    TVector3 beta;
    beta = (p2 + p4) * (1./(E2 + E4));
    lp2.Boost(-beta);
    lp4.Boost(-beta);
    double theta5_cm = Theta_K;
    double phi5_cm   = Phi_K;
    double E_cm      = lp2.E() + lp4.E();
    double E5_cm     = 0.5 * (E_cm + (M5*M5 - M6*M6)/E_cm);
    double E6_cm     = 0.5 * (E_cm - (M5*M5 - M6*M6)/E_cm);
    double p_cm      = sqrt(E5_cm*E5_cm - M5*M5);
    TVector3 p5_cm, p6_cm;
    p5_cm.SetMagThetaPhi(p_cm,theta5_cm,phi5_cm);
    p6_cm = -p5_cm;
    TLorentzVector lp5(p5_cm,E5_cm);
    TLorentzVector lp6(p6_cm,E6_cm);

// Lab system transfer
    lp2.Boost(beta);
    lp4.Boost(beta);
    lp5.Boost(beta);
    lp6.Boost(beta);

// Kaon
    double P5     = lp5.P();
    double E5     = lp5.E();
    double Theta5 = lp5.Theta();
    double Phi5   = lp5.Phi();
    TVector3 p5;
    p5 = lp5.Vect();

// Blur
    double Z = gRandom->Uniform(0,thick);
    double y, yy, dE = 0.;
//    double dE_max = 3.0;
    double y_max = 0.1;
#if 1
    while(1){
      y = gRandom->Uniform(0,y_max);
      dE = gRandom->Uniform(0,1.);//0-1 MeV
      yy = dEdx(dE,Z,Me,P1);
      if(y<yy){ 
	  //if(y_max<yy){cout<<"y_max is too small"<<endl;}
		//cout<<Form("e- %04d  dE [MeV], yy:  %.3lf  %.1lf",n,dE,yy)<<endl;
		 break; }
    }
#endif
    E1    -= dE;
    E1    += 0.311;//energy correction
    E1    += gRandom->Gaus(0,E_beam*Ereso_beam);
    P1     = E2p(E1,M1);
    Theta1 = 0.;
    Phi1   = 0.;
    p1.SetMagThetaPhi(P1,Theta1,Phi1);
    lp1.SetVect(p1);
    lp1.SetE(E1);
	h1_dEe->Fill(yy);

#if 1
    while(1){
      y = gRandom->Uniform(0,y_max);
      //dE = gRandom->Uniform(0,(thick-Z)*5.);
      dE = gRandom->Uniform(0,1.);
      yy = dEdx(dE,thick-Z,Me,P3);
      if(y<yy){
		 //cout<<Form("ep %04d  dE, yy:  %.3lf  %.1lf",n,dE,yy)<<endl;
		 break; }
	  if(y_max<yy){cout<<"y_max is too small"<<endl;}
    }
#endif
    P3     += gRandom->Gaus(0,P_scat*Preso_scat);
    E3      = p2E(M3,P3);
    E3     += dE;
    E3     -= 0.332;//energy correction
    P3      = E2p(E3,M3);
    Theta3 += gRandom->Gaus(1E-3,Thetareso_scat);
    Phi3   += gRandom->Gaus(0,Phireso_scat);
    p3.SetMagThetaPhi(P3,Theta3,Phi3);
    lp3.SetVect(p3);
    lp3.SetE(E3);
	h1_dEep->Fill(dE);

#if 1
    while(1){
      y = gRandom->Uniform(0,y_max);
      //dE = gRandom->Uniform(0,(thick-Z)*5.);
      dE = gRandom->Uniform(0,1.);
      yy = dEdx(dE,thick-Z,MK,P5);
      if(y<yy){ 
	//	cout<<Form("K+ %04d  dE, yy:  %.3lf  %.1lf",n,dE,yy)<<endl;
		 break; }
	  if(y_max<yy){cout<<"y_max is too small"<<endl;}
    }
#endif
    P5     += gRandom->Gaus(0,P5*Preso_K);
    E5      = p2E(M5,P5);
    E5     += dE;
    E5     -= 0.187;//energy correction
    P5      = E2p(E5,M5);
    Theta5 += gRandom->Gaus(0,Thetareso_K);
    Phi5   += gRandom->Gaus(0,Phireso_K);
    p5.SetMagThetaPhi(P5,Theta5,Phi5);
    lp5.SetVect(p5);
    lp5.SetE(E5);
	h1_dEk->Fill(dE);

    TLorentzVector lp;
    lp = lp1 + lp2 - lp3 - lp5;

    double mass = lp.M() - M_hyper;
	mass += -16.3+16.6807;//more energy correction

    h1_mm_p->Fill(mass);

//    cout<<lp.E()<<"  "<<lp.P()<<"  "<<lp.M() - M_hyper<<endl;
  }

#if 1
  for(int n=0;n<QF;n++){
    double mass = 0;
    while(1){
      double x = gRandom->Uniform(0,100);
      double y = -1./2500. * (x - 50) * (x - 50) + 1;
      double yy = gRandom->Uniform(0,1);
      if(yy<y){ mass = x; break; }
    }
    h1_mm_q->Fill(mass);
  }
#endif
#if 1
  for(int n=0;n<Acc;n++){
    //double mass = gRandom->Uniform(-100,100);
	double pep0  = 2740.*(1.+0.045*(1.-2.*gRandom->Uniform()));//MeV
	double pk0   = 1200.*(1.+0.125*(1.-2.*gRandom->Uniform()));//MeV
	double PI    = 4.*atan(1.);
//cout<<"PI="<<PI<<endl;
	double Ee0   = 4240.;//MeV
	double te0   = 6.5*PI/180.;//rad
	double tk0   = 11.5*PI/180.;//rad
	double tek0  = te0 + tk0;//rad
	double pe0  = sqrt(Ee0*Ee0-Me*Me);
	double Eep0 = sqrt(pep0*pep0+Me*Me);
	double Ek0  = sqrt(pk0*pk0+MK*MK);
	double MT = 25.133*1000.;// 27Al
	double MH = 25.319*1000.+14.;// 26Mg+L
	double MM0= sqrt(pow((Ee0-Eep0+MT-Ek0),2.)-(pe0*pe0+pep0*pep0+pk0*pk0-2.*pe0*pep0*cos(te0)-2.*pe0*pk0*cos(tk0)+2.*pep0*pk0*cos(tek0)));
    h1_mm_a->Fill(MM0-MH);
    h1_mm_a2->Fill(MM0-MH);
  }
#endif
  h1_mm->Add(h1_mm_p);
  h1_mm->Add(h1_mm_q);
  h1_mm->Add(h1_mm_a);
  h1_mm2->Add(h1_mm_p);
  h1_mm2->Add(h1_mm_q);
  h1_mm2->Add(h1_mm_a);
  h1_mm3->Add(h1_mm_p);
  h1_mm3->Add(h1_mm_q);
  h1_mm3->Add(h1_mm_a);
//#if 1
//  for(int n=0;n<Acc*99;n++){
//    double mass = gRandom->Uniform(-100,100);
//    //h1_mm_a->Fill(mass);
//  }
//#endif

// TF1 *f_mm = new TF1("f_mm","gaus",-10,10);
// SetTF1(f_mm);
 //h1_mm->Fit(f_mm,"","",-10,10);
  //cout<<"RMS: "<<f_mm->GetParameter(2)<<"  "<<f_mm->GetParError(2)<<endl;
  //cout<<"FWHM: "<<f_mm->GetParameter(2) * 2.35<<endl;

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->Divide(1,1,1E-4,1E-4);
  c1->cd(1);
  //c1->cd(1)->DrawFrame(-20.,0.,20.,1E+4);
    //h1_mm->Draw("E");
    h1_mm->Draw("E");
    //h1_mm->GetXaxis()->SetRangeUser(-25,10);
    //h1_mm->GetXaxis()->SetRangeUser(-30,10);
    //h1_mm_p->Draw("same");
    //h1_mm_q->Draw("same");
    //h1_mm_a->Scale(1./100.);
    h1_mm_a->Draw("same");
  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  c2->Divide(1,2,1E-4,1E-4);
  c2->cd(1);
	TF1* f_acc  = new TF1("f_acc" ,"pol0",-20.,10.);
	TF1* f_qf  = new TF1("f_qf" ,"pol2",-20,10.);
    h1_mm_a2->Fit(f_acc);
    h1_mm_q->Fit(f_qf,"","",0.,10.);
  c2->cd(2);
    h1_mm2->Draw("E");
    h1_mm_a->Draw("same");
	TF1* fs  = new TF1("fs" ,"gausn+pol0(3)",-20.,10.);
	TF1* fp1 = new TF1("fp1","gausn+pol0(3)",-20.,10.);
	TF1* fp2 = new TF1("fp2","gausn+pol0(3)",-20.,10.);
	TF1* fp3 = new TF1("fp3","gausn+pol0(3)",-20.,10.);
	fs ->SetParameter(0,100);
	fs ->SetParameter(1,-16.3);
	fs ->SetParameter(2,0.81);
	fs ->FixParameter(3,f_acc->GetParameter(0));
	fp1->SetParameter(0,50);
	fp1->SetParameter(1,-6.1);
	fp1->SetParameter(2,1.0);
	fp1->FixParameter(3,f_acc->GetParameter(0));
	fp2->SetParameter(0,50);
	fp2->SetParameter(1,-3.1);
	fp2->SetParameter(2,1.0);
	fp2->FixParameter(3,f_acc->GetParameter(0));
	fp3->SetParameter(0,50);
	fp3->SetParameter(1,-0.8);
	fp3->SetParameter(2,1.0);
	fp3->FixParameter(3,f_acc->GetParameter(0));
	h1_mm2->Fit(fs,"","",-18.,-15.);
	h1_mm2->Fit(fp1,"","",-8.0,-5.5);
	h1_mm2->Fit(fp2,"","",-5.0,-2.5);
	h1_mm2->Fit(fp3,"","",-1.5,0.);
    h1_mm2->Draw("E");
    h1_mm_a->Draw("same");
	fs->Draw("same");
	fp1->Draw("same");
	fp2->Draw("same");
	fp3->Draw("same");
cout<<"Energy Loss"<<endl;
  TCanvas *c3 = new TCanvas("c3","c3",800,600);
  c3->Divide(2,2,1E-4,1E-4);
  c3->cd(1);
  h1_dEe->Draw("");
  h1_dEe->Fit("gausn");
  c3->cd(2);
  h1_dEep->Draw("");
  h1_dEep->Fit("gausn");
  c3->cd(3);
  h1_dEk->Draw("");
  h1_dEk->Fit("gausn");
  TCanvas *c4 = new TCanvas("c4","c4",800,600);
  c4->Divide(1,1,1E-4,1E-4);
	TF1* fs_result  = new TF1( "fs_result","gausn",-20.,10.);
	TF1* fp1_result = new TF1("fp1_result","gausn",-20.,10.);
	TF1* fp2_result = new TF1("fp2_result","gausn",-20.,10.);
	TF1* fp3_result = new TF1("fp3_result","gausn",-20.,10.);
	TF1* ftot_result= new TF1("ftot_result",fmm_total,-20.,10.,16);
	fs_result->SetNpx(2000);
	fp1_result->SetNpx(2000);
	fp2_result->SetNpx(2000);
	fp3_result->SetNpx(2000);
	ftot_result->SetNpx(2000);
	fs_result->SetLineColor(kAzure);
	fp1_result->SetLineColor(kAzure);
	fp2_result->SetLineColor(kAzure);
	fp3_result->SetLineColor(kAzure);
	fs_result->SetFillColor(kAzure);
	fp1_result->SetFillColor(kAzure);
	fp2_result->SetFillColor(kAzure);
	fp3_result->SetFillColor(kAzure);
	fs_result->SetFillStyle(3004);
	fp1_result->SetFillStyle(3004);
	fp2_result->SetFillStyle(3004);
	fp3_result->SetFillStyle(3004);
	ftot_result->SetLineColor(kRed);
	fs_result ->SetParameter(0,fs->GetParameter(0));
	fs_result ->SetParameter(1,fs->GetParameter(1));
	fs_result ->SetParameter(2,fs->GetParameter(2));
	fp1_result->SetParameter(0,fp1->GetParameter(0));
	fp1_result->SetParameter(1,fp1->GetParameter(1));
	fp1_result->SetParameter(2,fp1->GetParameter(2));
	fp2_result->SetParameter(0,fp2->GetParameter(0));
	fp2_result->SetParameter(1,fp2->GetParameter(1));
	fp2_result->SetParameter(2,fp2->GetParameter(2));
	fp3_result->SetParameter(0,fp3->GetParameter(0));
	fp3_result->SetParameter(1,fp3->GetParameter(1));
	fp3_result->SetParameter(2,fp3->GetParameter(2));
	ftot_result ->SetParameter(0,fs->GetParameter(0));
	ftot_result ->SetParameter(1,fs->GetParameter(1));
	ftot_result ->SetParameter(2,fs->GetParameter(2));
	ftot_result ->SetParameter(3,fp1->GetParameter(0));
	ftot_result ->SetParameter(4,fp1->GetParameter(1));
	ftot_result ->SetParameter(5,fp1->GetParameter(2));
	ftot_result ->SetParameter(6,fp2->GetParameter(0));
	ftot_result ->SetParameter(7,fp2->GetParameter(1));
	ftot_result ->SetParameter(8,fp2->GetParameter(2));
	ftot_result ->SetParameter(9,fp3->GetParameter(0));
	ftot_result ->SetParameter(10,fp3->GetParameter(1));
	ftot_result ->SetParameter(11,fp3->GetParameter(2));
	ftot_result ->SetParameter(12,f_acc->GetParameter(0));
	ftot_result ->SetParameter(13,f_qf->GetParameter(0));
	ftot_result ->SetParameter(14,f_qf->GetParameter(1));
	ftot_result ->SetParameter(15,f_qf->GetParameter(2));
  h1_mm3->Draw("E");
    h1_mm_a->Draw("same");
	fs_result->Draw("same");
	fp1_result->Draw("same");
	fp2_result->Draw("same");
	fp3_result->Draw("same");
	ftot_result->Draw("same");
  //exit(1);
//  double Ex[19] = { -16.3, -15.2, -13.5, -12.0, Expval1, Expval1*7.5/8.3, Expval1*5.5/8.3, Expval1*4.0/8.3, Expval1*2.8/8.3, Expval2, (Expval2-Expval1)+Expval1*7.5/8.3, (Expval2-Expval1)+Expval1*5.5/8.3, (Expval2-Expval1)+Expval1*4.0/8.3, (Expval2-Expval1)+Expval1*2.8/8.3, Expval3, (Expval3-Expval1)+Expval1*7.5/8.3, (Expval3-Expval1)+Expval1*5.5/8.3, (Expval3-Expval1)+Expval1*4.0/8.3, (Expval3-Expval1)+Expval1*2.8/8.3};
//  double CS[19] = {  79.0,  26.3,  13.2,  35.1, CSpval1, CSpval1*100/320, CSpval1*60./320, CSpval1*100/320, CSpval1*60./320, CSpval2, CSpval2*100/320, CSpval2*60./320, CSpval2*100/320, CSpval2*60./320, CSpval3, CSpval3*100/320, CSpval3*60./320, CSpval3*100/320, CSpval3*60./320};
  
  c1->Print("pdf/Mg27L_scaled.pdf");
  c2->Print("pdf/Mg27L_scaled_fitprocess.pdf");
  c4->Print("pdf/Mg27L_scaled_fit.pdf");
  cout<<"Print complete!"<<endl;

  theApp.Run();
  return 0;

}

