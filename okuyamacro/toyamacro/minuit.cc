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

#include <TMinuit.h>
#include <TFitter.h>

#include "MathTools.h"
static int    lntrk_nHit;                 //   Number of hit on CDC/VDC in a cluster group
double lntrk_distance_to_wire[100];       //   Distance of trajectory to the wire
double lntrk_drift_length[200];           //   Drift length on the wire
double lntrk_position_sigma[200];         //   Position Resolution on the wire (defined layer by layer)
double lntrk_wire_vec[200][6];            //   The wire vector to be used in minuitFunction
double lntrk_par[6];            //   
//________________________________________________________________________________________________________________________
double calcDistanceFromWireToTrajectory( int ihit )
{
  // a point on line u and w
  //  x = u[0] + s u[3]             x = w[0] + t w[3]
  //  y = u[1] + s u[4]   and,      y = w[1] + t w[4]
  //  z = u[2] + s u[5]             z = w[2] + t w[5]

  // linear simultaneuous equation 
  //    _      _   _   _     _   _
  //   |  a  b  | |  s  |   |  e  |
  //   |        | |     | + |     |  = 0
  //   |_ c  d _| |_ t _|   |_ f _|
  //   
  //   
  if ( ihit < 0 || lntrk_nHit <= ihit ) {
    return  0.;
  }

  double dist = 1.0e+32;
  double wire_vec[6];
  for ( int i=0; i<6; i++ ) {
    wire_vec[i] = lntrk_wire_vec[ihit][i];
  }

  double dptl = MathTools::DistanceBetweenTwoLines( wire_vec, lntrk_par, false, false);
  if ( dptl < dist ) { 
    dist = dptl;
  }

  return dist;
}
//________________________________________________________________________________________________________________________
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg)
{


  double chi_square = 0;
  double delta;

  // set parameters in global variables
  lntrk_par[0]  =  par[4];//x
  lntrk_par[1]  =  par[5];//y
  lntrk_par[2]  =  par[2];//z
  lntrk_par[3]  =  par[3]*TMath::Sin(TMath::DegToRad()*par[0])*TMath::Cos(TMath::DegToRad()*par[1]);//x'
  lntrk_par[4]  =  par[3]*TMath::Sin(TMath::DegToRad()*par[0])*TMath::Sin(TMath::DegToRad()*par[1]);//y'
  lntrk_par[5]  =  par[3]*TMath::Cos(TMath::DegToRad()*par[0])                                     ;//z'
  //lntrk_par[0]  =  par[0];//x
  //lntrk_par[1]  =  par[1];//y
  //lntrk_par[2]  =  par[2];//z
  //lntrk_par[3]  =  par[3];
  //lntrk_par[4]  =  par[4];
  //lntrk_par[5]  =  par[5];

  for ( int ihit=0; ihit<lntrk_nHit; ihit++ ) {  // loop for each chamber hit

    lntrk_distance_to_wire[ihit] = calcDistanceFromWireToTrajectory( ihit );

    delta       = (lntrk_distance_to_wire[ihit] - lntrk_drift_length[ihit] ) / lntrk_position_sigma[ihit] ;
    chi_square += delta*delta;

  }

  result = chi_square;

}
//________________________________________________________________________________________________________________________
void SetVectorComponent(double x,double y,double z, double x1, double y1, double z1, double (&vec)[6]){
  vec[0]=x;
  vec[1]=y;
  vec[2]=z;
  vec[3]=x1 - x;
  vec[4]=y1 - y;
  vec[5]=z1 - z;
}


//////////////////////////// main //////////////////////////////////////////
int main(int argc, char** argv){
  TApplication *theApp = new TApplication("App", &argc, argv);
  gRandom -> SetSeed( time(NULL) ); //seed file set by time

  TH1F *h_theta = new TH1F("h_theta","#theta" ,1000,-200,200);
  TH1F *h_phi   = new TH1F("h_phi"  ,"#phi"   ,1000,-200,200);
  TH1F *h_offs  = new TH1F("h_offs" ,"offset" ,1000,-200,200);
  TCanvas *c1 = new TCanvas("c1","c1",800,800);

  //LNAnalysis::SearchTracks (Clear track) lntrk_nHit 0 nVDCCluster 4
  //LNAnalysis::SearchTracks VDCCluster 0 / 4
  //LNAnalysis::SearchTracks Elements 0 / 1
  lntrk_wire_vec[0][0]= 2.29505;
  lntrk_wire_vec[0][1]= -10.743;
  lntrk_wire_vec[0][2]= 0      ;
  lntrk_wire_vec[0][3]= 7.2076 ;
  lntrk_wire_vec[0][4]= 1.53978;
  lntrk_wire_vec[0][5]= 40     ;
  //LNAnalysis::SearchTracks VDCCluster 1 / 4
  //LNAnalysis::SearchTracks Elements 0 / 3
  lntrk_wire_vec[1][0]= -0.157108;
  lntrk_wire_vec[1][1]= -5.88203 ;
  lntrk_wire_vec[1][2]= 0        ;
  lntrk_wire_vec[1][3]= 4.44979  ;
  lntrk_wire_vec[1][4]= -0.118853;
  lntrk_wire_vec[1][5]= 40       ;
  //LNAnalysis::SearchTracks Elements 1 / 3
  lntrk_wire_vec[2][0]= 0.469035;
  lntrk_wire_vec[2][1]= -5.86541;
  lntrk_wire_vec[2][2]= 0       ;
  lntrk_wire_vec[2][3]= 4.43721 ;
  lntrk_wire_vec[2][4]= 0.354827;
  lntrk_wire_vec[2][5]= 40      ;
  //LNAnalysis::SearchTracks Elements 2 / 3
  lntrk_wire_vec[3][0]= 0.586238;
  lntrk_wire_vec[3][1]= -7.33107;
  lntrk_wire_vec[3][2]= 0       ;
  lntrk_wire_vec[3][3]= 4.91695 ;
  lntrk_wire_vec[3][4]= 0.39319 ;
  lntrk_wire_vec[3][5]= 40      ;
  //LNAnalysis::SearchTracks VDCCluster 2 / 4
  //LNAnalysis::SearchTracks Elements 0 / 2
  lntrk_wire_vec[4][0]= 1.74876 ;
  lntrk_wire_vec[4][1]= -7.88816;
  lntrk_wire_vec[4][2]= 0       ;
  lntrk_wire_vec[4][3]= 5.29737 ;
  lntrk_wire_vec[4][4]= 1.1744  ;
  lntrk_wire_vec[4][5]= 40      ;
  //LNAnalysis::SearchTracks Elements 1 / 2
  lntrk_wire_vec[5][0]= 1.90618 ;
  lntrk_wire_vec[5][1]= -8.59824;
  lntrk_wire_vec[5][2]= 0       ;
  lntrk_wire_vec[5][3]= 5.76667 ;
  lntrk_wire_vec[5][4]= 1.27844 ;
  lntrk_wire_vec[5][5]= 40      ;
  //LNAnalysis::SearchTracks VDCCluster 3 / 4
  //LNAnalysis::SearchTracks Elements 0 / 3
  lntrk_wire_vec[6][0]= 0.880499;
  lntrk_wire_vec[6][1]= -9.49149;
  lntrk_wire_vec[6][2]= 0       ;
  lntrk_wire_vec[6][3]= 6.37244 ;
  lntrk_wire_vec[6][4]= 0.591153;
  lntrk_wire_vec[6][5]= 40      ;
  //LNAnalysis::SearchTracks Elements 1 / 3
  lntrk_wire_vec[7][0]= 1.57902  ;
  lntrk_wire_vec[7][1]= -9.40055 ;
  lntrk_wire_vec[7][2]= 0        ;
  lntrk_wire_vec[7][3]= 6.31139  ;
  lntrk_wire_vec[7][4]= 1.06013  ;
  lntrk_wire_vec[7][5]= 40       ;
  //LNAnalysis::SearchTracks Elements 2 / 3
  lntrk_wire_vec[8][0]= 1.6993  ;
  lntrk_wire_vec[8][1]= -10.1166;
  lntrk_wire_vec[8][2]= 0       ;
  lntrk_wire_vec[8][3]= 6.79248 ;
  lntrk_wire_vec[8][4]= 1.14094 ;
  lntrk_wire_vec[8][5]= 40      ;
  //LNAnalysis::SearchTracks lntrk_nHit 9
  //DecayDisplay::DrawEventEvent Number = 417
  //Num. of track candidate : 1
  //init mean phi-81.1949 theta90
  //mean phi(after fit):-75.2106, theta:146.813, offs:7.16395e-322, chi2/ndf:632.867/6



  for(int n=0;n<1;n++){
    ///////////////////////
    //Hit information set//
    ///////////////////////
    lntrk_nHit = 9;                 //   Number of hit on CDC/VDC in a cluster group
    for(int i=0;i<lntrk_nHit;i++){
      lntrk_drift_length[i]     = 0.5;  //   Drift length on the wire
      lntrk_position_sigma[i]   = 0.05;  //   Position Resolution on the wire (defined layer by layer)
      //double w_vec[6];
      //double x  = 10. + gRandom->Gaus(0.,0.1);
      //double y  = 10. + gRandom->Gaus(0.,0.1);
      //double z  = 0.;
      //double x1 = 10. + gRandom->Gaus(0.,5.);
      //double y1 = 10. + gRandom->Gaus(0.,5.);
      //double z1 = 10.;
      //SetVectorComponent(x,y,z,x1,y1,z1,w_vec);
      //for(int k=0;k<6;k++){
      //  lntrk_wire_vec[i][k]    = w_vec[k];            //   The wire vector to be used in minuitFunction
      //}
    }


    //TMinuit *min = new TMinuit();
    TFitter *minimizer = new TFitter();
    //TFitter *min = new TFitter();
    //int p_level = min->SetPrintLevel(-1);
    //delete min;

    Double_t print_level[5];


    //bool                 is_verbose = false;                          //  if true, you will see so many message about fitting from TMinuit which is superclass of TFitter
    bool                 is_verbose = true;                          //  if true, you will see so many message about fitting from TMinuit which is superclass of TFitter
    print_level[0] = -1;
    if ( is_verbose == true ) {
      print_level[0] = 1;
    }
    minimizer -> ExecuteCommand("SET PRIntout", print_level, 1);            // Minuit command: set printout
    minimizer -> ExecuteCommand("SET NOWarnings", print_level, 0);            // Minuit command: set printout

    minimizer -> SetFCN(minuitFunction);


    double               theta_init=90., phi_init=-80., offs_init=0.;    // initial values for fitting
    // set initial values
    //           Define the parameters
    //                         arg1 - parameter number
    //                         |   arg2 - parameter name
    //                         |   |          arg3 - first guess at parameter value
    //                         |   |          |             arg4 - estimated distance to minimum
    //                         |   |          |             |      arg5 - lower limt
    //                         |   |          |             |      |         arg6 - higher limit
    //                         |   |          |             |      |         |
    minimizer -> SetParameter( 0, "theta"  , theta_init  , 0.01,  0.      , 180.       );//[degree]
    minimizer -> SetParameter( 1, "phi"    , phi_init    , 0.01,  -180.   , 180.       );//[degree]
    minimizer -> SetParameter( 2, "offs"   , offs_init   , 0.01,  0.      , 0.       );  //beam direction offset
    minimizer -> SetParameter( 3, "r"      , 1.          ,    0,  0.      , 0.       );
    minimizer -> SetParameter( 4, "x"      , 0.          , 0.01, -1.00    , 1.00     );   //(0.5 * VDC hight) / (CDC inner radius) = 19/20 = ~1
    minimizer -> SetParameter( 5, "y"      , 0.          , 0.01, -1.00    , 1.00     );
    //minimizer -> SetParameter( 0, "x"       , 0.          , 0.01,  0.      , 0.       );//[cm]
    //minimizer -> SetParameter( 1, "y"       , 0.          , 0.01,  0.      , 0.       );//[cm]
    //minimizer -> SetParameter( 2, "z"       , offs_init   , 0.01,  0.      , 0.       );//beam direction offset
    //minimizer -> SetParameter( 3, "x'"      , 1.          ,    0,  0.      , 0.       );
    //minimizer -> SetParameter( 4, "y'"      , 0.          , 0.01, -1.00    , 1.00     );   //(0.5 * VDC hight) / (CDC inner radius) = 19/20 = ~1
    //minimizer -> SetParameter( 5, "z'"      , 0.          , 0.01, -1.00    , 1.00     );

    //minimizer -> FixParameter(3);       // r is fixed
    //minimizer -> FixParameter(4);       // x is fixed
    //minimizer -> FixParameter(5);       // y is fixed
    //minimizer -> FixParameter(0);       // theta is fixed
    //minimizer -> FixParameter(1);       // phi is fixed
    //minimizer -> FixParameter(2);       // offset is fixed

    // minimization by "MIGRAD" method

    minimizer -> ExecuteCommand("SET PRI", print_level, 1);            // Minuit command: SET PRIntout, 1 means NO print of mssage
    minimizer -> ExecuteCommand("MIGRAD", 0, 0);

    // fit results
    double theta = minimizer -> GetParameter(0);
    double phi   = minimizer -> GetParameter(1);
    double offs  = minimizer -> GetParameter(2);
    //double theta = TMath::RadToDeg()*TMath::ACos(minimizer -> GetParameter(5));//cos(theta) = z'/r
    //double phi   = TMath::RadToDeg()*TMath::ATan(minimizer -> GetParameter(3)/minimizer -> GetParameter(4));//tan(phi) = x'/y' 
    //double offs  = minimizer -> GetParameter(2);

    cout<<Form("theta:%.2lf, phi:%.2lf, offs(z):%.2lf",theta,phi,offs)<<endl;
    h_theta->Fill(theta);
    h_phi  ->Fill(phi);
    h_offs ->Fill(offs);

    delete minimizer;
  }

  c1->Divide(1,3);
  c1->cd(1);h_theta->Draw();
  c1->cd(2);h_phi  ->Draw();
  c1->cd(3);h_offs ->Draw();
 // gSystem->Exit(0);
  theApp->Run();
  return 0;
}
