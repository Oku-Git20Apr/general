#include <iostream>
#include <sstream>
#include <iomanip>
#include <csignal>
#include <stdlib.h>
#include <climits>
#include <fstream>
#include <math.h>
#include <time.h>
#include <string>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <unistd.h>
#include <map>
using namespace std;

#include "TApplication.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
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
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h" //added 2016/07/15
#include "TPDF.h"
#include "TGaxis.h"

#include "TDatime.h"

#include "Settings.h"

void SetGr(TGraph *gr, TString hname, TString xname, TString yname, int LColor, int MColor, int MStyle, double Yoffset, double min, double max){
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
//____________________________________________________________________________________________
void FitGaus(TH1 *h, double &gamin, double &gamax, double range,int itr){
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


////////////////////////////////////////////////////////////////
/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */

int main(int argc, char** argv){
gStyle->SetOptStat(0);
gStyle->SetPadGridX(0);
gStyle->SetPadGridY(0);
gStyle->SetPadTopMargin(0.051);
gStyle->SetPadRightMargin(0.08);
gStyle->SetPadLeftMargin(0.10);
gStyle->SetPadBottomMargin(0.17);
Settings *set = new Settings();

string ifname = "../../data/date.txt";
  int s_time[3],e_time[3];//d,h,m
  s_time[0]=24;//d
  s_time[1]=13;//h
  s_time[2]=29;//m
  e_time[0]=29;//d
  e_time[1]=9;//h
  e_time[2]=0;//m

  TApplication theApp("App", &argc, argv);
  
  string line;
  string::size_type comment_start = 0;
  vector<double> timef, time,runnum;
  //vector<int> day,hour,min, sec, msec;
  double Timef, Time;
  int RunNum, Year, Mon, Day, Hour, Min, Sec, Msec;

  std::map<int, double> run2time;

  //Plot time region
  double UtimeOffset = 788918400.;
  double start,end;
  TDatime Ts(2017,4,s_time[0],s_time[1],s_time[0],0);
  TDatime Te(2017,4,e_time[0],e_time[1],e_time[0],0);
  start = Ts.Convert() - UtimeOffset;
  end   = Te.Convert() - UtimeOffset;

  TGraph *g_runnum;

  ifstream ifp( ifname.c_str() );
  if ( ifp.fail() ) {
    std::cout << "file open fail : " << ifname << std::endl;
    return 0;
  }
  while(!ifp.eof()){
    getline(ifp,line);
   if(line[0]=='#') continue; //skip if there is "#" in the 1st line.
  istringstream sline(line);
  sline >> RunNum;
  sline >> Year;
  sline >> Mon;
  sline >> Day;
  sline >> Hour;
  sline >> Min;

  TDatime T0(Year,Mon,Day,Hour,Min,0);
  Time = T0.Convert() - UtimeOffset;
   //Time  = 31.*24.*60.*60.*(Mon-13) + 24.*60.*60.*(Day-2) + 60.*60.*(Hour+6) + 60.*Min;

   //Fill
   if(Time >= start && Time <= end){
     time.push_back(Time); //timef.push_back(Timef);
     runnum.push_back(RunNum); //timef.push_back(Timef);
     run2time.insert(std::make_pair(RunNum,Time));
   }

 }
 ifp.close();
 
double runmin      =   10220;
double runmax      =   10360;
double ther_min     =   20.;
double ther_max     =   30.;
TH2F *h_frame = new TH2F("h_frame","h_frame",10 ,  start,  end,10,runmin,runmax );
set->SetTH2(h_frame , ""    , "", "runnum");
         h_frame->GetYaxis()->SetTitleOffset(0.9);
         h_frame->GetXaxis()->SetTitleOffset(0.9);
 g_runnum= new TGraph((int)runnum.size(),&time[0],&runnum[0]);
 SetGr(g_runnum, "Time vs. RunNum","time","RunNum",1,2,20,0.9,0,0);


/////////++++++++++++/////////////////
TString ifname1=Form("./root/all.root");//H3L
  TFile *ifp1  = new TFile(ifname1  );
      cout<<"input filename : "<<ifname1<<endl;
double run_beg = 10226;

TH2F *run_sigma_frame  = new TH2F("run_sigma_frame","run_sigma_frame",10,10225,10350,10,0,0.35);
         set->SetTH2(run_sigma_frame,"","Run Number","#sigma_{coin}[ns]");
         run_sigma_frame->GetYaxis()->SetTitleOffset(0.9);
         run_sigma_frame ->SetStats(0);

TH2F *time_sigma_frame  = new TH2F("time_sigma_frame","time_sigma_frame",10,start,end,10,0,0.35);
         set->SetTH2(time_sigma_frame,"","Date and Time","#sigma_{coin}[ns]");
         time_sigma_frame->GetYaxis()->SetTitleOffset(0.8);
         time_sigma_frame ->SetStats(0);

TH2F *time_mean_frame  = new TH2F("time_mean_frame","time_mean_frame",10,start,end,10,-0.1,0.1);
         set->SetTH2(time_mean_frame,"","Date and Time","mean[ns]");
         time_mean_frame->GetYaxis()->SetTitleOffset(0.8);
         time_mean_frame ->SetStats(0);


TH2F *run_m_frame  = new TH2F("run_m_frame","run_m_frame",10,10225,10350,10,-0.1,0.10);
         set->SetTH2(run_m_frame,"","Run Number","mean[ns]");
         run_m_frame->GetYaxis()->SetTitleOffset(0.9);
         run_m_frame ->SetStats(0);


  TH2F *h2_cointime_runnum;
  h2_cointime_runnum = (TH2F*)ifp1->Get("cointime/h2_cointime_runnum");

double min =-1.;
double max = 1.;
double merge=5.;
vector<double> runtime;
vector<double> sigma_run  ;
vector<double> e_sigma_run;
vector<double> mean_run, e_mean_run;
vector<double> run;
double par[3];
double e_par[3];

TF1 *f_run;
f_run = new TF1("f_run","gaus",min,max);
set->SetTF1(f_run, 2, 1, 2.5);
int aa=0;
TH1D *h_coin_run[200];
  for(run_beg;run_beg<10351;run_beg+=merge){
      TH1D *hpy;
      int bin = h2_cointime_runnum->GetYaxis()->FindBin(run_beg);
      cout<<"run:"<<run_beg<<"  bin:"<<bin<<endl;
      hpy=h2_cointime_runnum->ProjectionX(Form( "%d",run_beg),bin , bin+(int)merge-1);
      if(hpy->Integral()>10){
        h_coin_run[aa]=(TH1D*)hpy->Clone();
        FitGaus(h_coin_run[aa], min, max, 3.0, 5);
        h_coin_run[aa]->Fit(f_run,"QR","",min,max);
        mean_run.push_back(f_run->GetParameter(1));
        sigma_run.push_back(f_run->GetParameter(2));
        e_mean_run.push_back(f_run->GetParError(1));
        e_sigma_run.push_back(f_run->GetParError(2));
        run.push_back(run_beg+merge/2.);

        runtime.push_back(run2time.at((int)run_beg));
        aa++;
      }

  }

  int n = run.size();
cout<<n<<"  "<<aa<<endl;

  TGraphErrors *gr_run_sigma = new TGraphErrors(n,&run[0],&sigma_run[0],0,&e_sigma_run[0]);
  gr_run_sigma->SetMarkerStyle(22);
  gr_run_sigma->SetMarkerColor(2);
  gr_run_sigma->SetMarkerSize(1.5);

  TGraphErrors *gr_time_sigma = new TGraphErrors(n,&runtime[0],&sigma_run[0],0,&e_sigma_run[0]);
  gr_time_sigma->SetMarkerStyle(22);
  gr_time_sigma->SetMarkerColor(2);
  gr_time_sigma->SetMarkerSize(1.5);

  TGraphErrors *gr_run_mean = new TGraphErrors(n,&run[0],&mean_run[0],0,&e_mean_run[0]);
  gr_run_mean->SetMarkerStyle(24);
  gr_run_mean->SetMarkerColor(2);
  gr_run_mean->SetMarkerSize(1.5);

  TGraphErrors *gr_time_mean = new TGraphErrors(n,&runtime[0],&mean_run[0],0,&e_mean_run[0]);
  gr_time_mean->SetMarkerStyle(24);
  gr_time_mean->SetMarkerColor(2);
  gr_time_mean->SetMarkerSize(1.5);



 //Draw
 TCanvas *c[5];
for(int i=0;i<5;i++){ c[i] = new TCanvas(Form("c%d",i+1),Form("TimeDependence%d",i+1),1200,650);}
 c[0]->cd(1);//->SetMargin(0.1,0.1,0.1,0.07);
 c[0]->Clear();
 //g_runnum->GetXaxis()->SetTimeDisplay(1);
 //g_runnum->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M:%S}");
 //g_runnum->GetXaxis()->SetLabelOffset(0.02);
 //g_runnum->GetXaxis()->SetTitleOffset(1.1);
 //g_runnum->Draw("ALP");
 h_frame->GetXaxis()->SetTimeDisplay(1);
 h_frame->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");
 h_frame->GetXaxis()->SetLabelOffset(0.02);
 h_frame->GetXaxis()->SetTitleOffset(1.1);
 h_frame->Draw("");
 g_runnum->Draw("LP");


 c[1]->Clear();
 //gr_time_sigma->GetXaxis()->SetTimeDisplay(1);
 //gr_time_sigma->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M:%S}");
 //gr_time_sigma->GetXaxis()->SetLabelOffset(0.02);
 //gr_time_sigma->GetXaxis()->SetTitleOffset(1.1);
 time_sigma_frame->GetXaxis()->SetTimeDisplay(1);
 time_sigma_frame->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");
 time_sigma_frame->GetXaxis()->SetLabelOffset(0.02);
 time_sigma_frame->GetXaxis()->SetTitleOffset(1.1);
 time_sigma_frame->Draw("");
 gr_time_sigma->Draw("P");


 c[2]->Clear();
 run_sigma_frame       ->Draw("");
 gr_run_sigma->Draw("LP");
 
 c[3]->Clear();
 run_m_frame       ->Draw("");
 gr_run_mean->Draw("LP");

 c[4]->Clear();
 time_mean_frame->GetXaxis()->SetTimeDisplay(1);
 time_mean_frame->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");
 time_mean_frame->GetXaxis()->SetLabelOffset(0.02);
 time_mean_frame->GetXaxis()->SetTitleOffset(1.1);
 time_mean_frame       ->Draw("");
 gr_time_mean->Draw("LP");
 
// c[4]->Clear();
// gPad->SetLogz(1);
// h2_cointime_runnum->Draw("colz");

// c[3]->Clear();
// c[3]->Divide(5,4);
// for(int i=0;i<20;i++){
// c[3]->cd(i+1);h_coin_run[i]->Draw("");
// }
string ofname_pdf="coin_rundep";
  c[0]->Print(Form("./pdf/%s.pdf[",ofname_pdf.c_str())  );
  c[0]->Print(Form("./pdf/%s.pdf" ,ofname_pdf.c_str())  );
  c[1]->Print(Form("./pdf/%s.pdf" ,ofname_pdf.c_str())  );
  c[2]->Print(Form("./pdf/%s.pdf" ,ofname_pdf.c_str())  );
  c[3]->Print(Form("./pdf/%s.pdf" ,ofname_pdf.c_str())  );
  c[4]->Print(Form("./pdf/%s.pdf" ,ofname_pdf.c_str())  );
  c[3]->Print(Form("./pdf/%s.pdf]",ofname_pdf.c_str())  );
 
 theApp.Run();
  return 0;
}

