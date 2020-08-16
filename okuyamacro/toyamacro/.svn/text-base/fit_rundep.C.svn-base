#include "Settings.cc"
using namespace std;
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

//____________________________________________________________________________________________

void fit_rundep(){
gStyle->SetOptStat(0);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.06);
gStyle->SetPadLeftMargin(0.20);
gStyle->SetPadBottomMargin(0.15);
Settings *set = new Settings();
TString ifname1=Form("./root/all.root");//H3L
  TFile *ifp1  = new TFile(ifname1  );
      cout<<"input filename : "<<ifname1<<endl;
double run_beg = 10226;

TH2F *run_frame  = new TH2F("run_frame","run_frame",10,10225,10350,10,0,0.35);
         set->SetTH2(run_frame,"sigma","runnum","#sigma_{coin}[ns]");
         run_frame->GetYaxis()->SetTitleOffset(0.9);
         run_frame ->SetStats(0);

TH2F *run_m_frame  = new TH2F("run_m_frame","run_m_frame",10,10225,10350,10,-0.1,0.10);
         set->SetTH2(run_m_frame,"mean","runnum","#mean_{coin}[ns]");
         run_m_frame->GetYaxis()->SetTitleOffset(0.9);
         run_m_frame ->SetStats(0);


  TH2F *h2_cointime_runnum;
  h2_cointime_runnum = (TH2F*)ifp1->Get("cointime/h2_cointime_runnum");


double min =-1.;
double max = 1.;
double merge=5.;
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
        aa++;
      }

  }

  int n = run.size();
cout<<n<<endl;

  TGraphErrors *gr_run_sigma = new TGraphErrors(n,&run[0],&sigma_run[0],0,&e_sigma_run[0]);
  gr_run_sigma->SetMarkerStyle(22);
  gr_run_sigma->SetMarkerColor(2);
  gr_run_sigma->SetMarkerSize(1.5);

  TGraphErrors *gr_run_mean = new TGraphErrors(n,&run[0],&mean_run[0],0,&e_mean_run[0]);
  gr_run_mean->SetMarkerStyle(24);
  gr_run_mean->SetMarkerColor(2);
  gr_run_mean->SetMarkerSize(1.5);

TCanvas *c[6];
for(int i=0;i<6;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),2000/2,1500/2);}
c[0]->Clear();
gPad->SetLogz(1);
h2_cointime_runnum->Draw("colz");
c[1]->Clear();
   run_frame       ->Draw("");
gr_run_sigma->Draw("sameP");

c[2]->Clear();
   run_m_frame       ->Draw("");
gr_run_mean->Draw("sameP");

c[3]->Clear();
c[3]->Divide(5,4);
for(int i=0;i<20;i++){
c[3]->cd(i+1);h_coin_run[i]->Draw("");
}

#if 0
c[4]->Clear();
c[4]->Divide(5,5);
for(int i=25;i<50;i++){
c[4]->cd(i+1);h_coin_run[i]->Draw("");
}

c[5]->Clear();
c[5]->Divide(5,5);
for(int i=50;i<aa;i++){
c[5]->cd(i+1);h_coin_run[i]->Draw("");
}
#endif

string ofname_pdf="coin_rundep";
  c[0]->Print(Form("./pdf/fit_%s.pdf[",ofname_pdf.c_str())  );
  c[0]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  c[1]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  c[2]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  c[3]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  c[3]->Print(Form("./pdf/fit_%s.pdf]",ofname_pdf.c_str())  );




}
