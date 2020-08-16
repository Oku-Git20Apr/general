#include "Settings.cc"
using namespace std;
//____________________________________________________________________________________________
void FitGaus(TH1F *h, double &gamin, double &gamax, double range,int itr){
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

void fit_TDL_Tag(){
gStyle->SetOptStat(0);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.06);
gStyle->SetPadLeftMargin(0.16);
gStyle->SetPadBottomMargin(0.15);
TString ifname1=Form("./root/RKV_ana_acci/all.root");//H3L
  TFile *ifp1  = new TFile(ifname1  );
      cout<<"input filename : "<<ifname1<<endl;

TH1F *h_ctof_tdllr2_cut;
TH1F *h_ctof_tb_cut    ;//
Settings *set = new Settings();
h_ctof_tdllr2_cut = (TH1F*)ifp1->Get("h_ctof_tdllr2_cut");
h_ctof_tb_cut     = (TH1F*)ifp1->Get("h_ctof_tb_cut");
  set->SetTH1(h_ctof_tdllr2_cut    ,"TDLL2 - TDLR2"           ,"TDLL2 - TDLR2[ns]"      ,"Counts/25ps"       , 1, 3001, 7);
  set->SetTH1(h_ctof_tb_cut        ,"TagB23 - TagB24"         ,"TagB23 - TagB24[ns]"    ,"Counts/25ps"       , 1, 3001, 5);
h_ctof_tdllr2_cut->SetStats(0);
h_ctof_tb_cut    ->SetStats(0);
double min =-1.;
double max = 1.;
double sigma_tdl, sigma_tag;
TF1 *f_tdl, *f_tag;
f_tdl = new TF1("f_tdl","gaus",min,max);
f_tag = new TF1("f_tag","gaus",min,max);
set->SetTF1(f_tdl, 2, 1, 2.5);
set->SetTF1(f_tag, 2, 1, 2.5);
FitGaus(h_ctof_tdllr2_cut, min, max, 2.0, 5);
h_ctof_tdllr2_cut->Fit(f_tdl,"QR","",min,max);
FitGaus(h_ctof_tb_cut, min, max, 3.0, 5);
h_ctof_tb_cut    ->Fit(f_tag,"QR","",min,max);

sigma_tdl=f_tdl->GetParameter(2);
sigma_tag=f_tag->GetParameter(2);

cout<<"sigma TDL = "<<sigma_tdl*1000.<<"ps"<<endl;
cout<<"sigma Tag = "<<sigma_tag*1000.<<"ps"<<endl;

int p_flag=0;
cout<<"Do you want to save canvas?"<<endl;
cout<<"yes :1, no :0"<<endl;
cin>>p_flag;

TCanvas *c[2];
for(int i=0;i<2;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800/2,1500/2);}
c[0]->Clear();
c[0]->cd(1);
h_ctof_tdllr2_cut->Draw("");

c[1]->Clear();
c[1]->cd(1);
h_ctof_tb_cut->Draw("");

if(p_flag==1){
string ofname_pdf="TDL_Tag_ToF";
  c[0]->Print(Form("./pdf/fit_%s.pdf[",ofname_pdf.c_str())  );
  c[0]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  c[1]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  c[1]->Print(Form("./pdf/fit_%s.pdf]",ofname_pdf.c_str())  );
 }

} 
