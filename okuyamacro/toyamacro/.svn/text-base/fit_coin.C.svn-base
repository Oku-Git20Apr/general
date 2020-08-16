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

void fit_coin(){
gStyle->SetOptStat(0);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.06);
gStyle->SetPadLeftMargin(0.20);
gStyle->SetPadBottomMargin(0.15);
Settings *set = new Settings();
TString ifname1=Form("./root/RKV_ana_acci/true.root");
  TFile *ifp1  = new TFile(ifname1  );
      cout<<"input filename : "<<ifname1<<endl;

TH1F *h_cointime_tdll[10];
TH1F *h_cointime_tdlr[10];//
TH1F *h_cointime_tb[40];//


TH2F *tdll_frame  = new TH2F("tdll_frame","tdll_frame",10,1.5,9.5,10,0,0.35);
         set->SetTH2(tdll_frame,"sigma","TDLL seg","#sigma_{coin}[ns]");
         tdll_frame->GetYaxis()->SetTitleOffset(0.9);
         tdll_frame ->SetStats(0);

TH2F *tdlr_frame  = new TH2F("tdlr_frame","tdlr_frame",10,1.5,9.5,10,0,0.35);
         set->SetTH2(tdlr_frame,"sigma","TDLR seg","#sigma_{coin}[ns]");
         tdlr_frame->GetYaxis()->SetTitleOffset(0.9);
         tdlr_frame ->SetStats(0);

TH2F *tagb_frame  = new TH2F("tagb_frame","tagb_frame",10,0,26,10,0,0.35);
         set->SetTH2(tagb_frame,"sigma","TagB seg","#sigma_{coin}[ns]");
         tagb_frame->GetYaxis()->SetTitleOffset(0.9);
         tagb_frame ->SetStats(0);

for(int i=0;i<10;i++){
  h_cointime_tdll[i] = (TH1F*)ifp1->Get(Form("cointime/h_cointime_tdll%d",i+1));
  h_cointime_tdlr[i] = (TH1F*)ifp1->Get(Form("cointime/h_cointime_tdlr%d",i+1));
  //set->SetTH1(h_ctof_tdllr2_cut    ,"TDLL2 - TDLR2"           ,"TDLL2 - TDLR2[ns]"      ,"Counts/25ps"       , 1, 3001, 7);
  //set->SetTH1(h_ctof_tb_cut        ,"TagB23 - TagB24"         ,"TagB23 - TagB24[ns]"    ,"Counts/25ps"       , 1, 3001, 5);
  h_cointime_tdll[i] ->SetStats(0);
  h_cointime_tdlr[i] ->SetStats(0);
  h_cointime_tdll[i] ->GetXaxis()->SetRangeUser(-2,2);
  h_cointime_tdlr[i] ->GetXaxis()->SetRangeUser(-2,2);
}
for(int i=0;i<24;i++){
  h_cointime_tb[i] = (TH1F*)ifp1->Get(Form("cointime/h_cointime_tb%d",i+1));
  h_cointime_tb[i] ->SetStats(0);
  h_cointime_tb[i] ->GetXaxis()->SetRangeUser(-2,2);
}

double min =-1.;
double max = 1.;
double sigma_tdll[10], sigma_tdlr[10],sigma_tag[24];
double e_sigma_tdll[10], e_sigma_tdlr[10],e_sigma_tag[24];
double mean_tdll[10], mean_tdlr[10],mean_tag[24];
double seg[40];
sigma_tdll[0]=sigma_tdlr[0]=-1.;


TF1 *f_tdll, *f_tdlr, *f_tag;
f_tdll = new TF1("f_tdll","gaus",min,max);
f_tdlr = new TF1("f_tdlr","gaus",min,max);
f_tag = new TF1("f_tag","gaus",min,max);
set->SetTF1(f_tdll, 2, 1, 2.5);
set->SetTF1(f_tdlr, 2, 1, 2.5);
set->SetTF1(f_tag , 2, 1, 2.5);


  for(int i=1;i<9;i++){
    seg[i]=i+1;
  
    FitGaus(h_cointime_tdll[i], min, max, 3.0, 5);
    h_cointime_tdll[i]->Fit(f_tdll,"QR","",min,max);
    sigma_tdll[i]  =f_tdll->GetParameter(2);
    e_sigma_tdll[i]=f_tdll->GetParError(2);
  
    FitGaus(h_cointime_tdlr[i], min, max, 3.0, 5);
    h_cointime_tdlr[i]->Fit(f_tdlr,"QR","",min,max);
    sigma_tdlr[i]  =f_tdlr->GetParameter(2);
    e_sigma_tdlr[i]=f_tdlr->GetParError(2);
  }

  for(int i=0;i<24;i++){
    seg[i]=i+1;
    FitGaus(h_cointime_tb[i], min, max, 3.0, 5);
    h_cointime_tb[i]->Fit(f_tag,"QR","",min,max);
    sigma_tag[i]  =f_tag->GetParameter(2);
    e_sigma_tag[i]=f_tag->GetParError(2);
  }

  TGraphErrors *gr_tdll_sigma = new TGraphErrors(10,seg,sigma_tdll,0,e_sigma_tdll);
  gr_tdll_sigma->SetMarkerStyle(22);
  gr_tdll_sigma->SetMarkerColor(2);
  gr_tdll_sigma->SetMarkerSize(1.5);

  TGraphErrors *gr_tdlr_sigma = new TGraphErrors(10,seg,sigma_tdlr,0,e_sigma_tdlr);
  gr_tdlr_sigma->SetMarkerStyle(23);
  gr_tdlr_sigma->SetMarkerColor(4);
  gr_tdlr_sigma->SetMarkerSize(1.5);

  TGraphErrors *gr_tagb_sigma = new TGraphErrors(24,seg,sigma_tag,0,e_sigma_tag);
  gr_tagb_sigma->SetMarkerStyle(21);
  gr_tagb_sigma->SetMarkerColor(6);
  gr_tagb_sigma->SetMarkerSize(1.5);


int p_flag=0;
cout<<"Do you want to save canvas?"<<endl;
cout<<"yes :1, no :0"<<endl;
cin>>p_flag;

TCanvas *c[6];
for(int i=0;i<6;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),2000/2,1500/2);}

c[0]->Clear();
c[0]->Divide(3,4);
for(int i=0;i<10;i++){
  c[0]->cd(i+1);h_cointime_tdll[i]->Draw();
}

c[1]->Clear();
c[1]->cd(1);
   tdll_frame       ->Draw("");
gr_tdll_sigma->Draw("sameP");

c[2]->Clear();
c[2]->Divide(4,3);
for(int i=0;i<10;i++){
  c[2]->cd(i+1);h_cointime_tdlr[i]->Draw();
}

c[3]->Clear();
c[3]->cd(1);
   tdlr_frame       ->Draw("");
gr_tdlr_sigma->Draw("sameP");

c[4]->Clear();
c[4]->Divide(6,4);
for(int i=0;i<24;i++){
  c[4]->cd(i+1);h_cointime_tb[i]->Draw();
}

c[5]->Clear();
c[5]->cd(1);
   tagb_frame       ->Draw("");
   gr_tagb_sigma->Draw("sameP");

if(p_flag==1){
string ofname_pdf="TDL_coin";
  c[0]->Print(Form("./pdf/fit_%s.pdf[",ofname_pdf.c_str())  );
  c[0]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  c[1]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  c[2]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  c[3]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  c[4]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  c[5]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  c[5]->Print(Form("./pdf/fit_%s.pdf]",ofname_pdf.c_str())  );
 }

} 
