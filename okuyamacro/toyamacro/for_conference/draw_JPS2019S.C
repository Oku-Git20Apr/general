#include "Setting.cc"
const int NCanvas = 4;
void draw_JPS2019S(){
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.15);
  TString ifname1=Form("./r10634_tdl24.root");//
  TFile *ifp  = new TFile(ifname1  );
  cout<<"input filename1 : "<<ifname1<<endl;
  Setting *set=new Setting();
  TLatex *tex = new TLatex(0,0,"aaa");
  tex  -> SetTextSize(0.060);
  tex  -> SetTextFont(42);
  tex  -> SetTextAlign(21);

  TH1F *hc, *h_ctime_tdll_RF24,*h_ctime_tb_RF24;
  TH2F *h2_ctime_tb_RF24, *h2_ctime_tdllu_RF24, *h2_ctime_tdlld_RF24;

  TF1 *f_c, *f_tdl_RF, *f_tb_RF;
  f_c      = new TF1("f_c","gaus",-3,3);
  f_tdl_RF = new TF1("f_tdl_RF","gaus(0)+[3]",-3,3);
  f_tb_RF  = new TF1("f_tb_RF","gaus(0)+[3]",-3,3);
  set->SetTF1(f_c       , 2,1,2);
  set->SetTF1(f_tdl_RF  , 2,1,2);
  set->SetTF1(f_tb_RF   , 2,1,2);

  hc                = (TH1F*)ifp->Get("hc"               );
  hc                ->Rebin(4);
  //hc                ->GetXaxis()->SetRangeUser(0,2);
  hc                ->SetStats(0);
  set->SetTH1(hc,"TDL - TagB","ToF [ns]","Counts",1,3001,4);
  h_ctime_tdll_RF24 = (TH1F*)ifp->Get("h_ctime_tdll_RF24");
  h_ctime_tdll_RF24 ->SetStats(0);
  set->SetTH1(h_ctime_tdll_RF24,"TDL - RF","ToF [ns]","Counts",1,3001,3);
  h_ctime_tb_RF24   = (TH1F*)ifp->Get("h_ctime_tb_RF24"  );
  h_ctime_tb_RF24   ->SetStats(0);
  set->SetTH1(h_ctime_tb_RF24,"TagB - RF","ToF [ns]","Counts",1,3001,5);

  h2_ctime_tb_RF24    =(TH2F*)ifp->Get("h2_ctime_tb_RF24"   ); 
  h2_ctime_tdllu_RF24 =(TH2F*)ifp->Get("h2_ctime_tdllu_RF24"); 
  h2_ctime_tdlld_RF24 =(TH2F*)ifp->Get("h2_ctime_tdlld_RF24");

  double f_c_min=0.;
  double f_c_max=3.;
  double f_tdl_min=-2;
  double f_tdl_max= 2;
  double f_tb_min =-2;
  double f_tb_max = 2;

  set->FitGaus(hc               ,f_c_min  ,f_c_max  ,3.0,5);
  set->FitGaus(h_ctime_tdll_RF24,f_tdl_min,f_tdl_max,3.0,5);
  set->FitGaus(h_ctime_tb_RF24  ,f_tb_min ,f_tb_max ,3.0,5);

  //cout<<f_tdl_min<<" "<<f_tdl_max<<endl;
  f_tdl_RF ->SetParameter(1,0.5*(f_tdl_min+f_tdl_max));
  f_tdl_RF ->SetParameter(2,(f_tdl_max-f_tdl_min)/6.);
  f_tdl_RF ->SetParameter(3,h_ctime_tdll_RF24->GetBinContent(h_ctime_tdll_RF24->GetBin(f_tdl_min-0.3)));

  f_tb_RF ->SetParameter(1,0.5*(f_tb_min+f_tb_max));
  f_tb_RF ->SetParameter(2,(f_tb_max-f_tb_min)/6.);
  f_tb_RF ->SetParameter(3,h_ctime_tb_RF24->GetBinContent(h_ctime_tb_RF24->GetBin(f_tb_min-0.4)));


  hc                ->Fit(f_c     ,"","QR",f_c_min  ,f_c_max  );
  h_ctime_tdll_RF24 ->Fit(f_tdl_RF,"","QR",f_tdl_min-0.3,f_tdl_max+0.3);
  h_ctime_tb_RF24   ->Fit(f_tb_RF ,"","QR",f_tb_min-0.3 ,f_tb_max+0.3 );

  double cpar[3],tdlpar[4],tbpar[4];
  double er_c,er_tdl,er_tb;           //error of width parameter
  f_c      -> GetParameters(&cpar[0]);
  f_tdl_RF -> GetParameters(&tdlpar[0]);
  f_tb_RF  -> GetParameters(&tbpar[0]);

  er_c   = f_c      -> GetParError(2);
  er_tdl = f_tdl_RF -> GetParError(2);
  er_tb  = f_tb_RF  -> GetParError(2);
  hc                ->GetXaxis()->SetRangeUser(cpar[1]  -1.5,cpar[1]  +1.5);
  h_ctime_tdll_RF24 ->GetXaxis()->SetRangeUser(tdlpar[1]-1.5,tdlpar[1]+1.5);
  h_ctime_tb_RF24   ->GetXaxis()->SetRangeUser(tbpar[1] -1.5,tbpar[1] +1.5);

  //intrinsic time resolution//
  double reso_tdl = sqrt(0.5*(cpar[2]*cpar[2]+tdlpar[2]*tdlpar[2]-tbpar[2]*tbpar[2]));
  double reso_tb  = sqrt(0.5*(cpar[2]*cpar[2]-tdlpar[2]*tdlpar[2]+tbpar[2]*tbpar[2]));
  double reso_rf  = sqrt(0.5*(tdlpar[2]*tdlpar[2]+tbpar[2]*tbpar[2]-cpar[2]*cpar[2]));
  double er_prop = sqrt(er_c*er_c*cpar[2]*cpar[2]+er_tdl*er_tdl*tdlpar[2]*tdlpar[2]+er_tb*er_tb*tbpar[2]*tbpar[2]);//constant value
  double er_reso_tdl = 0.5*er_prop/reso_tdl;
  double er_reso_tb  = 0.5*er_prop/reso_tb;
  double er_reso_rf  = 0.5*er_prop/reso_rf;
  cout<<"-**-Time resolution-**-"<<endl;
  cout<<"TDL-TagB "<<1000.*cpar[2]  <<" +/- "<<1000.*er_c       <<endl; 
  cout<<"TDL-RF   "<<1000.*tdlpar[2]<<" +/- "<<1000.*er_tdl     <<endl; 
  cout<<"TagB-RF  "<<1000.*tbpar[2] <<" +/- "<<1000.*er_tb      <<endl; 
  cout<<"TDL      "<<1000.*reso_tdl <<" +/- "<<1000.*er_reso_tdl<<endl;
  cout<<"TagB     "<<1000.*reso_tb  <<" +/- "<<1000.*er_reso_tb <<endl;
  cout<<"RF       "<<1000.*reso_rf  <<" +/- "<<1000.*er_reso_rf <<endl;

  TCanvas *c[NCanvas];
  for(int i=0;i<NCanvas;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800/2,1200/2);
  }
  
  c[0]->Clear();
  c[0]->cd(1);hc               ->Draw();
  c[1]->Clear();
  c[1]->cd(1);h_ctime_tdll_RF24->Draw();
  c[2]->Clear();
  c[2]->cd(1);h_ctime_tb_RF24  ->Draw();

  c[3]->Clear();
  c[3]->Divide(2,2);
  c[3]->cd(1);h2_ctime_tb_RF24   ->Draw();
  c[3]->cd(2);h2_ctime_tdllu_RF24->Draw();
  c[3]->cd(3);h2_ctime_tdlld_RF24->Draw();

  c[0]->Print("pdf/JPS2019S.pdf[");
  for(int i=0;i<NCanvas;i++){
  c[i]->Print("pdf/JPS2019S.pdf");
  }
  c[NCanvas-1]->Print("pdf/JPS2019S.pdf]");
}
