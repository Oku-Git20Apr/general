#include "Settings.cc"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void draw_coin(){
gStyle->SetOptStat(0);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.13);
gStyle->SetPadLeftMargin(0.16);
gStyle->SetPadBottomMargin(0.15);
int rebin =1;

TString ifname1=Form("./root/all.root");//
  TFile *ifp1  = new TFile(ifname1  );
      cout<<"input filename1 : "<<ifname1<<endl;
Settings *set=new Settings();

  TH1F *h_inv_mass,*h_ctime_allpi_rf,*h_ctime_allpi;
  TH1F *h_ctime_tb_RF_all;
  TH2F *h2_pid_decut;
  TH1F *h_inv_coin[6];//cointime cut dependence

  TH2F *h2_cointime_beta_tdll_pipi[10], *h2_cointime_beta_tdlr_pipi[10], *h2_cointime_ftpi_tdll_pipi[10], *h2_cointime_ftpi_tdlr_pipi[10];
////////////
//Get hist//
////////////
  for(int i=0;i<10;i++){//TDL seg
    h2_cointime_beta_tdll_pipi[i]   = (TH2F*)ifp1->Get(Form("pipi/h2_cointime_beta_tdll_pipi%d",i+1) );
    h2_cointime_beta_tdlr_pipi[i]   = (TH2F*)ifp1->Get(Form("pipi/h2_cointime_beta_tdlr_pipi%d",i+1) );
    h2_cointime_ftpi_tdll_pipi[i]   = (TH2F*)ifp1->Get(Form("pipi/h2_cointime_ftpi_tdll_pipi%d",i+1) );
    h2_cointime_ftpi_tdlr_pipi[i]   = (TH2F*)ifp1->Get(Form("pipi/h2_cointime_ftpi_tdlr_pipi%d",i+1) );
    h2_cointime_beta_tdll_pipi[i]   ->Rebin2D(2,2);
    h2_cointime_beta_tdlr_pipi[i]   ->Rebin2D(2,2);
    h2_cointime_ftpi_tdll_pipi[i]   ->Rebin2D(2,2);
    h2_cointime_ftpi_tdlr_pipi[i]   ->Rebin2D(2,2);
  }


  TLatex *tex = new TLatex(0,0,"aaa");
  tex  -> SetTextSize(0.060);
  tex  -> SetTextFont(42);
  tex  -> SetTextAlign(21);

  TArrow *arrow = new TArrow(0,0,0,0,0.01,"");//x1,y1,x2,y2,arrow size,option
  arrow->SetAngle(45);
  arrow->SetLineColor(1);
  arrow->SetLineWidth(3);

  TLine *line = new TLine(0,0,0,0);//x1,y1,x2,y2,arrow size,option
  line->SetLineColor(1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);

  TF1 *f_guide = new TF1("f_guide","[0]*x+[1]",-2,2);
  f_guide ->SetParameter(0,1);
  f_guide ->SetParameter(1,0.45);
  set->SetTF1(f_guide,1,2,1.5);

TCanvas *c[5];
for(int i=0;i<4;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),2400/2,1800/3);}
c[0]->Clear();
c[0]->Divide(5,2);
for(int i=0;i<10;i++){
  c[0]->cd(i+1);gPad->SetLogz(0);h2_cointime_ftpi_tdll_pipi[i]->Draw("colz");
  f_guide->Draw("same");
}

c[1]->Clear();
c[1]->Divide(5,2);
for(int i=0;i<10;i++){
  c[1]->cd(i+1);gPad->SetLogz(0);h2_cointime_ftpi_tdlr_pipi[i]->Draw("colz");
  f_guide->Draw("same");
}

c[2]->Clear();
c[2]->Divide(5,2);
for(int i=0;i<10;i++){
  c[2]->cd(i+1);gPad->SetLogz(1);h2_cointime_beta_tdll_pipi[i]->Draw("colz");
}

c[3]->Clear();
c[3]->Divide(5,2);
for(int i=0;i<10;i++){
  c[3]->cd(i+1);gPad->SetLogz(1);h2_cointime_beta_tdlr_pipi[i]->Draw("colz");
}

#if 1
c[0]->Print("pdf/vertex/coin_all.pdf[");
c[0]->Print("pdf/vertex/coin_all.pdf");
c[1]->Print("pdf/vertex/coin_all.pdf");
c[2]->Print("pdf/vertex/coin_all.pdf");
c[3]->Print("pdf/vertex/coin_all.pdf");
//c[4]->Print("pdf/vertex/coin_all.pdf");
c[3]->Print("pdf/vertex/coin_all.pdf]");
#endif
}
