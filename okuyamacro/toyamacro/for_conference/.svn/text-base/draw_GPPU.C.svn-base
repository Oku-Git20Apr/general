#include "Settings.cc"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void draw_GPPU(){
gStyle->SetOptStat(0);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.13);
gStyle->SetPadLeftMargin(0.16);
gStyle->SetPadBottomMargin(0.15);
int rebin =1;

TString ifname1=Form("./root/ana10226_10350.root");//
  TFile *ifp1  = new TFile(ifname1  );
      cout<<"input filename1 : "<<ifname1<<endl;
Settings *set=new Settings();

  TH1F *h_inv_mass,*h_ctime_allpi_rf,*h_ctime_allpi;
  TH1F *h_ctime_tb_RF_all;
  TH2F *h2_pid_decut;

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

h_inv_mass   = (TH1F*)ifp1->Get("vertex/h_inv_mass");
h_inv_mass   ->Rebin(rebin);

h_ctime_tb_RF_all = (TH1F*)ifp1->Get("tag/h_ctime_tb_RF_all");
h_ctime_allpi     = (TH1F*)ifp1->Get("cointime/h_ctime_allpi");
h_ctime_allpi_rf  = (TH1F*)ifp1->Get("cointime/h_ctime_allpi_rf");
h2_pid_decut = (TH2F*)ifp1->Get("track/h2_pid_decut");
SetTH1(h_inv_mass,    "" ,"invariant mass [GeV/#it{c}^{2}]"      ,Form("counts/%.1lfMeV/#it{c}^{2}",2.5*rebin), 1, 3002, 2);
SetTH1(h_ctime_tb_RF_all,"TagB - RF" ,"TagB - RF [ns]","counts/6ps", 1, 3001, 9);
SetTH1(h_ctime_allpi    ,"" ,"coin. time [ns]","counts/25ps", 1, 3000, 0);
SetTH1(h_ctime_allpi_rf ,"" ,"coin. time [ns]","counts/25ps", 1, 3001, 4);
SetTH2(h2_pid_decut,"" ,"1/#beta","momentum[GeV/#it{c}]");
h_inv_mass       ->GetXaxis()->SetRangeUser(1.05,1.30);
h_ctime_tb_RF_all->GetXaxis()->SetRangeUser(-1,1);
h_ctime_allpi    ->GetXaxis()->SetRangeUser(-3,3);
h_ctime_allpi_rf ->GetXaxis()->SetRangeUser(-3,3);

TCanvas *c[4];
for(int i=0;i<4;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800/2,1500/3);}
c[0]->Clear();
gPad->SetLogz(1);h2_pid_decut->Draw("colz");
//  f_pi->Draw("same");
//  f_k ->Draw("same");
//  f_p ->Draw("same");

c[1]->Clear();
gPad->SetLogy(1);h_ctime_allpi->Draw("");h_ctime_allpi_rf->Draw("same");

c[2]->Clear();
h_inv_mass->Draw("");
  //tex  ->DrawLatex(0.001*ML,12,"#Lambda mass");
  //arrow->DrawArrow(0.001*ML, 9,0.001*ML,0.5,0.02,"|>");
  //line ->DrawLine(0.001*ML,0.5,0.001*ML,h_inv_mass->GetMaximum());

c[3]->Clear();
gPad->SetLogy(1);h_ctime_tb_RF_all->Draw("");
  line ->DrawLine( 0.5,0., 0.5, h_ctime_tb_RF_all->GetMaximum());
  line ->DrawLine(-0.5,0.,-0.5, h_ctime_tb_RF_all->GetMaximum());
  arrow->DrawArrow( 0.5, h_ctime_tb_RF_all->GetMaximum(), 0.3,h_ctime_tb_RF_all->GetMaximum(),0.02,"|>");
  arrow->DrawArrow(-0.5, h_ctime_tb_RF_all->GetMaximum(),-0.3,h_ctime_tb_RF_all->GetMaximum(),0.02,"|>");


c[0]->Print("pdf/GPPU2017A.pdf[");
c[0]->Print("pdf/GPPU2017A.pdf");
c[1]->Print("pdf/GPPU2017A.pdf");
c[2]->Print("pdf/GPPU2017A.pdf");
c[3]->Print("pdf/GPPU2017A.pdf");
c[3]->Print("pdf/GPPU2017A.pdf]");
}
