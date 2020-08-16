#include "Settings.cc"
void draw_decaytime_w_resfunc(){
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.12);
  Settings *set=new Settings();

  TString ifname1=Form("./root/RKV_ana_acci/true.root");
  TString ofname =Form("./root/dcaytime.root");//
  TFile *ifp1  = new TFile(ifname1  );
  cout<<"input filename1 : "<<ifname1<<endl;
  TFile *ofp   = new TFile(ofname,"recreate"  );

  TH1F *h_dtimeLam = (TH1F*)ifp1->Get("cointime/h_ctime_Lam_cor");
  TH1F *h_dtimeRes = (TH1F*)ifp1->Get("pipi/h_ctime_pipi_cor");
  TH1F *h_mm_pipi = (TH1F*)ifp1->Get("vertex/h_mm_pipi");

  TLatex *text[2];
  for(int i=0;i<2;i++){
    text[i] = new TLatex();
    set->SetTLatex(text[i],2-i,0.065,32);
  }

  int rebin = 4;
  h_dtimeLam ->Rebin(rebin);
  set->SetTH1(h_dtimeLam,"","",Form("Counts/%dps",25*rebin),2,3000,0);
  set->SetTH1(h_dtimeRes,"","decay time [ns]","Counts/25ps"                  ,1,3000,0);

  TH1F *h_dtimeRes_scale = (TH1F*)h_dtimeRes->Clone();
  double maxyLam = h_dtimeLam->GetBinContent(h_dtimeLam->GetMaximumBin());
  double maxyRes = h_dtimeRes->GetBinContent(h_dtimeRes->GetMaximumBin());
  h_dtimeRes_scale->Scale(maxyLam/maxyRes);
  h_dtimeRes_scale->SetMinimum(0.1);
  int NbinxLam = h_dtimeLam ->GetNbinsX();
  int NbinxRes = h_dtimeRes ->GetNbinsX();

  TCanvas *c[3];
  for(int i=0;i<3;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1000,1200);
  }

  c[0]->Divide(1,2);
  c[0]->cd(1);
  gPad->SetLogy(1);
  h_dtimeLam ->Draw();
  c[0]->cd(2);
  gPad->SetLogy(1);
  h_dtimeRes ->Draw();

  c[1]->Divide(1,2);
  c[1]->cd(1);
  gPad->SetLogy(1);
  h_dtimeRes_scale ->Draw();
  h_dtimeLam ->Draw("same");
  text[0]->DrawLatex(1.8,maxyLam    ,"#Lambda event" );
  text[1]->DrawLatex(1.8,maxyLam*0.5,"Response func.");

  c[1]->cd(2);
  gPad->SetLogy(1);
  h_dtimeRes_scale ->Draw();
  h_dtimeLam ->Draw("same");

  c[2]->Divide(1,2);
  c[2]->cd(1);
  h_mm_pipi->Draw();

  c[0]->Print("pdf/decaytime_w_resfunc.pdf[");
  c[0]->Print("pdf/decaytime_w_resfunc.pdf");
  c[1]->Print("pdf/decaytime_w_resfunc.pdf");
  c[2]->Print("pdf/decaytime_w_resfunc.pdf");
  c[2]->Print("pdf/decaytime_w_resfunc.pdf]");
}
