#include "Settings.cc"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void draw_DC2_2018(){
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.13);
  Settings *set=new Settings();

  /////////
  //Param//
  /////////
  int rebin =1;
  int rebin_time =2;
  double MassRange_Lam_min;
  double MassRange_Lam_max;
  MassRange_Lam_min   =  1.105;
  MassRange_Lam_max   =  1.122;       
  double scale = 4.0;

  TString ifname1=Form("../root/RKV_ana_acci/true.root");//
  TString ofname =Form("./root/dc2_2018.root");//
  TFile *ifp1  = new TFile(ifname1  );
      cout<<"input filename1 : "<<ifname1<<endl;

  TArrow *arrow = new TArrow(0,0,0,0,0.01,"");//x1,y1,x2,y2,arrow size,option
  arrow->SetAngle(45);
  arrow->SetLineColor(1);
  arrow->SetLineWidth(3);

  TLine *line = new TLine(0,0,0,0);//x1,y1,x2,y2,arrow size,option
  line->SetLineColor(1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);

  TLatex *text[2];
  for(int i=0;i<2;i++){
    text[i] = new TLatex();
    set->SetTLatex(text[i],2-i,0.065,32);
  }



  //Get histograms of true event
  TH1F *h_invm_lpi;
  h_invm_lpi   = (TH1F*)ifp1->Get("vertex/h_invm_lpi");
  h_invm_lpi   ->Rebin(rebin);
  h_invm_lpi   ->SetMinimum(0);

  TH1F *h_dtimeLam = (TH1F*)ifp1->Get("cointime/h_ctime_Lam_cor");
  TH1F *h_dtimeRes = (TH1F*)ifp1->Get("pipi/h_ctime_pipi_cor");
  h_dtimeLam ->Rebin(4);
  set->SetTH1(h_dtimeLam,"","","Counts",2,3004,2);
  set->SetTH1(h_dtimeRes,"","Decay time [ns]","Counts"                  ,1,3000,0);

  SetTH1(h_invm_lpi,    "" ,"Invariant mass [GeV/#it{c}^{2}]"      ,Form("Counts/%.1lfMeV/#it{c}^{2}",2.0*rebin), 1, 3004, 17);
  h_invm_lpi->GetXaxis()->SetTitleSize(0.070);
  h_invm_lpi->GetXaxis()->SetTitleOffset(0.85);
  h_invm_lpi->GetXaxis()->SetLabelSize(0.06);
  h_invm_lpi->GetXaxis()->SetLabelOffset(0.01);

  h_invm_lpi->GetYaxis()->SetTitleSize(0.070);
  h_invm_lpi->GetYaxis()->SetTitleOffset(0.7);
  h_invm_lpi->GetYaxis()->SetLabelSize(0.06);
  h_invm_lpi->GetYaxis()->SetLabelOffset(0.01);
  h_invm_lpi       ->GetXaxis()->SetRangeUser(1.05,1.30);

  h_dtimeRes->GetXaxis()->SetTitleSize(0.070);
  h_dtimeRes->GetXaxis()->SetTitleOffset(0.85);
  h_dtimeRes->GetXaxis()->SetLabelSize(0.06);
  h_dtimeRes->GetXaxis()->SetLabelOffset(0.01);

  h_dtimeRes->GetYaxis()->SetTitleSize(0.070);
  h_dtimeRes->GetYaxis()->SetTitleOffset(0.7);
  h_dtimeRes->GetYaxis()->SetLabelSize(0.06);
  h_dtimeRes->GetYaxis()->SetLabelOffset(0.01);

  TH1F *h_dtimeRes_scale = (TH1F*)h_dtimeRes->Clone();
  double maxyLam = h_dtimeLam->GetBinContent(h_dtimeLam->GetMaximumBin());
  double maxyRes = h_dtimeRes->GetBinContent(h_dtimeRes->GetMaximumBin());
  //cout<<maxyLam/maxyRes<<endl;
  h_dtimeRes_scale->Scale(0.9*maxyLam/maxyRes);
  h_dtimeRes_scale->SetMinimum(0.1);



  TCanvas *c[3];
  for(int i=0;i<3;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1200,1000);}
  c[0]->Clear();
  c[0]->Divide(1,2);
  c[0]->cd(1);
  gPad->SetLogy(0);
  h_invm_lpi->Draw("");
  line->DrawLine(MassRange_Lam_min,0,MassRange_Lam_min,h_invm_lpi->GetBinContent(h_invm_lpi->GetMaximumBin()));
  line->DrawLine(MassRange_Lam_max,0,MassRange_Lam_max,h_invm_lpi->GetBinContent(h_invm_lpi->GetMaximumBin()));
  c[0]->cd(2);
  gPad->SetLogy(1);
  h_dtimeRes_scale ->Draw();
  h_dtimeLam ->Draw("same");
  
  c[1]->Clear();
  c[1]->Divide(1,2);
  c[1]->cd(1);
  h_dtimeLam ->Draw("");
  c[1]->cd(2);
  gPad->SetLogy(1);
  h_dtimeRes_scale ->Draw();
  h_dtimeLam ->Draw("same");
  h_dtimeRes_scale ->Draw("same");
  text[0]->DrawLatex(1.8,maxyLam    ,"#Lambda decay" );
  text[1]->DrawLatex(1.8,maxyLam*0.5,"Response func.");

  c[2]->Clear();
  c[2]->Divide(1,2);
  c[2]->cd(1);
  gPad->SetLogy(1);
  h_dtimeRes_scale ->Draw();
  h_dtimeLam ->Draw("same");

  c[0]->Print("pdf/toyama_dc2_2018root.pdf[");
  c[0]->Print("pdf/toyama_dc2_2018root.pdf");
  c[1]->Print("pdf/toyama_dc2_2018root.pdf");
  c[2]->Print("pdf/toyama_dc2_2018root.pdf");
  c[2]->Print("pdf/toyama_dc2_2018root.pdf]");
}
