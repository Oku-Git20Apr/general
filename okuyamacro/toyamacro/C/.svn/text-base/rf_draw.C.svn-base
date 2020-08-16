//macro to draw histograms made by RF_ana
//TagB 1-24 vs RF
//TDL 21-24 vs RF
//2019.1.10 Y. Toyama
#include "../Setting.cc"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gaus_pol0bg(double *x, double *par) {
  //par[0]=area (gaus)
  //par[1]=mean (gaus)
  //par[2]=sigma (gaus)
  //par[3]=pol const
  //par[4]=pol tilt
  double val;
  double bg;
  double ga;

ga=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
bg =  par[3];
val = ga + bg;
  return val;
}

//____________________________________________________________________________________________

const int NCanvas =8;
void rf_draw(string root_file="../rf_10634.root"){
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetOptStat(0);

  cout<<root_file<<endl;
  TFile* ifp = new TFile( root_file.c_str(), "READONLY" );
  TLatex *latex = new TLatex();
  Setting *set = new Setting();
  set->SetTLatex(latex,1,0.05,33);

  TH1F *h_ctime_tb_RF[40], *h_ctime_tb_RF_fcut[40][4], *h_time_tb_RF[40];
  TH2F *h2_ctime_tb_RF[40],*h2_time_tb_RF[40];

  TGraphErrors *tg_tbrfreso = new TGraphErrors();
  set->SetGrErr(tg_tbrfreso,"","","",1,2,22);
  tg_tbrfreso->SetMarkerSize(1.5);

  TGraphErrors *tg_tbrfreso_tbtb = new TGraphErrors();
  set->SetGrErr(tg_tbrfreso_tbtb,"","","",1,3,33);
  tg_tbrfreso_tbtb->SetMarkerSize(1.8);

  TGraphErrors *tg_tbrfreso_tbtb_bar = new TGraphErrors();
  set->SetGrErr(tg_tbrfreso_tbtb_bar,"","","",1,2,25);
  tg_tbrfreso_tbtb_bar->SetMarkerSize(1.8);

  TGraphErrors *tg_tbrfreso_fbor = new TGraphErrors();
  set->SetGrErr(tg_tbrfreso_fbor,"","","",1,6,34);
  tg_tbrfreso_fbor->SetMarkerSize(1.5);

  TGraphErrors *tg_tbreso = new TGraphErrors();
  set->SetGrErr(tg_tbreso,"","","",1,2,23);
  tg_tbreso->SetMarkerSize(1.5);

  TH2F *h2_frame1 = new TH2F("h2_frame1","h2_frame1",10,0,25,10,0,200);
  set->SetTH2(h2_frame1,"TagB-RF resolution(#sigma_{t})","TagB segment","#sigma_{t} [ps]");

  TH2F *h2_frame2 = new TH2F("h2_frame2","h2_frame2",10,0,25,10,0,200);
  set->SetTH2(h2_frame2,"TagB resolution(#sigma_{tb})","TagB segment","#sigma_{tb} [ps]");

  TLegend *leg_reso = new TLegend( 0.7, 0.20, 0.85, 0.35);
  leg_reso -> SetBorderSize(1);
  leg_reso -> SetFillColor(0);
  leg_reso -> SetFillStyle(1);
  leg_reso -> SetTextFont(22);

  TH1F *h_ctime_tdll_RF[24];
  TH1F *h_ctime_tb_RF_tbtb[40], *h_ctime_tb_RF_fbor[40];
  TH1F *h_ctime_tb_RF_tbtb_bar[40];
  TH2F *h2_ctime_tdllu_RF[24], *h2_ctime_tdlld_RF[24];
  double reso_tbrf[40],er_reso_tbrf[40];
  double reso_tbrf_tbtb[40],er_reso_tbrf_tbtb[40];//tagb(x)tagb coin
  double reso_tbrf_tbtb_bar[40],er_reso_tbrf_tbtb_bar[40];//tagb(x)!tagb coin
  double reso_tbrf_fbor[40],er_reso_tbrf_fbor[40];//
  double reso_tb[40],er_reso_tb[40];
  double reso_rf=92.,er_reso_rf=3.;//[ps]
  double x[40],y[40];

  for(int i=0;i<24;i++){
    cout<<"\r"<<i<<flush;
    h_ctime_tb_RF[i]  = (TH1F*)ifp ->Get(Form("h_ctime_tb_RF%d" ,i+1));
    h_time_tb_RF[i]   = (TH1F*)ifp ->Get(Form("h_time_tb_RF%d"  ,i+1));
    h2_ctime_tb_RF[i] = (TH2F*)ifp ->Get(Form("h2_ctime_tb_RF%d",i+1));
    h2_time_tb_RF[i]  = (TH2F*)ifp ->Get(Form("h2_time_tb_RF%d" ,i+1));
    h_ctime_tb_RF_tbtb[i]     = (TH1F*)ifp ->Get(Form("h_ctime_tb_RF_tbtb%d" ,i+1));
    h_ctime_tb_RF_tbtb_bar[i] = (TH1F*)ifp ->Get(Form("h_ctime_tb_RF_tbtb_bar%d" ,i+1));
    h_ctime_tb_RF_fbor[i]     = (TH1F*)ifp ->Get(Form("h_ctime_tb_RF_fbor%d" ,i+1));

    double time_min  =h_time_tb_RF[i] ->GetXaxis()->GetBinCenter(h_time_tb_RF[i] ->GetMaximumBin()) -2.;
    double time_max  =h_time_tb_RF[i] ->GetXaxis()->GetBinCenter(h_time_tb_RF[i] ->GetMaximumBin()) +2.;
    double ctime_min =h_ctime_tb_RF[i]->GetXaxis()->GetBinCenter(h_ctime_tb_RF[i]->GetMaximumBin()) -2.;
    double ctime_max =h_ctime_tb_RF[i]->GetXaxis()->GetBinCenter(h_ctime_tb_RF[i]->GetMaximumBin()) +2.;

    h_ctime_tb_RF[i]  ->GetXaxis()->SetRangeUser(ctime_min,ctime_max);
    h_time_tb_RF[i]   ->GetXaxis()->SetRangeUser(time_min ,time_max );
    h2_ctime_tb_RF[i] ->GetXaxis()->SetRangeUser(ctime_min,ctime_max);
    h2_time_tb_RF[i]  ->GetXaxis()->SetRangeUser(time_min ,time_max );
    h_ctime_tb_RF_tbtb[i]     ->GetXaxis()->SetRangeUser(ctime_min,ctime_max);
    h_ctime_tb_RF_tbtb_bar[i] ->GetXaxis()->SetRangeUser(ctime_min,ctime_max);
    h_ctime_tb_RF_fbor[i]     ->GetXaxis()->SetRangeUser(ctime_min,ctime_max);

    h_ctime_tb_RF[i]  ->SetStats(0);
    h_time_tb_RF[i]   ->SetStats(0);
    h2_ctime_tb_RF[i] ->SetStats(0);
    h2_time_tb_RF[i]  ->SetStats(0);
    h_ctime_tb_RF_tbtb[i]     ->SetStats(0);
    h_ctime_tb_RF_tbtb_bar[i] ->SetStats(0);
    h_ctime_tb_RF_fbor[i]     ->SetStats(0);
    TF1 *fadd = new TF1("fadd",gaus_pol0bg,ctime_min,ctime_max,4);
    set->SetTF1(fadd,2,1,1);
    fadd -> SetParameter(0,h2_ctime_tb_RF[i] ->Integral()*0.7);
    fadd -> SetParameter(1,0.5*(ctime_min+ctime_max));
    fadd -> SetParameter(2,0.12);
    fadd -> SetParameter(3,h_ctime_tb_RF[i] ->GetBinContent(h2_ctime_tb_RF[i] ->FindBin(0.5*(ctime_min+ctime_max)-0.8)));

    h_ctime_tb_RF[i]  ->Fit(fadd,"QR","",ctime_min+1.,ctime_max-1);
    reso_tbrf[i]    = fadd->GetParameter(2) * 1000.;//[ps]
    er_reso_tbrf[i] = fadd->GetParError(2)  * 1000.;//[ps]
    x[i] = ctime_max;
    y[i] = fadd->GetMaximum();
    //cout<<y[i]<<endl;
    
    //Fitting fbor
    TF1 *fadd_fbor = new TF1("fadd_fbor",gaus_pol0bg,ctime_min,ctime_max,4);
    set->SetTF1(fadd_fbor,7,1,1);
    fadd_fbor -> SetParameter(0,h_ctime_tb_RF_fbor[i] ->Integral()*0.7);
    fadd_fbor -> SetParameter(1,0.5*(ctime_min+ctime_max));
    fadd_fbor -> SetParameter(2,0.12);
    fadd_fbor -> SetParameter(3,h_ctime_tb_RF_fbor[i] ->GetBinContent(h_ctime_tb_RF_fbor[i] ->FindBin(0.5*(ctime_min+ctime_max)-0.8)));
    for(int k=0;k<5;k++)h_ctime_tb_RF_fbor[i]  ->Fit(fadd_fbor,"QR","",ctime_min+0.1,ctime_max-0.1);
    reso_tbrf_fbor[i]    = fadd_fbor->GetParameter(2) * 1000.;//[ps]
    er_reso_tbrf_fbor[i] = fadd_fbor->GetParError(2)  * 1000.;//[ps]

    //Fitting tbtb
    TF1 *fadd_tbtb = new TF1("fadd_tbtb",gaus_pol0bg,ctime_min,ctime_max,4);
    set->SetTF1(fadd_tbtb,4,1,1);
    fadd_tbtb -> SetParameter(0,h_ctime_tb_RF_tbtb[i] ->Integral()*0.7);
    fadd_tbtb -> SetParameter(1,0.5*(ctime_min+ctime_max));
    fadd_tbtb -> SetParameter(2,0.12);
    fadd_tbtb -> SetParameter(3,h_ctime_tb_RF_tbtb[i] ->GetBinContent(h_ctime_tb_RF_tbtb[i] ->FindBin(0.5*(ctime_min+ctime_max)-0.8)));
    for(int k=0;k<5;k++)h_ctime_tb_RF_tbtb[i]  ->Fit(fadd_tbtb,"QR","",ctime_min+0.1,ctime_max-0.1);
    reso_tbrf_tbtb[i]    = fadd_tbtb->GetParameter(2) * 1000.;//[ps]
    er_reso_tbrf_tbtb[i] = fadd_tbtb->GetParError(2)  * 1000.;//[ps]


    //Fitting tbtb
    TF1 *fadd_tbtb_bar = new TF1("fadd_tbtb_bar",gaus_pol0bg,ctime_min,ctime_max,4);
    set->SetTF1(fadd_tbtb_bar,3,1,1);
    fadd_tbtb_bar -> SetParameter(0,h_ctime_tb_RF_tbtb_bar[i] ->Integral()*0.7);
    fadd_tbtb_bar -> SetParameter(1,0.5*(ctime_min+ctime_max));
    fadd_tbtb_bar -> SetParameter(2,0.12);
    fadd_tbtb_bar -> SetParameter(3,h_ctime_tb_RF_tbtb_bar[i] ->GetBinContent(h_ctime_tb_RF_tbtb_bar[i] ->FindBin(0.5*(ctime_min+ctime_max)-0.8)));
    for(int k=0;k<5;k++)h_ctime_tb_RF_tbtb_bar[i]  ->Fit(fadd_tbtb_bar,"QR","",ctime_min+0.1,ctime_max-0.1);
    reso_tbrf_tbtb_bar[i]    = fadd_tbtb_bar->GetParameter(2) * 1000.;//[ps]
    er_reso_tbrf_tbtb_bar[i] = fadd_tbtb_bar->GetParError(2)  * 1000.;//[ps]

  }

  for(int i=0;i<24;i++){
    h_ctime_tdll_RF[i]  =(TH1F*)ifp ->Get(Form("h_ctime_tdll_RF%d"  ,i+1));
    h2_ctime_tdllu_RF[i]=(TH2F*)ifp ->Get(Form("h2_ctime_tdllu_RF%d",i+1)); 
    h2_ctime_tdlld_RF[i]=(TH2F*)ifp ->Get(Form("h2_ctime_tdlld_RF%d",i+1));
  }

  TCanvas *c[NCanvas];
  for(int i=0;i<NCanvas;i++){
    c[i]=new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800,800);
  }

  for(int i=0;i<6;i++){ //Canvas
    c[i]->Divide(4,2);
    for(int k=0;k<4;k++){ //TagB seg
      seg = 4*i + k +1;
      c[i]->cd(1+k);gPad->SetLogz(1);h2_ctime_tb_RF[seg-1]->Draw("colz");
      //c[i]->cd(1+k);h_time_tb_RF[seg-1]->Draw();
      c[i]->cd(5+k);h_ctime_tb_RF[seg-1]->Draw();latex->DrawLatex(x[seg-1],y[seg-1],Form("#sigma = %.1lf #pm %.1lfps",reso_tbrf[seg-1],er_reso_tbrf[seg-1]));
                                                 h_ctime_tb_RF_tbtb_bar[seg-1]->Draw("same");
                                                 h_ctime_tb_RF_tbtb[seg-1]  ->Draw("same");
      tg_tbrfreso->SetPoint(seg-1,seg,reso_tbrf[seg-1]);
      tg_tbrfreso->SetPointError(seg-1,0.,er_reso_tbrf[seg-1]);

      tg_tbrfreso_fbor->SetPoint(seg-1,seg,reso_tbrf_fbor[seg-1]);
      tg_tbrfreso_fbor->SetPointError(seg-1,0.,er_reso_tbrf_fbor[seg-1]);

      tg_tbrfreso_tbtb->SetPoint(seg-1,seg,reso_tbrf_tbtb[seg-1]);
      tg_tbrfreso_tbtb->SetPointError(seg-1,0.,er_reso_tbrf_tbtb[seg-1]);

      tg_tbrfreso_tbtb_bar->SetPoint(seg-1,seg,reso_tbrf_tbtb_bar[seg-1]);
      tg_tbrfreso_tbtb_bar->SetPointError(seg-1,0.,er_reso_tbrf_tbtb_bar[seg-1]);

      reso_tb[seg-1]    = sqrt(reso_tbrf[seg-1]*reso_tbrf[seg-1]-reso_rf*reso_rf);
      er_reso_tb[seg-1] = sqrt(reso_tbrf[seg-1]*reso_tbrf[seg-1]*er_reso_tbrf[seg-1]*er_reso_tbrf[seg-1]+reso_rf*reso_rf*er_reso_rf*er_reso_rf)/reso_tb[seg-1];
      cout<<seg<<" "<<reso_tb[seg-1]<<" +/- "<<er_reso_tb[seg-1]<<endl;
      tg_tbreso->SetPoint(seg-1,seg,reso_tb[seg-1]);
      tg_tbreso->SetPointError(seg-1,0.,er_reso_tb[seg-1]);
    }
  }

  leg_reso -> AddEntry(tg_tbrfreso           ,"no cut","p");
  leg_reso -> AddEntry(tg_tbrfreso_fbor      ,"TagB#otimes(TagF or TagB)","p");
  leg_reso -> AddEntry(tg_tbrfreso_tbtb      ,"TagB#otimesTagB","p");
  leg_reso -> AddEntry(tg_tbrfreso_tbtb_bar  ,"TagB#otimes#bar{TagB}","p");
  c[6]->Clear();c[6]->Divide(1,1);
  c[6]->cd(1);gPad->SetGrid(1,1);
              h2_frame1 ->Draw();
              tg_tbrfreso         ->Draw("samePE");
              tg_tbrfreso_fbor    ->Draw("samePE");
              tg_tbrfreso_tbtb    ->Draw("samePE");
              tg_tbrfreso_tbtb_bar->Draw("samePE");
              leg_reso -> Draw("same");

  c[7]->Clear();c[7]->Divide(1,1);
  c[7]->cd(1);gPad->SetGrid(1,1);
              h2_frame2 ->Draw();
              tg_tbreso->Draw("samePE");

  string pdf_file = root_file;
  pdf_file.erase(pdf_file.size()-5);
  pdf_file.append("_reso.pdf");
  c[0]->Print(Form("%s[",pdf_file.c_str()));
  for(int i=0;i<NCanvas;i++)  c[i]->Print(Form("%s",pdf_file.c_str()));
  c[NCanvas-1]->Print(Form("%s]",pdf_file.c_str()));

  //TFile *ofp = new TFile("tdl24.root","recreate");
  //h_ctime_tb_RF[23]     ->Write();
  //h2_ctime_tb_RF[23]    ->Write();
  //h_ctime_tdll_RF[23]   ->Write();
  //h2_ctime_tdllu_RF[23] ->Write();
  //h2_ctime_tdlld_RF[23] ->Write();
  //ofp->Write();
  //ofp->Close();

}
