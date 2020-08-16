//macro to analyze histograms made by RF_ana
//Plz execute at toyamacro/C directory
//This macro shows TagB and TagF coincidence ratio etc...
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

const int NCanvas =6;
void rf_coinratio(string root_file="../rf_10634.root"){
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.13);

  cout<<root_file<<endl;
  TFile* ifp = new TFile( root_file.c_str(), "READONLY" );
  TLatex *latex = new TLatex();
  Setting *set = new Setting();
  set->SetTLatex(latex,1,0.05,33);

  TLegend *leg_coin = new TLegend( 0.7, 0.20, 0.85, 0.35);
  leg_coin -> SetBorderSize(1);
  leg_coin -> SetFillColor(0);
  leg_coin -> SetFillStyle(1);
  leg_coin -> SetTextFont(22);

  TLegend *leg_SN = new TLegend( 0.65, 0.20, 0.85, 0.4);
  leg_SN -> SetBorderSize(1);
  leg_SN -> SetFillColor(0);
  leg_SN -> SetFillStyle(1);
  leg_SN -> SetTextFont(22);
  
  TH1F *h_ctime_tb_RF[40],*h_ctime_tb_RF_cut[40], *h_ctime_tb_RF_f_or[40], *h_ctime_tb_RF_fcut[40][4], *h_time_tb_RF[40];
  TH1F *h_ctime_tb_RF_tbtb[40], *h_ctime_tb_RF_fbor[40];
  TH1F *h_ctime_tb_RF_f_or_scale[40];
  TH2F *h2_ctime_tb_RF[40],*h2_time_tb_RF[40];

  TGraphErrors *tg_tbcoinratio = new TGraphErrors();
  set->SetGrErr(tg_tbcoinratio,"","","",0,4,22);
  tg_tbcoinratio->SetMarkerSize(1.3);

  TGraphErrors *tg_fborcoinratio = new TGraphErrors();
  set->SetGrErr(tg_fborcoinratio,"","","",0,6,23);
  tg_fborcoinratio->SetMarkerSize(1.3);

  TGraphErrors *tg_tbSN = new TGraphErrors();
  set->SetGrErr(tg_tbSN,"","","",0,2,23);
  tg_tbSN->SetMarkerSize(1.3);

  TGraphErrors *tg_tbSNi = new TGraphErrors();
  set->SetGrErr(tg_tbSNi,"","","",0,2,21);
  tg_tbSNi->SetMarkerSize(1.3);

  TGraphErrors *tg_tbSN_tbtb = new TGraphErrors();
  set->SetGrErr(tg_tbSN_tbtb,"","","",0,3,34);
  tg_tbSN_tbtb->SetMarkerSize(1.5);

  TGraphErrors *tg_tbSN_fbor = new TGraphErrors();
  set->SetGrErr(tg_tbSN_fbor,"","","",0,6,29);
  tg_tbSN_fbor->SetMarkerSize(1.5);

  TGraphErrors *tg_tbSN_f_or = new TGraphErrors();
  set->SetGrErr(tg_tbSN_f_or,"","","",0,4,33);
  tg_tbSN_f_or->SetMarkerSize(1.5);

  TH2F *h2_frame1 = new TH2F("h2_frame1","h2_frame1",10,0,25,10,0,1);
  set->SetTH2(h2_frame1,"TagB-TagF coin. ratio","TagB segment","Coin. Ratio");

  TH2F *h2_frame2 = new TH2F("h2_frame2","h2_frame2",10,0,25,10,0,1.05);
  set->SetTH2(h2_frame2,"Signal to noise ratio","TagB segment","S/(S+N)");

  double reso_tbrf[40],er_reso_tbrf[40];
  double reso_tb[40],er_reso_tb[40];
  double reso_rf=92.,er_reso_rf=3.;//[ps]
  double x[40],y[40];
  double total_tb[40],total_fcoin[40];
  double tf_ratio[40],er_tf_ratio[40];

  double bin_w;//bin size [ns]
  double SN[40] ,area_S[40] ,area_BG[40];//area of Signal and BG obtained from fitting
  double SNi[40],area_Si[40],area_BGi[40];//area of Signal and BG obtained from integral
  double SN_f_or[40] ,area_S_f_or[40] ,area_BG_f_or[40];//area of Signal and BG obtained from fitting
  double SNi_f_or[40],area_Si_f_or[40],area_BGi_f_or[40];//area of Signal and BG obtained from integral
  double SN_tbtb[40] ,area_S_tbtb[40] ,area_BG_tbtb[40];//area of Signal and BG obtained from fitting
  double SNi_tbtb[40],area_Si_tbtb[40],area_BGi_tbtb[40];//area of Signal and BG obtained from integral
  double SN_fbor[40] ,area_S_fbor[40] ,area_BG_fbor[40];//area of Signal and BG obtained from fitting
  double SNi_fbor[40],area_Si_fbor[40],area_BGi_fbor[40];//area of Signal and BG obtained from integral
  

  for(int i=0;i<24;i++){
    cout<<"\r"<<i<<flush;
    h_ctime_tb_RF[i]      = (TH1F*)ifp ->Get(Form("h_ctime_tb_RF%d"      ,i+1));
    h_ctime_tb_RF_cut[i]  = (TH1F*)ifp ->Get(Form("h_ctime_tb_RF_cut%d"  ,i+1));
    h_time_tb_RF[i]       = (TH1F*)ifp ->Get(Form("h_time_tb_RF%d"       ,i+1));
    h2_ctime_tb_RF[i]     = (TH2F*)ifp ->Get(Form("h2_ctime_tb_RF%d"     ,i+1));
    h2_time_tb_RF[i]      = (TH2F*)ifp ->Get(Form("h2_time_tb_RF%d"      ,i+1));
    h_ctime_tb_RF_f_or[i] = (TH1F*)ifp ->Get(Form("h_ctime_tb_RF_f_or%d" ,i+1));
    h_ctime_tb_RF_tbtb[i] = (TH1F*)ifp ->Get(Form("h_ctime_tb_RF_tbtb%d" ,i+1));
    h_ctime_tb_RF_fbor[i] = (TH1F*)ifp ->Get(Form("h_ctime_tb_RF_fbor%d" ,i+1));
    h_ctime_tb_RF_f_or_scale[i] = (TH1F*)h_ctime_tb_RF_f_or[i]->Clone();
    for(int k=0;k<4;k++){
      h_ctime_tb_RF_fcut[i][k] = (TH1F*)ifp->Get(Form("h_ctime_tb_RF_f%dcut%d" ,k+1,i+1));
    }

           bin_w     =h_time_tb_RF[i] ->GetBinWidth(2);
    double time_min  =h_time_tb_RF[i] ->GetXaxis()->GetBinCenter(h_time_tb_RF[i] ->GetMaximumBin()) -2.;
    double time_max  =h_time_tb_RF[i] ->GetXaxis()->GetBinCenter(h_time_tb_RF[i] ->GetMaximumBin()) +2.;
    //double ctime_min =h_ctime_tb_RF_cut[i]->GetXaxis()->GetBinCenter(h_ctime_tb_RF_cut[i]->GetMaximumBin()) -2.;
    //double ctime_max =h_ctime_tb_RF_cut[i]->GetXaxis()->GetBinCenter(h_ctime_tb_RF_cut[i]->GetMaximumBin()) +2.;
    double ctime_min =-1.1;
    double ctime_max = 1.1;

    h_ctime_tb_RF_cut[i]  ->GetXaxis()->SetRangeUser(ctime_min,ctime_max);
    h_time_tb_RF[i]       ->GetXaxis()->SetRangeUser( time_min, time_max);
    h2_ctime_tb_RF[i]     ->GetXaxis()->SetRangeUser(ctime_min,ctime_max);
    h2_time_tb_RF[i]      ->GetXaxis()->SetRangeUser( time_min, time_max);
    h_ctime_tb_RF_tbtb[i] ->GetXaxis()->SetRangeUser(ctime_min,ctime_max);
    h_ctime_tb_RF_fbor[i] ->GetXaxis()->SetRangeUser(ctime_min,ctime_max);
    h_ctime_tb_RF_f_or_scale[i] = (TH1F*)h_ctime_tb_RF_f_or[i]->Clone();
    h_ctime_tb_RF_cut[i]  ->SetStats(0);
    h_ctime_tb_RF_cut[i]  ->SetMinimum(0.5);
    h_time_tb_RF[i]       ->SetStats(0);
    h2_ctime_tb_RF[i]     ->SetStats(0);
    h2_time_tb_RF[i]      ->SetStats(0);
    h_ctime_tb_RF_tbtb[i] ->SetStats(0);
    h_ctime_tb_RF_fbor[i] ->SetStats(0);

    //h_ctime_tb_RF_cut analysis
    //Fitting
    TF1 *fadd = new TF1("fadd",gaus_pol0bg,ctime_min,ctime_max,4);
    set->SetTF1(fadd,2,1,1);
    fadd -> SetParameter(0,h2_ctime_tb_RF[i] ->Integral()*0.7);
    fadd -> SetParameter(1,0.5*(ctime_min+ctime_max));
    fadd -> SetParameter(2,0.12);
    fadd -> SetParameter(3,h_ctime_tb_RF_cut[i] ->GetBinContent(h2_ctime_tb_RF[i] ->FindBin(0.5*(ctime_min+ctime_max)-0.8)));

    for(int k=0;k<5;k++)h_ctime_tb_RF_cut[i]  ->Fit(fadd,"QR","",ctime_min+0.1,ctime_max-0.1);
    reso_tbrf[i]    = fadd->GetParameter(2) * 1000.;//[ps]
    er_reso_tbrf[i] = fadd->GetParError(2)  * 1000.;//[ps]
    x[i] = ctime_max;

    total_tb[i]    = (double)h_ctime_tb_RF_cut[i] ->Integral();
    total_fcoin[i] = (double)h_ctime_tb_RF_f_or[i]->Integral();
    //tf_ratio[i]    = total_fcoin[i]/total_tb[i];
    //er_tf_ratio[i] = sqrt(total_tb[i]*tf_ratio*(1.-tf_ratio))/total_tb[i];

    //h_ctime_tb_RF_cut analysis
    //calc SN from fitting result
    double bg_min1 = -1.;
    double bg_max1 = fadd->GetParameter(1) - 3.*fadd->GetParameter(2);//mean - 3 * sigma
    double bg_min2 = fadd->GetParameter(1) + 3.*fadd->GetParameter(2);//mean + 3 * sigma
    double bg_max2 =  1.;
    double s_min   = bg_max1;
    double s_max   = bg_min2;
    double norm_bg = (bg_max2 - bg_min1)/((bg_max1 - bg_min1)+(bg_max2 - bg_min2));//normalize factor of SN
    double norm_s  = (bg_min2 - bg_max1)/((bg_max1 - bg_min1)+(bg_max2 - bg_min2));//normalize factor of SN

    area_BG[i] = fadd->GetParameter(3) * (bg_max2 - bg_min1) / bin_w;
    area_S[i]  = fadd->GetParameter(0) / bin_w;
    SN[i]      = area_S[i]/(area_BG[i]+area_S[i]);
    double total = fadd->Integral(bg_min1,bg_max2)/bin_w;

    //integral
    int bg_min1_bin = h_ctime_tb_RF_cut[i] ->FindBin(bg_min1);
    int bg_min2_bin = h_ctime_tb_RF_cut[i] ->FindBin(bg_min2);
    int bg_max1_bin = h_ctime_tb_RF_cut[i] ->FindBin(bg_max1);
    int bg_max2_bin = h_ctime_tb_RF_cut[i] ->FindBin(bg_max2);
    double BG_L = (double)h_ctime_tb_RF_cut[i] ->Integral(bg_min1_bin  ,bg_max1_bin  );
    double BG_R = (double)h_ctime_tb_RF_cut[i] ->Integral(bg_min2_bin  ,bg_max2_bin  );
    double S    = (double)h_ctime_tb_RF_cut[i] ->Integral(bg_max1_bin+1,bg_min2_bin-1);

    area_BGi[i] = (BG_L + BG_R) * norm_bg;
    area_Si[i]  = S - ((BG_L + BG_R) * norm_s);
    SNi[i]      = area_Si[i]/area_BGi[i];

    //cout<<Form("%.0lf %.0lf %.0lf %.0lf",total,area_BG[i]+area_S[i],area_BG[i],area_S[i])<<endl; 
    //cout<<Form("%.0lf %.0lf %.0lf %.0lf",area_BG[i],area_BGi[i],area_S[i],area_Si[i])<<endl; 

    tg_tbSN->SetPoint(i,i+1,SN[i]);
    tg_tbSN->SetPointError(i,0.,0.);

    tg_tbSNi->SetPoint(i,i+1,SNi[i]);
    tg_tbSNi->SetPointError(i,0.,0.);

    //h_ctime_tb_RF_f_or analysis
    //Fitting
    TF1 *fadd_f_or = new TF1("fadd_f_or",gaus_pol0bg,ctime_min,ctime_max,4);
    set->SetTF1(fadd_f_or,3,1,1);
    fadd_f_or -> SetParameter(0,h_ctime_tb_RF_f_or[i] ->Integral()*0.7);
    fadd_f_or -> SetParameter(1,0.5*(ctime_min+ctime_max));
    fadd_f_or -> SetParameter(2,0.12);
    fadd_f_or -> SetParameter(3,h_ctime_tb_RF_f_or[i] ->GetBinContent(h_ctime_tb_RF_f_or[i] ->FindBin(0.5*(ctime_min+ctime_max)-0.8)));

    for(int k=0;k<5;k++)h_ctime_tb_RF_f_or[i]  ->Fit(fadd_f_or,"QR","",ctime_min+0.1,ctime_max-0.1);
    //calc SN from fitting result
    bg_min1 = -1.;
    bg_max1 = fadd_f_or->GetParameter(1) - 3.*fadd_f_or->GetParameter(2);//mean - 3 * sigma
    bg_min2 = fadd_f_or->GetParameter(1) + 3.*fadd_f_or->GetParameter(2);//mean + 3 * sigma
    bg_max2 =  1.;
    s_min   = bg_max1;
    s_max   = bg_min2;
    norm_bg = (bg_max2 - bg_min1)/((bg_max1 - bg_min1)+(bg_max2 - bg_min2));//normalize factor of SN
    norm_s  = (bg_min2 - bg_max1)/((bg_max1 - bg_min1)+(bg_max2 - bg_min2));//normalize factor of SN

    area_BG_f_or[i] = fadd_f_or->GetParameter(3) * (bg_max2 - bg_min1) / bin_w;
    area_S_f_or[i]  = fadd_f_or->GetParameter(0) / bin_w;
    SN_f_or[i]      = area_S_f_or[i]/(area_BG_f_or[i]+area_S_f_or[i]);
    //cout<<"tagb "<<i+1<<" "<<SN_f_or[i]<<endl;
    tg_tbSN_f_or->SetPoint(i,i+1,SN_f_or[i]);
    tg_tbSN_f_or->SetPointError(i,0.,0.);
    
    //h_ctime_tb_RF_fbor analysis (TagB(x)TagB coin or TagB(x)TagF coin)
    //Fitting
    TF1 *fadd_fbor = new TF1("fadd_fbor",gaus_pol0bg,ctime_min,ctime_max,4);
    set->SetTF1(fadd_fbor,7,1,1);
    fadd_fbor -> SetParameter(0,h_ctime_tb_RF_fbor[i] ->Integral()*0.7);
    fadd_fbor -> SetParameter(1,0.5*(ctime_min+ctime_max));
    fadd_fbor -> SetParameter(2,0.12);
    fadd_fbor -> SetParameter(3,h_ctime_tb_RF_fbor[i] ->GetBinContent(h_ctime_tb_RF_fbor[i] ->FindBin(0.5*(ctime_min+ctime_max)-0.8)));

    for(int k=0;k<5;k++)h_ctime_tb_RF_fbor[i]  ->Fit(fadd_fbor,"QR","",ctime_min+0.1,ctime_max-0.1);
    //calc SN from fitting result
    bg_min1 = -1.;
    bg_max1 = fadd_fbor->GetParameter(1) - 3.*fadd_fbor->GetParameter(2);//mean - 3 * sigma
    bg_min2 = fadd_fbor->GetParameter(1) + 3.*fadd_fbor->GetParameter(2);//mean + 3 * sigma
    bg_max2 =  1.;
    s_min   = bg_max1;
    s_max   = bg_min2;
    norm_bg = (bg_max2 - bg_min1)/((bg_max1 - bg_min1)+(bg_max2 - bg_min2));//normalize factor of SN
    norm_s  = (bg_min2 - bg_max1)/((bg_max1 - bg_min1)+(bg_max2 - bg_min2));//normalize factor of SN

    area_BG_fbor[i] = fadd_fbor->GetParameter(3) * (bg_max2 - bg_min1) / bin_w;
    area_S_fbor[i]  = fadd_fbor->GetParameter(0) / bin_w;
    SN_fbor[i]      = area_S_fbor[i]/(area_BG_fbor[i]+area_S_fbor[i]);
    //cout<<"tagb "<<i+1<<" "<<SN_fbor[i]<<endl;
    
    //h_ctime_tb_RF_tbtb analysis
    //Fitting
    TF1 *fadd_tbtb = new TF1("fadd_tbtb",gaus_pol0bg,ctime_min,ctime_max,4);
    set->SetTF1(fadd_tbtb,9,1,1);
    fadd_tbtb -> SetParameter(0,h_ctime_tb_RF_tbtb[i] ->Integral()*0.7);
    fadd_tbtb -> SetParameter(1,0.5*(ctime_min+ctime_max));
    fadd_tbtb -> SetParameter(2,0.12);
    fadd_tbtb -> SetParameter(3,h_ctime_tb_RF_tbtb[i] ->GetBinContent(h_ctime_tb_RF_tbtb[i] ->FindBin(0.5*(ctime_min+ctime_max)-0.8)));

    for(int k=0;k<5;k++)h_ctime_tb_RF_tbtb[i]  ->Fit(fadd_tbtb,"QR","",ctime_min+0.1,ctime_max-0.1);
    //calc SN from fitting result
    bg_min1 = -1.;
    bg_max1 = fadd_tbtb->GetParameter(1) - 3.*fadd_tbtb->GetParameter(2);//mean - 3 * sigma
    bg_min2 = fadd_tbtb->GetParameter(1) + 3.*fadd_tbtb->GetParameter(2);//mean + 3 * sigma
    bg_max2 =  1.;
    s_min   = bg_max1;
    s_max   = bg_min2;
    norm_bg = (bg_max2 - bg_min1)/((bg_max1 - bg_min1)+(bg_max2 - bg_min2));//normalize factor of SN
    norm_s  = (bg_min2 - bg_max1)/((bg_max1 - bg_min1)+(bg_max2 - bg_min2));//normalize factor of SN

    area_BG_tbtb[i] = fadd_tbtb->GetParameter(3) * (bg_max2 - bg_min1) / bin_w;
    area_S_tbtb[i]  = fadd_tbtb->GetParameter(0) / bin_w;
    SN_tbtb[i]      = area_S_tbtb[i]/(area_BG_tbtb[i]+area_S_tbtb[i]);
    //cout<<"tagb "<<i+1<<" "<<SN_tbtb[i]<<endl;

    tg_tbSN_tbtb->SetPoint(i,i+1,SN_tbtb[i]);
    tg_tbSN_tbtb->SetPointError(i,0.,0.);

    tg_tbSN_fbor->SetPoint(i,i+1,SN_fbor[i]);
    tg_tbSN_fbor->SetPointError(i,0.,0.);

    tg_fborcoinratio->SetPoint(i,i+1,area_S_fbor[i]/area_S[i]);
    tg_fborcoinratio->SetPointError(i,0.,0.);

    tf_ratio[i]    = area_S_f_or[i]/area_S[i];
    er_tf_ratio[i] = 0.;
    tg_tbcoinratio->SetPoint(i,i+1,tf_ratio[i]);
    tg_tbcoinratio->SetPointError(i,0.,er_tf_ratio[i]);
  }
  cout<<"bin width "<<bin_w<<" ns"<<endl;

  TCanvas *c[NCanvas];
  for(int i=0;i<NCanvas;i++){
    c[i]=new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800,800);
  }


  for(int i=0;i<2;i++){ //Canvas
    c[i]->Divide(4,3);
    for(int k=0;k<12;k++){ //TagB seg
      int seg = 12*i + k +1;
      //h_ctime_tb_RF_f_or[i] = (TH1F*)ifp->Get(Form("h_ctime_tb_RF_f_or%d" ,i+1));
      c[i]->cd(1+k);gPad->SetLogy(1);h_ctime_tb_RF_cut[seg-1]->Draw("");
                                     h_ctime_tb_RF_f_or[seg-1]->Draw("same");
    }
  }

  leg_coin -> AddEntry(tg_tbcoinratio  ,"TagB#otimesTagF","p");
  leg_coin -> AddEntry(tg_fborcoinratio,"TagB#otimes(TagF or TagB)","p");
  c[2]->Divide(1,1);
  c[2]->cd(1);gPad->SetGrid(1,1);
              h2_frame1->Draw(); 
              tg_tbcoinratio->Draw("samePE");
              tg_fborcoinratio->Draw("samePE");
              leg_coin->Draw("same");

  for(int i=3;i<5;i++){ //Canvas
    c[i]->Divide(4,3);
    for(int k=0;k<12;k++){ //TagB seg
      int seg = 12*(i-3) + k +1;
      h_ctime_tb_RF_f_or_scale[seg-1] ->Scale(1./tf_ratio[seg-1]);
      c[i]->cd(1+k);gPad->SetLogy(1);h_ctime_tb_RF_cut[seg-1]->Draw("");
                                     //h_ctime_tb_RF_f_or_scale[seg-1]->Draw("same");
                                     h_ctime_tb_RF_fbor[seg-1]->Draw("same");
                                     h_ctime_tb_RF_tbtb[seg-1]->Draw("same");
    }
  }

  leg_SN -> AddEntry(tg_tbSN,"no cut","p");
  leg_SN -> AddEntry(tg_tbSN_tbtb,"TagB(odd)#otimesTagB(even)","p");
  leg_SN -> AddEntry(tg_tbSN_f_or,"TagB#otimesTagF","p");
  leg_SN -> AddEntry(tg_tbSN_fbor,"TagB#otimes(TagF or TagB)","p");
  c[5]->Divide(1,1);
  c[5]->cd(1);gPad->SetGrid(1,1);
              h2_frame2->Draw(); 
              tg_tbSN->Draw("samePE");//tg_tbSNi->Draw("samePE");
              tg_tbSN_tbtb->Draw("samePE");//tg_tbSNi_tbtb->Draw("samePE");
              tg_tbSN_f_or->Draw("samePE");//tg_tbSNi_tbtb->Draw("samePE");
              tg_tbSN_fbor->Draw("samePE");//tg_tbSNi_tbtb->Draw("samePE");
              leg_SN -> Draw("same");
 
  //save pdf
  string pdf_file = root_file;
  pdf_file.erase(pdf_file.size()-5);
  pdf_file.append("_tf_coin.pdf");
  c[0]->Print(Form("%s[",pdf_file.c_str()));
  for(int i=0;i<NCanvas;i++)  c[i]->Print(Form("%s",pdf_file.c_str()));
  c[NCanvas-1]->Print(Form("%s]",pdf_file.c_str()));
}
