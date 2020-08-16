using namespace std;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetTH1(TH1F *h1, TString hname, TString xname, TString yname, int LColor, int FStyle, int FColor){
  h1->SetTitle(hname);
  h1->GetXaxis()->SetTitle(xname);
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->SetTitle(yname);
  h1->GetYaxis()->CenterTitle();
  //h1->SetMinimum(0.0);
  h1->SetLineColor(LColor);
  h1->SetFillStyle(FStyle);
  h1->SetFillColor(FColor);
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->SetLineWidth(0);
  h1->SetTitleSize(0.05,"");
  h1->SetTitleSize(0.05,"x");
  h1->SetTitleSize(0.05,"y");
  h1->SetLabelSize(0.05,"x");
  h1->SetLabelSize(0.05,"y");
  ((TGaxis*)h1->GetYaxis())->SetMaxDigits(3);
}
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
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gaus_pol2bg(double *x, double *par) {
  //par[0]=area (gaus)
  //par[1]=mean (gaus)
  //par[2]=sigma (gaus)
  //par[3]=pol const
  //par[4]=pol tilt
  //par[5]=pol tilt
  double val;
  double bg;
  double ga;

  ga=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
  bg =  par[3]*pow(x[0]-par[4],2.) + par[5];
  val = ga + bg;
  return val;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gaus_pol3bg(double *x, double *par) {
  //par[0]=area (gaus)
  //par[1]=mean (gaus)
  //par[2]=sigma (gaus)
  //par[3]=pol tilt
  //par[4]=pol tilt
  //par[5]=pol tilt
  //par[6]=pol const
  double val;
  double bg;
  double ga;

  ga=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
  bg =  par[3]*pow(x[0]-par[4],3.) + par[5]*(x[0]-par[4])+par[6];
  val = ga + bg;
  return val;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double get_1st_mom(TH1 *h,double low_th = 0.,double high_th = 0.){
  double rel_freq;//soutai dosuu
  double xval,mom=0;
  int first_bin =h->FindFirstBinAbove(low_th); 
  int last_bin  =h->FindLastBinAbove(high_th);
  int total_event = h->Integral(first_bin,last_bin);
  for(int i=first_bin;i<last_bin;i++){
    rel_freq = (double)h->GetBinContent(i)/total_event;
    xval = h->GetBinCenter(i);
    mom += xval*rel_freq;
  }
  return mom;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          main
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void lifetime_ana(){
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.15);

  TString ifname1=Form("./root/RKV_ana/ppi_g1to9.root");//
  TFile *ifp1  = new TFile(ifname1  );
  cout<<"input filename : "<<ifname1<<endl;
  TH1F *h1_inv;
  h1_inv = (TH1F*)ifp1->Get("vertex/h_invm_lpi");
  h1_inv ->GetXaxis()->SetRangeUser(1.,1.3);
  SetTH1(h1_inv,"Invariant Mass","Invariant Mass[GeV/#it{c}^{2}]","Counts/2.0MeV/#it{c}^{2}",2,3000,0);
  h1_inv -> SetFillColor(0);
  h1_inv -> SetLineWidth(2);

  TH2F *h_frame = new TH2F("h_frame","h_frame",10,0,1,10,0,1);
  SetTH2(h_frame,"","","");
  h_frame->GetXaxis()->SetNdivisions(000);
  h_frame->GetYaxis()->SetNdivisions(000);
  h_frame->SetStats(0);

  double min =1.09;
  double max =1.15;
  double param_bg[4];
  double param[7];
  //TF1 *f_bg = new TF1("f_bg","pol0",1.14,1.18);
  //f_bg -> SetParameter(0,100);
  //TF1 *f_bg = new TF1("f_bg","[0]*(x-[1])*(x-[1])+[2]",1.00,1.28);
  TF1 *f_bg = new TF1("f_bg","[0]*(x-[1])*(x-[1])*(x-[1])+[2]*(x-[1])+[3]",1.00,1.28);
  f_bg -> SetParameter(0,100);
  f_bg -> SetParameter(1,1.3);
  //f_bg -> FixParameter(1,1.3);
  f_bg -> SetParameter(2,-10);
  f_bg -> SetParameter(3,1);
  h1_inv ->Fit(f_bg,"0QR","",1.06,1.30);
  f_bg -> GetParameters(&param_bg[0]);

  TF1 *fadd = new TF1("fadd",gaus_pol3bg,min,max,7);
  fadd -> SetParameter(0,150.);
  fadd -> SetParameter(1,1.115);
  fadd -> SetParameter(2,0.004);
  fadd -> SetParameter(3,param_bg[0]);
  fadd -> SetParameter(4,param_bg[1]);
  fadd -> SetParameter(5,param_bg[2]);
  fadd -> SetParameter(6,param_bg[3]);
  fadd ->SetLineWidth(2);  fadd ->SetLineColor(1);fadd ->SetLineStyle(1);fadd ->SetNpx(1000);
  h1_inv->Fit(fadd,"0QR","",min,max);
  fadd ->GetParameters(&param[0]);
  double chi2    = fadd->GetChisquare();
  double ndf     = fadd->GetNDF();
  double g_area  = param[0];
  double g_mean  = param[1];
  double g_sigma = param[2];
  double bg_p0   = param[3];
  double bg_p1   = param[4];
  double bg_p2   = param[5];
  double bg_p3   = param[6];
  double erg_area =fadd ->GetParError(0);
  double erg_mean =fadd ->GetParError(1);
  double erg_sigma=fadd ->GetParError(2);
  double erbg_p0  =fadd ->GetParError(3);
  TF1 *ga_inv = new TF1("ga_inv","gausn",min,max);
  ga_inv ->SetLineWidth(2);  ga_inv ->SetLineColor(2);ga_inv ->SetLineStyle(2);
  ga_inv -> SetParameter(0, g_area);  ga_inv -> SetParameter(1, g_mean);  ga_inv -> SetParameter(2, g_sigma);

  f_bg->SetParameter(0,bg_p0);
  f_bg->SetParameter(1,bg_p1);
  f_bg->SetParameter(2,bg_p2);
  f_bg->SetParameter(3,bg_p3);
  double x_min=g_mean-3.*g_sigma;
  double x_max=g_mean+3.*g_sigma;
  double bg_3s=(bg_p0*(x_max-x_min))/0.0020;


  //////////////
  //decay time//
  //////////////
  TH1F *h_cor, *h_BG_cor;//w/ flight time correction
  TH1F *h_BG_norm[5], *h_cor_sub[5];//w/ flight time correction after substruction B.G.
  h_cor = (TH1F*)ifp1->Get("cointime/h_ctime_Lam_cor");
  h_cor ->Rebin(10);
  h_cor->SetStats(kFALSE);
  h_cor->SetMinimum(0.1);
  h_cor->SetLineWidth(3);

  h_BG_cor = (TH1F*)ifp1->Get("pipi/h_ctime_pipi_cor");
  h_BG_cor ->Rebin(1);
  h_BG_cor->SetStats(kFALSE);
  //h_BG_cor->SetLineWidth(3);
  TF1 *f_res = new TF1("f_res","gausn",-1,2);
  h_BG_cor->Fit(f_res,"QR","",-1,1);
  double res_mean  = f_res->GetParameter(1);
  double res_sigma = f_res->GetParameter(2);
  f_res ->SetParameter(0,bg_3s);

  h_cor_sub[0] = (TH1F*)h_cor->Clone();
  h_BG_norm[0] = (TH1F*)h_cor->Clone();
  h_BG_norm[0] ->SetLineColor(4);
  for(int i=1;i<h_cor->FindLastBinAbove(-1.);i++){
    double low_edge  = h_cor->GetXaxis()->GetBinLowEdge(i);
    double up_edge   = h_cor->GetXaxis()->GetBinUpEdge(i);
    double content   = f_res->Integral(low_edge,up_edge);
    h_BG_norm[0] ->SetBinContent(i, content);
  }
  h_cor_sub[0]->Add(h_BG_norm[0],-1.);
  double mom1st= get_1st_mom(h_cor_sub[0],2.5, 0.5);
  cout<<"mom1st:"<<mom1st<<"  mean:"<<h_cor_sub[0]->GetMean()<<endl;
  

  SetTH1(h_cor_sub[0],"Decay time #Lambda (after B.G. subtraction)","decay time[ns]","Counts",4,3000,0);
  SetTH1(h_BG_norm[0],"Decay time B.G. "                           ,"decay time[ns]","Counts",1,3000,0);
  ////////
  //Draw//
  ////////
  TCanvas *c[2];
  for(int i=0;i<2;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800/2,1500/2);}
  c[0]->Clear();
  c[0]->Divide(2,2);
  c[0]->cd(1);
  h_frame->Draw("");
  TLatex *tex  = new TLatex(0.5,0.7,"");
  tex  -> SetTextSize(0.060);
  tex  -> SetTextAlign(22);
  tex->DrawLatex(0.5,0.9,Form("mean: %.01lf#pm%.01lf[MeV/#it{c}^{2}]"          ,g_mean*1000.,erg_mean*1000.            ) );
  tex->DrawLatex(0.5,0.7,Form("Num. of #Lambda : %.01lf#pm%.01lf"              ,g_area/0.0020,  erg_area/0.0020                      ) );
  tex->DrawLatex(0.5,0.5,Form("width(#sigma): %.01lf#pm%.01lf[MeV/#it{c}^{2}]" ,g_sigma*1000.                          ) );
  tex->DrawLatex(0.5,0.3,Form("B.G. (#pm 3#sigma): %.01lf"                     ,bg_3s                                  ) );
  //tex->DrawLatex(0.5,0.1,Form("S/#sqrt{S+N}: %.01lf"                           ,g_area/0.0025/sqrt(g_area/0.0025+bg_3s)) );
  tex->DrawLatex(0.5,0.1,Form("#chi^{2}/ndf= %.1lf / %d = %.3lf"                           ,fadd->GetChisquare(),fadd->GetNDF(), fadd->GetChisquare()/fadd->GetNDF()));
  c[0]->cd(2);
  gPad->SetLogy(0);h1_inv ->Draw("");f_bg->Draw("same");
  c[0]->cd(3);
  gPad->SetLogy(0);h1_inv ->Draw("");  fadd->Draw("same");ga_inv->Draw("same");f_bg->Draw("same");
  c[0]->cd(4);
  TH1F *h_zoom = (TH1F*)h1_inv->Clone("h_zoom");
  h_zoom->GetXaxis()->SetRangeUser(1.10,1.13);
  gPad->SetLogy(0);h_zoom->Draw("PE");fadd->Draw("same");ga_inv->Draw("same");f_bg->Draw("same");

  c[1]->Clear();
  c[1]->Divide(2,2);
  c[1]->cd(1);gPad->SetLogy(1);h_cor       ->Draw();h_BG_norm[0]->Draw("same");
  c[1]->cd(2);gPad->SetLogy(1);h_BG_cor    ->Draw();
  c[1]->cd(3);gPad->SetLogy(1);h_cor_sub[0]->Draw();
  //h_cor_sub[0]->Fit("expo","R","",0.1,1.0);
  c[1]->cd(4);gPad->SetLogy(1);h_BG_norm[0]->Draw();

  c[0]->Print("pdf/life_ana.pdf[");
  c[0]->Print("pdf/life_ana.pdf");
  c[1]->Print("pdf/life_ana.pdf");
  c[1]->Print("pdf/life_ana.pdf]");
}
