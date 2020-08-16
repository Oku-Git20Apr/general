using namespace std;
#include "Settings.cc"

void SetTH1org(TH1F *h1, TString hname, TString xname, TString yname, int LColor, int FStyle, int FColor){
  h1->SetTitle(hname);
  h1->GetXaxis()->SetTitle(xname);
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->SetTitle(yname);
  h1->GetYaxis()->CenterTitle();
  h1->SetMinimum(0.0);
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
  ((TGaxis*)h1->GetXaxis())->SetMaxDigits(3);
}
////////////////
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gaus_pol1bg(double *x, double *par) {
  //par[0]=area (gaus)
  //par[1]=mean (gaus)
  //par[2]=sigma (gaus)
  //par[3]=pol const
  //par[4]=pol tilt
  double val;
  double bg;
  double ga;

  ga=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
  bg =  par[3]+ x[0]*par[4] ;
  val = ga + bg;
  return val;
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
void fit_gaus_pol3bg(TH1F*h,double &area, double &area_error, double min, double max){
    TF1 *f = new TF1("f",gaus_pol3bg,min,max,7);
    f -> SetLineWidth(1);
    f -> SetLineColor(1);
    f -> SetLineStyle(1);
    f -> SetParameter(0,area);
    f -> SetParameter(1,1.115);
    f -> SetParameter(2,0.003);
    f -> SetParameter(3, 1.54848e+04);
    f -> SetParameter(4, 1.20643e+00);
    f -> SetParameter(5,-1.86645e+02);
    f -> SetParameter(6, 8.11975e+00);
    h -> Fit(f,"QR","",min,max);
    area    = f ->GetParameter(0);
    area_error = f ->GetParError(0);
}
//____________________________________________________________________________________________

void fit_inv_bgpol3(){
  gStyle->SetOptStat("ie");
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TString ifname1=Form("./root/RKV_ana_acci/true.root");//
  //TString ifname1=Form("./root/RKV_ana/ana_g1to9.root");//
  TFile *ifp1  = new TFile(ifname1  );
  Settings* set = new Settings();
  cout<<"input filename : "<<ifname1<<endl;

  double MassRange_Lam_min   =  1.105;
  double MassRange_Lam_max   =  1.125;       
  double param_bg[4];
  double param[7], para_error[7];
  double min =1.09;
  double max =1.25;
  double GeVperBin = 0.002;//Gev/bin
  TH1F *h1_inv_nocut;//
  h1_inv_nocut = (TH1F*)ifp1->Get("vertex/h_invm_lpi");
  h1_inv_nocut ->GetXaxis()->SetRangeUser(1.,1.3);
  SetTH1org(h1_inv_nocut,"Invariant Mass","Invariant Mass[GeV/#it{c}^{2}]",Form("Counts/%.1lfMeV/#it{c}^{2}",1000.*GeVperBin),1,3000,0);
  h1_inv_nocut -> SetFillColor(0); h1_inv_nocut -> SetLineWidth(1);

  TF1 *f_bg = new TF1("f_bg","[0]*(x-[1])*(x-[1])*(x-[1])+[2]*(x-[1])+[3]",1.00,1.28);
  f_bg -> SetParameter(0,100);
  f_bg -> SetParameter(1,1.3);
  //f_bg -> FixParameter(1,1.3);
  f_bg -> SetParameter(2,-10);
  f_bg -> SetParameter(3,1);
  h1_inv_nocut ->Fit(f_bg,"0QR","",1.06,1.30);
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
  h1_inv_nocut->Fit(fadd,"0QR","",min,max);
  fadd ->GetParameters(&param[0]);
  double chi2 = fadd->GetChisquare();
  double ndf  = fadd->GetNDF();
  double g_area  = param[0];
  double g_mean  = param[1];
  double g_sigma = param[2];
  double bg_p0   = param[3];
  double bg_p1   = param[4];
  double bg_p2   = param[5];
  double bg_p3   = param[6];
  cout<<g_area/GeVperBin<<endl;
  //cout<<"BG param(2nd step) ="<<bg_p0<<" "<<endl;//bg_p1<<endl;

  double erg_area =fadd ->GetParError(0);
  double erg_mean =fadd ->GetParError(1);
  double erg_sigma=fadd ->GetParError(2);
  double erbg_p0  =fadd ->GetParError(3);
  //double erbg_p1  =fadd ->GetParError(4);

  TF1 *ga = new TF1("ga","gausn",min,max);
  ga ->SetLineWidth(2);  ga ->SetLineColor(2);ga ->SetLineStyle(2);
  ga -> SetParameter(0, g_area);
  ga -> SetParameter(1, g_mean);
  ga -> SetParameter(2, g_sigma);
  double NLambda = (ga->Integral(MassRange_Lam_min,MassRange_Lam_max))/GeVperBin;
  

  f_bg->SetParameter(0,bg_p0);
  f_bg->SetParameter(1,bg_p1);
  f_bg->SetParameter(2,bg_p2);
  f_bg->SetParameter(3,bg_p3);

  double x_min=g_mean-3.*g_sigma;
  double x_max=g_mean+3.*g_sigma;
  double bg_3s=(f_bg->Integral(MassRange_Lam_min,MassRange_Lam_max))/GeVperBin;
  cout<<"Region of Lambda "<<MassRange_Lam_min<<" - "<<MassRange_Lam_max<<endl;
  cout<<"NLambda:"<<NLambda<<", BG:"<<bg_3s<<endl;
  //double bg_3s=(bg_p0*(x_max-x_min))/GeVperBin;

  int p_flag=0;
  cout<<"Do you want to save canvas?"<<endl;
  cout<<"yes :1, no :0"<<endl;
  cin>>p_flag;

////////////

  TH2F *h_frame = new TH2F("h_frame","h_frame",10,0,1,10,0,1);
  SetTH2(h_frame,"","","");
  h_frame->GetXaxis()->SetNdivisions(000);
  h_frame->GetYaxis()->SetNdivisions(000);
  h_frame->SetStats(0);

TCanvas *c[7];
for(int i=0;i<7;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800/2,1500/2);}
c[0]->Clear();
c[0]->Divide(2,2);
c[0]->cd(1);
gPad->SetLogy(0);h1_inv_nocut ->Draw("");
c[0]->cd(2);
gPad->SetLogy(0);h1_inv_nocut ->Draw("PE");f_bg->Draw("same");
c[0]->cd(3);
gPad->SetLogy(0);h1_inv_nocut ->Draw("");  fadd->Draw("same");ga->Draw("same");f_bg->Draw("same");
c[0]->cd(4);
gPad->SetLogy(0);h1_inv_nocut ->Draw("PE");fadd->Draw("same");ga->Draw("same");f_bg->Draw("same");

c[1]->Clear();
c[1]->Divide(2,2);
c[1]->cd(1);
TH1F *h_zoom = (TH1F*)h1_inv_nocut->Clone("h_zoom");
h_zoom->GetXaxis()->SetNdivisions(505);
h_zoom->GetXaxis()->SetRangeUser(x_min,x_max);
gPad->SetLogy(0);h_zoom->Draw("PE");fadd->Draw("same");ga->Draw("same");f_bg->Draw("same");
  c[1] -> cd(2);
h_frame->Draw("");
TLatex *tex  = new TLatex(0.5,0.7,"");
tex  -> SetTextSize(0.060);
tex  -> SetTextAlign(22);
                tex->DrawLatex(0.5,0.9,Form("mean: %.01lf#pm%.01lf[MeV/#it{c}^{2}]"          ,g_mean*1000.,erg_mean*1000.            ) );
                tex->DrawLatex(0.5,0.7,Form("Num. of #Lambda : %.01lf#pm%.01lf"              ,g_area/GeVperBin, erg_area/GeVperBin   ) );
                tex->DrawLatex(0.5,0.5,Form("width(#sigma): %.01lf#pm%.01lf[MeV/#it{c}^{2}]" ,g_sigma*1000.                          ) );
                tex->DrawLatex(0.5,0.3,Form("B.G. (#pm 3#sigma): %.01lf"                     ,bg_3s                                  ) );
                //tex->DrawLatex(0.5,0.1,Form("S/#sqrt{S+N}: %.01lf"                           ,g_area/0.0025/sqrt(g_area/0.0025+bg_3s)) );
				tex->DrawLatex(0.5,0.1,Form("#chi^{2}/ndf= %.1lf / %d = %.3lf"                           ,fadd->GetChisquare(),fadd->GetNDF(), fadd->GetChisquare()/fadd->GetNDF()));

  c[2]->Clear();
  c[2]->Divide(4,3);
  for(int i=0;i<12;i++){
   c[2]->cd(i+1);h_invm_chi[i]->Draw();
  }
  
  c[3]->Clear();
  c[3]->Divide(4,3);
  for(int i=0;i<12;i++){
   c[3]->cd(i+1);h_invm_dca[i]->Draw();
  }

  c[4]->Clear();
  c[4]->Divide(4,3);
  for(int i=0;i<12;i++){
   c[4]->cd(i+1);h_invm_mm[i]->Draw();
  }

  c[5]->Clear();
  c[5]->Divide(4,3);
  for(int i=0;i<12;i++){
   c[5]->cd(i+1);h_invm_p_min[i]->Draw();
  }


  TLine *tline = new TLine(0,0,0,0);
  tline ->SetLineWidth(1.5);
  tline ->SetLineColor(6);
  tline ->SetLineStyle(2);

if(p_flag==1){
string ofname_pdf="evalu_chi2";
  c[0]->Print(Form("./pdf/fit_%s.pdf[",ofname_pdf.c_str())  );
  for(int i=0;i<7;i++){
    c[i]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  }
  c[6]->Print(Form("./pdf/fit_%s.pdf]",ofname_pdf.c_str())  );
 }

} 
