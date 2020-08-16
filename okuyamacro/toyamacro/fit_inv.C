using namespace std;

void SetTH1(TH1F *h1, TString hname, TString xname, TString yname, int LColor, int FStyle, int FColor){
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

//____________________________________________________________________________________________

void fit_inv(){
gStyle->SetOptStat(0);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.06);
gStyle->SetPadLeftMargin(0.16);
gStyle->SetPadBottomMargin(0.15);
int runnum1 =1000;
TString ifname1=Form("./root/RKV_ana/ppi_g1to9.root");//
//TString ifname1=Form("./root/RKV_ana/ana_g1to9.root");//
//TString ifname1=Form("./root/RKV_ana/ana_g1to8.root");//
//TString ifname1=Form("./root/RKV_ana_pm/all.root");//
//TString ifname1=Form("./root/all_test.root");//
//TString ifname1=Form("./root/all.root");//
//TString ifname1=Form("./root/ana10226_10350.root");//
//TString ifname1=Form("./root/test.root");//
  TFile *ifp1  = new TFile(ifname1  );
      cout<<"input filename : "<<ifname1<<endl;

TH1F *h1_inv_nocut;//H3L mom
h1_inv_nocut = (TH1F*)ifp1->Get("vertex/h_invm_lpi");
//h1_inv_nocut = (TH1F*)ifp1->Get("vertex/h_inv_mass");
h1_inv_nocut ->GetXaxis()->SetRangeUser(1.,1.3);
SetTH1(h1_inv_nocut,"Invariant Mass","Invariant Mass[GeV/#it{c}^{2}]","Counts/2.5MeV/#it{c}^{2}",2,3000,0);
h1_inv_nocut -> SetFillColor(0);
h1_inv_nocut -> SetLineWidth(2);
TH1F* h1_inv_rebin2 = h1_inv_nocut -> Clone();
TH1F* h1_inv_rebin3 = h1_inv_nocut -> Clone();
h1_inv_rebin2 -> Rebin(2);
h1_inv_rebin3 -> Rebin(3);
h1_inv_rebin2 -> GetYaxis()->SetTitle("Counts/5.0MeV/#it{c}^{2}");
h1_inv_rebin3 -> GetYaxis()->SetTitle("Counts/7.5MeV/#it{c}^{2}");
double param_bg[2];
double param[5], para_error[5];
double min =1.09;
double max =1.15;
TF1 *f5 = new TF1("f5","pol0",min,max);
//TF1 *f5 = new TF1("f5","pol1",min,max);
f5->SetParameter(0,100);
//f5->SetParameter(1,-10);
    f5 ->SetLineWidth(2);  f5 ->SetLineColor(2);f5 ->SetLineStyle(3);
h1_inv_nocut ->Fit(f5,"0QR","",min,max+0.05);
f5 -> GetParameters(&param_bg[0]);
TLine *tline = new TLine(min,0,min,900);
tline->SetLineWidth(2.5);
tline->SetLineColor(1);
tline->SetLineStyle(1);
//TF1 *fadd = new TF1("fadd",gaus_pol1bg,min,max,5);
TF1 *fadd = new TF1("fadd",gaus_pol0bg,min,max,4);
fadd -> SetParameter(0,150.);
fadd -> SetParameter(1,1.115);
fadd -> SetParameter(2,0.004);
fadd -> SetParameter(3,param_bg[0]);
//fadd -> SetParameter(4,param_bg[1]);
//fadd -> FixParameter(3,param_bg[0]);
//fadd -> FixParameter(4,param_bg[1]);
    fadd ->SetLineWidth(2);  fadd ->SetLineColor(1);fadd ->SetLineStyle(1);fadd ->SetNpx(1000);
h1_inv_nocut->Fit(fadd,"0QR","",min,max);
fadd ->GetParameters(&param[0]);
double chi2 = fadd->GetChisquare();
double ndf  = fadd->GetNDF();
double g_area  = param[0];
double g_mean  = param[1];
double g_sigma = param[2];
double bg_p0   = param[3];
//double bg_p1   = param[4];
cout<<"BG param(2nd step) ="<<bg_p0<<" "<<endl;//bg_p1<<endl;

double erg_area =fadd ->GetParError(0);
double erg_mean =fadd ->GetParError(1);
double erg_sigma=fadd ->GetParError(2);
double erbg_p0  =fadd ->GetParError(3);
//double erbg_p1  =fadd ->GetParError(4);

TF1 *ga2 = new TF1("ga2","gausn",min,max);
    ga2 ->SetLineWidth(2);  ga2 ->SetLineColor(2);ga2 ->SetLineStyle(2);
ga2 -> SetParameter(0, g_area);
ga2 -> SetParameter(1, g_mean);
ga2 -> SetParameter(2, g_sigma);

TF1 *f6 = new TF1("f6","pol0",min,max);//to show fitting result
    f6 ->SetLineWidth(3);  f6 ->SetLineColor(4);f6 ->SetLineStyle(2);
f6->SetParameter(0,bg_p0);
//f6->SetParameter(1,bg_p1);

double x_min=g_mean-3.*g_sigma;
double x_max=g_mean+3.*g_sigma;
//double bg_3s=(bg_p1/2.*(x_max*x_max-x_min*x_min)+bg_p0*(x_max-x_min))/0.0025;
double bg_3s=(bg_p0*(x_max-x_min))/0.0020;

int p_flag=0;
cout<<"Do you want to save canvas?"<<endl;
cout<<"yes :1, no :0"<<endl;
cin>>p_flag;
  TH2F *h_frame = new TH2F("h_frame","h_frame",10,0,1,10,0,1);
  SetTH2(h_frame,"","","");
  h_frame->GetXaxis()->SetNdivisions(000);
  h_frame->GetYaxis()->SetNdivisions(000);
  h_frame->SetStats(0);

TCanvas *c[2];
for(int i=0;i<2;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800/2,1500/2);}
c[0]->Clear();
c[0]->Divide(2,2);
c[0]->cd(1);
gPad->SetLogy(0);h1_inv_nocut ->Draw("");
c[0]->cd(2);
gPad->SetLogy(0);h1_inv_nocut ->Draw("PE");f5->Draw("same");
c[0]->cd(3);
gPad->SetLogy(0);h1_inv_nocut ->Draw("");  fadd->Draw("same");ga2->Draw("same");f6->Draw("same");
c[0]->cd(4);
gPad->SetLogy(0);h1_inv_nocut ->Draw("PE");fadd->Draw("same");ga2->Draw("same");f6->Draw("same");

c[1]->Clear();
c[1]->Divide(2,2);
c[1]->cd(1);
TH1F *h_zoom = (TH1F*)h1_inv_nocut->Clone("h_zoom");
h_zoom->GetXaxis()->SetRangeUser(1.10,1.13);
gPad->SetLogy(0);h_zoom->Draw("PE");fadd->Draw("same");ga2->Draw("same");f6->Draw("same");
c[1]->cd(2);h1_inv_rebin2 -> Draw("");
c[1]->cd(4);h1_inv_rebin3 -> Draw("");
  c[1] -> cd(3);
h_frame->Draw("");
TLatex *tex  = new TLatex(0.5,0.7,"");
tex  -> SetTextSize(0.060);
tex  -> SetTextAlign(22);
                tex->DrawLatex(0.5,0.9,Form("mean: %.01lf#pm%.01lf[MeV/#it{c}^{2}]"          ,g_mean*1000.,erg_mean*1000.            ) );
                tex->DrawLatex(0.5,0.7,Form("Num. of #Lambda : %.01lf"                       ,g_area/0.0020                          ) );
                tex->DrawLatex(0.5,0.5,Form("width(#sigma): %.01lf#pm%.01lf[MeV/#it{c}^{2}]" ,g_sigma*1000.                          ) );
                tex->DrawLatex(0.5,0.3,Form("B.G. (#pm 3#sigma): %.01lf"                     ,bg_3s                                  ) );
                //tex->DrawLatex(0.5,0.1,Form("S/#sqrt{S+N}: %.01lf"                           ,g_area/0.0025/sqrt(g_area/0.0025+bg_3s)) );
				tex->DrawLatex(0.5,0.1,Form("#chi^{2}/ndf= %.1lf / %d = %.3lf"                           ,fadd->GetChisquare(),fadd->GetNDF(), fadd->GetChisquare()/fadd->GetNDF()));

if(p_flag==1){
string ofname_pdf="inv_30Jan";
  c[0]->Print(Form("./pdf/fit_%s.pdf[",ofname_pdf.c_str())  );
  c[0]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  c[1]->Print(Form("./pdf/fit_%s.pdf" ,ofname_pdf.c_str())  );
  c[1]->Print(Form("./pdf/fit_%s.pdf]",ofname_pdf.c_str())  );
 }

} 
