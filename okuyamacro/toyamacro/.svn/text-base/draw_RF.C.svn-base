#include "Settings.cc"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void draw_RF(){
gStyle->SetOptStat(0);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.13);
gStyle->SetPadLeftMargin(0.16);
gStyle->SetPadBottomMargin(0.15);

//TString ifname1=Form("./root/test.root");//
TString ifname1=Form("./root/RF_ana.root");//
  TFile *ifp1  = new TFile(ifname1  );
      cout<<"input filename1 : "<<ifname1<<endl;
Settings *set=new Settings();

TH1F *h_time_RFdiff[4], *h_frame;
TF1 *f_gaus[4],*f_lin;
double interval[4],resi[4],e_resi[4];
double mean[4],sigma[4];
double e_mean[4],e_sigma[4];
double max_y;
  for(int i=0;i<4;i++){
    interval[i]=165.950*(i+1);
    h_time_RFdiff[i]   = (TH1F*)ifp1->Get(Form("h_time_RFdiff%d",i+1) );
    set->SetTH1(h_time_RFdiff[i] ,Form("RF 1st - %d",i+2)            ,"diff. time[ns]"      ,"Counts/25ps"       , 1, 3000, 0);
    h_time_RFdiff[i]->SetStats(0);
    f_gaus[i] = new TF1(Form("f_gaus%d",i+1),"gaus", interval[i]-0.5, interval[i]+0.5);
    set->SetTF1(f_gaus[i],2,1,1.0 );
    h_time_RFdiff[i]->Fit(f_gaus[i],"QR","");

    mean[i]  = f_gaus[i]->GetParameter(1);
    sigma[i] = 1000.*f_gaus[i]->GetParameter(2);
    e_mean[i]  = f_gaus[i]->GetParError(1);
    e_sigma[i] = 1000.*f_gaus[i]->GetParError(2);

    resi[i]=1000.*(mean[i]-interval[i]);//ns->ps
    e_resi[i]=1000.*e_mean[i];//ns->ps
  }

TGraphErrors *gr_resi, *gr_mean, *gr_sigma;
gr_resi  = new TGraphErrors(4, interval,  resi, 0,  e_resi);
gr_mean  = new TGraphErrors(4, interval,  mean, 0,  e_mean);
gr_sigma = new TGraphErrors(4, interval, sigma, 0, e_sigma);
set->SetGrErr(gr_resi ,"residual","interval[ns]","residual[ps]", 1, 2, 34,1.0,0,0);
set->SetGrErr(gr_mean ,"mean"    ,"interval[ns]","mean[ns]"    , 1, 4, 22,1.0,0,0);
set->SetGrErr(gr_sigma,"sigma"   ,"interval[ns]","sigma[ps]"   , 1, 2, 33,1.0,0,0);
gr_resi ->SetMarkerSize(1.5);
gr_mean ->SetMarkerSize(1.5);
gr_sigma->SetMarkerSize(1.5);


 f_lin = new TF1("f_lin","[0]*x",0,800);
 set->SetTF1(f_lin,2,1,1.0 );
 gr_mean->Fit(f_lin,"R","");
double p0  = f_lin->GetParameter(0);
double e_p0= f_lin->GetParError(0);

TLatex *tex = new TLatex(0,0,"aaa");
tex  -> SetTextSize(0.060);
tex  -> SetTextFont(42);
tex  -> SetTextAlign(21);

  TArrow *arrow = new TArrow(0,0,0,0,0.01,"");//x1,y1,x2,y2,arrow size,option
  arrow->SetAngle(45);
  arrow->SetLineColor(1);
  arrow->SetLineWidth(3);

  TLine *line = new TLine(0,0,0,0);//x1,y1,x2,y2
  line->SetLineColor(1);
  line->SetLineWidth(1);
  line->SetLineStyle(2);


TCanvas *c[4];
for(int i=0;i<4;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800/2,1500/3);}
c[0]->Clear();
c[0]->Divide(2,2);
for(int i=0;i<4;i++){
  c[0]->cd(i+1);gPad->SetLogy(1);h_time_RFdiff[i]->Draw("");
  max_y = h_time_RFdiff[i]->GetBinContent(h_time_RFdiff[i]->GetMaximumBin());
                tex ->DrawLatex(interval[i]+0.2,0.90*max_y,Form("mean : %.4lf#pm%.4lf",mean[i],e_mean[i]));
}

c[1]->Clear();
gPad->SetLogy(0);
TH1F*frame_resi = gPad->DrawFrame(0,-11,700,11);
set->SetTH1(frame_resi,"","interval[ns]","residual[ps]", 1, 1000, 0);
gr_resi->Draw("P");
line->DrawLine(0,0,700,0);

c[2]->Clear();
TH1F*frame_mean = gPad->DrawFrame(0,0,700,700);
set->SetTH1(frame_mean,"","interval[ns]","mean[ns]", 1, 1000, 0);
gr_mean->Draw("P");
  tex  ->DrawLatex(300,600,Form("p0=%.5lf#pm%.5lf",p0,e_p0));
  //arrow->DrawArrow(0.001*ML,5,0.001*ML,0.5,0.02,"|>");
c[3]->Clear();
TH1F*frame_sigma = gPad->DrawFrame(0,0,700,50);
set->SetTH1(frame_sigma,"","interval[ns]","#sigma[ps]", 1, 1000, 0);
gr_sigma->Draw("P");

#if 1
  c[0]->Print("pdf/RF_diff.pdf[");
  c[0]->Print("pdf/RF_diff.pdf");
  c[1]->Print("pdf/RF_diff.pdf");
  c[2]->Print("pdf/RF_diff.pdf");
  c[3]->Print("pdf/RF_diff.pdf");
  c[3]->Print("pdf/RF_diff.pdf]");
#endif
}
