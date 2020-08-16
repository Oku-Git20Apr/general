#include "Settings.cc"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void draw_YONUPA(){
gStyle->SetOptStat(0);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.13);
gStyle->SetPadLeftMargin(0.16);
gStyle->SetPadBottomMargin(0.15);

//TString ifname1=Form("./root/test.root");//
TString ifname1=Form("./root/rkver_test17Aug.root");//
  TFile *ifp1  = new TFile(ifname1  );
      cout<<"input filename1 : "<<ifname1<<endl;
Settings *set=new Settings();

TH1F *h_inv_mass,*h_ctime_allpi;
TH2F *h2_pid_decut;
    TF1 *f_pi, *f_k, *f_p;
  f_pi = new TF1("f_pi","[0]/sqrt(x*x-1)", 1,7);
  f_k  = new TF1("f_k" ,"[0]/sqrt(x*x-1)", 1,7);
  f_p  = new TF1("f_p" ,"[0]/sqrt(x*x-1)", 1,7);
  f_pi->SetParameter(0, -0.001*Mpi);
  f_k ->SetParameter(0, 0.001*MK);
  f_p ->SetParameter(0, 0.001*Mp);
  set->SetTF1(f_pi ,6,1,2.5);
  set->SetTF1(f_k  ,3,1,2.5);
  set->SetTF1(f_p  ,2,1,2.5);
int SumDC[3], SumPt[3];
int SumDC2[3], SumPt2[3];
int TopDC[3], TopPt[3];
TLatex *tex = new TLatex(0,0,"aaa");
tex  -> SetTextSize(0.060);
tex  -> SetTextFont(42);
tex  -> SetTextAlign(21);

  TArrow *arrow = new TArrow(0,0,0,0,0.01,"");//x1,y1,x2,y2,arrow size,option
  arrow->SetAngle(45);
  arrow->SetLineColor(1);
  arrow->SetLineWidth(3);

h_inv_mass   = (TH1F*)ifp1->Get("vertex/h_inv_mass");
h_ctime_allpi= (TH1F*)ifp1->Get("cointime/h_ctime_allpi");
h2_pid_decut = (TH2F*)ifp1->Get("track/h2_pid_decut");
SetTH1(h_inv_mass,    "" ,"invariant mass [GeV/#it{c}^{2}]"      ,"counts/2.5MeV/#it{c}^{2}", 1, 3002, 2);
SetTH1(h_ctime_allpi ,"" ,"coin. time [ns]","counts/25ps", 1, 3001, 3);
SetTH2(h2_pid_decut,"" ,"1/#beta","momentum[GeV/#it{c}]");
h_inv_mass   ->GetXaxis()->SetRangeUser(1.05,1.30);
h_ctime_allpi->GetXaxis()->SetRangeUser(-3,3);

TCanvas *c[3];
for(int i=0;i<3;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800/2,1500/3);}
c[0]->Clear();
gPad->SetLogz(1);h2_pid_decut->Draw("colz");
  f_pi->Draw("same");
//  f_k ->Draw("same");
  f_p ->Draw("same");

c[1]->Clear();
gPad->SetLogy(0);h_ctime_allpi->Draw("");

c[2]->Clear();
h_inv_mass->Draw("");
  tex  ->DrawLatex(0.001*ML,6,"#Lambda mass");
  arrow->DrawArrow(0.001*ML,5,0.001*ML,0.5,0.02,"|>");
c[0]->Print("pdf/YONUPA2017.pdf[");
c[0]->Print("pdf/YONUPA2017.pdf");
c[1]->Print("pdf/YONUPA2017.pdf");
c[2]->Print("pdf/YONUPA2017.pdf");
c[2]->Print("pdf/YONUPA2017.pdf]");
}
