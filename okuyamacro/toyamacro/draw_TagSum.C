#include <iostream>
#include <sstream>
#include <iomanip>
#include <csignal>
#include <stdlib.h>
#include <climits>
#include <fstream>
#include <math.h>
#include <string>
#include <unistd.h>
using namespace std;

#include "TApplication.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include <TLegend.h>
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TString.h"
#include "TPaveText.h"
#include "TSystem.h"

#include "TRandom.h"


const int NIH   = 10;  // No. of IH
const int NTB   = 40;  // No. of TagB
////////////////
void SetTH1(TH1F *h1, TString hname, TString xname, TString yname, int LColor, int FStyle, int FColor){
  h1->SetTitle(hname);
  h1->GetXaxis()->SetTitle(xname);
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->SetTitle(yname);
  h1->GetYaxis()->CenterTitle();
  h1->SetMinimum(0.8);
  h1->SetLineColor(LColor);
  h1->SetFillStyle(FStyle);
  h1->SetFillColor(FColor);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->SetLineWidth(0);
  h1->SetTitleSize(0.05,"");
}
////////////////
void SetTH1(TH1D *h1, TString hname, TString xname, TString yname, int LColor, int FStyle, int FColor){
  h1->SetTitle(hname);
  h1->GetXaxis()->SetTitle(xname);
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->SetTitle(yname);
  h1->GetYaxis()->CenterTitle();
  h1->SetMinimum(0.8);
  h1->SetLineColor(LColor);
  h1->SetFillStyle(FStyle);
  h1->SetFillColor(FColor);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->SetLineWidth(0);
  h1->SetTitleSize(0.06,"x");
  h1->SetTitleSize(0.05,"y");
  h1->SetTitleSize(0.05,"");
  ((TGaxis*)h1->GetYaxis())->SetMaxDigits(3);
  //h1->GetYaxis()->SetNdivisions(505);
}
////////////////

void draw_TagSum(string root_file = "test.root"){
gStyle->SetOptStat("iem");
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.10);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadBottomMargin(0.15);
TString ifname1=;//H3L
  TFile *ifp1  = new TFile( Form("./root/%s",root_file.c_str()) );
      cout<<"input filename : "<<root_file<<endl;
TH1F *h_TagSumTDLL[NIH], *h_TagSumTDLR[NIH];
  for(int j=0;j<NIH;j++){
  h_TagSumTDLL[j] = new TH1F(Form("h_TagSumTDL L%d",j+1),Form("h_TagSumTDL L%d",j+1),20000,-1000,1000);
  h_TagSumTDLL[j]  = (TH1F*)ifp1    -> Get(Form("h_TagSumTDLL%d",j+1));
  //R
  h_TagSumTDLR[j] = new TH1F(Form("h_TagSumTDL R%d",j+1),Form("h_TagSumTDL R%d",j+1),20000,-1000,1000);
  h_TagSumTDLR[j]  = (TH1F*)ifp1    -> Get(Form("h_TagSumTDLR%d",j+1));

  SetTH1(h_TagSumTDLL[j], Form("ToF TagB Sum - TDL L%d",j+1), "coin. time[ns]", "count/0.1ns", 1, 0, 0);
  SetTH1(h_TagSumTDLR[j], Form("ToF TagB Sum - TDL R%d",j+1), "coin. time[ns]", "count/0.1ns", 1, 0, 0);
  }

int p_flag=1;
//cout<<"Do you want to save canvas?"<<endl;
//cout<<"yes :1, no :0"<<endl;
//cin>>p_flag;

TCanvas *c[6];
for(int i=0;i<6;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),900,1200);
                     }
c[0]->Clear();
c[0]->Divide(2,5);
for(int i=0;i<10;i++){
double xmin = h_TagSumTDLL[i] ->GetXaxis()->GetBinCenter(h_TagSumTDLL[i] ->GetMaximumBin() ) - 100;
double xmax = h_TagSumTDLL[i] ->GetXaxis()->GetBinCenter(h_TagSumTDLL[i] ->GetMaximumBin() ) + 100;
//h_TagSumTDLL[i] ->GetXaxis()->SetRangeUser(xmin,xmax);
h_TagSumTDLL[i] ->GetXaxis()->SetRangeUser(-80,80);
c[0]->cd(i+1);
gPad->SetLogy(0);
h_TagSumTDLL[i] ->Draw("");
}
////
c[1]->Clear();
c[1]->Divide(2,5);
for(int i=0;i<10;i++){
double xmin = h_TagSumTDLR[i] ->GetXaxis()->GetBinCenter(h_TagSumTDLR[i] ->GetMaximumBin() ) - 100;
double xmax = h_TagSumTDLR[i] ->GetXaxis()->GetBinCenter(h_TagSumTDLR[i] ->GetMaximumBin() ) + 100;
h_TagSumTDLR[i] ->GetXaxis()->SetRangeUser(-80,80);
//h_TagSumTDLR[i] ->GetXaxis()->SetRangeUser(xmin,xmax);
c[1]->cd(i+1);
gPad->SetLogy(0);
h_TagSumTDLR[i] ->Draw("");
}

if(p_flag==1){
  string ofname_pdf = root_file;
  ofname_pdf.erase(ofname_pdf.size()-5);
  ofname_pdf.append("_tagsum.pdf");
  c[0]->Print(Form("./pdf/%s[",ofname_pdf.c_str())  );
  c[0]->Print(Form("./pdf/%s" ,ofname_pdf.c_str())  );
  c[1]->Print(Form("./pdf/%s" ,ofname_pdf.c_str())  );
  c[1]->Print(Form("./pdf/%s]",ofname_pdf.c_str())  );
 }

} 
#if 0
c[4]->Clear();
c[4]->Divide(2,5);
for(int i=0;i<10;i++){
double xmin = h_TagTDLRU[24][i] ->GetXaxis()->GetBinCenter(h_TagTDLRU[24][i] ->GetMaximumBin() ) - 100;
double xmax = h_TagTDLRU[24][i] ->GetXaxis()->GetBinCenter(h_TagTDLRU[24][i] ->GetMaximumBin() ) + 100;
h_TagTDLRU[24][i] ->GetXaxis()->SetRangeUser(xmin,xmax);
c[4]->cd(i+1);
gPad->SetLogy(1);
h_TagTDLRU[24][i] ->Draw("");
}
c[5]->Clear();
c[5]->Divide(2,5);
for(int i=0;i<10;i++){
double xmin = h_TagTDLRD[24][i] ->GetXaxis()->GetBinCenter(h_TagTDLRD[24][i] ->GetMaximumBin() ) - 100;
double xmax = h_TagTDLRD[24][i] ->GetXaxis()->GetBinCenter(h_TagTDLRD[24][i] ->GetMaximumBin() ) + 100;
h_TagTDLRD[24][i] ->GetXaxis()->SetRangeUser(xmin,xmax);
c[5]->cd(i+1);
gPad->SetLogy(1);
h_TagTDLRD[24][i] ->Draw("");
c[1]->Clear();
c[1]->Divide(2,5);
for(int i=0;i<10;i++){
double xmin = h_TagTDLLU[24][i] ->GetXaxis()->GetBinCenter(h_TagTDLLU[24][i] ->GetMaximumBin() ) - 100;
double xmax = h_TagTDLLU[24][i] ->GetXaxis()->GetBinCenter(h_TagTDLLU[24][i] ->GetMaximumBin() ) + 100;
h_TagTDLLU[24][i] ->GetXaxis()->SetRangeUser(xmin,xmax);
c[1]->cd(i+1);
gPad->SetLogy(1);
h_TagTDLLU[24][i] ->Draw("");
}
c[2]->Clear();
c[2]->Divide(2,5);
for(int i=0;i<10;i++){
double xmin = h_TagTDLLD[24][i] ->GetXaxis()->GetBinCenter(h_TagTDLLD[24][i] ->GetMaximumBin() ) - 100;
double xmax = h_TagTDLLD[24][i] ->GetXaxis()->GetBinCenter(h_TagTDLLD[24][i] ->GetMaximumBin() ) + 100;
h_TagTDLLD[24][i] ->GetXaxis()->SetRangeUser(xmin,xmax);
c[2]->cd(i+1);
gPad->SetLogy(1);
h_TagTDLLD[24][i] ->Draw("");
}
}
#endif
