#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TAxis.h"
#include "TLatex.h"
#include "Math/SpecFunc.h"

void init() {
  gROOT->SetStyle("Plain");
  //gStyle->SetPaperSize(10,10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetFuncWidth(2);
  gStyle->SetHistLineWidth(2);
  gStyle->SetGridWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetTickLength(0.04,"x");
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetPadTickX(0); // 1: top tick on 
  gStyle->SetPadTickY(0); // 1: left tick on
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(0.02); // graph title y position
  gStyle->SetTitleY(0.02); // graph title x position
  gStyle->SetTitleXOffset(0.85); // length : y title <--> y axis
  gStyle->SetTitleYOffset(0.573); // length : x title <--> x axis
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleColor(1);
  gStyle->SetLabelFont(132,"xyz");
  gStyle->SetTitleFont(132,"xyz");
  gStyle->SetTitleFont(132,"");
  gStyle->SetTextFont(132);
  gStyle->SetStatFont(132);
  gStyle->SetNdivisions(505,"xyz");
  gStyle->SetLabelSize(0.15,"xyz");
  gStyle->SetTitleSize(0.15,"xyz");
  gStyle->SetTitleSize(0.06,"");
  gStyle->SetPadTopMargin(0.083);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.22);
  gStyle->SetPadBottomMargin(0.25);
  gStyle->SetLabelOffset(0.015,"xy");
  gStyle->SetLineStyleString(5,"");
  gStyle->SetLineStyleString(6,"48 32");
  gStyle->SetLineStyleString(7,"16 32");
  gStyle->SetLineStyleString(8,"48 32 16 32");
  gStyle->SetLineStyleString(9,"96 32");
}

Double_t redef_sph_bessel(Double_t *x, Double_t *par) {
  return ROOT::Math::sph_bessel((Int_t)par[0],x[0]);
}
Double_t redef_sph_neumann(Double_t *x, Double_t *par) {
  return ROOT::Math::sph_neumann((Int_t)par[0],x[0]);
}

void plot_bes() {
//  init();
  TF1 *f1[5];
  TF1 *f2[5];
  Int_t colarr[5] = {1,2,kGreen+1,4,kMagenta+1};
  Int_t styarr[5] = {1,2,3,4,7};
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  for (Int_t i =0; i < 5; i++) {
    f1[i] = new TF1("f1",redef_sph_bessel,0,10,1);
    f1[i]->SetParameter(0,i);
    f1[i]->SetNpx(200);
    f1[i]->SetLineColor(colarr[i]);
    f1[i]->SetLineStyle(i+5);
  }
//  for (Int_t i =0; i < 5; i++) {
//    f2[i] = new TF1("f2",redef_sph_neumann,0.0,10,1);
//    f2[i]->SetParameter(0,i);
//    f2[i]->SetNpx(200);
//    f2[i]->SetLineColor(colarr[i]);
//    f2[i]->SetLineStyle(i+5);
//  }
//
//  TPad *p1 = new TPad("p1", "p1", 0.0, 0.55, 1.0, 1.0);
//  p1->SetBottomMargin(0);
//  p1->Draw();

//  TPad *p2 = new TPad("p2", "p2", 0.0, 0.0, 1.0, 0.55);
//  p2->SetTopMargin(0);
//  p2->Draw();
  
//  p1->cd();
//  TH1* frame = p1->DrawFrame(0,-0.499,10,1.5);
//  frame->GetXaxis()->SetTitle("#it{x}");
//  frame->GetYaxis()->SetTitle("#it{j_{#kern[-1.0]{l}}}(#it{x})");   
//  frame->GetXaxis()->CenterTitle();
//  frame->GetYaxis()->CenterTitle();
	f1[0]->Draw();
  for (Int_t i =1; i < 5; i++) {
    f1[i]->Draw("same");
  }
//
//  gStyle->SetLabelSize(0.123,"xyz");
//  gStyle->SetTitleSize(0.123,"xyz");
//  gStyle->SetTitleYOffset(0.75);

//  p2->cd();
//  frame = p2->DrawFrame(0,-1.0,10,0.999);
//  frame->GetXaxis()->SetTitle("#it{x}");
//  frame->GetYaxis()->SetTitle("#it{n_{l}}(#it{x})");   
//  frame->GetXaxis()->CenterTitle();
//  frame->GetYaxis()->CenterTitle();
//  for (Int_t i =0; i < 5; i++) {
//    f2[i]->Draw("SAME");
//  }
//   
//  p1->cd();
//  TLatex Tl;
//  Tl.SetTextSize(0.08*0.55/0.45);
//  Tl.SetTextAngle(-40);
//  Tl.SetTextColor(colarr[0]);
//  Tl.DrawLatex(1.0,1.05,"#it{j}_{#kern[-1.0]{0}}(#it{x})");
//  Tl.SetTextAngle(-10);
//  Tl.SetTextColor(colarr[1]);
//  Tl.DrawLatex(2.2,0.60,"#it{j}_{#kern[-1.0]{1}}(#it{x})");
//  Tl.SetTextAngle(0);
//  Tl.SetTextColor(colarr[2]);
//  Tl.DrawLatex(3.3,0.45,"#it{j}_{#kern[-1.0]{2}}(#it{x})");
//  Tl.SetTextAngle(0);
//  Tl.SetTextColor(colarr[3]);
//  Tl.DrawLatex(4.7,0.38,"#it{j}_{#kern[-1.0]{3}}(#it{x})");
//  Tl.SetTextAngle(0);
//  Tl.SetTextColor(colarr[4]);
//  Tl.DrawLatex(6.0,0.32,"#it{j}_{#kern[-1.0]{4}}(#it{x})");

//  p2->cd();
//  Tl.SetTextSize(0.08);
//  Tl.SetTextAngle(0);
//  Tl.SetTextColor(colarr[0]);
//  Tl.DrawLatex(2.3,0.45,"#it{n}_{0}(#it{x})");
//  Tl.SetTextColor(colarr[1]);
//  Tl.DrawLatex(4.0,0.35,"#it{n}_{1}(#it{x})");
//  Tl.SetTextColor(colarr[2]);
//  Tl.DrawLatex(5.3,0.30,"#it{n}_{2}(#it{x})");
//  Tl.SetTextColor(colarr[3]);
//  Tl.DrawLatex(6.6,0.25,"#it{n}_{3}(#it{x})");
//  Tl.SetTextColor(colarr[4]);
//  Tl.DrawLatex(8.0,0.22,"#it{n}_{4}(#it{x})");
//  c1->Update();
}

