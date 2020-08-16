#include "Settings.cc"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void draw_vertex(){
gStyle->SetOptStat(0);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.13);
gStyle->SetPadLeftMargin(0.16);
gStyle->SetPadBottomMargin(0.15);

TString ifname1=Form("./root/RKV_ana/all.root");//
//TString ifname1=Form("./root/RKV_ana/ana10246_10256.root");//
  TFile *ifp1  = new TFile(ifname1  );
      cout<<"input filename1 : "<<ifname1<<endl;
Settings *set=new Settings();

    TH1F *h_oa_ee, *h_oa_ppi, *h_dca, *h_dca_cut, *h_dca_ppi;
    TH1F *h_inv_mass, *h_invm_lpi, *h_mis_mass;
    TH1F *h_invm_eg1, *h_invm_eg2, *h_invm_eg3;
    TH1F *h_invm_mm[5];//missing mass cut dependence
    TH1F *h_invm_chi[6], *h_invm_chi_pi[6], *h_invm_chi_p[6];//chi square cut dependence
    TH1F *h_inv_coin[6];//cointime cut dependence

    TH2F *h2_dvol, *h2_dvol_oa, *h2_dvol_ppi, *h2_dvol_cut, *h2_dvol_Lam, *h2_dvol_bg;
    TH2F *h2_im_vs_mm;
    TH2F *h2_im_pimom;

   h_inv_mass = (TH1F*)ifp1->Get("vertex/h_inv_mass" );
   h_invm_lpi = (TH1F*)ifp1->Get("vertex/h_invm_lpi" );
   h_invm_eg1 = (TH1F*)ifp1->Get("vertex/h_invm_eg1" );
   h_invm_eg2 = (TH1F*)ifp1->Get("vertex/h_invm_eg2" );
   h_invm_eg3 = (TH1F*)ifp1->Get("vertex/h_invm_eg3" );
   h_mis_mass = (TH1F*)ifp1->Get("vertex/h_mis_mass" );
   h_oa_ee    = (TH1F*)ifp1->Get("vertex/h_oa_ee"    );
   h_oa_ppi   = (TH1F*)ifp1->Get("vertex/h_oa_ppi"   );
   h_dca      = (TH1F*)ifp1->Get("vertex/h_dca"      );
   h_dca_ppi  = (TH1F*)ifp1->Get("vertex/h_dca_ppi"  );
   h_dca_cut  = (TH1F*)ifp1->Get("vertex/h_dca_cut"  );
  for(int i=0;i<5;i++){
    h_invm_chi[i]   = (TH1F*)ifp1->Get(Form("vertex/h_invm_chi"   ,i+1) );
    h_invm_chi_pi[i]= (TH1F*)ifp1->Get(Form("vertex/h_invm_chi_pi",i+1) );
    h_invm_chi_p[i] = (TH1F*)ifp1->Get(Form("vertex/h_invm_chi_p" ,i+1) );
    h_inv_coin[i]   = (TH1F*)ifp1->Get(Form("vertex/h_inv_coin"   ,i+1) );
    h_invm_mm[i]    = (TH1F*)ifp1->Get(Form("vertex/h_invm_mm"    ,i+1) );
  }

   h2_dvol    = (TH2F*)ifp1->Get("vertex/h2_dvol"     );
   h2_dvol_oa = (TH2F*)ifp1->Get("vertex/h2_dvol_oa"  );
   h2_dvol_ppi= (TH2F*)ifp1->Get("vertex/h2_dvol_ppi" );
   h2_dvol_cut= (TH2F*)ifp1->Get("vertex/h2_dvol_cut" );
   h2_dvol_Lam= (TH2F*)ifp1->Get("vertex/h2_dvol_Lam" );
   h2_dvol_bg = (TH2F*)ifp1->Get("vertex/h2_dvol_bg"  );
   h2_im_vs_mm= (TH2F*)ifp1->Get("vertex/h2_im_vs_mm" );
   h2_im_pimom= (TH2F*)ifp1->Get("vertex/h2_im_pimom" );


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
c[0]->Divide(3,2);
c[0]->cd(1);gPad->SetLogz(0);h2_dvol    ->Draw("colz");
c[0]->cd(2);gPad->SetLogz(1);h2_dvol_oa ->Draw("colz");
c[0]->cd(3);gPad->SetLogz(1);h2_dvol_ppi->Draw("colz");
c[0]->cd(4);gPad->SetLogz(0);h2_dvol_cut->Draw("colz");
c[0]->cd(5);gPad->SetLogz(0);h2_dvol_Lam->Draw("colz");
c[0]->cd(6);gPad->SetLogz(0);h2_dvol_bg ->Draw("colz");

c[1]->Clear();
c[1]->Divide(3,2);
c[1]->cd(1);h_inv_mass->Draw();
c[1]->cd(2);h_invm_lpi->Draw();
c[1]->cd(3);h_invm_eg1->Draw();
c[1]->cd(4);h_invm_eg2->Draw();
c[1]->cd(5);h_invm_eg3->Draw();
c[1]->cd(6);h_mis_mass->Draw();

c[2]->Clear();
c[2]->Divide(3,2);
c[2]->cd(1);gPad->SetLogy(1);h_oa_ee  ->Draw();
c[2]->cd(2);gPad->SetLogy(1);h_oa_ppi ->Draw();
c[2]->cd(3);gPad->SetLogy(1);h_dca    ->Draw();
c[2]->cd(4);gPad->SetLogy(1);h_dca_ppi->Draw();
c[2]->cd(5);gPad->SetLogy(1);h_dca_cut->Draw();

#if 1
  c[0]->Print("pdf/dvol.pdf[");
  c[0]->Print("pdf/dvol.pdf");
  c[1]->Print("pdf/dvol.pdf");
  c[2]->Print("pdf/dvol.pdf");
  c[3]->Print("pdf/dvol.pdf");
  c[3]->Print("pdf/dvol.pdf]");
#endif
}
