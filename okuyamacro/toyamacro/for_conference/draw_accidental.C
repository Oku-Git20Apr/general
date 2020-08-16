#include "Settings.cc"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void draw_accidental(){
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.15);
  Settings *set=new Settings();

  /////////
  //Param//
  /////////
  int rebin =1;
  int rebin_time =2;
  double MassRange_Lam_min;
  double MassRange_Lam_max;
  MassRange_Lam_min   =  1.105;
  MassRange_Lam_max   =  1.125;       
  double scale = 4.0;

  //TString ifname1=Form("../root/RKV_ana/all.root");//
  //TString ifname2=Form("../root/RKV_ana_acci/all_4ns.root");//
  TString ifname1=Form("../root/RKV_ana_acci/true.root");//
  TString ifname2=Form("../root/RKV_ana_acci/all.root");//
  //TString ifname3=Form("../root/RKV_ana_acci/all_RF.root");//
  TString ofname =Form("./root/accidental.root");//
  TFile *ifp1  = new TFile(ifname1  );
  TFile *ifp2  = new TFile(ifname2  );
      cout<<"input filename1 : "<<ifname1<<endl;
      cout<<"input filename2 : "<<ifname2<<endl;

  TArrow *arrow = new TArrow(0,0,0,0,0.01,"");//x1,y1,x2,y2,arrow size,option
  arrow->SetAngle(45);
  arrow->SetLineColor(1);
  arrow->SetLineWidth(3);

  TLine *line = new TLine(0,0,0,0);//x1,y1,x2,y2,arrow size,option
  line->SetLineColor(1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);



  TFile *ofp   = new TFile(ofname,"recreate"  );
  //Get histograms of true event
  TH1F *h_invm_lpi, *h_ctime_allpi_rf, *h_ctime_tb_RF_all, *h_cor;
  h_invm_lpi   = (TH1F*)ifp1->Get("vertex/h_invm_lpi");
  h_invm_lpi   ->Rebin(rebin);
  h_invm_lpi   ->SetMinimum(0);
  h_ctime_allpi_rf  = (TH1F*)ifp1->Get("cointime/h_ctime_allpi_rf");
  h_ctime_allpi_rf->SetMinimum(0.8);
  h_ctime_tb_RF_all = (TH1F*)ifp1->Get("tag/h_ctime_tb_RF_all");
  h_cor = (TH1F*)ifp1->Get("cointime/h_ctime_Lam_cor");
  h_cor ->Rebin(rebin_time);
  h_cor ->SetMinimum(0.1);
  TH2 *h2_pid = (TH2F*)ifp1->Get("track/h2_pid_decut");
  TH2 *h2_pid_ppi = (TH2F*)ifp1->Get("track/h2_pid_ppi");
  SetTH2(h2_pid     ,"" ,"1/#beta"        ,"momentum[GeV/#it{c}]");
  SetTH2(h2_pid_ppi ,"" ,"1/#beta"        ,"momentum[GeV/#it{c}]");

  SetTH1(h_invm_lpi,    "" ,"invariant mass [GeV/#it{c}^{2}]"      ,Form("counts/%.1lfMeV/#it{c}^{2}",2.0*rebin), 1, 3002, 2);
  SetTH1(h_ctime_allpi_rf ,"" ,"coin. time [ns]","counts/25ps", 1, 3001, 4);
  SetTH1(h_ctime_tb_RF_all,"TagB - RF" ,"TagB - RF [ns]","counts/6ps", 1, 3001, 7);
  SetTH1(h_cor            ,"decay time of #Lambda" ,"decay time[ns]","counts", 2, 3000, 0);
  h_invm_lpi       ->GetXaxis()->SetRangeUser(1.05,1.30);
  h_ctime_allpi_rf ->GetXaxis()->SetRangeUser(-6,6);
  h_ctime_tb_RF_all->GetXaxis()->SetRangeUser(-1,1);


  TH1F *h_invm_acc;//pion low momentum event only
  h_invm_acc  = (TH1F*)ifp2 ->Get("vertex/h_invm_lpi");
  h_invm_acc   ->Rebin(rebin);
  h_invm_acc   ->Scale(1./scale);
  h_cor_acc = (TH1F*)ifp2->Get("cointime/h_ctime_Lam_cor");
  h_cor_acc    ->Scale(1./scale);
  h_cor_acc ->Rebin(rebin_time);
  SetTH1(h_invm_acc,    "" ,"invariant mass [GeV/#it{c}^{2}]"      ,Form("counts/%.1lfMeV/#it{c}^{2}",2.0*rebin), 1, 3002, 4);
  SetTH1(h_cor_acc        ,"decay time of #Lambda" ,"decay time[ns]","counts", 1, 3000, 0);
  ofp->cd();
  //h_invm_lpi   ->Write();
  ofp->Write();


  //Fit//
  TF1 *f_ga1 = new TF1("f_ga1","gausn",-1,1);
  TF1 *f_ga1_pol0 = new TF1("f_ga1_pol0","gausn(0)+pol0(3)",-1,1);

  TF1 *f_ga2 = new TF1("f_ga2","gausn",-1,1);
  TF1 *f_ga3 = new TF1("f_ga3","gausn", 1,3);
  //f_ga1->SetParameter(0,100000);
  //f_ga1->SetParameter(1,0.);
  //f_ga1->SetParameter(2,0.2);
  f_ga1_pol0->SetNpx(2000);  f_ga1_pol0->SetLineColor(2);   f_ga1_pol0->SetLineWidth(2.7);  f_ga1_pol0->SetLineStyle(1);
  f_ga2     ->SetNpx(2000);  f_ga2     ->SetLineColor(2);   f_ga2     ->SetLineWidth(2.7);  f_ga2->SetLineStyle(1);
  f_ga3     ->SetNpx(2000);  f_ga3     ->SetLineColor(6);   f_ga3     ->SetLineWidth(2.7);  f_ga3->SetLineStyle(1);

  f_ga1_pol0->SetParameter(0,100000);
  f_ga1_pol0->SetParameter(1,0.);
  f_ga1_pol0->SetParameter(2,0.2);
  f_ga1_pol0->SetParameter(3,h_ctime_tb_RF_all->GetBinContent(h_ctime_tb_RF_all->FindBin(-0.7)));
  h_ctime_tb_RF_all ->Fit(f_ga1_pol0);

  f_ga2->SetParameter(0,10000);
  f_ga2->SetParameter(1,0.);
  f_ga2->SetParameter(2,0.2);
  h_ctime_allpi_rf->Fit(f_ga2,"QR","");

  f_ga3->SetParameter(0,10000);
  f_ga3->SetParameter(1,0.);
  f_ga3->SetParameter(2,0.2);
  h_ctime_allpi_rf->Fit(f_ga3,"QR","");

  cout<<"True:acc= "<<f_ga2->GetParameter(0)<<" : "<<f_ga3->GetParameter(0)<<endl;

  TCanvas *c[6];
  for(int i=0;i<6;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800/2,1500/3);}
  c[0]->Clear();
  gPad->SetLogy(0);
  h_invm_lpi->Draw("");
  h_invm_acc->Draw("same");
  //h_invm_acc->Draw("sameE");
  line->DrawLine(MassRange_Lam_min,0,MassRange_Lam_min,20);
  line->DrawLine(MassRange_Lam_max,0,MassRange_Lam_max,20);
  c[1]->Clear();
  c[1]->cd(1);
  gPad->SetLogy(0);
  h_ctime_allpi_rf->Draw("");
  f_ga2->Draw("same");
  f_ga3->Draw("same");
  
  c[2]->Clear();
  c[2]->cd(1);
  gPad->SetLogy(1);h_ctime_tb_RF_all->Draw("");
  line ->DrawLine( 0.5,0., 0.5, h_ctime_tb_RF_all->GetMaximum());
  line ->DrawLine(-0.5,0.,-0.5, h_ctime_tb_RF_all->GetMaximum());
  arrow->DrawArrow( 0.5, h_ctime_tb_RF_all->GetMaximum(), 0.3,h_ctime_tb_RF_all->GetMaximum(),0.02,"|>");
  arrow->DrawArrow(-0.5, h_ctime_tb_RF_all->GetMaximum(),-0.3,h_ctime_tb_RF_all->GetMaximum(),0.02,"|>");
  //f_ga1->Draw("same");
  f_ga1_pol0->Draw("same");

  c[3]->Clear();
  c[3]->cd(1);
  gPad->SetLogy(1);h_cor->Draw("");h_cor_acc->Draw("same");

  c[4]->Clear();
  c[4]->cd(1);
  TF1 *f_pi, *f_k, *f_p;
  f_pi = new TF1("f_pi","[0]/sqrt(x*x-1)", 1,7);
  f_pi->SetParameter(0, -0.001*Mpi);
  set->SetTF1(f_pi ,6,1,2.5);
  f_p  = new TF1("f_p" ,"[0]/sqrt(x*x-1)", 1,7);
  f_p ->SetParameter(0, 0.001*Mp);
  set->SetTF1(f_p  ,2,1,2.5);
  gPad->SetLogz(1);h2_pid->Draw("colz");
  f_pi->Draw("same");
  f_p ->Draw("same");

  c[5]->Clear();
  c[5]->cd(1);
  gPad->SetLogz(1);h2_pid_ppi->Draw("colz");

  c[0]->Print("pdf/accidental.pdf[");
  c[0]->Print("pdf/accidental.pdf");
  c[1]->Print("pdf/accidental.pdf");
  c[2]->Print("pdf/accidental.pdf");
  c[3]->Print("pdf/accidental.pdf");
  c[4]->Print("pdf/accidental.pdf");
  c[5]->Print("pdf/accidental.pdf");
  c[5]->Print("pdf/accidental.pdf]");
  ofp->Close();
}
