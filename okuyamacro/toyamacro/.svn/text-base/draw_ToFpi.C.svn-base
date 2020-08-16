void draw_ToFpi(){
  gStyle->SetOptStat("");

  TFile *file = new TFile("./root/param/group1to9_tof.root");
  TH1F* h_oa           = (TH1F*)file->Get("vertex/h_oa");
  TH2F* h_pid          = (TH2F*)file->Get("vertex/h2_pid_all");
  TH2F* h_ft_mom_tdll6 = (TH2F*)file->Get("cointime/h2_ftime_mom_tdll6");
  TH2F* h_ct_mom_tdll6 = (TH2F*)file->Get("cointime/h2_cointime_mom_tdll6");
  TH2F* h_vt_mom_tdll6 = (TH2F*)file->Get("cointime/h2_time_at_ver_mom_tdll6");
  SetTH1(h_oa,"opening angle","cos#theta","counts",1,3000,0);
  SetTH2(h_pid,"","1/#beta","momentum[GeV/#it{c}]",0);
  SetTH2(h_ft_mom_tdll6,"ToF vs mom(TDLL6)"           ,"ToF_{#pi}"     ,"momentum[GeV/#it{c}]",0);
  SetTH2(h_ct_mom_tdll6,"coin. time vs mom(TDLL6)"    ,"coin time"     ,"momentum[GeV/#it{c}]",0);
  SetTH2(h_vt_mom_tdll6,"time at vertex vs mom(TDLL6)","time at vertex","momentum[GeV/#it{c}]",0);
  h_ft_mom_tdll6 ->GetXaxis()->SetRangeUser(0,1);

  h_pid->Rebin2D(4,4);

  TCanvas *c[5];
  for(int i=0;i<5;i++){
    c[i]= new TCanvas(Form("c_%d",i+1),Form("c_%d",i+1),900,800 );
  }
  c[0]->Clear();
  c[0]->cd(1);
  gPad->SetLogy(1);
  h_oa->Draw();

  c[1]->Clear();
  c[1]->cd(1);
  gPad->SetLogz(0);
  h_pid->Draw("colz");

  c[2]->Clear();
  c[2]->cd(1);
  gPad->SetLogz(1);
  h_ft_mom_tdll6->Draw("colz");

  c[3]->Clear();
  c[3]->cd(1);
  gPad->SetLogz(0);
  h_ct_mom_tdll6->Draw("colz");

  c[4]->Clear();
  c[4]->cd(1);
  gPad->SetLogz(0);
  h_vt_mom_tdll6->Draw("colz");

  c[0]->Print("pdf/ToFpi.pdf[");
  c[0]->Print("pdf/ToFpi.pdf");
  c[1]->Print("pdf/ToFpi.pdf");
  c[2]->Print("pdf/ToFpi.pdf");
  c[3]->Print("pdf/ToFpi.pdf");
  c[4]->Print("pdf/ToFpi.pdf");
  c[4]->Print("pdf/ToFpi.pdf]");

}
