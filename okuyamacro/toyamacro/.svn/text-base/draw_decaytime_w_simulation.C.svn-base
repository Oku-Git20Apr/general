#include "Settings.cc"
const int color[]={2,3,4,6,7,8,9};
const int lifetime[21]={150, 160, 170, 180, 190 ,
                        200, 210, 220, 230, 240 ,
                        250, 260, 270, 280, 290 ,
                        300, 310, 320, 330, 340 ,350};
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
void draw_decaytime_w_simulation(double NLambda=150,double inv_min=1.105, double inv_max=1.125,int res_width = 180){
  if(res_width != 180 && res_width != 182 &&  res_width != 185 &&  res_width != 188 &&  res_width != 190){
    cout<<"invalid width of responese function :"<<res_width<<endl;
    return;
  }
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.12);
  Settings *set=new Settings();

  TString ifname1=Form("./root/RKV_ana_acci/true.root");
  TString ifname2=Form("./root/pdf.root");
  TFile *ifp1  = new TFile(ifname1  );
  TFile *ifp2  = new TFile(ifname2  );
  cout<<"input filename1 : "<<ifname1<<endl;
  cout<<"input filename2 : "<<ifname2<<endl;

  //TH1F *h_dtimeLam = (TH1F*)ifp1->Get("cointime/h_ctime_Lam_cor");
  #if 1
  TTree *tree = (TTree*)ifp1->Get("tree");
  TH1F *h_dtimeLam = new TH1F("h_dtimeLam","h_dtimeLam",160,-2,2);
  tree->Project("h_dtimeLam","dtime",Form("inv_mass>%lf && inv_mass<%lf",inv_min,inv_max));
  #endif
  TH1F *h_dtimeRes = (TH1F*)ifp1->Get("pipi/h_ctime_pipi_cor");
  TH1F *h_dtimeBG  = (TH1F*)h_dtimeRes->Clone();
  int rebin = 4;
  h_dtimeLam ->Rebin(rebin);
  set->SetTH1(h_dtimeLam,"","",Form("Counts/%dps",25*rebin),2,3000,0);
  set->SetTH1(h_dtimeRes,"","decay time [ns]","Counts/25ps"                  ,1,3000,0);

  TH1F *h_pdf[21];
  TH1F *h_pdf_added[21];
  TGraphErrors *tg_stand=new TGraphErrors();
  TGraph *tg_resi[21];
  for(int i=0;i<21;i++){
    tg_resi[i]=new TGraph();
    set->SetGr(tg_resi[i],"","","",1,color[i%7],20+i,0);
  }


  //int lifetime[6] = {100,150,200,263,300,500};
  //double NLambda =150.;//num of Lambda
  double NBG     = (double)h_dtimeLam->Integral() - NLambda;
  h_dtimeBG ->Scale( NBG/((double)h_dtimeBG ->Integral()) );
  cout<<"NL:NBG= "<<NLambda<<" "<<NBG<<endl;
  for(int i=0;i<21;i++){
    cout<<i<<endl;
    h_pdf[i]       = (TH1F*)ifp2->Get(Form("h_pdf_L%dps_R%dps",lifetime[i],res_width));
    h_pdf_added[i] = (TH1F*)ifp2->Get(Form("h_pdf_L%dps_R%dps",lifetime[i],res_width));
    double nevent =h_pdf_added[i]->Integral();
    //cout<<"NBin : "<<h_dtimeBG->GetNbinsX()<<endl;
    h_pdf_added[i] ->Scale( NLambda/nevent ); 
    h_pdf_added[i] ->Add( h_dtimeBG  ); 
    h_pdf_added[i] ->Rebin( rebin  ); 
    h_pdf_added[i] ->SetMarkerStyle( 20+i  ); 
    h_pdf_added[i] ->SetMarkerColor( color[i%7]  ); 
  }

  TLatex *text[2];
  for(int i=0;i<2;i++){
    text[i] = new TLatex();
    set->SetTLatex(text[i],2-i,0.065,32);
  }
  TLegend *legend = new TLegend(0.70,0.50,0.87,0.90);
  legend -> SetBorderSize(1);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1);
  legend -> SetTextFont(22);
  legend -> AddEntry(h_dtimeLam,"data","LE");


  TH1F *h_dtimeRes_scale = (TH1F*)h_dtimeRes->Clone();
  double maxyLam = h_dtimeLam->GetBinContent(h_dtimeLam->GetMaximumBin());
  double maxyRes = h_dtimeRes->GetBinContent(h_dtimeRes->GetMaximumBin());
  h_dtimeRes_scale->Scale(maxyLam/maxyRes);
  h_dtimeRes_scale->SetMinimum(0.1);
  int NbinxLam = h_dtimeLam ->GetNbinsX();
  int NbinxRes = h_dtimeRes ->GetNbinsX();

  //////////////////
  //calc. residual//
  //////////////////
  for(int i=0;i<NbinxLam;i++){
    double x = h_dtimeLam->GetXaxis()->GetBinCenter(i);
    tg_stand->SetPoint(i,x,0);
    tg_stand->SetPointError(i,0,h_dtimeLam->GetBinError(i));
  }

  double chi_sq[21];
  for(int i=0;i<21;i++){
    chi_sq[i]=0;
    for(int k=0;k<NbinxLam;k++){
      double x = h_dtimeLam->GetXaxis()->GetBinCenter(k);
      double y = h_pdf_added[i]->GetBinContent(k) - h_dtimeLam->GetBinContent(k);
      tg_resi[i]->SetPoint(k,x,y);
      double binerror = h_dtimeLam->GetBinError(k);
      if(binerror< 0.05)binerror=1.;
      if(x>0.5&&x<1.0)chi_sq[i] += pow( y/binerror,2);
    }
    //cout<<"tau "<<lifetime[i]<<" "<<chi_sq[i]<<endl;
  }

  ////////
  //draw//
  ////////
  TCanvas *c[3];
  for(int i=0;i<3;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1000,1200);
  }

  c[0]->Divide(1,2);
  c[0]->cd(1);
  gPad->SetLogy(0);
  tg_stand->Draw("AP");
  for(int i=0;i<21;i++){
    tg_resi[i]->Draw("sameP");
  }
  c[0]->cd(2);
  gPad->SetLogy(1);
  h_dtimeRes ->Draw();

  c[1]->Divide(1,2);
  c[1]->cd(1);
  gPad->SetLogy(1);
  h_dtimeRes_scale ->Draw();
  h_dtimeLam ->Draw("same");
  text[0]->DrawLatex(1.8,maxyLam    ,"#Lambda event" );
  text[1]->DrawLatex(1.8,maxyLam*0.5,"Response func.");

  c[1]->cd(2);
  gPad->SetLogy(1);
  h_dtimeRes_scale ->Draw();
  h_dtimeLam ->Draw("same");

  c[2]->Divide(1,2);
  c[2]->cd(1);
  gPad->SetLogy(1);
  h_dtimeLam ->Draw("PE");
  for(int i=0;i<21;i++){
    h_pdf_added[i]->Draw("sameP");
  }
  c[2]->cd(2);
  h_dtimeLam ->Draw("PE");
  for(int i=0;i<21;i++){
    h_pdf_added[i]->Draw("sameP");
    legend -> AddEntry(h_pdf_added[i],Form("#tau=%d ps",lifetime[i]),"P");
  }
  legend->Draw("same");

  c[0]->Print(Form("pdf/decaytime/decaytime_w_simulation_NL%d_inv%.0lf_%.0lf_w%dps.pdf[",(int)NLambda,inv_min*1000.,inv_max*1000.,res_width));
  c[0]->Print(Form("pdf/decaytime/decaytime_w_simulation_NL%d_inv%.0lf_%.0lf_w%dps.pdf" ,(int)NLambda,inv_min*1000.,inv_max*1000.,res_width));
  c[1]->Print(Form("pdf/decaytime/decaytime_w_simulation_NL%d_inv%.0lf_%.0lf_w%dps.pdf" ,(int)NLambda,inv_min*1000.,inv_max*1000.,res_width));
  c[2]->Print(Form("pdf/decaytime/decaytime_w_simulation_NL%d_inv%.0lf_%.0lf_w%dps.pdf" ,(int)NLambda,inv_min*1000.,inv_max*1000.,res_width));
  c[2]->Print(Form("pdf/decaytime/decaytime_w_simulation_NL%d_inv%.0lf_%.0lf_w%dps.pdf]",(int)NLambda,inv_min*1000.,inv_max*1000.,res_width));
  ofstream fout;
  fout.open(Form("Life_chi2_list/decaytime_chisq_inv%.0lf_%.0lf_w%dps.dat",inv_min*1000.,inv_max*1000.,res_width), ios::out|ios::app);
  for(int i=0;i<21;i++){
    fout<<NLambda<<" "<<lifetime[i]<<" "<<chi_sq[i]<<endl;
  }
}
////////////++++++++++++++++++++///////////////
void draw_result(double inv_min=1.105, double inv_max=1.125,int res_width = 180){
  if(res_width != 180 && res_width != 182 &&  res_width != 185 &&  res_width != 188 &&  res_width != 190){
    cout<<"invalid width of responese function :"<<res_width<<endl;
    return;
  }
  Settings *set = new Settings();
  TH2F *h_frame[2];
  h_frame[0] = new TH2F("h_frame1","h_frame1",100,125, 205, 100, 0, 3);
  h_frame[1] = new TH2F("h_frame2","h_frame2",100,105, 355, 100, 0, 3);
  set->SetTH2(h_frame[0] ,"" ,"Num. of #Lambda","#chi^{2}");
  set->SetTH2(h_frame[1] ,"" ,"#tau_{#Lambda}","#chi^{2}");
  TGraph *tg_NL_chi[21];
  TGraph *tg_tau_chi[10];
  for(int i=0;i<21;i++){
    tg_NL_chi[i]=new TGraph();
    set->SetGr(tg_NL_chi[i],"","","",1,color[i%7],20+i,0);
    tg_NL_chi[i]->SetMarkerSize(1.5);
  }
  for(int i=0;i<10;i++){
    tg_tau_chi[i]=new TGraph();
    set->SetGr(tg_tau_chi[i],"","","",1,color[i%7],20+i,0);
    tg_tau_chi[i]->SetMarkerSize(1.5);
  }

  TLegend *legend_NL = new TLegend(0.80,0.50,0.87,0.90);
  legend_NL -> SetBorderSize(1);
  legend_NL -> SetFillColor(0);
  legend_NL -> SetFillStyle(1);
  legend_NL -> SetTextFont(22);

  TLegend *legend_tau = new TLegend(0.80,0.50,0.87,0.90);
  legend_tau -> SetBorderSize(1);
  legend_tau -> SetFillColor(0);
  legend_tau -> SetFillStyle(1);
  legend_tau -> SetTextFont(22);

  TLegend *legend_tau_dec = new TLegend(0.80,0.50,0.87,0.90);
  legend_tau_dec -> SetBorderSize(1);
  legend_tau_dec -> SetFillColor(0);
  legend_tau_dec -> SetFillStyle(1);
  legend_tau_dec -> SetTextFont(22);

  int NLam[10]={110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
  int point_NL[21], point_tau[10];
  for(int i=0;i<21;i++){
    point_NL[i] = 0;
  }
  for(int i=0;i<10;i++){
    if(inv_min>1.109&&inv_min<1.112)     NLam[i]=NLam[i]- 40;//better S/N
    else if(inv_max<1.118)NLam[i]=NLam[i]-100;//best S/N
    point_tau[i] = 0;
  }

  string ifname = Form("Life_chi2_list/decaytime_chisq_inv%.0lf_%.0lf_w%dps.dat",inv_min*1000,inv_max*1000,res_width);
  std::ifstream fp(ifname.c_str());
  if(fp.fail()){ cout<<"file open fail!  "<<ifname<< endl; exit(1); }
  string line, aa;
  int NL,tau;
  double chisq;
  while(1){
    getline(fp,line);
    if(line[0]=='#') continue;
    if( fp.eof() ) break;
    istringstream sline(line);
    sline >> NL >> tau >> chisq;
    //cout<< NL <<" "<<tau<<" "<<chisq<<endl;
    for(int i=0;i<21;i++){
      if(lifetime[i]==tau){
        tg_NL_chi[i]->SetPoint(point_NL[i],NL,chisq);
        point_NL[i]++;
      }
    }
    /////////
    for(int i=0;i<10;i++){
      if(NLam[i]==NL){
        tg_tau_chi[i]->SetPoint(point_tau[i],tau,chisq);
        point_tau[i]++;
      }
    }
  }
  ////////
  //draw//
  ////////
  TCanvas *c[3];
  for(int i=0;i<3;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),2000,1200);
  }
  c[0]->Clear();
  c[0]->cd(1);
  h_frame[0]->Draw();
  for(int i=0;i<21;i++){
    legend_NL -> AddEntry(tg_NL_chi[i],Form("#tau=%d ps",lifetime[i]),"P");
    tg_NL_chi[i]->Draw("sameP");
  }
  legend_NL->Draw("same");

  c[1]->Clear();
  c[1]->cd(1);
  h_frame[1]->Draw();
  for(int i=0;i<10;i++){
    legend_tau -> AddEntry(tg_tau_chi[i],Form("N#Lambda =%d",NLam[i]),"P");
    tg_tau_chi[i]->Draw("sameP");
  }
  legend_tau->Draw("same");

  c[2]->Clear();
  c[2]->cd(1);
  h_frame[1]->Draw();
  legend_tau_dec -> AddEntry(tg_tau_chi[2],Form("N#Lambda =%d",NLam[2]),"P");
  legend_tau_dec -> AddEntry(tg_tau_chi[5],Form("N#Lambda =%d",NLam[5]),"P");
  legend_tau_dec -> AddEntry(tg_tau_chi[8],Form("N#Lambda =%d",NLam[8]),"P");
  tg_tau_chi[2]->Draw("samePL");
  tg_tau_chi[5]->Draw("samePL");
  tg_tau_chi[8]->Draw("samePL");
  legend_tau_dec ->Draw("same");
  
  c[0]->Print(Form("pdf/decaytime/decaytime_w_simulation_result_invm%.0lf_%.0lf_w%dps.pdf[",inv_min*1000,inv_max*1000,res_width));
  c[0]->Print(Form("pdf/decaytime/decaytime_w_simulation_result_invm%.0lf_%.0lf_w%dps.pdf" ,inv_min*1000,inv_max*1000,res_width));
  c[1]->Print(Form("pdf/decaytime/decaytime_w_simulation_result_invm%.0lf_%.0lf_w%dps.pdf" ,inv_min*1000,inv_max*1000,res_width));
  c[2]->Print(Form("pdf/decaytime/decaytime_w_simulation_result_invm%.0lf_%.0lf_w%dps.pdf" ,inv_min*1000,inv_max*1000,res_width));
  c[2]->Print(Form("pdf/decaytime/decaytime_w_simulation_result_invm%.0lf_%.0lf_w%dps.pdf]",inv_min*1000,inv_max*1000,res_width));
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void draw_mom_of_pdf(){
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.12);
  Settings *set=new Settings();
  int res_width = 180;

  TString ifname2=Form("./root/pdf.root");
  TFile *ifp2  = new TFile(ifname2  );
  cout<<"input filename2 : "<<ifname2<<endl;

  TH2F *h_frame[2];
  h_frame[0] = new TH2F("h_frame1","h_frame1",100,105, 355, 100, 105, 300);
  h_frame[1] = new TH2F("h_frame2","h_frame2",100,105, 355, 100, 0, 3);
  set->SetTH2(h_frame[0] ,"" ,"#tau_{#Lambda}","#mu");
  set->SetTH2(h_frame[1] ,"" ,"#tau_{#Lambda}","#chi^{2}");

  TH1F *h_pdf[21];
  double pdf_mom[21];
  
  TGraph *tg_mom=new TGraph();
  set->SetGr(tg_mom,"","","",1,color[2],28,0);

  for(int i=0;i<21;i++){
    cout<<i<<endl;
    h_pdf[i]       = (TH1F*)ifp2->Get(Form("h_pdf_L%dps_R%dps",lifetime[i],res_width));
    //cout<<"NBin : "<<h_dtimeBG->GetNbinsX()<<endl;
    h_pdf[i] ->SetMarkerStyle( 20+i  ); 
    h_pdf[i] ->SetMarkerColor( color[i%7]  ); 
    h_pdf[i] ->SetLineColor( color[i%7]  ); 
    pdf_mom[i] = get_1st_mom(h_pdf[i],0.5,0.5);
    cout<<pdf_mom[i]<<endl;
    tg_mom->SetPoint(i,lifetime[i],1000.*pdf_mom[i]);
  }

  ////////
  //draw//
  ////////
  TCanvas *c[3];
  for(int i=0;i<3;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),2000,1200);
  }
  c[0]->Clear();
  c[0]->cd(1);
  h_frame[0]->Draw();
  tg_mom->Draw("samePL");

  c[1]->Clear();
  c[1]->cd(1);
  gPad->SetLogy(1);
    h_pdf[0]->Draw("");
  for(int i=0;i<21;i++){
    h_pdf[i]->Draw("same");
  }

  c[2]->Clear();
  c[2]->cd(1);

}
