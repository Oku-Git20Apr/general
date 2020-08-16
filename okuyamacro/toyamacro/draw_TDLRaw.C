
#include "Settings.h"
////////////////

void draw_TDLRaw(string root_file = "test.root"){
gStyle->SetOptStat("iem");
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.10);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadBottomMargin(0.15);
TString ifname1=;//H3L
  TFile *ifp1  = new TFile( Form("./root/%s",root_file.c_str()) );
      cout<<"input filename : "<<root_file<<endl;
TH1F *h1_luTDC_l[NIH], *h1_ldTDC_l[NIH], *h1_ruTDC_l[NIH], *h1_rdTDC_l[NIH];
TH1F *h1_luTDC_t[NIH], *h1_ldTDC_t[NIH], *h1_ruTDC_t[NIH], *h1_rdTDC_t[NIH];
TH1F *h1_lutime_w[NIH], *h1_ldtime_w[NIH], *h1_rutime_w[NIH], *h1_rdtime_w[NIH];
TH2F *h2_luWvsToF[NIH],*h2_luWvsToF_cut[NIH];
TH2F *h2_ldWvsToF[NIH],*h2_ldWvsToF_cut[NIH];
TH2F *h2_ruWvsToF[NIH],*h2_ruWvsToF_cut[NIH];
TH2F *h2_rdWvsToF[NIH],*h2_rdWvsToF_cut[NIH];
TH2F *h2_TagBWvsToF[NTB];
TH1F *h_TagBW[NTB];
TH1F *h_TagBToF[NTB];
TH1F *h_TagFTDC[NTF];

TH2D *h2_daq0trig_evnum = (TH2D*)ifp1 ->Get("h2_daq0trig_evnum");
TH2D *h2_daq1trig_evnum = (TH2D*)ifp1 ->Get("h2_daq1trig_evnum");;

   for(int i=0;i<NIH;i++){
   
     h1_luTDC_l[i] = (TH1F*)ifp1    -> Get(Form("h1_luTDC_l%d",i+1));
     h1_ldTDC_l[i] = (TH1F*)ifp1    -> Get(Form("h1_ldTDC_l%d",i+1));
     h1_ruTDC_l[i] = (TH1F*)ifp1    -> Get(Form("h1_ruTDC_l%d",i+1));
     h1_rdTDC_l[i] = (TH1F*)ifp1    -> Get(Form("h1_rdTDC_l%d",i+1));
     
     h1_luTDC_t[i] = (TH1F*)ifp1    -> Get(Form("h1_luTDC_t%d",i+1));
     h1_ldTDC_t[i] = (TH1F*)ifp1    -> Get(Form("h1_ldTDC_t%d",i+1));
     h1_ruTDC_t[i] = (TH1F*)ifp1    -> Get(Form("h1_ruTDC_t%d",i+1));
     h1_rdTDC_t[i] = (TH1F*)ifp1    -> Get(Form("h1_rdTDC_t%d",i+1));
   
     h1_lutime_w[i] = (TH1F*)ifp1    -> Get(Form("h1_lutime_w%d",i+1));
     h1_ldtime_w[i] = (TH1F*)ifp1    -> Get(Form("h1_ldtime_w%d",i+1));
     h1_rutime_w[i] = (TH1F*)ifp1    -> Get(Form("h1_rutime_w%d",i+1));
     h1_rdtime_w[i] = (TH1F*)ifp1    -> Get(Form("h1_rdtime_w%d",i+1));
     h1_lutime_w[i] -> GetXaxis()->SetRangeUser(1,60);
     h1_ldtime_w[i] -> GetXaxis()->SetRangeUser(1,60);
     h1_rutime_w[i] -> GetXaxis()->SetRangeUser(1,60);
     h1_rdtime_w[i] -> GetXaxis()->SetRangeUser(1,60);

     h2_luWvsToF[i] = (TH2F*)ifp1   -> Get(Form("h2_luWvsToF%d",i+1));
     h2_ldWvsToF[i] = (TH2F*)ifp1   -> Get(Form("h2_ldWvsToF%d",i+1));
     h2_ruWvsToF[i] = (TH2F*)ifp1   -> Get(Form("h2_ruWvsToF%d",i+1));
     h2_rdWvsToF[i] = (TH2F*)ifp1   -> Get(Form("h2_rdWvsToF%d",i+1));
     h2_luWvsToF[i] -> GetXaxis()->SetRangeUser(810,830);
     h2_ldWvsToF[i] -> GetXaxis()->SetRangeUser(810,830);
     h2_ruWvsToF[i] -> GetXaxis()->SetRangeUser(810,830);
     h2_rdWvsToF[i] -> GetXaxis()->SetRangeUser(810,830);
   }

   for(int i=0;i<NTB;i++){
     h2_TagBWvsToF[i] = (TH2F*)ifp1 ->Get(Form("h2_TagB%dWvsToF",i+1));
     h_TagBW[i]       = (TH1F*)ifp1 ->Get(Form("h_TagB%dW",i+1)     ); 
     h_TagBToF[i]     = (TH1F*)ifp1 ->Get(Form("h_TagB%dToF",i+1)   ); 
   }

   for(int i=0;i<NTF;i++){
      h_TagFTDC[i] = (TH1F*)ifp1 -> Get(Form("h_TagF%dTDC",i+1));
   }

int p_flag=1;
//cout<<"Do you want to save canvas?"<<endl;
//cout<<"yes :1, no :0"<<endl;
//cin>>p_flag;

TCanvas *c[10];
for(int i=0;i<10;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),900,1200);
                     }
c[0]->Divide(3,10);
//c[0]->cd(1) ;gPad->SetLogz(0); h2_daq0trig_evnum->Draw("colz");
//c[0]->cd(2) ;gPad->SetLogz(0); h2_daq1trig_evnum->Draw("colz");
for(int i=0;i<10;i++){
c[0]->cd(3*i+1) ;gPad->SetLogz(1);h2_luWvsToF[i] ->Draw("colz");
c[0]->cd(3*i+2) ;gPad->SetLogy(1);h1_luTDC_l[i]  ->Draw("");
c[0]->cd(3*i+3) ;gPad->SetLogy(1);h1_lutime_w[i] ->Draw("");
}

c[1]->Divide(3,10);
for(int i=0;i<10;i++){
c[1]->cd(3*i+1) ;gPad->SetLogz(1);h2_ldWvsToF[i] ->Draw("colz");
c[1]->cd(3*i+2) ;gPad->SetLogy(1);h1_ldTDC_l[i]  ->Draw("");
c[1]->cd(3*i+3) ;gPad->SetLogy(1);h1_ldtime_w[i] ->Draw("");
}

c[2]->Divide(3,10);
for(int i=0;i<10;i++){
c[2]->cd(3*i+1) ;gPad->SetLogz(1);h2_ruWvsToF[i] ->Draw("colz");
c[2]->cd(3*i+2) ;gPad->SetLogy(1);h1_ruTDC_l[i]  ->Draw("");
c[2]->cd(3*i+3) ;gPad->SetLogy(1);h1_rutime_w[i] ->Draw("colz");
}

c[3]->Divide(3,10);
for(int i=0;i<10;i++){
c[3]->cd(3*i+1) ;gPad->SetLogz(1);h2_rdWvsToF[i] ->Draw("colz");
c[3]->cd(3*i+2) ;gPad->SetLogy(1);h1_rdTDC_l[i]  ->Draw("");
c[3]->cd(3*i+3) ;gPad->SetLogy(1);h1_rdtime_w[i] ->Draw("colz");
}

c[4]->Divide(3,11);
c[4]->cd(1);h_TagBhit->Draw();
c[4]->cd(2);h_TagFhit->Draw();
for(int i=0;i<10;i++){
 c[4]->cd(3*i+1+3) ;gPad->SetLogz(1); h2_TagBWvsToF[i]  ->Draw("colz");
 c[4]->cd(3*i+2+3) ;gPad->SetLogy(0); h_TagBToF[i]      ->Draw("colz");
 c[4]->cd(3*i+3+3) ;gPad->SetLogy(1); h_TagBW[i]        ->Draw("");
}

c[5]->Divide(3,10);
for(int i=0;i<10;i++){
 c[5]->cd(3*i+1) ;gPad->SetLogz(1); h2_TagBWvsToF[i+10]  ->Draw("colz");
 c[5]->cd(3*i+2) ;gPad->SetLogy(0); h_TagBToF[i+10]      ->Draw("colz");
 c[5]->cd(3*i+3) ;gPad->SetLogy(1); h_TagBW[i+10]        ->Draw("");
}

c[6]->Divide(3,10);
for(int i=0;i<10;i++){
 c[6]->cd(3*i+1) ;gPad->SetLogz(1); h2_TagBWvsToF[i+20]  ->Draw("colz");
 c[6]->cd(3*i+2) ;gPad->SetLogy(0); h_TagBToF[i+20]      ->Draw("colz");
 c[6]->cd(3*i+3) ;gPad->SetLogy(1); h_TagBW[i+20]        ->Draw("");
}

c[7]->Divide(3,10);
for(int i=0;i<10;i++){
 c[7]->cd(3*i+1) ;gPad->SetLogz(1); h2_TagBWvsToF[i+30]  ->Draw("colz");
 c[7]->cd(3*i+2) ;gPad->SetLogy(0); h_TagBToF[i+30]      ->Draw("colz");
 c[7]->cd(3*i+3) ;gPad->SetLogy(1); h_TagBW[i+30]        ->Draw("");
}
     
c[8]->Divide(6,20);
for(int i=0;i<120;i++){
 c[8]->cd(i+1) ;gPad->SetLogz(1); h_TagFTDC[i]  ->Draw("");
}
     

if(p_flag==1){
  string ofname_pdf = root_file;
  ofname_pdf.erase(ofname_pdf.size()-5);
  ofname_pdf.append("_TDLRaw.pdf");
  c[0]->Print(Form("./pdf/%s[",ofname_pdf.c_str())  );
  c[0]->Print(Form("./pdf/%s" ,ofname_pdf.c_str())  );
  c[1]->Print(Form("./pdf/%s" ,ofname_pdf.c_str())  );
  c[2]->Print(Form("./pdf/%s" ,ofname_pdf.c_str())  );
  c[3]->Print(Form("./pdf/%s" ,ofname_pdf.c_str())  );
  c[4]->Print(Form("./pdf/%s" ,ofname_pdf.c_str())  );
  c[5]->Print(Form("./pdf/%s" ,ofname_pdf.c_str())  );
  c[6]->Print(Form("./pdf/%s" ,ofname_pdf.c_str())  );
  c[7]->Print(Form("./pdf/%s" ,ofname_pdf.c_str())  );
  c[8]->Print(Form("./pdf/%s" ,ofname_pdf.c_str())  );
  c[8]->Print(Form("./pdf/%s]",ofname_pdf.c_str())  );
 }

} 
