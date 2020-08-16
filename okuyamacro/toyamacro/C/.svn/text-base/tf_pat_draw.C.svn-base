#include "../Setting.cc"
const int NCanvas =3;

void tf_pat_draw(string root_file="../rf_10634.root"){
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetOptStat(0);

  cout<<root_file<<endl;
  TFile* ifp = new TFile( root_file.c_str(), "READONLY" );
  TLatex *latex = new TLatex();
  Setting *set = new Setting();

  TCanvas *c[NCanvas];
  for(int i=0;i<NCanvas;i++){
    c[i]=new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800,800);
  }

  TH1F *h_tagfhitpat_tb[40];
  for(int i=0;i<24;i++){
    h_tagfhitpat_tb[i] = (TH1F*)ifp->Get(Form("h_tagfhitpat_tb%d",i+1));
    h_tagfhitpat_tb[i] -> SetTitle(Form("TagF hit pattern(w/ TagB%d hit)"          ,i+1));

    int tf_min = i*4 -8;
    int tf_max = i*4 +12;
    h_tagfhitpat_tb[i] ->GetXaxis()->SetRangeUser(tf_min,tf_max);
  }

  for(int i=0;i<3;i++){ //Canvas
    c[i]->Divide(4,2);
    for(int k=0;k<8;k++){ //TagB seg
      int seg = 8*i + k +1;
      c[i]->cd(1+k);gPad->SetLogy(0);h_tagfhitpat_tb[seg-1]->Draw("");
      if(seg%2==0)h_tagfhitpat_tb[seg-2]->Draw("same");
    }
  }

  //save pdf
  string pdf_file = root_file;
  pdf_file.erase(pdf_file.size()-5);
  pdf_file.append("_tf_hitpat.pdf");
  c[0]->Print(Form("%s[",pdf_file.c_str()));
  for(int i=0;i<NCanvas;i++)  c[i]->Print(Form("%s",pdf_file.c_str()));
  c[NCanvas-1]->Print(Form("%s]",pdf_file.c_str()));

}
