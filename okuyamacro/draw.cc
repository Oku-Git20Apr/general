//--------------------------------------------------------------------
void SetTH1(TH1 *h, TString name, TString xname, TString yname, int LColor=1, int FStyle=0, int FColor=0){
  h->SetTitle(name);
  h->SetLineColor(LColor);
  h->SetLineWidth(0);
  h->SetTitleSize(0.04,"");
  h->SetTitleFont(42,"");
  h->SetFillStyle(FStyle);
  h->SetFillColor(FColor);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(0.90);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.00);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
}

//--------------------------------------------------------------------
void init(){
  gROOT->SetStyle("Plain");
  gStyle->SetOptDate(0);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetLineWidth(0);
  gStyle->SetOptDate(0);
  gStyle->SetTextFont(42);
  gStyle->SetGridWidth(0);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetNdivisions(510); // tertiary*10000 + secondary*100 + first
// Stat box
  gStyle->SetOptStat("ei");
  gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFont(42);
// Pad
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadBottomMargin(0.13);
// Title
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleFontSize(0.04);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFont(42,"XYZ");
// Label
  gStyle->SetStripDecimals(kFALSE);
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetLabelOffset(0.012,"X");

//  const Int_t NRGBs = 5;
//  const Int_t NCont = 99;
//  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
//  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
//  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
//  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
//  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
//  gStyle->SetNumberContours(NCont);
}

void draw(){
  init();

  TChain *tree = new TChain("tree");
  tree->Add("r0000138.root");


  TH1D *h[6];
  for(int i=0;i<6;i++){
    h[i] = new TH1D(Form("h%d",i),Form("h%d",i),1000,0,4000);
    SetTH1(h[i],Form("ch%0",i),"QDC [ch] (0.1pC/ch)","counts / 2ch");
    tree->Project(Form("h%d",i),Form("qdc[%d]",i),"","");
  }

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(3,2);
  for(int i=0;i<6;i++){
    c1->cd(i+1); gPad->SetLogy(); h[i]->Draw();
  }
}
