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
void gain_fitting_Left(){

  TChain *tree = new TChain("tree");
  tree->Add("elph0179.root");
  TCanvas *c = new TCanvas("c","c",800,600);

  double range = 100.;//Hist range, [peak-range,peak+range]
  double nbins = 200;//nbins of Hist
  double view_range = 4000.;//View range, [peak-view_range,peak+view_range]
  double nbins_of_view = 4000;//nbins of View
  double fit_range = 5.;//Fit range, [peak-fit_range,peak+fit_range]

  TH1D *View,*Hist;     
  View = new TH1D("View","View",nbins_of_view,0.,view_range);
  tree->Project("View","qdc[13]","","");
  double main_part = (View->GetBinCenter(View->GetMaximumBin()))-view_range/nbins_of_view/2;
  Hist = new TH1D("Hist","Hist",nbins,main_part-range,main_part+range);
  tree->Project("Hist","qdc[13]","","");
  SetTH1(Hist,"OKACL Pedestal","QDC [ch]","Counts"); 

  char condi[1000];

  TF1 *f1 = new TF1("f1","gausn",main_part-range,main_part+range);
  double center = (Hist->GetBinCenter(Hist->GetMaximumBin()))-range/(double)nbins;
  cout << "Maximum at " << center << " [ch]" <<endl; 

  c->SetLogy(0);//Linear
  Hist->Draw("");
  Hist->Fit("f1","","",center-fit_range,center+fit_range);
  f1->SetLineStyle(2);
  f1->Draw("same");

  Hist -> SetStats(0); //statistic box does not appear
	

}
