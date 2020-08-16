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
