const int nCanvas=2;

//void SetTH1(TH1 *h, TString name, TString xname, TString yname, int LC    olor=1, int FStyle=0, int FColor=0){
//   h->SetTitle(name);
//   h->SetLineColor(LColor);
//   h->SetLineWidth(0);
//   h->SetTitleSize(0.04,"");
//   h->SetTitleFont(42,"");
//   h->SetFillStyle(FStyle);
//   h->SetFillColor(FColor);
// 
//   h->GetXaxis()->SetTitle(xname);
//   h->GetXaxis()->CenterTitle();
//   h->GetXaxis()->SetTitleFont(42);
//   h->GetXaxis()->SetTitleOffset(0.90);
//   h->GetXaxis()->SetTitleSize(0.06);
//   h->GetXaxis()->SetLabelFont(42);
//   h->GetXaxis()->SetLabelOffset(0.01);
// 
//   h->GetYaxis()->SetTitle(yname);
//   h->GetYaxis()->CenterTitle();
//   h->GetYaxis()->SetTitleFont(42);
//   h->GetYaxis()->SetTitleOffset(1.00);
//   h->GetYaxis()->SetTitleSize(0.06);
//   h->GetYaxis()->SetLabelFont(42);
//   h->GetYaxis()->SetLabelOffset(0.01);
//   ((TGaxis*)h->GetYaxis())->SetMaxDigits(4);
//}

 double sqrt2pi = sqrt(2.*3.141592);
 double F_Pois( double *x, double *par )
 {
   /*
     par[0] : lambda, average of #photon
     par[1] : energy resolution factor
     par[2] : menseki
   */
   double val = 0;
   for( int np=1; np<200; np++ ){
     double pois;
     double sigma = par[1]*sqrt(np);
     if(np<50){
       pois = par[2]*pow( par[0], np )*exp(-par[0])/TMath::Gamma(np+1);
     }
     else{ // stirling's approximation
       pois = par[2]*pow( par[0]/np, np )*exp(-par[0]+np)/(sqrt2pi*pow(    np,0.5));
     }
     val += pois/(sqrt2pi*sigma)*exp( -pow(x[0]-np,2)/2./sigma/sigma );     // adding gauss    ian distribution
   }
   return val;
}



void gain_npe() {
 
  cout<<"before file open"<<endl;

  TFile *ifp = new TFile("r0000073.root","READONLY");
  TTree *tree = (TTree*)ifp->Get("tree");
 // Setting *set = new Setting();

  int tdc[6],qdc[6];
  tree->SetBranchAddress("tdc", tdc); 
  tree->SetBranchAddress("qdc", qdc); 

//  gROOT -> Reset();
//  gROOT -> SetStyle( "Plain" );// Canvas -> clean

  cout<<"after file open"<<endl;


  int ENum =tree-> GetEntries();
  cout << "Total Event Number is " << ENum << endl;


	double lambda = 4.; //average of #photon
	double res = 0.68; //energy resolution factor
	int roop = 200.; //menseki

	float pedestal[2] = {340.,650.};
	float peak_pe[2]  = {500.,740.};

    TH1D *h_qdc[4], *h_tdc[4], *h_pe[4];
	TH2D  *h2_qdc_tdc[4];
    for(int i=0;i<4;i++){
      h_qdc[i] = new TH1D(Form("h_qdc%d",i),Form("h_qdc%d",i),500,0.,2000.);
      h_tdc[i] = new TH1D(Form("h_tdc%d",i),Form("h_tdc%d",i),125,500.,1000.);
      h_pe[i] = new TH1D(Form("h_pe%d",i),Form("h_pe%d",i),125,-5.,20.);
	  h2_qdc_tdc[i] = new TH2D(Form("h2_qdc_tdc%d",i),Form("h2_qdc_tdc%d",i),500,0.,2000.,500,0.,2000.);
}
 
  for(int n=0;n<ENum;n++){
	if(n%1000==0)cout<<n<<" / "<<ENum<<endl;
     tree->GetEntry(n);

       for(int i=0;i<4;i++){
         h_qdc[i] ->Fill(qdc[i]);
	
//         h_tdc[i] ->Fill(tdc[i]); //tdc numbering is little bit different
       }
	h_tdc[0] -> Fill(tdc[0]);
	h_tdc[1] -> Fill(tdc[1]);
	h_tdc[2] -> Fill(tdc[2]);
	h_tdc[3] -> Fill(tdc[3]);
	h2_qdc_tdc[0] -> Fill(qdc[1],tdc[0]);//for AC1
	h2_qdc_tdc[1] -> Fill(qdc[2],tdc[1]);//for AC1
     
	h_pe[0] -> Fill((qdc[1]-pedestal[0])/(peak_pe[0]-pedestal[0])+(qdc[2]-pedestal[1])/(peak_pe[1]-pedestal[1])); 
	h2_qdc_tdc[2] -> Fill(qdc[1],qdc[2]);//for AC1
}

TCanvas *c[nCanvas];
    for(int i=0;i<nCanvas;i++){
      c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),800,800);
    }
    c[0]->Clear();
    c[0]->Divide(2,2);
	c[0]->cd(1);h_qdc[0]->Draw();
	c[0]->cd(2);h_qdc[1]->Draw();
	c[0]->cd(3);h_qdc[3]->Draw();
	c[0]->cd(4);h_qdc[2]->Draw();
//    for(int i=0;i<4;i++){
//      c[0]->cd(i+1);h_qdc[i]->Draw();
//    }
    c[1]->Clear();
    c[1]->Divide(2,2);
	c[1]->cd(1);h_pe[0]->Draw();
	c[1]->cd(2);h2_qdc_tdc[2]->Draw("colz");//qdc_qdc_corl
//	c[1]->cd(1);h_tdc[0]->Draw();
//	c[1]->cd(2);h_tdc[1]->Draw();
	c[1]->cd(3);h2_qdc_tdc[0]->Draw("colz");
	c[1]->cd(4);h2_qdc_tdc[1]->Draw("colz");
//    for(int i=0;i<4;i++){
//      c[1]->cd(i+1);h_tdc[i]->Draw();
//    }
}
