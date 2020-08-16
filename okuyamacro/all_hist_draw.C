const int nCanvas=1;
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
         pois = par[2]*pow( par[0]/np, np )*exp(-par[0]+np)/(sqrt2pi*pow (np,0.5));
       }
       val += pois/(sqrt2pi*sigma)*exp( -pow(x[0]-np,2)/2./sigma/sigma ); // adding gaussian distribution
  }
    return val;
  }


void all_hist_draw() {
 
    cout<<"before file open"<<endl;

//  TFile *ifp = new TFile("../r0000114.root","READONLY");
//  TTree *tree = (TTree*)ifp->Get("tree");
//
//  //Definition for Fill command
//  int tdc[6],qdc[6];
//  tree->SetBranchAddress("tdc", tdc); 
//  tree->SetBranchAddress("qdc", qdc); 

	TChain *tree = new TChain("tree");
	tree->Add("elph0142.root");

    cout<<"after file open"<<endl;

	double lambda = 1.; //average of #photon
	double res = 0.68; //energy resolution factor
	int roop = 8000.; //menseki
	float pedestal[2] = { 1197., 1228.};
	float ope[2]      = { 1217., 1248.};

    int ENum =tree-> GetEntries();
    cout << "Total Event Number is " << ENum << endl;

	char condition1[1000];
	sprintf(condition1,"(qdc[13]-%f)/%f",pedestal[0],ope[0]-pedestal[0]);
	char condition2[1000];
	sprintf (condition2,"(qdc[14]-%f)/%f",pedestal[1],ope[1]-pedestal[1]);
	char condition3[1000];
	sprintf (condition3,"(qdc[13]-%f)/%f+(qdc[14]-%f)/%f",pedestal[0],ope[0]-pedestal[0],pedestal[1],ope[1]-pedestal[1]);
    
	TH1D *h_qdc[2],*h_qdc_cut[2],*h_pe[3];
	TH2D  *h2_qdc[2];
      h_qdc[0] = new TH1D("h_qdc[0]","h_qdc[0]",200,1150.,1350.);
	  tree->Project("h_qdc[0]","qdc[13]","","");
	  SetTH1(h_qdc[0],"OKACL QDC","QDC [ch]","Counts");
      h_qdc[1] = new TH1D("h_qdc[1]","h_qdc[1]",200,1150.,1350.);
	  tree->Project("h_qdc[1]","qdc[14]","","");
	  SetTH1(h_qdc[1],"OKACR QDC","QDC [ch]","Counts");
      
	  h_qdc_cut[0] = new TH1D("h_qdc_cut[0]","h_qdc_cut[0]",200,1150.,1350.);
	  tree->Project("h_qdc_cut[0]","qdc[13]","tdc[6]>2850&&tdc[6]<2900","");
      h_qdc_cut[1] = new TH1D("h_qdc_cut[1]","h_qdc_cut[1]",200,1150.,1350.);
	  tree->Project("h_qdc_cut[1]","qdc[14]","tdc[6]>2850&&tdc[6]<2900","");
     
	  h_pe[0] = new TH1D("h_pe[0]","h_pe[0]",100,-5.,15.);
	  tree->Project("h_pe[0]",condition1,"tdc[6]>2850&&tdc[6]<2900","");
	  SetTH1(h_pe[0],"Number of Photoelectrons","N.P.E.","Counts");
      h_pe[1] = new TH1D("h_pe[1]","h_pe[1]",100,-5.,15.);
	  tree->Project("h_pe[1]",condition2,"tdc[6]>2850&&tdc[6]<2900","");
	  SetTH1(h_pe[1],"Number of Photoelectrons","N.P.E.","Counts");
      h_pe[2] = new TH1D("h_pe[2]","h_pe[2]",20,-5.,15.);
	  tree->Project("h_pe[2]",condition3,"tdc[6]>2850&&tdc[6]<2900","");
	  SetTH1(h_pe[2],"Number of Photoelectrons","N.P.E.(SUM)","Counts");
	  

	  char correlation[1000];
	  sprintf(correlation,"%s:%s",condition1,condition2);
	  h2_qdc[0] = new TH2D("h2_qdc[0]","h2_qdc[0]",20,-5.,20.,20,-5.,20.);
	  tree->Project("h2_qdc[0]",correlation,"","");
	  SetTH2(h2_qdc[0],"N.P.E. correlation","N.P.E.2","N.P.E.1");

	  TF1 *f1 = new TF1("f1","F_Pois",0.,8.,3);
	  f1->SetNpx(500);
	  f1->SetParameters(lambda,res,roop);

	TCanvas *c[nCanvas];
    for(int i=0;i<nCanvas;i++){
      c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),800,800);
    }
    c[0]->Clear();
    c[0]->Divide(2,3);
    c[0]->cd(1);h_qdc[0]->Draw("");h_qdc_cut[0]->Draw("same");
    c[0]->cd(2);h_qdc[1]->Draw("");h_qdc_cut[1]->Draw("same");
    c[0]->cd(3);h_pe[0]->Draw("");
    c[0]->cd(4);h_pe[1]->Draw("");
    c[0]->cd(5);h_pe[2]->Draw("");h_pe[2]->Fit("f1","","",1.,5.);f1->Draw("same");
    c[0]->cd(6);h2_qdc[0]->Draw("colz");
//    c[1]->Clear();
//    c[1]->Divide(2,2);
//	c[1]->cd(1);h_tdc[0]->Draw("");
//	c[1]->cd(2);h_tdc[1]->Draw("");
//	c[1]->cd(3);h_tdc[2]->Draw("");
//	c[1]->cd(4);h_tdc[3]->Draw("");
//    for(int i=0;i<4;i++){
//      c[1]->cd(i+1);h_tdc[i]->Draw();
//    }
}
