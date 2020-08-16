const int nCanvas=1;
void sum_ratio() {
 
    cout<<"before file open"<<endl;

//  TFile *ifp = new TFile("../r0000114.root","READONLY");
//  TTree *tree = (TTree*)ifp->Get("tree");
//
//  //Definition for Fill command
//  int tdc[6],qdc[6];
//  tree->SetBranchAddress("tdc", tdc); 
//  tree->SetBranchAddress("qdc", qdc); 

	TChain *tree = new TChain("tree");
	tree->Add("elph0179.root");

    cout<<"after file open"<<endl;

	float pedestal = 200.;
	float ope = 220.;

    int ENum =tree-> GetEntries();
    cout << "Total Event Number is " << ENum << endl;

	char condition[1000];
	sprintf(condition,"(qdc[12]-%f)/%f",pedestal,ope-pedestal);
	
    TH1D *h_qdc[2], *h_pe[2];
	TH2D *h2_qdc_tdc[2];
      h_qdc[0] = new TH1D("h_qdc[0]","h_qdc[0]",150,1050.,1200.);
	  tree->Project("h_qdc[0]","qdc[12]","","");
	  SetTH1(h_qdc[0],"SUM QDC RATIO","QDC [ch]","Counts");
      h_qdc[1] = new TH1D("h_qdc[1]","h_qdc[1]",150,1050.,1200.);
	  tree->Project("h_qdc[1]","qdc[12]","qdc[10]>500&&qdc[10]<800","");
	  SetTH1(h_qdc[1],"SUM QDC RATIO","QDC [ch]","Counts");

//	  SetTH1(h_qdc[0],"Lead Glass QDC","QDC [ch]","Counts");
//      h_pe[0] = new TH1D("h_pe[0]","h_pe[0]",20,-5.,15.);
//	  tree->Project("h_pe[0]",condition1,"","");
//	  SetTH1(h_pe[0],"Number of Photoelectrons","N.P.E.","Counts");
//      h_pe[1] = new TH1D("h_pe[1]","h_pe[1]",20,-5.,15.);
//	  tree->Project("h_pe[1]",condition2,"","");
//	  SetTH1(h_pe[1],"Number of Photoelectrons","N.P.E.","Counts");
//      h_pe[2] = new TH1D("h_pe[2]","h_pe[2]",20,-5.,15.);
//	  tree->Project("h_pe[2]",condition3,"","");
//	  SetTH1(h_pe[2],"Number of Photoelectrons","N.P.E.(SUM)","Counts");
	  

//	  char correlation[1000];
//	  sprintf(correlation,"%s:%s",condition1,condition2);
//	  h2_qdc_tdc[0] = new TH2D("h2_qdc[0]","h2_qdc[0]",20,-5.,20.,20,-5.,20.);
//	  tree->Project("h2_qdc[0]",correlation,"","");
//	  SetTH2(h2_qdc[0],"N.P.E. correlation","N.P.E.2","N.P.E.1");



TCanvas *c[nCanvas];
    for(int i=0;i<nCanvas;i++){
      c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),800,800);
    }
    c[0]->Clear();
	c[0]->cd(1);h_qdc[0]->Draw("");h_qdc[0]->Divide(h_qdc[1]);h_qdc[0]->Draw("");

//    c[0]->Divide(2,3);
//    c[0]->cd(1);h_qdc[0]->Draw("");
//    c[0]->cd(2);h_qdc[1]->Draw("");
//    c[0]->cd(3);h_pe[0]->Draw("");
//    c[0]->cd(4);h_pe[1]->Draw("");
//    c[0]->cd(5);h_pe[2]->Draw("");
//    c[0]->cd(6);h2_qdc[0]->Draw("colz");
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
