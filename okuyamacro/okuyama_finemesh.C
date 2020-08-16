const int nCanvas=2;
void okuyama_finemesh() {
 
  cout<<"before file open"<<endl;

  TFile *ifp = new TFile("r0000124.root","READONLY");
  TTree *tree = (TTree*)ifp->Get("tree");

  int tdc[6],qdc[6];
  tree->SetBranchAddress("tdc", tdc); 
  tree->SetBranchAddress("qdc", qdc); 
  cout<<"after file open"<<endl;

  int ENum =tree-> GetEntries();
  cout << "Total Event Number is " << ENum << endl;

	float pedestal[2] = {340.,650.};
	float peak_pe[2]  = {500.,740.};

    TH1D *h_qdc[4], *h_tdc[4], *h_pe[4];
	TH2D  *h2_qdc_tdc[4];
    for(int i=0;i<4;i++){
      h_qdc[i] = new TH1D(Form("h_qdc%d",i),Form("h_qdc%d",i),2000,0.,2000.);
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
	c[0]->cd(1);h_qdc[0]->Draw();h_qdc[0]->Fit("gausn");
	c[0]->cd(2);h_qdc[1]->Draw();h_qdc[1]->Fit("gausn");
	c[0]->cd(3);h_qdc[3]->Draw();h_qdc[2]->Fit("gausn");
	c[0]->cd(4);h_qdc[2]->Draw();h_qdc[3]->Fit("gausn");

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
