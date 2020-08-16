 double Noise ( double *x, double *par)  //Fitting function
{
  double val = 0;
  val = par[0];//func.
  return val;
}

 const int NCanvas =2;
 const int NumOfDCLayers = 18;
 const int nwires[18] = {
                    59, 59, 72, 72,  85,  85,  97,  97,          //VDC
                    67, 68, 89, 90, 110, 120, 143, 144, 175, 176 //CDC
                 };//num of wire of each layer
 const double phi0[NumOfDCLayers] = {
-183.06, -183.06, -182.50, -182.50, -182.105, -182.105, -181.849, -159.5, -161.5, -163.875, -157.5, -159.25, -162.0, -163.35, -158.4, -159.5, -162.8, -163.725
};//phi0 parameter from ../../param/dcgeo.param

 const double dphi[NumOfDCLayers] = {
6.102, 6.102, 5.000, 5.000, 4.235, 4.235, 3.711, 3.711, 4.75, 4.75, 3.5, 3.5, 2.7, 2.7, 2.2, 2.2, 1.85, 1.85
};//dphi parameter from ../../param/dcgeo.param

void sn_cdc (string run_number ) {
  
  char root_file[200] = "";  sprintf( root_file, "../../root/all/%s.root" , run_number.c_str() );
cout<< run_number << "  "<<  root_file << endl;
cout<<"input file name is " << root_file <<endl;

//  gROOT -> Reset();
  gROOT -> SetStyle( "Plain" );// Canvas -> clean

//  gStyle -> SetTitleH( 0.30 );  //height of Title box
//  gStyle -> SetTitleX( 0.30 );  //x-coordinate of Title box
  gStyle -> SetStatH( 0.35 );  //height of Stat box
  gStyle -> SetStatW( 0.3 );  //width of Stat box
  gStyle -> SetStatY( 1.0 );  //y-coordinate of Statistical box
  gStyle -> SetOptStat( 1000000001 ); //display title-name only
cout<<"before file open"<<endl;
  TFile* f = new TFile( root_file, "READONLY" );
  TTree *tree = (TTree*)f->Get("tree");
cout<<"after file open"<<endl;

  ///set tree info.//
  int runnum, evnum;
  int dc_tdc[18][176];
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("runnum",1);     tree->SetBranchAddress("runnum"     ,&runnum    );
  tree->SetBranchStatus("evnum",1);      tree->SetBranchAddress("evnum"      ,&evnum     );
  tree->SetBranchStatus("dc_tdc",1);     tree->SetBranchAddress("dc_tdc"     , dc_tdc     );


  const double angle_pattern[] = {  //chosen angle
	-45, 10, 45, 90
};

  const int NumOfAngle = 4;
//  const int NumOfAngle = sizeof(angle_pattern) / sizeof(angle_pattern[0]);

  double wire[NumOfDCLayers][NumOfAngle] = {0};  //Def(WireNum correnspond to the angle)
  TH1F *h_cdc_wire[NumOfDCLayers*NumOfAngle];  //Def(hist) , Title
  for(int i=0;i<NumOfDCLayers;i++){
	for(int k=0;k<NumOfAngle;k++){
    wire[i][k] = (angle_pattern[k]-phi0[i])/dphi[i];
    int w = (int)wire[i][k];
    h_cdc_wire[4*i+k+1] = new TH1F(Form("(Layer,Wire)=(%d,%d)",i+1, w+1),"",1000,0,2000);
  }
}

  int ENum =tree-> GetEntries();
  cout << "Total Event Number is " << ENum << endl;
//  if(ENum>5000)ENum=5000;  //for test
  
  for(int n=0;n<ENum;n++){  //Fill
	if(n%10000==0)cout<<n<<" / "<<ENum<<endl;
  tree->GetEntry(n);
  for(int i=8;i<NumOfDCLayers;i++){ //CDC only
    for(int k=0;k<NumOfAngle;k++){
    wire[i][k] = (angle_pattern[k]-phi0[i])/dphi[i];
    int w = (int)wire[i][k];
    h_cdc_wire[4*i+k+1]->Fill(dc_tdc[i][w]);
    }
} 
}
  cout << "Total Event Number is " << ENum << endl;

  TCanvas *c[NCanvas];  //Def(Canvas)
  for(int i=0;i<NCanvas;i++){
    c[i] =  new TCanvas(Form("c%d",i+1), Form("c%d",i+1), 1400, 900);
  }


  TF1 *fit[NumOfDCLayers][NumOfAngle]; //Def(fit)
  for(int i=0;i<NumOfDCLayers;i++){
	for(int k=0;k<NumOfAngle;k++){
	fit[i][k] = new TF1(Form("fit%d-%d",i,k), Noise, 400., 1200., 1);
	fit[i][k]->SetLineColor(kRed);
	fit[i][k]->SetParameter(0,10.0);
	}
  }

/////////////////////////////////////////////////////////////////Display
  c1->Divide(4,5);
  for(int i=8;i<13;i++){
	for(int k=0;k<NumOfAngle;k++){
	h_cdc_wire[4*i+k+1]->Fit(fit[i][k],"0","",400.,600.);
	double ymax = (h_cdc_wire[4*i+k+1]->GetBinContent(h_cdc_wire[4*i+k+1]->GetMaximumBin()));
	c1->cd(4*(i-8)+k+1);h_cdc_wire[4*i+k+1]->Draw();fit[i][k]->Draw("same");
	TLine *ly1 = new TLine(400.,0.,400.,0.9*ymax);
	TLine *ly2 = new TLine(600.,0.,600.,0.9*ymax);
	TLine *ly3 = new TLine(1050.,0.,1050.,0.9*ymax);
	TLine *ly4 = new TLine(1200.,0.,1200.,0.9*ymax);
	ly1->SetLineColor(kBlue);ly1->Draw();
	ly2->SetLineColor(kBlue);ly2->Draw();
	ly3->SetLineColor(kBlue);ly3->Draw();
	ly4->SetLineColor(kBlue);ly4->Draw();
	double p0 = fit[i][k]->GetParameter(0);
	double yhmax = (h_cdc_wire[4*i+k+1]->GetMaximum());
	double yc = yhmax/2;
	if(yc == 0) yc = 0.5; //when data is empty
	double Noise = p0*450.;
	double DummySignal = 0.;
	double Signal = 0.;
	for(int b=300;b<525;b++){
	DummySignal += 2*(h_cdc_wire[4*i+k+1]->GetBinContent(b));
	}
	double Signal = DummySignal - Noise;
	double SN = Signal/Noise;
	cout << "S/N = " << SN <<endl;
	gPad->Update();
	//TPaveStats *ps = (TPaveStats*)h_cdc_wire[4*i+k+1]->FindObject("stats");
	//TList *listOfLines = ps->GetListOfLines();
	if (yhmax == 0){
	TLatex *myt = new TLatex(1000,yc,"S/N is Not Defined");
	}
	else{
	TLatex *myt = new TLatex(1400,yc,Form("S/N=%.2f",SN));
	}
	//listOfLines->Add(myt);
	myt->SetTextSize(0.1);
	myt->Draw();
	}
}
  
  c2->Divide(4,5);
  for(int i=13;i<18;i++){
	for(int k=0;k<NumOfAngle;k++){
    h_cdc_wire[4*i+k+1]->Fit(fit[i][k],"0","",400, 600);
	double ymax = (h_cdc_wire[4*i+k+1]->GetBinContent(h_cdc_wire[4*i+k+1]->GetMaximumBin()));
	c2->cd(4*(i-13)+k+1);h_cdc_wire[4*i+k+1]->Draw();fit[i][k]->Draw("same");
	TLine *ly1 = new TLine(400.,0.,400.,0.9*ymax);
	TLine *ly2 = new TLine(600.,0.,600.,0.9*ymax);
	TLine *ly3 = new TLine(1050.,0.,1050.,0.9*ymax);
	TLine *ly4 = new TLine(1200.,0.,1200.,0.9*ymax);
	ly1->SetLineColor(kBlue);ly1->Draw();
	ly2->SetLineColor(kBlue);ly2->Draw();
	ly3->SetLineColor(kBlue);ly3->Draw();
	ly4->SetLineColor(kBlue);ly4->Draw();
	double p0 = fit[i][k]->GetParameter(0);
	double yhmax = (h_cdc_wire[4*i+k+1]->GetMaximum());
	double yc = yhmax/2;
	if(yc == 0) yc = 0.5; //when data is empty
	double Noise = p0*500.;
	double DummySignal = 0.;
	double Signal = 0.;
	for(int b=300;b<525;b++){
	DummySignal += 2*(h_cdc_wire[4*i+k+1]->GetBinContent(b));
	}
	double Signal = DummySignal - Noise;
	double SN = Signal/Noise;
	cout << "S/N = " << SN <<endl;
	gPad->Update();
	//TPaveStats *ps = (TPaveStats*)h_cdc_wire[4*i+k+1]->FindObject("stats");
	//TList *listOfLines = ps->GetListOfLines();
	if (yhmax == 0){
	TLatex *myt = new TLatex(1000,yc,"S/N is Not Defined");
	}
	else{
	TLatex *myt = new TLatex(1400,yc,Form("S/N=%.2f",SN));
	}
	//listOfLines->Add(myt);
	myt->SetTextSize(0.1);
	myt->Draw();
	}
}

/////////////////////////////////////////////////display all on one Canvas
//  c1->Divide(NumOfAngle, NumOfDCLayers);
//  for(int i=8;i<NumOfDCLayers;i++){ //CDC only
//	for(int k=0;k<NumOfAngle;k++){
//    c1->cd(4*(i-8)+k+1);h_cdc_wire[4*i+k+1]->Draw();f1->Draw("same");
//  }
//}

}
