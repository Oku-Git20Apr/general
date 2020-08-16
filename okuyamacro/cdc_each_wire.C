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

void cdc_each_wire (string run_number ) {
  
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


  const double angle_pattern[] = {
	-45, 10, 45, 90
};

  const int NumOfAngle = 4;
//  const int NumOfAngle = sizeof(angle_pattern) / sizeof(angle_pattern[0]);

  double wire[NumOfDCLayers][NumOfAngle] = {0};

  TH1F *h_cdc_wire[NumOfDCLayers*NumOfAngle];
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
  
  for(int n=0;n<ENum;n++){
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

  TCanvas *c[NCanvas];
  for(int i=0;i<NCanvas;i++){
    c[i] =  new TCanvas(Form("c%d",i+1), Form("c%d",i+1), 1400, 900);
  }

  c[0]->Divide(4,5);
  for(int i=8;i<13;i++){
	for(int k=0;k<NumOfAngle;k++){
	c1->cd(4*(i-8)+k+1);h_cdc_wire[4*i+k+1]->Draw();
	}
}
  
  c[1]->Divide(4,5);
  for(int i=13;i<18;i++){
	for(int k=0;k<NumOfAngle;k++){
	c2->cd(4*(i-13)+k+1);h_cdc_wire[4*i+k+1]->Draw();
	}
}

/////////////////////////////////////////////////display all on one Canvas
//  c1->Divide(NumOfAngle, NumOfDCLayers);
//  for(int i=8;i<NumOfDCLayers;i++){ //CDC only
//	for(int k=0;k<NumOfAngle;k++){
//    c1->cd(4*(i-8)+k+1);h_cdc_wire[4*i+k+1]->Draw();
//  }
//}


}
