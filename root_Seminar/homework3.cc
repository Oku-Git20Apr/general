void homework3(){
  TH3D *h = new TH3D("h","h",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
//  SetTH1(h,"homework3","x","y");

//	double r = 1.;
  gRandom->SetSeed(0);
  for(int i=0;i<1E+6;i++){
	double r = pow(3*(gRandom->Uniform(0.,1.)),1/3);
	double theta = acos( gRandom->Uniform(-1.,1.));
	double phi = gRandom->Uniform(0.,2.*TMath::Pi());
	double x = r*sin(theta)*cos(phi);
	double y = r*sin(theta)*sin(phi);
	double z = r*cos(theta);
    h->Fill(x,y,z);
  }

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  h->Draw();

}
