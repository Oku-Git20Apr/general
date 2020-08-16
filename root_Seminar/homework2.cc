void homework2(){
  TH1D *h = new TH1D("h","h",100,0.,10);
  SetTH1(h,"homework2","x","y");

  double e = 2.718; //Napier's constant
  gRandom->SetSeed(0);
  for(int i=0;i<1E+6;i++){
    double x = gRandom->Uniform(0.,10.);
	double y = gRandom->Uniform(0.,0.3); //greater than Max[Poisson(x,2.)]
	double fx = TMath::Poisson(x,2.);
    if(y < fx) h->Fill(x);
  }

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  h->Draw();

}
