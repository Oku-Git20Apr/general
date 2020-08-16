void homework1(){
  TH1D *h = new TH1D("h","h",100,-1.5,1.5);
  SetTH1(h,"homework1","x","y");

  gRandom->SetSeed(0);
  for(int i=0;i<1E+5;i++){
    double x = -1+TMath::Sqrt(4*gRandom->Uniform(-1.,1.));
    h->Fill(x);
  }

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  h->Draw();

}
