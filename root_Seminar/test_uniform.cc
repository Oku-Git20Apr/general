#include <TMath.h>
void test_uniform(){
  TH3D *h = new TH3D("h","h",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
  TH3D *h_cm = new TH3D("h_cm","h_cm",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
  TH1D *h_test = new TH1D("h_test","h_test",1000,-4.,4.);
  TH1D *h_test2 = new TH1D("h_test2","h_test2",1000,-4.,4.);
//  SetTH1(h,"homework3","x","y");

  gRandom->SetSeed(0);
  int N_cm = 0;
  int N_lab = 0;
  for(int i=0;i<1E+6;i++){

	double theta_gen = acos( (2.*gRandom->Uniform()-1.)*gRandom->Uniform());
	double theta_gen2 = acos(gRandom->Uniform());
    h_test->Fill(theta_gen);
    h_test2->Fill(theta_gen2);

	
  }

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  h_test->Draw();
  h_test2->Draw("same");

}
