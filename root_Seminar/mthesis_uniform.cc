#include <TMath.h>
void mthesis_uniform(){
  double PI = 4.*atan(1.);

  TH3D *h3_test1 = new TH3D("h3_test1","h3_test1",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
  //TH3D *h3_test2 = new TH3D("h3_test2","h3_test2",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
  TH3D *h3_test2 = new TH3D("h3_test2","",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
  h3_test2->SetNdivisions(505,"XYZ");
  h3_test2->GetXaxis()->SetTitle("x");
  h3_test2->GetXaxis()->SetTitleSize(0.06);
  h3_test2->GetXaxis()->CenterTitle();
  h3_test2->GetYaxis()->SetTitle("y");
  h3_test2->GetYaxis()->SetTitleSize(0.06);
  h3_test2->GetYaxis()->CenterTitle();
  h3_test2->GetZaxis()->SetTitle("z");
  h3_test2->GetZaxis()->SetTitleSize(0.06);
  h3_test2->GetZaxis()->SetTitleOffset(0.8);
  h3_test2->GetZaxis()->CenterTitle();
  TH3D *h3_test3 = new TH3D("h3_test3","h3_test3",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
  TH1D *h_test1  = new TH1D("h_test1","h_test1",1000,0.,PI);
  TH1D *h_test2  = new TH1D("h_test2","h_test2",1000,0.,PI);
  TH1D *h_test3  = new TH1D("h_test3","h_test3",1000,0.,PI);
//  SetTH1(h,"homework3","x","y");
  gStyle->SetOptStat(0);

cout<<"PI="<<PI<<endl;
  gRandom->SetSeed(0);
  int N1 = 0;
  int N2 = 0;
  int N3 = 0;
  for(int i=0;i<1E+6;i++){

	double theta_gen1 = acos( -(2.0*gRandom->Uniform()*-1.0)*gRandom->Uniform());//SIMC
	double theta_gen2 = acos(gRandom->Uniform(-1.,1.));//UNIFORM
	double theta_gen3 = acos(1.-gRandom->Uniform()*(1.-cos(0.1)));//uniform in 0--pi/2 rad
	double phi_gen = 2.*PI*gRandom->Uniform(0.,1.);
	double x1 = sin(theta_gen1)*cos(phi_gen);
	double y1 = sin(theta_gen1)*sin(phi_gen);
	double z1 = cos(theta_gen1);
	double x2 = sin(theta_gen2)*cos(phi_gen);
	double y2 = sin(theta_gen2)*sin(phi_gen);
	double z2 = cos(theta_gen2);
	double x3 = sin(theta_gen3)*cos(phi_gen);
	double y3 = sin(theta_gen3)*sin(phi_gen);
	double z3 = cos(theta_gen3);
	h3_test1->Fill(x1,y1,z1);
	h3_test2->Fill(x2,y2,z2);
	h3_test3->Fill(x3,y3,z3);
    h_test1->Fill(theta_gen1);
    h_test2->Fill(theta_gen2);
    h_test3->Fill(theta_gen3);
	if(theta_gen1<=0.1)N1++;
	if(theta_gen2<=0.1)N2++;
	if(theta_gen3<=0.1)N3++;

	
  }

cout<<"N1="<<N1<<endl;
cout<<"N2="<<N2<<endl;
cout<<"N3="<<N3<<endl;
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  h_test1->Scale((double)N2/(double)N1);
  h_test3->Scale((double)N2/(double)N3);
  //h_test1->Draw();
  h_test2->SetLineColor(kRed);
  h_test3->SetLineColor(kGreen);
  h_test2->Draw("same");
  h_test3->Draw("same");

  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  h3_test2->Draw("");
//  c3->Print("uniform_3d.pdf");

}
