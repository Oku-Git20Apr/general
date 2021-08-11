#include <TMath.h>
#include <TAttPad.h>

//spherical harmonics
//Aug. 7, 2021

void harmonic(){
  double PI = 4.*atan(1.);
cout<<"PI="<<PI<<endl;

  TH3D *Y00 = new TH3D("Y00","Y_{0}^{0}",100,-3.,3.,100.,-3.,3.,100.,-3.,3.);
  TH3D *Y10 = new TH3D("Y10","Y_{1}^{0}",100,-3.,3.,100.,-3.,3.,100.,-3.,3.);
  TH3D *Y11 = new TH3D("Y11","Y_{1}^{1}",100,-3.,3.,100.,-3.,3.,100.,-3.,3.);
  TH3D *Y1_1 = new TH3D("Y1_1","Y_{1}^{-1}",100,-3.,3.,100.,-3.,3.,100.,-3.,3.);
  TF3 *Y10_f = new TF3("Y10_f","(x*x+y*y+z*z)-(3.*4.*atan(1.)/4.)*z*z/(x*x+y*y+z*z)",-3,3,-3,3,-3,3);
  //TF3 *Y10_f = new TF3("Y10_f","sin(x*x+y*y*z*z-36)",-3,3,-3,3,-3,3);
  //TH3D *h3_test2 = new TH3D("h3_test2","h3_test2",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
  TH3D *h3_test2 = new TH3D("h3_test2","",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
  TH3D *h3_test3 = new TH3D("h3_test3","h3_test3",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
  gStyle->SetOptStat(0);





  //TF1* Y00 = new TF1("Y00","1/sqrt(4.*PI)",-1,1);
  gRandom->SetSeed(0);
  int N1 = 0;
  int N2 = 0;
  int N3 = 0;
  for(int i=0;i<1E+6;i++){

	double theta = acos(gRandom->Uniform(-1.,1.));//UNIFORM
	double phi = 2.*PI*gRandom->Uniform(0.,1.);//UNIFORM
	//double r = 1./sqrt(4.*PI);
	double r = (3.*PI/4.);
	double x1 = r*sin(theta)*cos(phi);
	double y1 = r*sin(theta)*sin(phi);
	double z1 = r*cos(theta);
	Y00->Fill(x1,y1,z1);
	r = (3.*PI/4.)*cos(theta)*cos(theta);
	x1 = r*sin(theta)*cos(phi);
	y1 = r*sin(theta)*sin(phi);
	z1 = r*cos(theta);
	Y10->Fill(x1,y1,z1);
	r = 4.*(3.*PI/8.)*sin(theta)*sin(theta)*sin(phi)*sin(phi);
	x1 = r*sin(theta)*cos(phi);
	y1 = r*sin(theta)*sin(phi);
	z1 = r*cos(theta);
	Y11->Fill(x1,y1,z1);
	r = 4.*(3.*PI/8.)*sin(theta)*sin(theta)*cos(phi)*cos(phi);
	x1 = r*sin(theta)*cos(phi);
	y1 = r*sin(theta)*sin(phi);
	z1 = r*cos(theta);
	Y1_1->Fill(x1,y1,z1);
	//double x1 = gRandom->Uniform(-1.,1.);
	//double y1 = gRandom->Uniform(-1.,1.);
	//double z1 = gRandom->Uniform(-1.,1.);
	//double x2 = sin(theta)*cos(phi);
	//double y2 = sin(theta)*sin(phi);
	//double z2 = cos(theta);
	//double x3 = sin(theta)*cos(phi);
	//double y3 = sin(theta)*sin(phi);
	//double z3 = cos(theta);
	//h3_test2->Fill(x2,y2,z2);
	//h3_test3->Fill(x3,y3,z3);
    //h_test1->Fill(theta_gen1);
    //h_test2->Fill(theta_gen2);
    //h_test3->Fill(theta_gen3);
	//if(theta_gen1<=0.1)N1++;
	//if(theta_gen2<=0.1)N2++;
	//if(theta_gen3<=0.1)N3++;

	
  }

cout<<"N1="<<N1<<endl;
cout<<"N2="<<N2<<endl;
cout<<"N3="<<N3<<endl;
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->Divide(3,2);
  c1->cd(1);
  Y00->SetFillColor(kAzure);
  Y00->SetLineColor(kAzure);
  Y00->SetAxisColor(10,"XYZ");
  Y00->SetLabelColor(10,"XYZ");
  Y00->SetLineColor(10);
  Y00->Draw("isoFBBB");
  c1->cd(4);
  Y10->SetFillColor(kAzure);
  Y10->SetLineColor(kAzure);
  Y10->SetAxisColor(10,"XYZ");
  Y10->SetLabelColor(10,"XYZ");
  Y10->SetLineColor(10);
  Y10->Draw("isoFBBB");
  c1->cd(5);
  Y11->SetFillColor(kAzure);
  Y11->SetLineColor(kAzure);
  Y11->SetAxisColor(10,"XYZ");
  Y11->SetLabelColor(10,"XYZ");
  Y11->SetLineColor(10);
  Y11->Draw("isoFBBB");
  c1->cd(6);
  Y1_1->SetFillColor(kAzure);
  Y1_1->SetLineColor(kAzure);
  Y1_1->SetAxisColor(10,"XYZ");
  Y1_1->SetLabelColor(10,"XYZ");
  Y1_1->SetLineColor(10);
  Y1_1->Draw("isoFBBB");

  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  TH3D *frame = new TH3D("frame","frame",100,-3.,3.,100.,-3.,3.,100.,-3.,3.);
  Y10_f->SetFillColor(10);
  Y10_f->SetLineColor(kAzure);
  frame->SetAxisColor(10,"XYZ");
  frame->SetLabelColor(10,"XYZ");
  frame->Draw("FBBB");
  Y10_f->Draw("FBBBsame");

}
