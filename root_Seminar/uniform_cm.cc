#include <TMath.h>
void uniform_cm(){
  TH3D *h = new TH3D("h","h",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
  TH3D *h_cm = new TH3D("h_cm","h_cm",100,-1.,1.,100.,-1.,1.,100.,-1.,1.);
//  SetTH1(h,"homework3","x","y");

  gRandom->SetSeed(0);
  int N_cm = 0;
  int N_lab = 0;
  for(int i=0;i<1E+6;i++){
//	double r = pow(3*(gRandom->Uniform(0.,1.)),1/3);
	double r = 1.;
	double beta = 0.74;
	double gamma = 1./sqrt(1.-beta*beta);
	double theta_gen_cm = acos( gRandom->Uniform(-1.,1.));
	double phi_gen_cm = gRandom->Uniform(0.,2.*TMath::Pi());
 	double theta_gk = 1.20*3.14/180.0;
 	double phi_gen_n = sin(theta_gen_cm)*sin(phi_gen_cm)*cos(theta_gk)-cos(theta_gen_cm)*sin(theta_gk);
 	double phi_gen_d = sin(theta_gen_cm)*cos(theta_gen_cm);
 	double phi = atan(phi_gen_n/phi_gen_d);
 	double theta_gen_n1 = pow(sin(theta_gen_cm)*cos(phi_gen_cm),2.);
 	double theta_gen_n2 = pow((sin(theta_gen_cm)*cos(phi_gen_cm)*cos(theta_gk)-cos(theta_gen_cm)*sin(theta_gk)),2.);
 	double theta_gen_d1 = sin(theta_gen_cm)*sin(phi_gen_cm)*sin(theta_gk);
 	double theta_gen_d2 = cos(theta_gen_cm)*cos(theta_gk);
 	double theta = atan((theta_gen_n1*theta_gen_n2)/(theta_gen_d1+theta_gen_d2));
	double x_cm = r*sin(theta_gen_cm)*cos(phi_gen_cm);
	double y_cm = r*sin(theta_gen_cm)*sin(phi_gen_cm);
	double z_cm = r*cos(theta_gen_cm);
	//double x = r*sin(theta)*cos(phi);
	//double y = r*sin(theta)*sin(phi);
	//double z = r*cos(theta);
	double E = 494./1800.;
	double x = r*sin(theta_gen_cm)*cos(phi_gen_cm);
	double y = r*sin(theta_gen_cm)*sin(phi_gen_cm)*cos(theta_gk)-r*(gamma*beta*E+gamma*cos(theta_gk))*sin(theta_gk);
	double z = r*sin(theta_gen_cm)*sin(phi_gen_cm)*sin(theta_gk)+r*(gamma*beta*E+gamma*cos(theta_gk))*cos(theta_gk);
	double r_cm = sqrt(x*x+y*y+z*z);
	x /= r_cm;
	y /= r_cm;
	z /= r_cm;

	double theta_test = acos(z);
	//double theta_test = acos(z/r_cm);
	double phi_test = TMath::ATan2(y,x);
	double x_test = r_cm*sin(theta_test)*cos(phi_test);
	double y_test = r_cm*sin(theta_test)*sin(phi_test);
	double z_test = r_cm*cos(theta_test);
	//x -= x_test;
	//y -= y_test;
	//z -= z_test;

    h->Fill(x,y,z);
    //h->Fill(x_test,y_test,z_test);
    h_cm->Fill(x_cm,y_cm,z_cm);

	// 6msr --> theta=0.043
	if(theta_gen_cm<0.043)N_cm++;
	if(theta_test<0.043)N_lab++;
	
  }
	cout<<"N_cm="<<N_cm<<endl;
	cout<<"N_lab="<<N_lab<<endl;

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  h->Draw();
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  c2->cd();
  h_cm->Draw();

}
