#include <TMath.h>
#include <TAttPad.h>

//Description of triaxial deformed nuclei
//with spherical harmonics
//Aug. 11, 2021 

double fY00(double *var, double *par){
	double x = var[0];
	double y = var[1];
	double z = var[2];
	return (x*x+y*y+z*z)-(1./(4.*4.*atan(1.))); 
}
double fY10(double *var, double *par){
	double x = var[0];
	double y = var[1];
	double z = var[2];
	return (x*x+y*y+z*z)-(3./(4.*4.*atan(1.)))*z*z/(x*x+y*y+z*z); 
}
double fY11p(double *var, double *par){
	double x = var[0];
	double y = var[1];
	double z = var[2];
	return (x*x+y*y+z*z)-(3./(4.*4.*atan(1.)))*(y*y/(x*x+y*y+z*z)); 
	//return (x*x+y*y+z*z)-(3./(4.*4.*atan(1.)))*(y*y/(x*x+y*y+z*z));//other expression
}
double fY11m(double *var, double *par){
	double x = var[0];
	double y = var[1];
	double z = var[2];
	return (x*x+y*y+z*z)-(3./(4.*4.*atan(1.)))*(x*x/(x*x+y*y+z*z)); 
}
double fY20(double *var, double *par){
	double x = var[0];
	double y = var[1];
	double z = var[2];
	return (x*x+y*y+z*z)-(5./(16.*4.*atan(1.)))*(-x*x-y*y+2.*z*z)*(-x*x-y*y+2.*z*z)/(x*x+y*y+z*z)/(x*x+y*y+z*z); 
}
double fY21p(double *var, double *par){
	double x = var[0];
	double y = var[1];
	double z = var[2];
	return (x*x+y*y+z*z)-(15./(8.*4.*atan(1.)))*y*y*z*z/(x*x+y*y+z*z)/(x*x+y*y+z*z); 
}
double fY21m(double *var, double *par){
	double x = var[0];
	double y = var[1];
	double z = var[2];
	return (x*x+y*y+z*z)-(15./(8.*4.*atan(1.)))*x*x*z*z/(x*x+y*y+z*z)/(x*x+y*y+z*z); 
}
double fY22p(double *var, double *par){
	double x = var[0];
	double y = var[1];
	double z = var[2];
	return (x*x+y*y+z*z)-(15./(32.*4.*atan(1.)))*4.*x*x*y*y/(x*x+y*y+z*z)/(x*x+y*y+z*z); 
}
double fY22m(double *var, double *par){
	double x = var[0];
	double y = var[1];
	double z = var[2];
	return (x*x+y*y+z*z)-(15./(32.*4.*atan(1.)))*(x*x-y*y)*(x*x-y*y)/(x*x+y*y+z*z)/(x*x+y*y+z*z); 
}
double fR(double *var, double *par){
	double x = var[0];
	double y = var[1];
	double z = var[2];
	double beta = par[0];
	double gamma = par[1];
	double R0 = par[2];
	double pi = 4.*atan(1.);
	return (x*x+y*y+z*z)-sqrt(R0)*pow(1.+beta*cos(gamma)*sqrt(5./16./pi)*(-x*x-y*y+2.*z*z)/(x*x+y*y+z*z)+(1./sqrt(2.))*beta*sin(gamma)*sqrt(15./32./pi)*2.*x*y/(x*x+y*y+z*z),2.); 
}

void triaxial(){
  double PI = 4.*atan(1.);
cout<<"PI="<<PI<<endl;

  TF3 *Y00 = new TF3("Y00",fY00,-0.5,0.5,-0.5,0.5,-0.5,0.5);
  TF3 *Y10 = new TF3("Y10",fY10,-0.5,0.5,-0.5,0.5,-0.5,0.5);
  TF3 *Y11p = new TF3("Y11p",fY11p,-0.5,0.5,-0.5,0.5,-0.5,0.5);
  TF3 *Y11m = new TF3("Y11m",fY11m,-0.5,0.5,-0.5,0.5,-0.5,0.5);
  TF3 *Y20 = new TF3("Y20",fY20,-0.5,0.5,-0.5,0.5,-1.0,1.0);
  TF3 *Y21p = new TF3("Y21p",fY21p,-0.5,0.5,-0.5,0.5,-0.5,0.5);
  TF3 *Y21m = new TF3("Y21m",fY21m,-0.5,0.5,-0.5,0.5,-0.5,0.5);
  TF3 *Y22p = new TF3("Y22p",fY22p,-0.5,0.5,-0.5,0.5,-0.5,0.5);
  TF3 *Y22m = new TF3("Y22m",fY22m,-0.5,0.5,-0.5,0.5,-0.5,0.5);
  TF3 *nuclei = new TF3("nuclei",fR,-0.5,0.5,-0.5,0.5,-0.5,0.5,3);
  nuclei->SetParameter(0,1.);//beta
  nuclei->SetParameter(1,PI/3.);//gamma
  nuclei->SetParameter(2,0.01);//scale, R0
  double x[2]={-0.75,0.75};
  double y[2]={-0.75,0.75};
  double z[2]={-0.75,0.75};
  double z_long[2]={-1.5,1.5};
  double null[2]={0,0};
  TPolyLine3D *xaxis = new TPolyLine3D(2,x,null,null);
  TPolyLine3D *yaxis = new TPolyLine3D(2,null,y,null);
  TPolyLine3D *zaxis = new TPolyLine3D(2,null,null,z);
  TPolyLine3D *zaxis_long = new TPolyLine3D(2,null,null,z_long);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.1);



  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  TH3D *frame_y00 = new TH3D("frame_y00","Y_{0}^{0}",100,-0.5,0.5,100.,-0.5,0.5,100.,-0.5,0.5);
  frame_y00->SetAxisColor(10,"XYZ");
  frame_y00->SetLabelColor(10,"XYZ");
  frame_y00->Draw("FBBB");
  Y00->SetFillColor(10);
  Y00->SetLineColor(kAzure);
  Y00->Draw("FBBBsame");
  xaxis->Draw("FBBBsame");
  yaxis->Draw("FBBBsame");
  zaxis->Draw("FBBBsame");
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(3,1);
  c2->cd(1);
  TH3D *frame_y10 = new TH3D("frame_y10","Y_{1}^{ 0}",100,-0.5,0.5,100.,-0.5,0.5,100.,-0.5,0.5);
  frame_y10->SetAxisColor(10,"XYZ");
  frame_y10->SetLabelColor(10,"XYZ");
  frame_y10->Draw("FBBB");
  Y10->SetFillColor(10);
  Y10->SetLineColor(kAzure);
  Y10->Draw("FBBBsame");
  xaxis->Draw("FBBBsame");
  yaxis->Draw("FBBBsame");
  zaxis->Draw("FBBBsame");
  c2->cd(2);
  TH3D *frame_y11 = new TH3D("frame_y11","Y_{1}^{ 1}+Y_{1}^{ -1}",100,-0.5,0.5,100.,-0.5,0.5,100.,-0.5,0.5);
  frame_y11->SetAxisColor(10,"XYZ");
  frame_y11->SetLabelColor(10,"XYZ");
  frame_y11->Draw("FBBB");
  Y11p->SetFillColor(10);
  Y11p->SetLineColor(kAzure);
  Y11p->Draw("FBBBsame");
  xaxis->Draw("FBBBsame");
  yaxis->Draw("FBBBsame");
  zaxis->Draw("FBBBsame");
  c2->cd(3);
  TH3D *frame_y11_ = new TH3D("frame_y11_","Y_{1}^{ 1}-Y_{1}^{ 1}",100,-0.5,0.5,100.,-0.5,0.5,100.,-0.5,0.5);
  frame_y11_->SetAxisColor(10,"XYZ");
  frame_y11_->SetLabelColor(10,"XYZ");
  frame_y11_->Draw("FBBB");
  Y11m->SetFillColor(10);
  Y11m->SetLineColor(kAzure);
  Y11m->Draw("FBBBsame");
  xaxis->Draw("FBBBsame");
  yaxis->Draw("FBBBsame");
  zaxis->Draw("FBBBsame");
  TCanvas *c3 = new TCanvas("c3","c3",2000,600);
  c3->Divide(5,1);
  c3->cd(1);
  TH3D *frame_y20 = new TH3D("frame_y20","Y_{2}^{ 0}",100,-0.5,0.5,100.,-0.5,0.5,100.,-1.0,1.0);
  frame_y20->SetAxisColor(10,"XYZ");
  frame_y20->SetLabelColor(10,"XYZ");
  frame_y20->Draw("FBBB");
  Y20->SetFillColor(10);
  Y20->SetLineColor(kAzure);
  Y20->Draw("FBBBsame");
  xaxis->Draw("FBBBsame");
  yaxis->Draw("FBBBsame");
  zaxis_long->Draw("FBBBsame");
  c3->cd(2);
  TH3D *frame_y21p = new TH3D("frame_y21p","Y_{2}^{ 1}+Y_{2}^{ -1}",100,-0.5,0.5,100.,-0.5,0.5,100.,-0.5,0.5);
  frame_y21p->SetAxisColor(10,"XYZ");
  frame_y21p->SetLabelColor(10,"XYZ");
  frame_y21p->Draw("FBBB");
  Y21p->SetFillColor(10);
  Y21p->SetLineColor(kAzure);
  Y21p->Draw("FBBBsame");
  xaxis->Draw("FBBBsame");
  yaxis->Draw("FBBBsame");
  zaxis->Draw("FBBBsame");
  c3->cd(3);
  TH3D *frame_y21m = new TH3D("frame_y21m","Y_{2}^{ 1}-Y_{2}^{ -1}",100,-0.5,0.5,100.,-0.5,0.5,100.,-0.5,0.5);
  frame_y21m->SetAxisColor(10,"XYZ");
  frame_y21m->SetLabelColor(10,"XYZ");
  frame_y21m->Draw("FBBB");
  Y21m->SetFillColor(10);
  Y21m->SetLineColor(kAzure);
  Y21m->Draw("FBBBsame");
  xaxis->Draw("FBBBsame");
  yaxis->Draw("FBBBsame");
  zaxis->Draw("FBBBsame");
  c3->cd(4);
  TH3D *frame_y22p = new TH3D("frame_y22p","Y_{2}^{ 2}+Y_{2}^{ -2}",100,-0.5,0.5,100.,-0.5,0.5,100.,-0.5,0.5);
  frame_y22p->SetAxisColor(10,"XYZ");
  frame_y22p->SetLabelColor(10,"XYZ");
  frame_y22p->Draw("FBBB");
  Y22p->SetFillColor(10);
  Y22p->SetLineColor(kAzure);
  Y22p->Draw("FBBBsame");
  xaxis->Draw("FBBBsame");
  yaxis->Draw("FBBBsame");
  zaxis->Draw("FBBBsame");
  c3->cd(5);
  TH3D *frame_y22m = new TH3D("frame_y22m","Y_{2}^{ 2}-Y_{2}^{ -2}",100,-0.5,0.5,100.,-0.5,0.5,100.,-0.5,0.5);
  frame_y22m->SetAxisColor(10,"XYZ");
  frame_y22m->SetLabelColor(10,"XYZ");
  frame_y22m->Draw("FBBB");
  Y22m->SetFillColor(10);
  Y22m->SetLineColor(kAzure);
  Y22m->Draw("FBBBsame");
  xaxis->Draw("FBBBsame");
  yaxis->Draw("FBBBsame");
  zaxis->Draw("FBBBsame");

  TCanvas *c4 = new TCanvas("c4","c4",800,800);
  TH3D *frame_nuclei = new TH3D("frame_nuclei","#beta=0.5, #gamma=#pi/4",100,-0.5,0.5,100.,-0.5,0.5,100.,-0.5,0.5);
  frame_nuclei->SetAxisColor(10,"XYZ");
  frame_nuclei->SetLabelColor(10,"XYZ");
  frame_nuclei->Draw("FBBB");
  nuclei->SetFillColor(10);
  nuclei->SetLineColor(kAzure);
  nuclei->Draw("FBBBsame");
  xaxis->Draw("FBBBsame");
  yaxis->Draw("FBBBsame");
  zaxis->Draw("FBBBsame");

//c1->Print("pdf/Y_0.pdf");
//c2->Print("pdf/Y_1.pdf");
//c3->Print("pdf/Y_2.pdf");

}
