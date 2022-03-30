//K. Okuyama (Mar. 30, 2022)
//Drawing B-field of Q-magnet

const double PI = 4.*atan(1.);
const double mu0 = 1.25663706212*1e-6;//N/A^2
const double Qmin = -0.5;
const double Qmax =  0.5;

double phi(double *x, double *par){//magnetic field scalar potential
double qm = 224000./0.12;//Wb --KQ1
double a  = par[0];//half-gap
double xx = x[0];
double yy = x[1];
double b1 = sqrt((a-xx)*(a-xx)+(a-yy)*(a-yy));
double b2 = sqrt((a-xx)*(a-xx)+(a+yy)*(a+yy));
double b3 = sqrt((a+xx)*(a+xx)+(a+yy)*(a+yy));
double b4 = sqrt((a+xx)*(a+xx)+(a-yy)*(a-yy));
double val = (qm*mu0/4./PI)*(1./b1-1./b2+1./b3-1./b4);
if(b1==0||b3==0){val=Qmax;}
if(b2==0||b4==0){val=Qmin;}
return val;
}

double F_Bxfield(double *x, double *par){//magnetic field function
double qm = 224000./0.12;//Wb --KQ1
double a  = par[0];//half-gap
double xx = x[0];
double yy = x[1];
double b1 = sqrt((a-xx)*(a-xx)+(a-yy)*(a-yy));
double b2 = sqrt((a-xx)*(a-xx)+(a+yy)*(a+yy));
double b3 = sqrt((a+xx)*(a+xx)+(a+yy)*(a+yy));
double b4 = sqrt((a+xx)*(a+xx)+(a-yy)*(a-yy));
double db1x = (a-xx)/b1/b1/b1;
double db2x = (a-xx)/b2/b2/b2;
double db3x = -(a+xx)/b3/b3/b3;
double db4x = -(a+xx)/b4/b4/b4;
double db1y = (a-yy)/b1/b1/b1;
double db2y = -(a+yy)/b2/b2/b2;
double db3y = -(a+yy)/b3/b3/b3;
double db4y = (a-yy)/b4/b4/b4;
double val = (qm*mu0/4./PI)*(1./b1-1./b2+1./b3-1./b4);
double Bx = (qm*mu0/4./PI)*(db1x-db2x+db3x-db4x);
double By = (qm*mu0/4./PI)*(db1y-db2y+db3y-db4y);
return Bx;
}
double F_Byfield(double *x, double *par){//magnetic field function
double qm = 224000./0.12;//Wb --KQ1
double a  = par[0];//half-gap
double xx = x[0];
double yy = x[1];
double b1 = sqrt((a-xx)*(a-xx)+(a-yy)*(a-yy));
double b2 = sqrt((a-xx)*(a-xx)+(a+yy)*(a+yy));
double b3 = sqrt((a+xx)*(a+xx)+(a+yy)*(a+yy));
double b4 = sqrt((a+xx)*(a+xx)+(a-yy)*(a-yy));
double db1x = (a-xx)/b1/b1/b1;
double db2x = (a-xx)/b2/b2/b2;
double db3x = -(a+xx)/b3/b3/b3;
double db4x = -(a+xx)/b4/b4/b4;
double db1y = (a-yy)/b1/b1/b1;
double db2y = -(a+yy)/b2/b2/b2;
double db3y = -(a+yy)/b3/b3/b3;
double db4y = (a-yy)/b4/b4/b4;
double val = (qm*mu0/4./PI)*(1./b1-1./b2+1./b3-1./b4);
double Bx = (qm*mu0/4./PI)*(db1x-db2x+db3x-db4x);
double By = (qm*mu0/4./PI)*(db1y-db2y+db3y-db4y);
return By;
}

double F_test(double *x, double *par){//magnetic field function
double a  = par[0];//half-gap
double xx = x[0];
return a*(1./pow(sqrt((a-xx)*(a-xx)+a*a),3)-1./pow(sqrt((a+xx)*(a+xx)+a*a),3));
}
double F_taylor1(double *x, double *par){//Taylor expansion
double a  = par[0];//half-gap
double xx = x[0];
return 3.*xx/2./sqrt(2.)/a/a/a;
}
double F_taylor3(double *x, double *par){//Taylor expansion
double a  = par[0];//half-gap
double xx = x[0];
return 3.*xx/2./sqrt(2.)/a/a/a+5.*xx*xx*xx/16./sqrt(2.)/a/a/a/a/a;
}
double F_taylor5(double *x, double *par){//Taylor expansion
double a  = par[0];//half-gap
double xx = x[0];
return 3.*xx/2./sqrt(2.)/a/a/a+5.*xx*xx*xx/16./sqrt(2.)/a/a/a/a/a-147.*a*xx*xx*xx*xx*xx/256./sqrt(2.)/a/a/a/a/a/a/a/a;
}
double F_taylor7(double *x, double *par){//Taylor expansion
double a  = par[0];//half-gap
double xx = x[0];
return 3.*xx/2./sqrt(2.)/a/a/a+5.*xx*xx*xx/16./sqrt(2.)/a/a/a/a/a-147.*a*xx*xx*xx*xx*xx/256./sqrt(2.)/a/a/a/a/a/a/a/a-243.*xx*xx*xx*xx*xx*xx*xx/2048./sqrt(2.)/a/a/a/a/a/a/a/a/a;
}

void Qmagnet_field(){
gROOT->Reset();
gROOT->SetStyle("Plain");

TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);

TF2* func = new TF2("Scalar potential",phi,-10,10,-10,10,1);
//func->SetNpx(20000);
//func->SetNpy(20000);
func->SetParameter(0,5.);
func->Draw("surf1");

TH2D* Bfield = new TH2D("Bfield","Bfield",20,-10,10,20,-10,10);
double x, y, k, val;
for(int i=0;i<1000000;i++){
x = gRandom->Uniform(-10.,10.);
y = gRandom->Uniform(-10.,10.);
k = gRandom->Uniform(Qmin,Qmax);
val = func->Eval(x,y);
if(k<val){Bfield->Fill(x,y);}
}
TCanvas* c2 = new TCanvas("c2", "c2", 800, 800);
Bfield->Draw("ARRcolz");

TCanvas* c3 = new TCanvas("c3", "c3", 800, 800);
TF2* f_Bx = new TF2("Bx(x,y)",F_Bxfield,-4,4,-4,4,1);
//func->SetNpx(20000);
//func->SetNpy(20000);
f_Bx->SetParameter(0,5.);
f_Bx->Draw("surf1");

TCanvas* c4 = new TCanvas("c4", "c4", 800, 800);
TF2* f_By = new TF2("By(x,y)",F_Byfield,-4,4,-4,4,1);
//func->SetNpx(20000);
//func->SetNpy(20000);
f_By->SetParameter(0,5.);
f_By->Draw("surf1");

TCanvas* c5 = new TCanvas("c5", "c5", 800, 800);
TF1* f_test = new TF1("By(x,y=0)",F_test,-8,8,1);
TF1* f_taylor1 = new TF1("f_taylor1",F_taylor1,-8,8,1);
TF1* f_taylor3 = new TF1("f_taylor3",F_taylor3,-8,8,1);
TF1* f_taylor5 = new TF1("f_taylor5",F_taylor5,-8,8,1);
TF1* f_taylor7 = new TF1("f_taylor7",F_taylor7,-8,8,1);
f_test->SetParameter(0,5.);
f_taylor1->SetParameter(0,5.);
f_taylor3->SetParameter(0,5.);
f_taylor5->SetParameter(0,5.);
f_taylor7->SetParameter(0,5.);
f_taylor1->SetLineColor(kRed);
f_taylor3->SetLineColor(kRed);
f_taylor5->SetLineColor(kRed);
f_taylor7->SetLineColor(kRed);
f_taylor1->SetLineStyle(4);
f_taylor3->SetLineStyle(2);
f_taylor5->SetLineStyle(7);
f_taylor7->SetLineStyle(9);
f_test->Draw("");
f_taylor1->Draw("same");
f_taylor3->Draw("same");
f_taylor5->Draw("same");
f_taylor7->Draw("same");

}
