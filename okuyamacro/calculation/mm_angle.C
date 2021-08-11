#include <TMath.h>


double F_theta(double *x, double *par){
	double xp = x[0];
	double yp = x[1];
	double theta_0 = par[0];
		double s0 = sin(theta_0);
		double c0 = cos(theta_0);
		double R = sqrt(1.+xp*xp+yp*yp);
		double ccc = (-1.*yp*s0+c0)/R;//cos(theta)
		double sss = sqrt(1.-ccc*ccc);//sin(theta)
return acos(ccc);
}
double F_delxp(double *x, double *par){
	double xp = x[0];
	double yp = x[1];
	double theta_0 = par[0];
		double s0 = sin(theta_0);
		double c0 = cos(theta_0);
		double R = sqrt(1.+xp*xp+yp*yp);
		double ccc = (-1.*yp*s0+c0)/R;//cos(theta)
		double sss = sqrt(1.-ccc*ccc);//sin(theta)
		double delxp = (-1.*yp*s0+c0)*xp/(R*R*R*sss);//d(theta)/dX'
		double delyp = (-1.*yp*s0+c0)*yp/(R*R*R*sss)+s0/(R*sss);//d(theta)/dY'
return delxp;
}
double F_delyp(double *x, double *par){
	double xp = x[0];
	double yp = x[1];
	double theta_0 = par[0];
		double s0 = sin(theta_0);
		double c0 = cos(theta_0);
		double R = sqrt(1.+xp*xp+yp*yp);
		double ccc = (-1.*yp*s0+c0)/R;//cos(theta)
		double sss = sqrt(1.-ccc*ccc);//sin(theta)
		double delxp = (-1.*yp*s0+c0)*xp/(R*R*R*sss);//d(theta)/dX'
		double delyp = (-1.*yp*s0+c0)*yp/(R*R*R*sss)+s0/(R*sss);//d(theta)/dY'
return delyp;
}

void mm_angle(){

	TH2D *h2_xpyp = new TH2D("h2_xpyp","h2_xpyp",100,-0.2,0.2,100,-0.2,0.2);
	TH1D *h_ratio = new TH1D("h_ratio","h_ratio",1000,-10.,10.);
	TH2D *h2_xpyp_ypdominance = new TH2D("h2_xpyp_ypdominance","h2_xpyp_ypdominance",100,-0.2,0.2,100,-0.2,0.2);
	TH2D *h2_delxp = new TH2D("h2_delxp","h2_delxp",80,-0.2,0.2,80,-0.2,0.2);
	TH2D *h2_delyp = new TH2D("h2_delyp","h2_delyp",80,-0.2,0.2,80,-0.2,0.2);

	gStyle->SetOptStat(0);

	gRandom->SetSeed(0);

cout<<"PI="<<PI<<endl;

	double xs_theta=0.;

	for(int i=0;i<40;i++){
	for(int j=0;j<40;j++){
		//double xp = 0.1*(2.*gRandom->Uniform()-1.);//X'
		//double yp = 0.1*(2.*gRandom->Uniform()-1.);//Y'
		double xp = -0.1 + (double)i*0.005;
		double yp = -0.1 + (double)j*0.005;
		//double xp = 0.;
		//double yp = 13.2*PI/180.;
		double theta_0 = 13.2*PI/180.;// \Theta_0
		double s0 = sin(theta_0);
		double c0 = cos(theta_0);
		double R = sqrt(1.+xp*xp+yp*yp);
		double ccc = (-1.*yp*s0+c0)/R;//cos(theta)
		double sss = sqrt(1.-ccc*ccc);//sin(theta)
		
		double delxp = (-1.*yp*s0+c0)*xp/(R*R*R*sss);//d(theta)/dX'
		double delyp = (-1.*yp*s0+c0)*yp/(R*R*R*sss)+s0/(R*sss);//d(theta)/dY'
	
		h2_xpyp->Fill(xp,yp);
		h2_delxp->SetBinContent(i+21,j+21,abs(delxp));
		h2_delyp->SetBinContent(i+21,j+21,abs(delyp));
	
		double ratio_yx = yp/xp;
		ratio_yx += R*R*sin(theta_0)/(-yp*sin(theta_0)+cos(theta_0))/xp;
	
		h_ratio->Fill(ratio_yx);
		if(ratio_yx>1.)h2_xpyp_ypdominance->Fill(xp,yp);//Y' is dominant for "theta (polar angle)" determination
		xs_theta += delxp*delxp*0.0025*0.0025+delyp*delyp*0.0012*0.0012;
//cout<<"delxp="<<delxp<<endl;
//cout<<"delyp="<<delyp<<endl;
//cout<<"R="<<R<<endl;
	}
	}

	TF2  *f_theta = new TF2("f_theta",F_theta,-0.4,0.4,-0.4,0.4,1);	
	f_theta->SetNpx(2000);
	f_theta->GetXaxis()->SetNdivisions(000);
	f_theta->GetYaxis()->SetNdivisions(000);
	f_theta->SetParameter(0,13.2*PI/180.);
	f_theta->SetLineColor(kAzure);
	f_theta->SetLineWidth(3);
	TF2  *f_delxp = new TF2("f_delxp",F_delxp,-0.2,0.2,-0.2,0.2,1);	
	f_delxp->SetNpx(2000);
	f_delxp->SetParameter(0,13.2*PI/180.);
	f_delxp->SetLineColor(kAzure);
	f_delxp->SetLineWidth(3);
	TF2  *f_delyp = new TF2("f_delyp",F_delyp,-0.2,0.2,-0.2,0.2,1);	
	f_delyp->SetNpx(2000);
	f_delyp->SetParameter(0,13.2*PI/180.);
	f_delyp->SetLineColor(kAzure);
	f_delyp->SetLineWidth(3);
	
	TCanvas *c0 = new TCanvas("c0","c0",800,800);
	TF1 *fx = new TF1("fx","-x",-0.4,0.4);
	TF1 *fy = new TF1("fy","x",-0.4,0.4);
	TGaxis *ax = new TGaxis(-0.4,-0.4,0.4,-0.4,"fx");
	TGaxis *ay = new TGaxis(-0.4,-0.4,-0.4,0.4,"fy");
	ax->SetLineWidth(1);
	ay->SetLineWidth(1);
	f_theta->Draw("colz");
	f_theta->GetXaxis()->SetLabelOffset(100);
	f_theta->GetYaxis()->SetLabelOffset(100);
	f_theta->GetXaxis()->SetTickLength(0);
	f_theta->GetYaxis()->SetTickLength(0);
	ax->Draw("same");
	ay->Draw("same");

	TCanvas *c0x = new TCanvas("c0x","c0x",800,800);
	f_delxp->Draw("colz");
	TCanvas *c0y = new TCanvas("c0y","c0y",800,800);
	f_delyp->Draw("colz");
	
	TCanvas *c1 = new TCanvas("c1","c1",800,800);
	h2_xpyp->Draw("colz");

	TCanvas *c2 = new TCanvas("c2","c2",800,800);
	h_ratio->Draw("");

	TCanvas *c3 = new TCanvas("c3","c3",800,800);
	h2_xpyp_ypdominance->Draw("colz");

	TCanvas *c4 = new TCanvas("c4","c4",800,800);
	h2_delxp->Draw("colz");

	TCanvas *c5 = new TCanvas("c5","c5",800,800);
	h2_delyp->Draw("colz");

cout<<"xs_theta [mrad]="<<sqrt(xs_theta/40./40)/1000.<<endl;
cout<<"Well done!"<<endl;
}
	
