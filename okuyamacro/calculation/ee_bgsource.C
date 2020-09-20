/**********************************************/
/*  Background source of Scattered Electrons  */
/**********************************************/
//
// K. Okuyama (Sep. 7, 2020)
//
//VP Flux
//Moller scattering
//Bremsstrahlung

double PI = 4.*atan(1.);

//E. Amaldi, S. Fubini, and G. Furlan, Pion Electroproduction. Electroproduction at Low- Energy and Hadron Form-Factors, Vol. 83 of Springer Tracts in Mod. Phys., Springer, Berlin, 1979.
double vpflux_lab(double *x, double *par){
	double Einc  = 4.3;//[GeV]
	//double Escat = 2.1;//[GeV]
	double Escat = par[0];//[GeV]
	double theta = x[0];	

	double Me=pow(511,-6.);//[GeV/c^2]
	double Mp=0.9382720;//[GeV/c^2]
	double Qsq=2*Einc*Escat*(1-cos(theta));
	double omega = Einc - Escat;
	double q2=Qsq+omega*omega;
	double kg=omega-Qsq/(2*Mp);
	double eps=1/(1+2*(q2/Qsq)*tan(theta/2)*tan(theta/2));
//if(theta>0.230&&theta<0.231){
//cout<<"Qsq="<<Qsq<<endl;
//cout<<"q2="<<q2<<endl;
//cout<<"kg="<<kg<<endl;
//cout<<"eps="<<eps<<endl;
//}

	double vpflux=Escat*kg/(137*2*PI*PI*Einc*Qsq*(1-eps));
if(vpflux<0.)vpflux=0.;	
	return vpflux;
}


double moller_lab(double *x, double *par){
	double Einc  = 2.344*1000000.;//[GeV]
	//double Escat = 2.1*1000000.;//[GeV]
	double Escat = par[0];//[GeV]
	double theta = x[0];	
	//double Me=pow(511,-3.);//[GeV/c^2]
	double Me=511.;//[GeV/c^2]
	double Mp=0.9382720*1000000.;//[GeV/c^2]
	double pinc = sqrt(Einc*Einc-Me*Me);
	double pscat = sqrt(Escat*Escat-Me*Me);
//cout<<"pscat="<<pscat<<endl;
//Lorentz transformation
	double beta = pinc/(Einc+Me);//-->gamma diverges
	//double beta=0.99;
	double gamma=1./sqrt(1.-beta*beta);
	//double gamma=(Einc+Me)/sqrt((Einc+Me)*(Einc+Me)-pinc*pinc);
	if(isfinite(gamma));
	else{
	cout<<"gamma infinite"<<endl;
	beta = 0.9999988;//-->gamma diverges
	gamma= 649.;
	}
//cout<<"beta="<<beta<<endl;
//cout<<"gamma="<<gamma<<endl;

	TLorentzVector B_4vec;
	TLorentzVector T_4vec;
	TLorentzVector L_4vec;
	B_4vec.SetPxPyPzE(0,0.,pinc,Einc);
//cout<<"B.Mass="<<B_4vec.M()<<endl;
	T_4vec.SetPxPyPzE(0.,0.,0.,Me);
	L_4vec.SetPxPyPzE(pscat*sin(theta),0.,pscat*cos(theta),Escat);
//cout<<"L.Mass="<<L_4vec.M()<<endl;
//cout<<"L.Mom="<<L_4vec.Rho()<<endl;
//cout<<"L.px="<<L_4vec.Px()<<endl;
//cout<<"L.py="<<L_4vec.Py()<<endl;
//cout<<"L.pz="<<L_4vec.Pz()<<endl;
	TVector3 boost;
	TLorentzVector BT_4vec;
	BT_4vec=B_4vec+T_4vec;
	boost=BT_4vec.BoostVector();
	//boost.SetXYZ(0.,0.,beta);
//cout<<"boost.x="<<boost.X()<<endl;
//cout<<"boost.y="<<boost.Y()<<endl;
//cout<<"boost.z="<<boost.Z()<<endl;
	L_4vec.Boost(-boost);
//cout<<"after"<<endl;
//cout<<"L.Mass="<<L_4vec.M()<<endl;
//cout<<"L.Mom="<<L_4vec.Rho()<<endl;
//cout<<"L.px="<<L_4vec.Px()<<endl;
//cout<<"L.py="<<L_4vec.Py()<<endl;
//cout<<"L.pz="<<L_4vec.Pz()<<endl;
	B_4vec.Boost(-boost);
	//theta = L_4vec.Angle(B_4vec.Vect());// -->180 deg (Back-to-back in CM)
	theta=L_4vec.Theta();
//cout<<"theta="<<theta<<endl;
//if(fabs(x[0]-0.8)<0.01)cout<<"theta="<<theta*180./PI<<endl;
	Einc=B_4vec.E();
//cout<<"Einc="<<Einc<<endl;
	Escat=L_4vec.E();
//cout<<"Escat="<<Escat<<endl;
//cout<<"Escat="<<gamma*(sqrt(2.3*2.3+Me*Me)-beta*cos(0.23)*2.3)<<endl;
	pinc=B_4vec.Rho();
//cout<<"pinc="<<pinc<<endl;
	pscat=L_4vec.Rho();
//cout<<"pscat="<<pscat<<endl;

	double Qsq=2*Einc*Escat*(1-cos(theta));
	double omega = Einc - Escat;
	double q2=Qsq+omega*omega;
	double kg=omega-Qsq/(2*Mp);
	double eps=1/(1+2*(q2/Qsq)*tan(theta/2)*tan(theta/2));
	double alpha = 1./137.;

	double factor1 = 2*Einc*Einc-Me*Me;
	double factor2 = Einc*Einc-Me*Me;
	double dif_cm = ((alpha*alpha*factor1*factor1)/(4*Einc*Einc*factor2*factor2))*(4./pow(sin(theta),4.)-3./pow(sin(theta),2.)+(1.+4./(sin(theta)*sin(theta)))*(factor2*factor2)/(factor1*factor1));
	//double dif_cm = ((alpha*alpha)/(4*Einc*Einc*pow(pinc*sin(theta),4.)))*(4.*(Me*Me+2.*pinc*pinc)*(Me*Me+2.*pinc*pinc)+(4.*pow(pinc,4.)-3.*(Me*Me+2.*pinc*pinc)*(Me*Me+2.*pinc*pinc))*sin(theta)*sin(theta)+pow(pinc,4.)*pow(sin(theta),4.));

	double labtocm = (gamma*pscat*pscat*(pscat*cos(theta)+beta*Escat))/(pow(sqrt(pscat*pscat*sin(theta)*sin(theta)+gamma*gamma*(pscat*cos(theta)+beta*Escat)*(pscat*cos(theta)+beta*Escat)),3.));
//cout<<"labtocm="<<labtocm<<endl;
//labtocm=1.;
	double dif_lab =dif_cm/labtocm;
if(dif_lab<0.)dif_lab=0.;
	//return dif_lab*0.0004*pow(10.,12.);//[b/sr]
	return dif_lab*400.*pow(10.,6.);//[b/sr]

}

//Y. Tsai, "Pair production and bremsstrahlung of charged meson", Rev. Mod. Phys., Vol 46, No4 (1974)
double brems_lab(double *x, double *par){
	//double Einc  = 4.523*1000.;//[MeV]
	//double Escat = 3.030*1000.;//[MeV]
	double Einc  = 2.344*1000.;//[MeV]
	double Escat = 0.844*1000.;//[MeV]
	double Z=par[0];
//	double Escat = par[0];//[MeV]
	double thetaee = x[0];	
	double Me=pow(511,-3.);//[MeV/c^2]
	double Mp=0.9382720*1000.;//[MeV/c^2]
	double pinc = sqrt(Einc*Einc-Me*Me);
	double pscat = sqrt(Escat*Escat-Me*Me);
	double theta = asin(pscat*sin(thetaee)/(Einc-Escat));
//cout<<"theta_ee:theta_k="<<thetaee<<":"<<theta<<endl;
//cout<<"pscat="<<pscat<<endl;
//Lorentz transformation
	double beta = pinc/(Einc+Me);//-->gamma diverges
	double gamma=1./sqrt(1-beta*beta);
	double omega=Einc-Escat;

	double alpha = 1./137.;
	double y=omega/Einc;
//if(fabs(theta-0.2)<0.0001)cout<<"y="<<y<<endl;
	double l=theta*theta*Einc*Einc/Me/Me;
//if(fabs(theta-0.2)<0.0001)cout<<"l="<<l<<endl;
	double b1=((2.*y-2.)/(l+1.)/(l+1.))+(12.*l*(1.-y)/(pow((l+1.),4.)));
//if(fabs(theta-0.2)<0.0001)cout<<"b1="<<b1<<endl;
	double b2=((2.-2.*y+y*y)/(l+1.)/(l+1.))-(4.*l*(1.-y)/(pow((l+1.),4.)));
//b2=-b2;
//if(fabs(theta-0.2)<0.0001)cout<<"b2="<<b2<<endl;
	double G2_infty=Z*Z+Z;
//if(fabs(theta-0.2)<0.0001)cout<<"G2_infty="<<G2_infty<<endl;
	double a=184.15/(sqrt(2.718)*pow(Z,1./3.)*Me);
//if(fabs(theta-0.2)<0.0001)cout<<"a="<<a<<endl;
	double ap=1194./(sqrt(2.718)*pow(Z,2./3.)*Me);
//if(fabs(theta-0.2)<0.0001)cout<<"ap="<<ap<<endl;
	//double tmin=(omega*Me*Me*(l+1.)*(l+1.)/(2.*Einc*(Einc-omega)))*(omega*Me*Me*(l+1.)*(l+1.)/(2.*Einc*(Einc-omega)));
	double tmin=(omega*Me*Me*(l+1.)/(2.*Einc*(Einc-omega)));
	tmin=tmin*tmin;
//if(fabs(theta-0.2)<0.0001)cout<<"tmin="<<tmin<<endl;
	double X=Z*Z*(log(a*a*Me*Me*(l+1.)*(l+1.)/(a*a*tmin+1.))-1.)+Z*(log(ap*ap*Me*Me*(l+1.)*(l+1.)/(ap*ap*tmin+1.))-1.);
//if(fabs(theta-0.2)<0.0001)cout<<"X="<<X<<endl;
//X=-X;
	double var=alpha*Z;
	double fun=1.202*var-1.0369*var*var+1.008*var*var*var/(1.+var);
//if(fabs(theta-0.2)<0.0001)cout<<"fun="<<fun<<endl;
	double dif_lab=(2.*alpha*alpha*alpha*Einc*Einc/(PI*omega*Me*Me*Me*Me))*(b1*G2_infty+b2*(X-2.*Z*Z*fun*fun));
	//return dif_lab*400.*0.001*0.001;//[mb/MeV/sr]
	return dif_lab*400.*1000.;//[mb/MeV/sr]

}

void ee_bgsource(){
gROOT->Reset();
gROOT->SetStyle("Plain");

TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
c1->SetLeftMargin(0.15);
c1->SetBottomMargin(0.15);
c1->SetTopMargin(0.15);
c1->SetRightMargin(0.15);
c1->SetLogy(1);
c1->SetGrid(1);
c1->Draw();


//TH1F* h1 = c1 -> DrawFrame(0.15,0.,0.3,0.01);
//h1->GetXaxis()->SetTitleFont(12);
//h1->GetXaxis()->SetTitleSize(0.05);
//h1->GetXaxis()->SetTitleOffset(1.05);
//h1->GetXaxis()->SetTitle("theta [rad]");
//h1->GetXaxis()->SetLabelOffset(0.02);
////h1->GetYaxis()->SetTitleFont(22);
//h1->GetYaxis()->SetTitleOffset(1.1);
//h1->GetYaxis()->SetTitleSize(0.05);
//h1->GetYaxis()->SetTitle("vpflux [/GeV/sr]");

cout<<"f(theta=13.2 degree) is shown:"<<endl;
TF1* func1 = new TF1("func1",brems_lab, 0.001, 0.3,1);
TF1* func11 = new TF1("func11",brems_lab, 0.001, 0.3,1);
TF1* func111 = new TF1("func111",brems_lab, 0.001, 0.3,1);
func1->SetNpx(600000);
func1->SetParameter(0,3.);
func1->SetLineColor(kViolet);
func1->SetLineWidth(4);
func1->SetLineStyle(1);
func11->SetNpx(600000);
func11->SetParameter(0,6.);
func11->SetLineColor(kViolet);
func11->SetLineWidth(4);
func11->SetLineStyle(9);
func111->SetNpx(600000);
func111->SetParameter(0,24.);
func111->SetLineColor(kViolet);
func111->SetLineWidth(4);
func111->SetLineStyle(7);
func1->Draw("");
func11->Draw("same");
func111->Draw("same");
cout<<"val_lab"<<func1->Eval(13.2*PI/180)<<endl;
cout<<"val_lab"<<func1->Eval(0.2)<<endl;

TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
c2->SetLogy(1);
c2->SetGrid(1);
TF1* func2 = new TF1("func2",moller_lab, 0.001, 0.3,1);
func2->SetNpx(600000);
func2->SetParameter(0,0.844*1000000);
//func2->SetParameter(0,13.);
//func2->SetParameter(0,100.);
func2->SetLineColor(kAzure);
func2->SetLineWidth(4);
func2->Draw("");
cout<<"val_lab"<<func2->Eval(13.2*PI/180)<<endl;
cout<<"val_lab"<<func2->Eval(0.2)<<endl;

TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);
c3->SetLogy(1);
c3->SetGrid(1);
TF1* func3 = new TF1("func3",vpflux_lab, 0.001, 0.3,1);
func3->SetNpx(600000);
func3->SetParameter(0,2.1);
func3->SetLineColor(kGreen);
func3->SetLineWidth(4);
func3->Draw("");
cout<<"val_lab"<<func3->Eval(13.2*PI/180)<<endl;

TCanvas* c4 = new TCanvas("c4", "c4", 800, 600);
c4->SetLogy(1);
c4->SetGrid(1);
func1->Draw("");
func2->Draw("same");
func3->Draw("same");
cout<<"log(5)"<<log(5)<<endl;
//c1->Print("montecarlo_data0826_1kai.png");

}
