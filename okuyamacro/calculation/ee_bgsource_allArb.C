/**********************************************/
/*  Background source of Scattered Electrons  */
/**********************************************/
//
// K. Okuyama (Sep. 7, 2020)
// K. Okuyama (Nov. 1, 2020)
// -revised with Katayama-kun's help
// -M-thesis (ee_bgsource.pdf)
// -Arbitrary Unit (scaled at theta_ee=0.001)
//
//VP Flux
//Moller scattering
//Bremsstrahlung

double PI=4.*atan(1.);

//E. Amaldi, S. Fubini, and G. Furlan, Pion Electroproduction. Electroproduction at Low- Energy and Hadron Form-Factors, Vol. 83 of Springer Tracts in Mod. Phys., Springer, Berlin, 1979.
double vpflux_lab(double *x, double *par){
	double Einc  = 4.319;//[GeV]
	//double Escat = 2.1;//[GeV]
	double Escat = par[0];//[GeV]
	double theta = x[0];	

	double Me=511.*pow(10,-6.);//[GeV/c^2]
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


double moller_lab_nnlHRS(double *x, double *par){
	//double Einc  = 2.344*1000000.;//[keV]
	double Einc  = 4.319*1000000.;//[keV]
	//double Escat = 2.1*1000000.;//[keV]
	double Escat = par[0];//[keV]
	double theta = x[0];	
	//double Me=511.*pow(10,-6.);//[keV/c^2]
	double Me=511.;//[GeV/c^2]
	double Mp=0.9382720*1000000.;//[keV/c^2]
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
//if(dif_lab<0.)dif_lab=0.;
if(fabs(Escat-2100000.)>100000.)dif_lab=0.;
	return dif_lab*400.*pow(10.,6.)*20624/(4.5*pow(10.,7.));//[Arb.]
	//return dif_cm*400.*pow(10.,6.);//[b/sr] (CM)

}
double moller_lab(double *x, double *par){
	//double Einc  = 2.344*1000000.;//[keV]
	double Einc  = 4.319*1000000.;//[keV]
	//double Escat = 2.1*1000000.;//[keV]
	double Escat = par[0];//[keV]
	double theta = x[0];	
	//double Me=511.*pow(10,-6.);//[keV/c^2]
	double Me=511.;//[GeV/c^2]
	double Mp=0.9382720*1000000.;//[keV/c^2]
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
//if(dif_lab<0.)dif_lab=0.;
	return dif_lab*400.*pow(10.,6.)*20624/(4.5*pow(10.,7.));//[Arb.]
	//return dif_cm*400.*pow(10.,6.);//[b/sr] (CM)

}

//Y. Tsai, "Pair production and bremsstrahlung of charged meson", Rev. Mod. Phys., Vol 46, No4 (1974)
double brems_lab(double *x, double *par){
	//double Einc  = 4.523*1000.;//[MeV]
	//double pscat = 3.030*1000.;//[MeV]
	//double Einc  = 2.344*1000.;//[MeV]
	//double pscat = 0.844*1000.;//[MeV]
	double Einc  = 4.319*1000.;//[MeV]
	double pscat = 2.1*1000.;//[MeV]
	double Z=par[0];
//	double Escat = par[0];//[MeV]
	double thetaee = x[0];	
	double Me=511.*pow(10,-3.);//[MeV/c^2]
	double Mp=0.9382720*1000.;//[MeV/c^2]
	double pinc = sqrt(Einc*Einc-Me*Me);
	double Escat = sqrt(pscat*pscat+Me*Me);
	double sine = pscat*sin(thetaee)/(Einc-Escat);
	double theta = asin(sine);
	//theta=thetaee;
	double omega = Einc - Escat;
	double alpha = 1./137.;


	double y=omega/Einc;
	double l=theta*theta*Einc*Einc/Me/Me;
	double b1=((2.*y-2.)/(l+1.)/(l+1.))+(12.*l*(1.-y)/(pow((l+1.),4.)));
	double b2=((2.-2.*y+y*y)/(l+1.)/(l+1.))-(4.*l*(1.-y)/(pow((l+1.),4.)));
	double G2_infty=Z*Z+Z;
	double a, ap;
	if(Z==1){a=122.8/Me;ap=282.4/Me;}
	else if(Z==2){a=90.8/(pow(Z,1./3.)*Me);ap=265.8/(pow(Z,2./3.)*Me);}
	else if(Z==3){a=100./(pow(Z,1./3.)*Me);ap=418.6/(pow(Z,2./3.)*Me);}
	else if(Z==4){a=106./(pow(Z,1./3.)*Me);ap=571.4/(pow(Z,2./3.)*Me);}
	else if(Z==5){a=111.7/(pow(Z,1./3.)*Me);ap=724.2/(pow(Z,2./3.)*Me);}
	else{a=184.15/(sqrt(2.718)*pow(Z,1./3.)*Me);ap=1194./(sqrt(2.718)*pow(Z,2./3.)*Me);}
	//double tmin=(omega*Me*Me*(l+1.)*(l+1.)/(2.*Einc*(Einc-omega)))*(omega*Me*Me*(l+1.)*(l+1.)/(2.*Einc*(Einc-omega)));
	double tmin=(omega*Me*Me*(l+1.)/(2.*Einc*(Einc-omega)));
	tmin=tmin*tmin;
	double X=Z*Z*(log(a*a*Me*Me*(l+1.)*(l+1.)/(a*a*tmin+1.))-1.)+Z*(log(ap*ap*Me*Me*(l+1.)*(l+1.)/(ap*ap*tmin+1.))-1.);
	double var=alpha*Z*alpha*Z;
	double fun=1.202*var-1.0369*var*var+1.008*var*var*var/(1.+var);
	double dif_lab=(2.*alpha*alpha*alpha*Einc*Einc/(PI*omega*Me*Me*Me*Me))*(b1*G2_infty+b2*(X-2.*Z*Z*fun));
	//return dif_lab*197.*197.*1000./100.;//[mb/MeV/sr]
	return dif_lab*197.*197.*1000./100.*20624/74846;//[Arb.]
}

void ee_bgsource_allArb(){
gROOT->Reset();
gROOT->SetStyle("Plain");

//TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
//c1->SetLeftMargin(0.15);
//c1->SetBottomMargin(0.15);
//c1->SetTopMargin(0.15);
//c1->SetRightMargin(0.15);
//c1->SetLogy(1);
//c1->SetGrid(1);
//c1->Draw();


//TH1F* h1 = c1 -> DrawFrame(0.,pow(10.,-8.),0.3,pow(10.,4.));
//h1->GetXaxis()->SetTitleSize(0.05);
//h1->GetXaxis()->SetTitleOffset(1.05);
//h1->GetXaxis()->SetTitle("Emitted Photon Angle [rad]");
//h1->GetXaxis()->SetLabelOffset(0.02);
//h1->GetYaxis()->SetTitleOffset(1.1);
//h1->GetYaxis()->SetTitleSize(0.05);
//h1->GetYaxis()->SetTitle("Differential C.S. [mb/MeV/sr]");

cout<<"f(theta=13.2 degree) is shown:"<<endl;
TF1* func1 = new TF1("func1",brems_lab, 0.0001, 0.35,1);
TF1* func11 = new TF1("func11",brems_lab, 0.001, 0.35,1);
TF1* func111 = new TF1("func111",brems_lab, 0.001, 0.35,1);
TF1* func1111 = new TF1("func1111",brems_lab, 0.001, 0.35,1);
func1->SetNpx(600);
func1->SetParameter(0,1.);
func1->SetLineColor(kCyan);
func1->SetLineWidth(3);
func1->SetLineStyle(1);
//func11->SetNpx(600000);
//func11->SetParameter(0,6.);
//func11->SetLineColor(kOrange);
//func11->SetLineWidth(3);
////func11->SetLineStyle(9);
//func11->SetLineStyle(1);
//func111->SetNpx(600000);
//func111->SetParameter(0,13.);
//func111->SetLineColor(kGray);
//func111->SetLineWidth(3);
////func111->SetLineStyle(7);
//func111->SetLineStyle(1);
//func1111->SetNpx(600000);
//func1111->SetParameter(0,82.);
//func1111->SetLineColor(kBlack);
//func1111->SetLineWidth(3);
////func1111->SetLineStyle(10);
//func1111->SetLineStyle(1);
//func1->Draw("same");
//func11->Draw("same");
//func111->Draw("same");
//func1111->Draw("same");
// TLegend *tl = new TLegend(0.65,0.50,0.80,0.80);
// tl->AddEntry(func1111,"Pb","l");;
// tl->AddEntry(func111,"Al","l");;
// tl->AddEntry(func11,"C","l");
// tl->AddEntry(func1,"H","l");
// tl->Draw("same");
//cout<<"val_lab"<<func1->Eval(13.2*PI/180)<<endl;
//cout<<"val_lab"<<func1->Eval(0.2)<<endl;

//TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
//TH1F* h2 = c2 -> DrawFrame(0.,pow(10.,-5.),0.3,pow(10.,3.));
//h2->GetXaxis()->SetTitleSize(0.05);
//h2->GetXaxis()->SetTitleOffset(0.95);
//h2->GetXaxis()->SetTitle("Scattered Electron Angle [rad]");
//h2->GetXaxis()->SetLabelOffset(0.02);
//h2->GetYaxis()->SetTitleOffset(1.00);
//h2->GetYaxis()->SetTitleSize(0.05);
//h2->GetYaxis()->SetTitle("Differential C.S. (CM) [b/sr]");
//c2->SetLogy(1);
//c2->SetGrid(1);
//h2->Draw("");
TF1* func2_nnl = new TF1("func2_nnl",moller_lab_nnlHRS, 0.0001, 0.35,1);
func2_nnl->SetNpx(600);
//func2_nnl->SetParameter(0,0.844*1000000);
func2_nnl->SetParameter(0,2.1*1000000);
//func2_nnl->SetParameter(0,13.);
//func2_nnl->SetParameter(0,100.);
func2_nnl->SetLineColor(kAzure);
func2_nnl->SetLineWidth(3);
TF1* func2 = new TF1("func2",moller_lab, 0.0001, 0.35,1);
func2->SetNpx(600);
//func2->SetParameter(0,0.844*1000000);
func2->SetParameter(0,2.1*1000000);
//func2->SetParameter(0,13.);
//func2->SetParameter(0,100.);
func2->SetLineColor(kAzure);
func2->SetLineWidth(3);
func2->SetLineStyle(2);
//func2->Draw("same");
//cout<<"val_lab"<<func2->Eval(13.2*PI/180)<<endl;
//cout<<"val_lab"<<func2->Eval(0.2)<<endl;

TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);
TH1F* h3 = c3 -> DrawFrame(0.,pow(10.,-3.),0.3,pow(10.,3.));
h3->GetXaxis()->SetTitleSize(0.05);
h3->GetXaxis()->SetTitleOffset(0.95);
h3->GetXaxis()->SetTitle("Scattered Electron Angle [rad]");
h3->GetXaxis()->SetLabelOffset(0.02);
h3->GetYaxis()->SetTitleOffset(0.90);
h3->GetYaxis()->SetTitleSize(0.05);
h3->GetYaxis()->SetTitle("VP Flux [#gamma/MeV/sr/electron]");
c3->SetLogy(1);
c3->SetGrid(1);
h3->Draw("");
TF1* func3 = new TF1("func3",vpflux_lab, 0.0001, 0.35,1);
func3->SetNpx(600);
func3->SetParameter(0,2.1);
func3->SetLineColor(kGreen);
func3->SetLineWidth(3);
func3->Draw("same");
cout<<"val_lab"<<func3->Eval(13.2*PI/180)<<endl;

TCanvas* c4 = new TCanvas("c4", "c4", 800, 600);
TH1F* h4 = c4 -> DrawFrame(0.,pow(10.,-9.),0.3,pow(10.,3.));
h4->GetXaxis()->SetTitleSize(0.05);
h4->GetXaxis()->SetTitleOffset(0.95);
h4->GetXaxis()->SetTitle("Scattered Electron Angle [rad]");
h4->GetXaxis()->SetLabelOffset(0.02);
h4->GetYaxis()->SetTitleOffset(1.00);
h4->GetYaxis()->SetTitleSize(0.05);
h4->GetYaxis()->SetTitle("Arbitrary Unit");
c4->SetLogy(1);
c4->SetGrid(1);
h4->Draw("");
func1->Draw("same");
func2_nnl->Draw("same");
func2->Draw("same");
func3->Draw("same");
cout<<"log(5)="<<log(5)<<endl;
cout<<"(8)^(1/3)="<<pow(8.,1./3.)<<endl;
//c1->Print("brems_calc.pdf");
//c2->Print("moller_calc.pdf");
//c3->Print("vpflux_calc.pdf");
c4->Update();
c4->Print("ee_bgsource.pdf");

cout<<"Values at 0.001"<<endl;
cout<<"Bremsstrahlung:"<<func1->Eval(0.0001)<<endl;
cout<<"Moller scat.:"<<func2->Eval(0.0001)<<endl;
cout<<"VP Flux:"<<func3->Eval(0.0001)<<endl;
}
