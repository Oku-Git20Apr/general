/////Mass[GeV]/////
const double mp = Mp*0.001;//proton
const double mn = Mn*0.001;//neutron
const double mpi =Mpi*0.001;//pion
const double mL = ML*0.001;//Lambda
const double mK = MK*0.001;//K+ 
double pik_mom(double *x, double *par){

double momK, momL;
double ppi = x[0];
double theta_lab = par[0];

double beta = ppi/(sqrt(ppi*ppi+mpi*mpi)+mn);
double gamma = 1./sqrt(1.-beta*beta);

double s = sqrt(abs( pow(sqrt(ppi*ppi+mpi*mpi)+mn,2.)-ppi*ppi ));//s^2=E^2-P^2
if (pow(s*s - mL*mL - mK*mK,2.) -4.*mL*mL*mK*mK < 0) return -1.;
double pKcm =sqrt( pow(s*s - mL*mL - mK*mK,2.) -4.*mL*mL*mK*mK  )/(2.*s); //kaon momentum (CM)


double C = gamma*gamma-gamma*gamma*beta*beta*pow(cos(theta_lab * deg_to_rad),2.);
double D = 2.*gamma*beta*cos(theta_lab*deg_to_rad)*sqrt(abs(pKcm*pKcm+mK*mK));
double F = gamma*gamma*mK*mK-pKcm*pKcm-mK*mK;

momK = (D + sqrt(abs(D*D-4.*C*F)))/(2.*C);//kaon momentum (Lab)
if (ppi*ppi+momK*momK-2.*ppi*momK*cos(theta_lab*deg_to_rad) < 0 ) return -1.;
momL = sqrt(ppi*ppi+momK*momK-2.*ppi*momK*cos(theta_lab*deg_to_rad));//Lambda momentum (Lab)

return momL;
}

double eek_mom(double *x, double *par){

double momK, momL;
double Eg = x[0];
double theta_lab = par[0];

double beta = Eg/(Eg+mp);
double gamma = 1./sqrt(1.-beta*beta);

double s = sqrt(abs( pow(Eg+mp,2.)-Eg*Eg) );//s^2=E^2-P^2
double pKcm =sqrt( abs(pow(s*s - mL*mL - mK*mK,2.) -4.*mL*mL*mK*mK ) )/(2.*s); //pion momentum (CM)


double C = gamma*gamma-gamma*gamma*beta*beta*pow(cos(theta_lab * deg_to_rad),2.);
double D = 2.*gamma*beta*cos(theta_lab*deg_to_rad)*sqrt(pKcm*pKcm+mK*mK);
double F = gamma*gamma*mK*mK-pKcm*pKcm-mK*mK;

momK = (D + sqrt(abs(D*D-4.*C*F)))/(2.*C);//kaon momentum (Lab)
momL = sqrt(abs(Eg*Eg+momK*momK-2.*Eg*momK*cos(theta_lab*deg_to_rad)));//Lambda momentum (Lab)

return momL;
}

double kpi_mom(double *x, double *par){

double mompi, momL;
double pk = x[0];
double theta_lab = par[0];

double beta = pk/(sqrt(pk*pk+mK*mK)+mn);
double gamma = 1./sqrt(1.-beta*beta);

double s = sqrt( pow(sqrt(pk*pk+mK*mK)+mn,2.)-pk*pk );//s^2=E^2-P^2
double ppicm =sqrt( pow(s*s - mL*mL - mpi*mpi,2.) -4.*mL*mL*mpi*mpi  )/(2.*s); //pion momentum (CM)


double C = gamma*gamma-gamma*gamma*beta*beta*pow(cos(theta_lab * deg_to_rad),2.);
double D = 2.*gamma*beta*cos(theta_lab*deg_to_rad)*sqrt(ppicm*ppicm+mpi*mpi);
double F = gamma*gamma*mpi*mpi-ppicm*ppicm-mpi*mpi;

mompi = (D + sqrt(D*D-4.*C*F))/(2.*C);//pion momentum (Lab)
momL = sqrt(pk*pk+mompi*mompi-2.*pk*mompi*cos(theta_lab*deg_to_rad));//Lambda momentum (Lab)

return momL;
}

/////////////////////////////////////////////
////////   Reaction Threshold  //////////////
/////////////////////////////////////////////
double pik_lim_func(double theta_lab) {

double ppi = 0.01;
double momK, momL;

while(1){
double beta = ppi/(sqrt(ppi*ppi+mpi*mpi)+mn);
double gamma = 1./sqrt(1.-beta*beta);
double s = sqrt(abs( pow(sqrt(ppi*ppi+mpi*mpi)+mn,2.)-ppi*ppi ));//s^2=E^2-P^2
//if (pow(s*s - mL*mL - mK*mK,2.) -4.*mL*mL*mK*mK < 0 )ppi += 0.01;

double pKcm =sqrt( pow(s*s - mL*mL - mK*mK,2.) -4.*mL*mL*mK*mK  )/(2.*s); //kaon momentum (CM)


double C = gamma*gamma-gamma*gamma*beta*beta*pow(cos(theta_lab * deg_to_rad),2.);
double D = 2.*gamma*beta*cos(theta_lab*deg_to_rad)*sqrt(abs(pKcm*pKcm+mK*mK));
double F = gamma*gamma*mK*mK-pKcm*pKcm-mK*mK;

//momK = (D + sqrt(abs(D*D-4.*C*F)))/(2.*C);//kaon momentum (Lab)
if (D*D-4.*C*F >= 0. ) break;
//if (ppi*ppi+momK*momK-2.*ppi*momK*cos(theta_lab*deg_to_rad) > 0.01 ) break;
	ppi += 0.01;
}
	return ppi;
}

/////////////////////////////////////////////////
double eek_lim_func(double theta_lab) {

double momK, momL;
double Eg = 0.01;

while(1){
double beta = Eg/(Eg+mp);
double gamma = 1./sqrt(1.-beta*beta);

double s = sqrt( pow(Eg+mp,2.)-Eg*Eg );//s^2=E^2-P^2
double pKcm =sqrt( pow(s*s - mL*mL - mK*mK,2.) -4.*mL*mL*mK*mK  )/(2.*s); //kaon momentum (CM)

double C = gamma*gamma-gamma*gamma*beta*beta*pow(cos(theta_lab * deg_to_rad),2.);
double D = 2.*gamma*beta*cos(theta_lab*deg_to_rad)*sqrt(pKcm*pKcm+mK*mK);
double F = gamma*gamma*mK*mK-pKcm*pKcm-mK*mK;
if (D*D-4.*C*F >= 0. ) break;
	Eg += 0.01;
}
	return Eg;
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//                   main                                  //
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
void recoil_2body(){



/////////////////////////////////////////////kpi reaction
	TF1 *f0 = new TF1("f0",kpi_mom, 0.01, 9.0, 1);
	f0->SetParameters(0,0);
	f0->SetLineColor(kAzure);
	f0->SetLineStyle(1);
	f0->SetLineWidth(3);

	TF1 *f10 = new TF1("f10",kpi_mom, 0.01, 8.0, 1);
	f10->SetParameters(10,0);
	f10->SetLineColor(kAzure);
	f10->SetLineStyle(9);
	f10->SetLineWidth(3);

	TF1 *f20 = new TF1("f20",kpi_mom, 0.01, 7.0, 1);
	f20->SetParameters(20,0);
	f20->SetLineColor(kAzure);
	f20->SetLineStyle(7);
	f20->SetLineWidth(3);

	TF1 *f30 = new TF1("f30",kpi_mom, 0.01, 6.0, 1);
	f30->SetParameters(30,0);
	f30->SetLineColor(kAzure);
	f30->SetLineStyle(2);
	f30->SetLineWidth(3);

/////////////////////////////////////////////pik reaction
	double pik_lim0 = pik_lim_func(0);
	double pik_lim10 = pik_lim_func(10);
	double pik_lim20 = pik_lim_func(20);
	double pik_lim30 = pik_lim_func(30);
	cout << "limit0 is " << pik_lim0 << endl;
	cout << "limit10 is " << pik_lim10 << endl;
	cout << "limit20 is " << pik_lim20 << endl;
	cout << "limit30 is " << pik_lim30 << endl;
	
	TF1 *ff0 = new TF1("ff0",pik_mom, 0.01, 9, 1);
	ff0->SetRange(pik_lim0,6.);
	ff0->SetParameters(0,0);
	ff0->SetLineColor(kGreen);
	ff0->SetLineStyle(1);
	ff0->SetLineWidth(3);

	TF1 *ff10 = new TF1("ff10",pik_mom, 0.01, 8.0, 1);
	ff10->SetRange(pik_lim10,6.);
	ff10->SetParameters(10,0);
	ff10->SetLineColor(kGreen);
	ff10->SetLineStyle(9);
	ff10->SetLineWidth(3);

	TF1 *ff20 = new TF1("ff20",pik_mom, 0.01, 7.0, 1);
	ff20->SetRange(pik_lim20,6.);
	ff20->SetParameters(20,0);
	ff20->SetLineColor(kGreen);
	ff20->SetLineStyle(7);
	ff20->SetLineWidth(3);

	TF1 *ff30 = new TF1("ff30",pik_mom, 0.01, 6.0, 1);
	ff30->SetRange(pik_lim30,6.);
	ff30->SetParameters(30,0);
	ff30->SetLineColor(kGreen);
	ff30->SetLineStyle(2);
	ff30->SetLineWidth(3);
	
/////////////////////////////////////////////eek reaction
	double eek_lim0  = eek_lim_func(0);
	double eek_lim10 = eek_lim_func(10);
	double eek_lim20 = eek_lim_func(20);
	double eek_lim30 = eek_lim_func(30);
	cout << "limit0 is "  << eek_lim0 << endl;
	cout << "limit10 is " << eek_lim10 << endl;
	cout << "limit20 is " << eek_lim20 << endl;
	cout << "limit30 is " << eek_lim30 << endl;
	
	TF1 *fff0 = new TF1("fff0",eek_mom, 0.01, 9, 1);
	fff0->SetRange(eek_lim0,6.);
	fff0->SetParameters(0,0);
	fff0->SetLineColor(kRed);
	fff0->SetLineStyle(1);
	fff0->SetLineWidth(3);

	TF1 *fff10 = new TF1("fff10",eek_mom, 0.01, 8.0, 1);
	fff10->SetRange(eek_lim10,6.);
	fff10->SetParameters(10,0);
	fff10->SetLineColor(kRed);
	fff10->SetLineStyle(9);
	fff10->SetLineWidth(3);

	TF1 *fff20 = new TF1("fff20",eek_mom, 0.01, 7.0, 1);
	fff20->SetRange(eek_lim20,6.);
	fff20->SetParameters(20,0);
	fff20->SetLineColor(kRed);
	fff20->SetLineStyle(7);
	fff20->SetLineWidth(3);

	TF1 *fff30 = new TF1("fff30",eek_mom, 0.01, 6.0, 1);
	fff30->SetRange(eek_lim30,6.);
	fff30->SetParameters(30,0);
	fff30->SetLineColor(kRed);
	fff30->SetLineStyle(2);
	fff30->SetLineWidth(3);


///////////////////////////////////////////////////drawing hist
	TH2D *h = new TH2D("h","h",10,0.,2.,10,0,0.71);
    SetTH2(h,"Momentum tansfer","beam momentum (GeV/c)","recoil momentum (GeV/c)",0.4);
	h->SetStats(kFALSE);
	//h1->GetXaxis()->SetTitleFont(12);
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetTitleOffset(1.02);
	h->GetXaxis()->SetTitle("Incident Beam Momentum [GeV/#it{c}]");
	h->GetXaxis()->SetLabelOffset(0.02);
	//h1->GetYaxis()->SetTitleFont(22);
	h->GetYaxis()->SetTitleOffset(0.8);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitle("Produced #Lambda Momentum [GeV/#it{c}]");

	TLegend *tl = new TLegend(0.20,0.65,0.40,0.85);
    tl->SetTextFont(42);
	tl->AddEntry(fff0,"(#gamma,K^{+})","l");;
    tl->AddEntry(ff0,"(#pi^{+},K^{+})","l");
    tl->AddEntry(f0,"(K^{-},#pi^{-})","l");



//TCanvas *c = new TCanvas("c","canvas",800,600);
//	c->cd();
	h   -> Draw();
	f0  -> Draw("same");
	f10 -> Draw("same");
	f20 -> Draw("same");
	f30 -> Draw("same");
	ff0  -> Draw("same");
	ff10 -> Draw("same");
	ff20 -> Draw("same");
	ff30 -> Draw("same");
	fff0  -> Draw("same");
	fff10 -> Draw("same");
	fff20 -> Draw("same");
	fff30 -> Draw("same");
    tl->Draw("same");

}
