//------------------------//
//--  Differential C.S. --//
//------------------------//
//
//K. Okuyama (Sep. 18, 2020)
//
//Exp.
//CLAS2010
//CLAS2006
//
//Theory
//Kaon-Maid
//SLA
//RPR


//arXiv.1601.03840v2 [nucl-th] 22 Frb 2016
//Photoproduction of KL on the proton
//D. Skoupil and P. Bydzovsky
double RPR( double *x, double *par){

 const double PI=3.14159265359;
	double theta = PI*x[0]/180.;
	//double theta = acos(x[0]);
	double W = par[0];
	double Qsq = par[1];
	double phiK = 0.;
	
	
//---Physics Constant---//
 
 const double Mpi = 0.13957018;         // charged pion mass (GeV/c2)
 const double MP = 0.938272046*1000.;         // proton       mass (GeV/c2)
 const double MK = 0.493677*1000.;            // charged Kaon mass (GeV/c2)
 const double Me = 0.510998928e-3;      // electron     mass (GeV/c2)
 const double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
 const double ML = 1.115683*1000.;            // Lambda       mass (GeV/c2)
 const double MS0 = 1.192642*1000.;           // Sigma Zero   mass (GeV/c2)

	double alpha = 1./137.;
	double pK = (sqrt((W*W-ML*ML-MK*MK)*(W*W-ML*ML-MK*MK)-4.*ML*ML*MK*MK))/(2.*W);//~750MeV/c
//	double pg = sqrt((W*W+Qsq-MP*MP)*(W*W+Qsq-MP*MP)+4.*Qsq*MP*MP)/2./W;//~750MeV/c
	double pg = sqrt((W-MP)*(W-MP)+Qsq);//~1470MeV/c
cout<<"pg="<<pg<<endl;
	
	double pL = pK;
	double EL = sqrt(pL*pL+ML*ML);
	double pP = pg;
	double EP = sqrt(pP*pP+MP*MP);
	double Eg = sqrt(pg*pg-Qsq);
cout<<"Eg="<<Eg<<endl;
	double EK = sqrt(pK*pK+MK*MK);
	double k_k = -Qsq;
	double k_pL = Eg*EL-pg*pL;
	double k_p = Eg*EP-pg*pP;

	double C = (197.*197.*alpha/4./PI)*ML*pK/pg/W;
	double N = sqrt((EL+ML)*(EP+MP)/4./ML/MP);
	
//Mandelstam variables
	double s = W*W;//~5*10^6
	double t = (Eg-EK)*(Eg-EK)-pg*pg-pK*pK+2.*pg*pK*cos(theta);
	//double t = (Eg-EL)*(Eg-EL)-pg*pg-pL*pL+2.*pg*pL*cos(PI-theta);
	double u = (EP-EK)*(EP-EK)-pP*pP-pK*pK+2.*pP*pK*cos(PI-theta);
cout<<"s+t+u="<<s+t+u<<endl;
	double pgpg = -1.*Qsq;
	double pgpL = 0.5*(ML*ML+pgpg-u);
	double pgpp = 0.5*(s-pgpg-MP*MP);
	double pgpk = Eg*EK-pg*pK*cos(theta);//???
cout<<"s="<<s<<endl;
cout<<"t="<<t<<endl;
cout<<"u="<<u<<endl;
	k_k=pgpg;
	k_pL=pgpL;
	k_p=pgpp;
	double psi=0.23;
	double eps=1/(1+2*(pg*pg/Qsq)*tan(psi/2)*tan(psi/2));
	//double eps=1.-2.*(pg*pg/pgpg)*tan(theta/2.)*tan(theta/2.);
	//eps=0.7;
cout<<"eps="<<eps<<endl;
	double epsL=Qsq*eps/Eg/Eg;
cout<<"epsL="<<epsL<<endl;

//Form factors (SLA)
	double kappa_iv = 3.706;
	double kappa_is = -0.12;
	double g_rho=0.631;
	double kappa_rho=3.3;
	double g_omega=0.658;
	double kappa_omega=0.4;
	double Lambda_1_rho=0.863*1000.;//MeV
	double Lambda_1_D=1.21*1000.;//MeV
	double Lambda_2=2.1*1000.;//MeV
	double Lambda_QCD=0.33*1000.;//MeV
	double Mrho=775.8;//MeV
	double Momega=782.65;//MeV
	double kappa_p = 1.79;//anomalous proton magnetic moment??
	double kappa_n = -2.04;//anomalous proton magnetic moment??
	double kappa_L = -0.73;//from SLA
	double kappa_S = 1.02;//from SLA
	double Q2_tilde=Qsq*log((Lambda_2*Lambda_2+Qsq)/(Lambda_QCD*Lambda_QCD))/log(Lambda_2*Lambda_2/Lambda_QCD/Lambda_QCD);
	double F1_rho=Lambda_1_rho*Lambda_1_rho*Lambda_2*Lambda_2/((Lambda_1_rho*Lambda_1_rho+Q2_tilde)*(Lambda_2*Lambda_2+Q2_tilde));
	double F1_omega=F1_rho;
	double F1_D=Lambda_1_D*Lambda_1_D*Lambda_2*Lambda_2/((Lambda_1_D*Lambda_1_D+Q2_tilde)*(Lambda_2*Lambda_2+Q2_tilde));
	double F2_rho=(Lambda_1_rho*Lambda_1_rho/(Lambda_1_rho*Lambda_1_rho+Q2_tilde))*(Lambda_1_rho*Lambda_1_rho/(Lambda_1_rho*Lambda_1_rho+Q2_tilde))*(Lambda_2*Lambda_2/(Lambda_2*Lambda_2+Q2_tilde));
	double F2_omega=F2_rho;
	double F2_D=(Lambda_1_D*Lambda_1_D/(Lambda_1_D*Lambda_1_D+Q2_tilde))*(Lambda_1_D*Lambda_1_D/(Lambda_1_D*Lambda_1_D+Q2_tilde))*(Lambda_2*Lambda_2/(Lambda_2*Lambda_2+Q2_tilde));
	double F1_iv = g_rho*(Mrho*Mrho/(Mrho*Mrho+Qsq))*F1_rho+(1.-g_rho)*F1_D;
	double F2_iv = (kappa_rho*g_rho*(Mrho*Mrho/(Mrho*Mrho+Qsq))*F2_rho+(kappa_iv-kappa_rho*g_rho)*F2_D)/kappa_iv;
	double F1_is = g_omega*(Momega*Momega/(Momega*Momega+Qsq))*F1_omega+(1.-g_omega)*F1_D;
	double F2_is = (kappa_omega*g_omega*(Momega*Momega/(Momega*Momega+Qsq))*F2_omega+(kappa_is-kappa_omega*g_omega)*F2_D)/kappa_is;
	double F1_p = 0.5*(F1_is+F1_iv);
	double F2_p = (kappa_is*F2_is+kappa_iv*F2_iv)/2./kappa_p;
	double F1_L = 0.5*(F1_is-F1_iv);
	double F2_L = (kappa_is*F2_is-kappa_iv*F2_iv)/2./kappa_n;
	double F1_S = (F1_is-F1_iv)/2.;
	double F2_S = (kappa_is*F2_is-kappa_iv*F2_iv)/2./kappa_n;
	double a = 0.398;
	double Lambda1_k=0.642*1000.;
	double Lambda2_k=1.386*1000.;
	double F_k  = a/(1.+Qsq/Lambda1_k/Lambda1_k)+(1.-a)/((1.+Qsq/Lambda2_k/Lambda2_k)*(1.+Qsq/Lambda2_k/Lambda2_k));
	double Lambda_ks=0.95*1000.;
	double Lambda_k1=0.55*1000.;
	double F_ks=1./(1.+Qsq/Lambda_ks/Lambda_ks);
	double F_k1=1./(1.+Qsq/Lambda_k1/Lambda_k1);
cout<<"F1_is="<<F1_iv<<endl;
cout<<"F1_is="<<F1_is<<endl;
cout<<"F2_iv="<<F2_iv<<endl;
cout<<"F2_iv="<<F2_iv<<endl;
cout<<"F1_p="<<F1_p<<endl;
cout<<"F1_L="<<F1_L<<endl;
cout<<"F_k="<<F_k<<endl;
cout<<"F_ks="<<F_ks<<endl;
cout<<"F_k1="<<F_k1<<endl;
	
	double F1 = 1.;
	F1=F1_p;
	double F2 = kappa_p;
	F2=F2_p;
	//double g_KLP = -3.00;//BS1
	double g_KLP = -3.16;//SLA
	g_KLP = g_KLP*sqrt(4.*PI*alpha);
	double g_KSP = 0.91;
	g_KSP = g_KSP*sqrt(4.*PI*alpha);
	double kappa_SL=1.61;

	//SLA
	double A1_born = (g_KLP/(s-MP*MP))*(F1_p+kappa_p*F2_p)+(g_KLP/(u-ML*ML))*(F1_L+kappa_L*F2_L)+(g_KSP/(u-MS0*MS0))*(((MS0+ML)*kappa_SL)/2./MP)*F2_S;
	double A2_born = (g_KLP/(s-MP*MP)/(t-MK*MK))*((F_k+F1_p)+(F_k-F1_p)*((pgpk+pgpp)/(pgpL)));
	double A3_born = (g_KLP/(s-MP*MP))*(kappa_p/MP)*F2_p;
	double A4_born = (g_KLP/(u-ML*ML))*(kappa_L/ML)*F2_L+(g_KSP/(u-MS0*MS0))*(kappa_SL/MP)*F2_S;
	double A5_born = (g_KLP/(s-MP*MP))*(kappa_p/2./MP)*F2_p-(g_KLP/(u-ML*ML))*(kappa_L/2./ML)*F2_L-(g_KSP/(u-MS0*MS0))*(kappa_SL/2./MP)*F2_S;
	double A6_born = ((-2.*g_KLP)/(s-MP*MP)/(t-MK*MK))*(F_k+(F_k-F1_p)*(pgpg-2.*pgpk)*(pgpp/pgpg/pgpL))+((g_KLP)/(u-ML*ML))*2.*F1_L/pgpg;
	//A1_born += (g_KSP/(s-MP*MP))*(F1_p+kappa_p*F2_p)+g_KSP*(F1_S+kappa_L*F2_S)+(g_KLP/(u-ML*ML))*(((ML+MS0)*kappa_SL)/2./MP)*F2_L;
	//A2_born += (g_KSP/(s-MP*MP)/(t-MK*MK))*((F_k+F1_p)+(F_k-F1_p)*((pgpk+pgpp)/(pgpL)));
	//A3_born += (g_KSP/(s-MP*MP))*(kappa_p/MP)*F2_p;
	//A4_born += (g_KSP/(u-MS0*MS0))*(kappa_L/MS0)*F2_S+(g_KLP/(u-ML*ML))*(kappa_SL/MP)*F2_L;
	//A5_born += (g_KSP/(s-MP*MP))*(kappa_p/2./MP)*F2_p-(g_KSP/(u-MS0*MS0))*(kappa_L/2./MS0)*F2_S-(g_KLP/(u-ML*ML))*(kappa_SL/2./MP)*F2_L;
    //A6_born += ((-2.*g_KSP)/(s-MP*MP)/(t-MK*MK))*(F_k+(F_k-F1_p)*(pgpg-2.*pgpk)*(pgpp/pgpg/pgpL))+((g_KSP)/(u-MS0*MS0))*2.*F1_S/pgpg;
cout<<"A1="<<A1_born<<endl;
cout<<"A2="<<A2_born<<endl;
cout<<"A3="<<A3_born<<endl;
cout<<"A4="<<A4_born<<endl;
cout<<"A5="<<A5_born<<endl;
cout<<"A6="<<A6_born<<endl;
	double MKS=896.1;//K*0
	double Gamma_KS=50.5;
	double MK1=1273.0;
	double Gamma_K1=90.0;
	double MN1=1440.0;
	double Gamma_N1=350.0;
	double MN7=1720.0;
	double Gamma_N7=150.0;
	double MN8=1680.0;
	double Gamma_N8=150.0;
	double ML1=1407.0;
	double Gamma_L1=50.0;
	double ML3=1670.0;
	double Gamma_L3=35.0;
	double ML5=1810.0;
	double Gamma_L5=150.0;
	double MS1=1660.0;
	double Gamma_S1=100.0;
	double GV=-0.05*4.*PI;
	double GT=0.16*4.*PI;
	double GV1=-0.19*4.*PI;
	double GT1=-0.35*4.*PI;
	GV=GV/4./PI;
	GT=GT/4./PI;
	GV1=GV1/4./PI;
	GT1=GT1/4./PI;
	double GN1=-0.01*sqrt(4.*PI);
	double GN7_a=-0.04*4.*PI;
	double GN7_b=-0.14*4.*PI;
	double GN8_a=-0.63*4.*PI;
	double GN8_b=-0.05*4.*PI;
	double GL1=-0.31*sqrt(4.*PI*alpha);
	double GL3=1.18*sqrt(4.*PI*alpha);
	double GL5=-1.25*sqrt(4.*PI*alpha);
	double GS1=-4.96*sqrt(4.*PI*alpha);

	double M=1000.;//[MeV/c^2] renormalization

	double re;
	double im;
	double com;
	double real_part;
	double imaginary_part;

	double re2;
	double im2;
	double com2;
	double real_part2;
	double imaginary_part2;

	
	/*A1_kstar*/
	com=(GV/M)*(ML+MP)*F_ks;
	re=t-MKS*MKS;
	im=MKS*Gamma_KS;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	com2=(GT/M)*(t/(ML+MP))*F_ks;
	re2=t-MKS*MKS;
	im2=MKS*Gamma_KS;
	real_part2=com2*(re2/(re2*re2+im2*im2));
	imaginary_part2=com2*(-1.*im2/(re2*re2+im2*im2));
	TComplex A1_kstar(real_part+real_part2,imaginary_part+imaginary_part2);
	/*A2_kstar*/
	com=(GT/M)*(1./(ML+MP))*F_ks;
	re=t-MKS*MKS;
	im=MKS*Gamma_KS;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A2_kstar(real_part,imaginary_part);
	/*A3_kstar*/
	com=(GV/M)*F_ks;
	re=t-MKS*MKS;
	im=MKS*Gamma_KS;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	com2=-1.*(GT/M)*((ML-MP)/(ML+MP))*F_ks;
	re2=t-MKS*MKS;
	im2=MKS*Gamma_KS;
	real_part2=com2*(re2/(re2*re2+im2*im2));
	imaginary_part2=com2*(-1.*im2/(re2*re2+im2*im2));
	TComplex A3_kstar(real_part+real_part2,imaginary_part+imaginary_part2);
	/*A4_kstar*/
	com=(GV/M)*F_ks;
	re=t-MKS*MKS;
	im=MKS*Gamma_KS;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	com2=(GT/M)*((ML-MP)/(ML+MP))*F_ks;
	re2=t-MKS*MKS;
	im2=MKS*Gamma_KS;
	real_part2=com2*(re2/(re2*re2+im2*im2));
	imaginary_part2=com2*(-1.*im2/(re2*re2+im2*im2));
	TComplex A4_kstar(real_part+real_part2,imaginary_part+imaginary_part2);
	/*A5_kstar*/
	TComplex A5_kstar(0.,0.);
	/*A6_kstar*/
	TComplex A6_kstar(0.,0.);
	
cout<<"A1_kstar"<<A1_kstar<<endl;
cout<<"A2_kstar"<<A2_kstar<<endl;
cout<<"A3_kstar"<<A3_kstar<<endl;
cout<<"A4_kstar"<<A4_kstar<<endl;
cout<<"A5_kstar"<<A5_kstar<<endl;
cout<<"A6_kstar"<<A6_kstar<<endl;

	/*A1_k1*/
	TComplex A1_k1(0.,0.);
	/*A2_k1*/
	com=-1.*(GT1/M)*(1./(ML+MP))*F_k1;
	re=t-MK1*MK1;
	im=MK1*Gamma_K1;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A2_k1(real_part,imaginary_part);
	/*A3_k1*/
	com=(GV1/M)*F_k1;
	re=t-MK1*MK1;
	im=MK1*Gamma_K1;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	com2=(GT1/M)*((ML-MP)/(ML+MP))*F_k1;
	re2=t-MK1*MK1;
	im2=MK1*Gamma_K1;
	real_part2=com2*(re2/(re2*re2+im2*im2));
	imaginary_part2=com2*(-1.*im2/(re2*re2+im2*im2));
	TComplex A3_k1(real_part+real_part2,imaginary_part+imaginary_part2);
	/*A4_k1*/
	com=-1.*(GV1/M)*F_k1;
	re=t-MK1*MK1;
	im=MK1*Gamma_K1;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	com2=-1.*(GT1/M)*((ML-MP)/(ML+MP))*F_k1;
	re2=t-MK1*MK1;
	im2=MK1*Gamma_K1;
	real_part2=com2*(re2/(re2*re2+im2*im2));
	imaginary_part2=com2*(-1.*im2/(re2*re2+im2*im2));
	TComplex A4_k1(real_part+real_part2,imaginary_part+imaginary_part2);
	/*A5_k1*/
	TComplex A5_k1(0.,0.);
	/*A6_k1*/
	TComplex A6_k1(0.,0.);

cout<<"A1_k1"<<A1_k1<<endl;
cout<<"A2_k1"<<A2_k1<<endl;
cout<<"A3_k1"<<A3_k1<<endl;
cout<<"A4_k1"<<A4_k1<<endl;
cout<<"A5_k1"<<A5_k1<<endl;
cout<<"A6_k1"<<A6_k1<<endl;

//L1: Lambda(1405), 1/2-
	/*A1_L1*/
	com=(GL1)*(ML1+ML)*kappa_L*(ML1-ML)*F2_L/(2.*MP*(ML1+ML));
	re=u-ML1*ML1;
	im=ML1*Gamma_L1;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A1_L1(real_part,imaginary_part);
	/*A2_L1*/
	TComplex A2_L1(0.,0.);
	/*A3_L1*/
	TComplex A3_L1(0.,0.);
	/*A4_L1*/
	com=(-1.*GL1)*kappa_L*F2_L/(MP);
	re=u-ML1*ML1;
	im=ML1*Gamma_L1;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A4_L1(real_part,imaginary_part);
	/*A5_L1*/
	com=(GL1)*kappa_L*F2_L/(2.*2.*MP);
	re=u-ML1*ML1;
	im=ML1*Gamma_L1;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A5_L1(real_part,imaginary_part);
	/*A6_L1*/
	TComplex A6_L1(0.,0.);
cout<<"A1_L1"<<A1_L1<<endl;
cout<<"A2_L1"<<A2_L1<<endl;
cout<<"A3_L1"<<A3_L1<<endl;
cout<<"A4_L1"<<A4_L1<<endl;
cout<<"A5_L1"<<A5_L1<<endl;
cout<<"A6_L1"<<A6_L1<<endl;

//L3: Lambda(1670), 1/2-
	/*A1_L3*/
	com=(GL3)*(ML3+ML)*kappa_L*(ML3-ML)*F2_L/(2.*MP*(ML3+ML));
	re=u-ML3*ML3;
	im=ML3*Gamma_L3;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A1_L3(real_part,imaginary_part);
	/*A2_L3*/
	TComplex A2_L3(0.,0.);
	/*A3_L3*/
	TComplex A3_L3(0.,0.);
	/*A4_L3*/
	com=(-1.*GL3)*kappa_L*F2_L/(MP);
	re=u-ML3*ML3;
	im=ML3*Gamma_L3;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A4_L3(real_part,imaginary_part);
	/*A5_L3*/
	com=(GL3)*kappa_L*F2_L/(2.*2.*MP);
	re=u-ML3*ML3;
	im=ML3*Gamma_L3;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A5_L3(real_part,imaginary_part);
	/*A6_L3*/
	TComplex A6_L3(0.,0.);
cout<<"A1_L3"<<A1_L3<<endl;
cout<<"A2_L3"<<A2_L3<<endl;
cout<<"A3_L3"<<A3_L3<<endl;
cout<<"A4_L3"<<A4_L3<<endl;
cout<<"A5_L3"<<A5_L3<<endl;
cout<<"A6_L3"<<A6_L3<<endl;

//L5: Lambda(1810), 1/2+
	/*A1_L5*/
	com=(GL5)*(ML5+ML)*kappa_L*(ML5-ML)*F2_L/(2.*MP*(ML5-ML));
	re=u-ML5*ML5;
	im=ML5*Gamma_L5;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A1_L5(real_part,imaginary_part);
	/*A2_L5*/
	TComplex A2_L5(0.,0.);
	/*A3_L5*/
	TComplex A3_L5(0.,0.);
	/*A4_L5*/
	com=(GL5)*kappa_L*F2_L/(MP);
	re=u-ML5*ML5;
	im=ML5*Gamma_L5;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A4_L5(real_part,imaginary_part);
	/*A5_L5*/
	com=(-1.*GL5)*kappa_L*F2_L/(2.*2.*MP);
	re=u-ML5*ML5;
	im=ML5*Gamma_L5;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A5_L5(real_part,imaginary_part);
	/*A6_L5*/
	TComplex A6_L5(0.,0.);
cout<<"A1_L5"<<A1_L5<<endl;
cout<<"A2_L5"<<A2_L5<<endl;
cout<<"A3_L5"<<A3_L5<<endl;
cout<<"A4_L5"<<A4_L5<<endl;
cout<<"A5_L5"<<A5_L5<<endl;
cout<<"A6_L5"<<A6_L5<<endl;

//S1: Sigma(1660), 1/2+
	/*A1_S1*/
	com=(GS1)*(MS1+ML)*kappa_S*(MS1-ML)*F2_S/(2.*MP*(MS1+ML));
	re=u-MS1*MS1;
	im=MS1*Gamma_S1;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A1_S1(real_part,imaginary_part);
	/*A2_S1*/
	TComplex A2_S1(0.,0.);
	/*A3_S1*/
	TComplex A3_S1(0.,0.);
	/*A4_S1*/
	com=(GS1)*kappa_S*F2_S/(MP);
	re=u-MS1*MS1;
	im=MS1*Gamma_S1;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A4_S1(real_part,imaginary_part);
	/*A5_S1*/
	com=(-1.*GS1)*kappa_S*F2_S/(2.*2.*MP);
	re=u-MS1*MS1;
	im=MS1*Gamma_S1;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A5_S1(real_part,imaginary_part);
	/*A6_S1*/
	TComplex A6_S1(0.,0.);
cout<<"A1_S1"<<A1_S1<<endl;
cout<<"A2_S1"<<A2_S1<<endl;
cout<<"A3_S1"<<A3_S1<<endl;
cout<<"A4_S1"<<A4_S1<<endl;
cout<<"A5_S1"<<A5_S1<<endl;
cout<<"A6_S1"<<A6_S1<<endl;

//N7: N(1720), 3/2+
//okuyama
	double G1=GN7_a*F2_p/4./PI/(-3.16*sqrt(4.*PI));
	double G2=GN7_b*F2_p/4./PI/(-3.16*sqrt(4.*PI));
	double pLq = 0.5*(s+ML*ML-MK*MK);
	double A32 = 2.*pLq-W*ML;
	double B32 = pLq+W*ML;
	double P11 = (3.*W-MP)*A32/6./s+ML-pgpL/(W+MP);
	double P12 = -1./(W+MP);
	double P13 = ML/(W+MP)-MP*A32/(3.*s*(W+MP));
	double P14 = 1.;
	double P21 = -B32*(W-MP)/(6.*W*(W+MP));
	double P22 = (W-MP)/2./(W+MP)/(W+MP);
	double P23 = pgpL/(W+MP)/(W+MP)-B32*(W-MP)/(3.*(W+MP)*(W+MP)*W);
	double P24 = -(W-MP)/(2.*(W+MP));
	double R11 = A32/(6.*s*(W+MP)); 
	double R12 = A32/(6.*s*(W+MP)*pgpL);
	double R13 = 0.;
	double R14 = 0.;
	double R21 = B32/(6.*(W+MP)*(W+MP)*W);
	double R22 = -(W*ML*MP-2.*MP*pLq+3.*W*pLq)/(6.*pgpL*(W+MP)*(W+MP)*s);
	double R23 = -A32/(3.*(W+MP)*(W+MP)*s);
	double R24 = 1./(2.*(W+MP)*(W+MP));
	double E11 = P11+pgpg*R11;
	double E12 = P12+pgpg*R12;
	double E13 = P13+pgpg*R13;
	double E14 = P14+pgpg*R14;
	double E21 = P21+pgpg*R21;
	double E22 = P22+pgpg*R22;
	double E23 = P23+pgpg*R23;
	double E24 = P24+pgpg*R24;
	double E15 = -A32/3./s;
	double E16 = pgpp*A32/(3.*s*(W+MP)*pgpL)+1./(W+MP);
	double E25 = pgpp*A32/(3.*(W+MP)*(W+MP)*s);
	double E26 = (-pgpp*(W*ML*MP-2.*MP*pLq+3.*W*pLq))/(3.*pgpL*(W+MP)*(W+MP)*s);
cout<<"P11="<<P11<<endl;
cout<<"P12="<<P12<<endl;
cout<<"P13="<<P13<<endl;
cout<<"P14="<<P14<<endl;
cout<<"P21="<<P21<<endl;
cout<<"P22="<<P22<<endl;
cout<<"P23="<<P23<<endl;
cout<<"P24="<<P24<<endl;

	com=G1*E11+G2*E21;
	re=s-MN7*MN7;
	im=MN7*Gamma_N7;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A1_N7(real_part,imaginary_part);

	com=G1*E12+G2*E22;
	re=s-MN7*MN7;
	im=MN7*Gamma_N7;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A2_N7(real_part,imaginary_part);
	
	com=G1*E13+G2*E23;
	re=s-MN7*MN7;
	im=MN7*Gamma_N7;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A3_N7(real_part,imaginary_part);

	com=G1*E14+G2*E24;
	re=s-MN7*MN7;
	im=MN7*Gamma_N7;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A4_N7(real_part,imaginary_part);

	com=G1*E15+G2*E25;
	re=s-MN7*MN7;
	im=MN7*Gamma_N7;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A5_N7(real_part,imaginary_part);

	com=G1*E16+G2*E26;
	re=s-MN7*MN7;
	im=MN7*Gamma_N7;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A6_N7(real_part,imaginary_part);
cout<<"A1_N7"<<A1_N7<<endl;
cout<<"A2_N7"<<A2_N7<<endl;
cout<<"A3_N7"<<A3_N7<<endl;
cout<<"A4_N7"<<A4_N7<<endl;
cout<<"A5_N7"<<A5_N7<<endl;
cout<<"A6_N7"<<A6_N7<<endl;


//N8: N(1675), 5/2-
//A(-)=-A(+)[W --> -W]
	double Wp = W+MP;
	double Wm = W-MP;
//P11 = ((2.*ML*W*(2.*W-MP)*pLq-2.*(4.*W+MP)*pLq*pLq)*Wp-ML*ML*W*W*Wp*Wp)*Wm*Wm/(2.*Wm*W*W*W*W);
//P11 += ((4.*W*W*(7.*W+4.*MP)*pLq*pgpL-4.*ML*W*W*W*(2.*W-MP)*pgpL)*Wm-20.*pgpL*pgpL*W*W*W*W)/(2.*Wm*W*W*W*W);
//P12 = ((4.*pLq+ML*W)*Wp*Wm-10.*pgpL*W*W)/(Wm*W*W);
//P13 = ((2.*pLq*pLq*MP-2.*ML*W*(2.*W-MP)*pLq)*Wp*Wm-ML*ML*W*W*Wp*Wm*Wm-8.*W*W*pgpL*MP*pLq)/(Wm*W*W*W*W);
//P13 += 2.*ML*W*W*W*(5.*W-MP)*pgpL/(Wm*W*W*W*W);
//P14 = (2.*Wp*(3.*W-2.*MP)*pLq-10.*pgpL*W*W-ML*W*Wp*Wp)/W/W;
//P21 = ((pLq*pLq-ML*W*pLq)*Wp*Wp*Wm+(2.*ML*W*W*W*pgpL-2.*pgpL*W*W*pLq)*Wp)/W/W;
//P22 = ((-4.*pLq*-ML*W)*Wp*Wp*Wm+10.*W*W*Wp*pgpL)/(2.*W*W*Wm*Wm);
//P23 = ((ML*W*pLq-pLq*pLq)*Wp*Wp*Wm+(4.*W*(2.*W-MP)*pgpL*pLq-ML*W*W*(3.*W+MP)*pgpL)*Wp)/(W*W*W*Wm*Wm)-10.*pgpL*pgpL*W*W*W/(W*W*W*Wm*Wm);
//P24 = (10.*W*W*Wp*pgpL-2.*(3.*W-2.*MP)*Wp*Wp*pLq+ML*W*Wp*Wp*Wp)/(2.*W*W*Wm);
	Wp = -Wp;
	Wm = -Wm;
P11 = ((-2.*ML*W*(-2.*W+MP)*pLq-2.*(-4.*W-MP)*pLq*pLq)*Wp-ML*ML*W*W*Wp*Wp)*Wm*Wm/(2.*Wm*W*W*W*W);
P11 = P11+((4.*W*W*(-7.*W-4.*MP)*pLq*pgpL+4.*ML*W*W*W*(-2.*W+MP)*pgpL)*Wm-20.*pgpL*pgpL*W*W*W*W)/(2.*Wm*W*W*W*W);
P12 = ((4.*pLq-ML*W)*Wp*Wm-10.*pgpL*W*W)/(Wm*W*W);
P13 = ((-2.*pLq*pLq*MP+2.*ML*W*(-2.*W+MP)*pLq)*Wp*Wm-ML*ML*W*W*Wp*Wm*Wm+8.*W*W*pgpL*MP*pLq)/(Wm*W*W*W*W);
P13 = P13-2.*ML*W*W*W*(-5.*W+MP)*pgpL/(Wm*W*W*W*W);
P14 = (2.*Wp*(-3.*W+2.*MP)*pLq-10.*pgpL*W*W+ML*W*Wp*Wp)/W/W;
P21 = ((pLq*pLq+ML*W*pLq)*Wp*Wp*Wm+(-2.*ML*W*W*W*pgpL-2.*pgpL*W*W*pLq)*Wp)/W/W;
P22 = ((-4.*pLq*+ML*W)*Wp*Wp*Wm+10.*W*W*Wp*pgpL)/(2.*W*W*Wm*Wm);
P23 = ((-ML*W*pLq-pLq*pLq)*Wp*Wp*Wm+(-4.*W*(-2.*W+MP)*pgpL*pLq-ML*W*W*(-3.*W-MP)*pgpL)*Wp)/(-W*W*W*Wm*Wm)-10.*pgpL*pgpL*W*W*W/(W*W*W*Wm*Wm);
P24 = (10.*W*W*Wp*pgpL-2.*(-3.*W+2.*MP)*Wp*Wp*pLq-ML*W*Wp*Wp*Wp)/(2.*W*W*Wm);
cout<<"P11="<<P11<<endl;
cout<<"P12="<<P12<<endl;
cout<<"P13="<<P13<<endl;
cout<<"P14="<<P14<<endl;
cout<<"P21="<<P21<<endl;
cout<<"P22="<<P22<<endl;
cout<<"P23="<<P23<<endl;
cout<<"P24="<<P24<<endl;

	G1=GN8_a*F2_p/4./PI/(-3.16*sqrt(4.*PI))/MN8/MN8/MN8/MN8;
	G2=GN8_b*F2_p/4./PI/(-3.16*sqrt(4.*PI))/MN8/MN8/MN8/MN8;
	
	com=0.1*(G1*P11+G2*P21);
	re=s-MN8*MN8;
	im=MN8*Gamma_N8;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A1_N8(real_part,imaginary_part);
	
	com=0.1*(G1*P12+G2*P22);
	re=s-MN8*MN8;
	im=MN8*Gamma_N8;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A2_N8(real_part,imaginary_part);
	
	com=0.1*(G1*P13+G2*P23);
	re=s-MN8*MN8;
	im=MN8*Gamma_N8;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A3_N8(real_part,imaginary_part);

	com=0.1*(G1*P14+G2*P24);
	re=s-MN8*MN8;
	im=MN8*Gamma_N8;
	real_part=com*(re/(re*re+im*im));
	imaginary_part=com*(-1.*im/(re*re+im*im));
	TComplex A4_N8(real_part,imaginary_part);
	TComplex A5_N8(0.,0.);
	TComplex A6_N8(0.,0.);
cout<<"A1_N8"<<A1_N8<<endl;
cout<<"A2_N8"<<A2_N8<<endl;
cout<<"A3_N8"<<A3_N8<<endl;
cout<<"A4_N8"<<A4_N8<<endl;
cout<<"A5_N8"<<A5_N8<<endl;
cout<<"A6_N8"<<A6_N8<<endl;


	TComplex A1, A2, A3, A4, A5, A6;
	A1=A1_born+A1_L1+A1_L5+A1_kstar+A1_k1+A1_L3+A1_S1+A1_N7;// +A1_N8
	A2=A2_born+A2_L1+A2_L5+A2_kstar+A2_k1+A2_L3+A2_S1+A2_N7;// +A2_N8
	A3=A3_born+A3_L1+A3_L5+A3_kstar+A3_k1+A3_L3+A3_S1+A3_N7;// +A3_N8
	A4=A4_born+A4_L1+A4_L5+A4_kstar+A4_k1+A4_L3+A4_S1+A4_N7;// +A4_N8
	A5=A5_born+A5_L1+A5_L5+A5_kstar+A5_k1+A5_L3+A5_S1+A5_N7;// +A5_N8
	A6=A6_born+A6_L1+A6_L5+A6_kstar+A6_k1+A6_L3+A6_S1+A6_N7;// +A6_N8
	A1=A1_born;
	A2=A2_born;
	A3=A3_born;
	A4=A4_born;
	A5=A5_born;
	A6=A6_born;
cout<<"A1_all"<<A1<<endl;
cout<<"A2_all"<<A2<<endl;
cout<<"A3_all"<<A3<<endl;
cout<<"A4_all"<<A4<<endl;
cout<<"A5_all"<<A5<<endl;
cout<<"A6_all"<<A6<<endl;


	TComplex f1 = ((W-MP)*A1-pgpp*A3-pgpL*A4-pgpg*A5);
	TComplex f2 = ((pg*pK)/(EL+ML)/(EP+MP))*((W+MP)*A1+pgpp*A3+pgpL*A4+pgpg*A5);
	TComplex f3 = (pg*pK/(EP+MP))*(-2.*pgpp*A2+(W+MP)*A4+pgpg*A6);
	TComplex f4 = (pK*pK/(EL+ML))*(2.*pgpp*A2+(W-MP)*A4-pgpg*A6);
	TComplex f5 = (pg*pg/(EP+MP))*(-1.*A1+2.*pgpL*A2+(W+MP)*(A3-A5)+pgpL*A6);
	TComplex f6 = (pg*pK/(EL+ML))*(-2.*pgpL*A2+(W-MP)*A3-pgpL*A6-(Eg*A1+pgpp*A3+pgpL*A4+Eg*(W+MP)*A5)/(EP+MP));
	TComplex f7 = f1+f3*cos(theta)+f5;
	TComplex f8 = f4*cos(theta)+f6;
	double AA = (pK/2./W)*(EP+MP)*(EL+ML)/(s-MP*MP);
cout<<"A3factor="<<(pgpp)<<endl;//*(-2.*pgpp)<<endl;
cout<<"f1="<<f1<<endl;
cout<<"f2="<<f2<<endl;
cout<<"f3="<<f3<<endl;
cout<<"f4="<<f4<<endl;
cout<<"f5="<<f5<<endl;
cout<<"f6="<<f6<<endl;
cout<<"f7="<<f7<<endl;
cout<<"f8="<<f8<<endl;
cout<<"AA="<<AA<<endl;
	
	TComplex test1(3,5);
	TComplex test2(1,7);
cout<<"test"<<endl;
cout<<test1<<endl;
cout<<test1.Rho2()<<endl;
cout<<test2<<endl;
cout<<(test1.Conjugate(test1)*test2).Re()<<endl;
	double sigma_U = AA*(f1.Rho2()+f2.Rho2()+2.*(f1.Conjugate(f1)*f2).Re()*cos(theta)+0.5*sin(theta)*sin(theta)*(f3.Rho2()+f4.Rho2()+2.*(f1.Conjugate(f1)*f4-f2.Conjugate(f2)*f3+f3.Conjugate(f3)*f4*cos(theta)).Re()));
	double sigma_L = AA*(f7.Rho2()+f8.Rho2()+2.*(f7.Conjugate(f7)*f8).Re()*cos(theta));
	double sigma_P = AA*(0.5*f3.Rho2()+0.5*f4.Rho2()+(f1.Conjugate(f1)*f4-f2.Conjugate(f2)*f3+f3.Conjugate(f3)*f4*cos(theta)).Re());
	double sigma_I = AA*(f7*(-f2.Conjugate(f2)+f3.Conjugate(f3)+f4.Conjugate(f4)*cos(theta))+f8*(f1.Conjugate(f1)+f3.Conjugate(f3)*cos(theta)+f4.Conjugate(f4))).Re();



//	double A1 = (g_KLP/(s-MP*MP))*(F1+F2);
//	double A2 = 2.*(g_KLP/(s-MP*MP))*F1;
//	double A3 = (g_KLP/(s-MP*MP))*(F2/MP);
//	double A4 = 0.;
//	double A5 = 0.;
//	double A6 = -1.*A3/2.;
	
//	double f1 = N*(-(W-MP)*A1+k_p*A4+k_pL*A5-k_k*A6);
//	double f2 = N*((Eg*pK)/(EL+ML)/(EP+MP))*((W+MP)*A1+k_p*A4+k_pL*A5-k_k*A6);
//	double f3 = -N*(Eg*pK/(EP+MP))*(A3+(W+MP)*A5);
//	double f4 = N*(pK*pK/(EL+ML))*(A3-(W+MP)*A5);
//	double f5 = N*(Eg*Eg/(EP+MP))*(A1-((k_k+k_p)*A2+k_pL*A3)/k_k - (W+MP)*(A4+A6));
//	double f6 = N*(Eg*Eg*pK/(EL+ML)/(EP+MP))*(A1-MP*A4+(k_pL/Eg)*A5+((EP+MP)/(Eg*k_k))*((k_k+k_p)*A2+k_pL*A3)-(W+MP)*A6);
	double f5t = f1+f3*cos(theta)+f5;
	double f6t = f4*cos(theta)+f6;

	double sigma_T = C*(f1*f1+f2*f2-2.*f1*f2*cos(theta)+sin(theta)*sin(theta)*((f3*f3+f4*f4)/2.+f1*f4+f2*f3+f3*f4*cos(theta)));
//	double sigma_L = C*(f5t*f5t+f6t*f6t+2.*f5t*f6t*cos(theta));
	double sigma_TT = C*((f3*f3+f4*f4)/2.+f1*f4+f2*f3+f4*f3*cos(theta))*sin(theta)*sin(theta);
	double sigma_LT = -C*((f1+f4)*f6t+(f2+f3)*f5t+(f3*f6t+f4*f5t)*cos(theta));

	//return sigma_T+eps*sigma_L+eps*sigma_TT*cos(2.*phiK)+sqrt(2.*epsL*(eps+1))*sigma_LT*cos(phiK);
	//return (sigma_U+epsL*sigma_L)*400.*1000000.;//[ub/sr]
	return (sigma_U+epsL*sigma_L+eps*sigma_P*sin(theta)*sin(theta)*cos(2.*phiK)+sqrt(2.*epsL*(1.+eps))*sigma_I*sin(theta)*phiK)*1000000.;//[ub/sr]
}

void difcs(){

	string clas2010 = "./CLAS_2010.dat";//CLAS(2010) McCracken D-Thesis 
	string clas2006 = "./CLAS_2006.dat";//CLAS(2006) Bradform Phys.Rev.C73, 035202 
	string kaonmaid = "./KaonMaid.dat";//Kaon-Maid web page 

 const double PI=3.14159265359;

	string buf;
	int nofdata = 100;
	int npoint = 0;
	int npoint2 = 0;
	int npoint3 = 0;
	double cos_gk, dif_cs, stat_err, syst_err, theta_gk;
	double x[nofdata], y[nofdata], xe[nofdata], ye[nofdata];
	double x2[nofdata], y2[nofdata], xe2[nofdata], ye2[nofdata];
	double x3[nofdata], y3[nofdata], xe3[nofdata], ye3[nofdata];



/*----- CLAS 2010 -----*/
	ifstream ifp(clas2010.c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << clas2010.c_str() << endl;
	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		sbuf >> cos_gk >> dif_cs >> stat_err >> syst_err;
		cout << cos_gk << ", " << dif_cs << ", " << stat_err << ", " << syst_err  <<endl;

		x[npoint] = 180.*acos(cos_gk)/PI;
		y[npoint] = dif_cs/2./PI;
		xe[npoint]= 0.; 
		ye[npoint]= stat_err; 
		npoint++;

	}

/*----- CLAS 2006 -----*/
	ifstream ifp2(clas2006.c_str(),ios::in);
	if (ifp2.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << clas2010.c_str() << endl;
	while(1){
		getline(ifp2,buf);
		if(buf[0]=='#'){continue;}
		if(ifp2.eof())break;
		stringstream sbuf(buf);
		sbuf >> cos_gk >> dif_cs >> stat_err;
		cout << cos_gk << ", " << dif_cs << ", " << stat_err  <<endl;

		x2[npoint2] = 180.*acos(cos_gk)/PI;
		y2[npoint2] = dif_cs/2./PI;
		xe2[npoint2]= 0.; 
		ye2[npoint2]= stat_err; 
		npoint2++;
	}

/*----- Kaon-Maid -----*/
	ifstream ifp3(kaonmaid.c_str(),ios::in);
	if (ifp3.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << kaonmaid.c_str() << endl;
	while(1){
		getline(ifp3,buf);
		if(buf[0]=='#'){continue;}
		if(ifp3.eof())break;
		stringstream sbuf(buf);
		sbuf >> theta_gk >> dif_cs;
		cout << theta_gk << ", " << dif_cs << endl;

		x3[npoint3] = theta_gk;
		y3[npoint3] = dif_cs;
		xe3[npoint3]= 0.; 
		ye3[npoint3]= 0.; 
		npoint3++;
	}
		TGraphErrors *g = new TGraphErrors( npoint, x, y, xe, ye);
		g->SetMarkerStyle(21);
		g->SetMarkerColor(kAzure);
		g->SetMarkerSize(2.0);
		TGraphErrors *g2 = new TGraphErrors( npoint2, x2, y2, xe2, ye2);
		g2->SetMarkerStyle(21);
		g2->SetMarkerColor(kRed);
		g2->SetMarkerSize(2.0);
		TGraphErrors *g3 = new TGraphErrors( npoint3, x3, y3, xe3, ye3);
		g3->SetMarkerStyle(3);
		g3->SetMarkerColor(kGreen);
		g3->SetMarkerSize(1.0);
		TF1 *f1 = new TF1("f1",RPR,0.,180.,2);
		//f1->SetParameter(0,2250);
		//f1->SetParameter(1,450000);
		f1->SetParameter(0,2200);
		f1->SetParameter(1,70000);
		//f1->SetParameter(1,0.01);
		f1->SetLineColor(kViolet);
	TCanvas *c = new TCanvas("c","plot",800.,600.);
	c->cd()->SetGrid();
	TH1 *frame = c->DrawFrame(0.,0.,180.,0.65);
frame->GetXaxis()->SetTitleSize(0.05);
frame->GetXaxis()->SetTitleOffset(1.05);
frame->GetXaxis()->SetTitle("#theta_{#gammaK}^{CM} [degree]");
frame->GetXaxis()->SetLabelOffset(0.02);
frame->GetYaxis()->SetTitleOffset(1.1);
frame->GetYaxis()->SetTitleSize(0.05);
frame->GetYaxis()->SetTitle("d#sigma/d#Omega [#mub/sr] (CM)");
		g3->Draw("psame");
		g2->Draw("psame");
		g->Draw("psame");
		f1->Draw("");

}
