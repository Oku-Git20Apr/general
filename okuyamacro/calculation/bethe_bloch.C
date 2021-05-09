//Bethe-Bloch Calc.
//
//K.Okuyama (Apr. 2, 2021)

double gauss1(double nsigmax){
//**************************************************************************
//* This subroutine generates a random number distributed about a Gaussian
//* centered at zero with a standrd deviation of 1.  The gaussian is
//* truncated at nsigmax.  Use as a function like: VAL' = VAL + WIDTH*X
//*
//* Algorithm is from PHYSICS LETTERS B V.204, P.83 ''Review of particle
//* properties", Particle Data Group. (with nsigmax cutoff added).
//**************************************************************************
	double u1,u2,v1,v2,s;
	double grnd, val;

	s = 0.;
	val = 0.;
	while (s>1. || s==0. || abs(val)>nsigmax){
    gRandom->SetSeed(0);
	u1 = gRandom->Uniform();
	u2 = gRandom->Uniform();
	v1 = 2.0*u1-1.0;
	v2 = 2.0*u2-1.0;
	s  = v1*v1+v2*v2;
	val = v1*sqrt(-2.*log(s)/s);  //! <--want a natural log here
	}

	return val;
}

void bethe_bloch(){

int typeflag = 4;//
	//1=normal eloss (picked from distribution)
	//2=min eloss
	//3=max eloss
	//4=most probable eloss
	
double logterm;
double c = 2.998*pow(10.0,8.0);
double me = 0.511; //mass of electron[MeV/c^2]
double z = 1.0; //charge of electron
double M, Z, A, I, K, rho, delta, C, Kend, step, interval, Wmax, energyloss, Range, gamma, betasq, s; 
double beta, Eloss_temp, Eloss_mp, Eloss, denscorr, hnup, log10bg, CO, chsi, lambda, x;



cout<<"入射粒子が物質中で減速するときのenergy depositとRangeをBethe-Blochの式から計算します。\n"; 
cout<<"入射粒子の質量[MeV/c^2]を入力してください:"; cin>>M;
cout<<"物質の原子番号(Z):"; cin>>Z;
cout<<"物質の原子量(A):"; cin>>A;
//cout<<"イオン化エネルギー(I)[eV]:"; cin>>I;
//I = I * pow(10.0,-6.0);
cout<<"入射粒子の入射エネルギー(Ein)[MeV]:"; cin>>K;
cout<<"物質の密度[g/cm^3]:"; cin>>rho;
//cout<<"density correction:"; cin>>delta;
//cout<<"shell correction:"; cin>>C;
cout<<"入射粒子の最終エネルギー(Eout)[MeV]:"; cin>>Kend;
cout<<"step数:"; cin>>step;

//M=105.66;Z=82;A=208.2;K=750;rho=11.35;delta=0;C=0;Kend=0;step=15; //for example (incident particle:muon, material:Pb, 750->0[MeV])

//Ionization potential
	if(Z==1) I = 21.8 * pow(10.,-6.);
	else I = (16.*pow(Z,0.9)) * pow(10.,-6.);



s = me/M;
interval = (K-Kend)/step;
cout<<"interval="<<interval<<"[MeV]\n";
Range = 0.0;

cout<<setw(12)<<"range[MeV]:"<<setw(11)<<"dE/dρx"<<setw(4)<<"/"<<setw(12)<<"Δx=ΔE(dE/dρx)^(-1)"<<endl;
for( ; M-Kend<=K; K -= (double)(int)interval){ //表示のときだけ少数点以下を表示しない

	gamma = K/M;
	betasq = 1 - 1/(gamma*gamma);
	beta = sqrt(betasq);

	hnup = 28.816*pow(10.,-6.)*sqrt(rho*Z/A); //plasma frequency
	log10bg = log(beta*gamma)/log(10.);
	CO=log(hnup)-log(I)+0.5;

//Density correction
	if(log10bg<0.) denscorr=0.;
	else if(log10bg<3.) denscorr=CO+log(10.)*log10bg+abs(CO/27.)*pow((3.-log10bg),3.);
	else if(log10bg<4.7) denscorr=CO+log(10.)*log10bg;
	else denscorr=CO+log(10.)*4.7;
denscorr = 0.;

	Wmax = (2*me*betasq*pow(gamma,2.0))/(1+2*s*sqrt(1+betasq*pow(gamma,2.0))+pow(s,2.0));
	logterm = (2*me*pow(gamma,2.0)*betasq*Wmax)/pow(I,2.0);
	
	energyloss = 0.1535*(Z/A)*(pow(z,2.0)/betasq)*(log(logterm)-2*betasq-delta-2*(C/Z));

	  Eloss_temp = 0.1536*pow(10.,-3.) * Z/A/betasq* (log(me/I/I) + 1.063 + 2.*log(gamma*beta) +	log(0.1536*Z/A/betasq)-betasq-denscorr);
	  Eloss_mp = Eloss_temp*1000.;
	  chsi = 0.307075/2.*Z/A/betasq;
	  if(typeflag==1) x=abs(gauss1(10.));
	  else if(typeflag==2) x=3.;
	  else if(typeflag==3) x=0.0067;
	  else if(typeflag==4) x=1.;

	  if(x>0.) lambda = -2.0*log(x);
	  else lambda = 100000.;

	  Eloss = lambda*chsi+Eloss_mp;
	  energyloss = Eloss;

	if (energyloss<0) break;
	cout<<setw(5)<<K<<setw(1)<<"-"<<setw(5)<<K-(double)(int)interval<<setw(1)<<":"<<setw(10)<<energyloss<<setw(4)<<"/"<<setw(5)<<interval/energyloss<<endl;
	Range += interval/energyloss;

}
cout<<"Range:"<<Range<<"[g/cm^2]="<<Range/rho<<"[cm]"<<endl;

	

} 
