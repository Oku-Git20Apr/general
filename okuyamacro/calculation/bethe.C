void bethe()
{

double logterm;
double c = 2.998*pow(10.0,8.0);
double me = 0.511; //mass of electron[MeV/c^2]
double z = 1.0; //charge of electron
double M, Z, A, I, K, rho, delta, C, Kend, step, interval, Wmax, energyloss, Range, gamma, betasq, s; 



//cout<<"入射粒子が物質中で減速するときのenergy depositとRangeをBethe-Blochの式から計算します。\n"; 
//cout<<"入射粒子の質量[MeV/c^2]を入力してください:"; cin>>M;
//cout<<"物質の原子番号:"; cin>>Z;
//cout<<"物質の原子量:"; cin>>A;
//cout<<"イオン化エネルギー[eV]:"; cin>>I;
//I = I * pow(10.0,-6.0);
//cout<<"入射粒子の入射エネルギー[MeV]:"; cin>>K;
//cout<<"物質の密度[g/cm^3]:"; cin>>rho;
//cout<<"density correction:"; cin>>delta;
//cout<<"shell correction:"; cin>>C;
//cout<<"入射粒子の最終エネルギー[MeV]:"; cin>>Kend;
//cout<<"step数:"; cin>>step;

//M=105.66;Z=82;A=208.2;I=700E-6;K=750;rho=11.35;delta=0;C=0;Kend=0;step=15; //for example (incident particle:muon, material:Pb, 750->0[MeV])
//M=me;Z=13.;A=27.;I=166E-6;K=4240;rho=2.7;delta=0;C=0;Kend=4230;step=10; //for example (incident particle:electron, material:Al, 4240->4239[MeV])
M=me;Z=2.;A=3.;I=41.8E-6;K=4240;rho=0.0095;delta=0;C=0;Kend=4230;step=10; //for example (incident particle:electron, material:Al, 4240->4239[MeV])

s = me/M;
interval = (K-Kend)/step;
cout<<"interval="<<interval<<"[MeV]\n";
Range = 0.0;

cout<<setw(12)<<"range[MeV]:"<<setw(11)<<"dE/dρx"<<setw(4)<<"/"<<setw(12)<<"Δx=ΔE(dE/dρx)^(-1)"<<endl;
for( ; Kend<=K; K -= (double)(int)interval){ //表示のときだけ少数点以下を表示しない

	gamma = K/M;
	betasq = 1 - 1/pow(gamma,2.0);

	Wmax = (2*me*betasq*pow(gamma,2.0))/(1+2*s*sqrt(1+betasq*pow(gamma,2.0))+pow(s,2.0));
	logterm = (2*me*pow(gamma,2.0)*betasq*Wmax)/pow(I,2.0);
	
	energyloss = 0.1535*(Z/A)*(pow(z,2.0)/betasq)*(log(logterm)-2*betasq-delta-2*(C/Z));
	if (energyloss<0) break;
	cout<<setw(5)<<K<<setw(1)<<"-"<<setw(5)<<K-(double)(int)interval<<setw(1)<<":"<<setw(10)<<energyloss<<setw(4)<<"/"<<setw(5)<<interval/energyloss<<endl;
	Range += interval/energyloss;

}
cout<<"Range:"<<Range<<"[g/cm^2]="<<Range/rho<<"[cm]"<<endl;

	

} 
