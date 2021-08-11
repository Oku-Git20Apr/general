//-- Comparison (SIMC vs. Geant4)  --//
//SIMC: Bremsstrahlung(Internal + External), peaking approx.
//Geant4: Bremsstrahlung(External only?), Ionization, Multiple Scat.
//
//	K. Okuyama, May 8, 2021
//	K. Okuyama, June 1, 2021 (SIMC code is directly copied)
//
//taken over from SIMC/rootfiles/okuyamacro/radiative_Al.C

const double alpha = 1./137.;
const double PI=3.14159265359;
const double euler=0.577215665;//Euler's constant
const double Me = 0.510998928e-3;      // electron     mass (GeV/c2)

double Internal( double *x, double *par){
		//par[0]: Mom. of scattered electrons
		//par[1]: Zenith angles of scattered electrons
	//double Qsq = 2.*4.318*0.001*par[0]*(1.-cos(par[1]));//(GeV/c)^2
	double Qsq = 2.*0.001*par[0]*2.1*(1.-cos(par[1]));//(GeV/c)^2
	double aa = (alpha/PI)*log((Qsq/Me/Me)-1.);
//cout<<"a="<<aa<<endl;
	//double val = (aa/x[0])*pow((x[0]/par[0]),aa);
	double val = (aa/x[0])*pow((x[0]/par[0]),aa);
	return val+par[2];
	}


void compare(){

  
  //TFile *file_G4 = new TFile("/data/41a/ELS/okuyama/g4-gene/root/kout13.root","read");
  //TFile *file_G4 = new TFile("/data/41a/ELS/okuyama/g4-gene/root/kout17.root","read");
  TFile *file_G4 = new TFile("/data/41a/ELS/okuyama/g4-gene/root/kout10009.root","read");
  TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/rad_Alupdate.root","read");

  TTree *tree_G4 = (TTree*)file_G4->Get("tree");
  TTree *tree_simc = (TTree*)file_simc->Get("SNT");

//Output file should be loaded after loading input files.
  TFile *file_out = new TFile("ana_root/temp.root","recreate");
  //TTree *tree_out = new TTree("tree_out","tree_out");


	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);



//---Physics Constant---//
 
 const double Mpi = 0.13957018;         // charged pion mass (GeV/c2)
 const double Mp = 0.938272046;         // proton       mass (GeV/c2)
 const double MK = 0.493677;            // charged Kaon mass (GeV/c2)
 const double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
 const double ML = 1.115683;            // Lambda       mass (GeV/c2)
 const double MS0 = 1.192642;           // Sigma Zero   mass (GeV/c2)
 const double def_n_L=250.; 
 const double def_sig_L=0.003; 
 const double def_mean_L=0.0;
 const double def_n_S=70.; 
 const double def_sig_S=0.004;
 const double def_mean_S=MS0-ML;
 const double def_sig_p=0.852;
 const double def_mean_p=-8.0;
 const double def_sig_pi=0.443;
 const double def_mean_pi=3.0;
 const double def_sig_k=0.644;
 const double def_mean_k=0.0;
 const double def_acc=27.7;
 const double min_coin_c=-20.0;
 const double max_coin_c=20.0;
 const double tdc_time=0.056;//ns
 int bin_coin_c=(int)((max_coin_c-min_coin_c)/tdc_time);
 const double min_mm=-0.1;//GeV/c^2
 const double max_mm=0.2;//GeV/c^2
 int bin_mm=(max_mm-min_mm)/0.001; //Counts/2 MeV
 bin_mm=(int)bin_mm;
 //const double fit_min_mm=-0.006;
 const double fit_min_mm=-0.01;
 const double fit_max_mm=0.095;
 const int fit_bin_mm = (fit_max_mm-fit_min_mm)/0.001;
 const double fit_bin_width = (fit_max_mm-fit_min_mm)/fit_bin_mm;



//---------------------------------------//
//               Branch                  //
//---------------------------------------//

//GEANT4
  int evid;
  int pid;
  int NofHit;
  double Mom_mag_i;
  double Zenith_i;
  double Azimuth_i;//Initial State
  double Mom_mag_f[10];
  double Zenith_f[10];
  double Azimuth_f[10];//Final State
  double Pos_x;
  double Pos_y;
  double Pos_z;//Position
  double edep_test;//Energy deposit


	tree_G4->SetBranchStatus("*",0);
	tree_G4->SetBranchStatus("eventid",1);tree_G4->SetBranchAddress("eventid",&evid);
	tree_G4->SetBranchStatus("vdpid",1);tree_G4->SetBranchAddress("vdpid",&pid);
	tree_G4->SetBranchStatus("vdnhit",1);tree_G4->SetBranchAddress("vdnhit",&NofHit);
	tree_G4->SetBranchStatus("mom_i",1);tree_G4->SetBranchAddress("mom_i",&Mom_mag_i);
	tree_G4->SetBranchStatus("theta_i",1);tree_G4->SetBranchAddress("theta_i",&Zenith_i);
	tree_G4->SetBranchStatus("phi_i",1);tree_G4->SetBranchAddress("phi_i",&Azimuth_i);
	tree_G4->SetBranchStatus("x_i",1);tree_G4->SetBranchAddress("x_i",&Pos_x);
	tree_G4->SetBranchStatus("y_i",1);tree_G4->SetBranchAddress("y_i",&Pos_y);
	tree_G4->SetBranchStatus("z_i",1);tree_G4->SetBranchAddress("z_i",&Pos_z);
	tree_G4->SetBranchStatus("vdp",1);tree_G4->SetBranchAddress("vdp",Mom_mag_f);
	tree_G4->SetBranchStatus("vdtheta",1);tree_G4->SetBranchAddress("vdtheta",Zenith_f);
	tree_G4->SetBranchStatus("vdphi",1);tree_G4->SetBranchAddress("vdphi",Azimuth_f);
	//tree_G4->SetBranchStatus("edep",1);tree_G4->SetBranchAddress("edep",&edep_test);//for test
	
//SIMC
	float mm_simc;
	float brems_beam;
	float brems_scat;
	tree_simc->SetBranchStatus("*",0);
	tree_simc->SetBranchStatus("missmass",1);tree_simc->SetBranchAddress("missmass",&mm_simc);
	tree_simc->SetBranchStatus("brems_e" ,1);tree_simc->SetBranchAddress("brems_e" ,&brems_beam);
	tree_simc->SetBranchStatus("brems_ep",1);tree_simc->SetBranchAddress("brems_ep",&brems_scat);



//---------------------------------------//
//             Histrograms               //
//---------------------------------------//

//GEANT4
  TH1D* h_mom_i  = new TH1D("h_mom_i","Initial Momentum;Momentum [MeV/c];Counts/200keV",1000,1800.,2200.);
  h_mom_i->SetLineColor(kAzure);
  TH1D* h_mom_i2  = new TH1D("h_mom_i2","Initial Momentum;Momentum [MeV/c];Counts/200keV",1000,4100.,4500.);
  h_mom_i2->SetLineColor(kAzure);
  //h_mom_i->GetXaxis()->SetTitle("Momentum [MeV/c]");
  //h_mom_i->GetYaxis()->SetTitle("Counts/200keV");
  //h_mom_i->GetXaxis()->SetRangeUser(-14.0,17.);
  //double xmin = -100., xmax = 200.; int xbin = 300; // 1 MeV / bin
  TH1D* h_mom_f  = new TH1D("h_mom_f","Final Momentum;Momentum [MeV/c];Counts",1000,1800.,2200.);
  h_mom_f->SetLineColor(kRed);
  TH1D* h_mom_f2  = new TH1D("h_mom_f2","Final Momentum;Momentum [MeV/c];Counts",1000,4100.,4500.);
  h_mom_f2->SetLineColor(kRed);
  TH1D* h_mom_dif_G4  = new TH1D("h_mom_dif_G4","Momentum difference;Momentum [MeV/c];Counts",100000,0.,2200.);
  h_mom_dif_G4->SetLineColor(kGreen);
  TH1D* h_mom_dif2_G4  = new TH1D("h_mom_dif2_G4","Momentum difference;Momentum [MeV/c];Counts",100000,0.,4500.);
  h_mom_dif2_G4->SetLineColor(kGreen);
  TH1D* h_mom_dif3_G4  = new TH1D("h_mom_dif3_G4","Momentum difference;Momentum [MeV/c];Counts",100000,0.,4500.);
  h_mom_dif3_G4->SetLineColor(kRed);

//SIMC
  double xmin = -100., xmax = 200.; int xbin = 300; // 1 MeV / bin
  TH1F* hmm_simc= new TH1F("hmm_simc","hmm_simc",xbin,xmin,xmax);
  TH1D* h_mom_dif_simc  = new TH1D("h_mom_dif_simc","Momentum difference;Momentum [MeV/c];Counts",100000,0.,2200.);
  h_mom_dif_simc->SetLineColor(kAzure);
  TH1D* h_mom_dif2_simc  = new TH1D("h_mom_dif2_simc","Momentum difference;Momentum [MeV/c];Counts",100000,0.,4500.);
  h_mom_dif2_simc->SetLineColor(kAzure);

//--------------------------------//
//      SIMC Brems. Calc.         //
//--------------------------------//

	gRandom->SetSeed(0);
	for(int i=0;i<1200000;i++){
	double theta_temp = 0.;
	while(1){// Solid Angle = 6msr
		theta_temp = acos(1.-gRandom->Uniform()*(1.-cos(0.5)));
		if(abs(theta_temp-13.2*PI/180.)<0.044)break;
	}
	double theta_scat = theta_temp;
	double Escat = 2.*94.5*gRandom->Uniform()+2005.5;
//cout<<"theta_scat="<<theta_scat*180./PI<<endl;
//cout<<"Escat="<<Escat<<endl;
	double Xrad = 89.;//mm; Al radiation length
	double thick = 0.357/sin(theta_scat)/Xrad;//in radiation length 
//cout<<"thick="<<thick<<endl;
	//double Z=1.;
	double Z=13.;//Al
	//double L1 = 5.31;//if Z=1
	double L1 = log(184.15)-log(Z)/3.;
	//double L2 = 6.144;//if Z=1
	double L2 = log(1194.)-2.*log(Z)/3.;
	double etatzai = (12.0+(Z+1.)/(Z*L1+L2))/9.0;
	double bt = etatzai * thick;
	double alpi = alpha/PI;
	double lambda = alpi*(2.*log(2.*Escat*0.001/Me) -1.+log((1-cos(theta_scat))/2.));
	double ggg = lambda + bt;
	//double ggg = bt;//external only

	double Egmin = 0.;
	double Egmax = Escat;
	double ymin = pow(Egmin/Egmax,ggg);
	double x = pow(gRandom->Uniform()*(1.-ymin)+ymin,1./ggg)*Egmax;
	//if(Escat-x<2005.5){i--;continue;}
	h_mom_dif3_simc->Fill(x);
	}

//--------------------------------//
//             Fill               //
//--------------------------------//

//GEANT4
  int ENum_G4 = tree_G4->GetEntries(); 
cout<<"Entries(G4): "<<ENum_G4<<endl;

  for(int i=0;i<ENum_G4;i++){
	if(i%10000==0)cout<<i<<"/"<<ENum_G4<<endl;
	tree_G4->GetEntry(i);
	h_mom_i->Fill(Mom_mag_i);
	double temp = Mom_mag_i-Mom_mag_f[0];
	h_mom_dif_G4->Fill(Mom_mag_i-Mom_mag_f[0]);
	h_mom_dif2_G4->Fill(Mom_mag_i-Mom_mag_f[0]);
	}

	TF1* f_internal = new TF1("f_internal",Internal,0.,200.,3);
if(ENum_G4>2000){ENum_G4=2000;cout<<"loading only "<<ENum_G4<<" events..."<<endl;}
  for(int i=0;i<ENum_G4;i++){
	if(i%10==0)cout<<i<<"/"<<ENum_G4<<endl;
	tree_G4->GetEntry(i);
	if(Mom_mag_f[0]>0.){
	f_internal->SetParameter(0,Mom_mag_f[0]);
	f_internal->SetParameter(1,13.2*PI/180.-Zenith_i);
	f_internal->SetParameter(2,Mom_mag_i-Mom_mag_f[0]);
	}else{
	f_internal->SetParameter(0,4318.);
	f_internal->SetParameter(1,13.2*PI/180.);
	f_internal->SetParameter(2,0.);
	}
	//f_internal->SetParameter(0,2100.);
	//f_internal->SetParameter(1,0.23);
	//h_mom_dif_G4->FillRandom("f_internal",10000);
	//h_mom_dif2_G4->FillRandom("f_internal",10000);
	h_mom_dif3_G4->FillRandom("f_internal",100);
	}

//SIMC
  int ENum_simc = tree_simc->GetEntries(); 
cout<<"Entries(SIMC): "<<ENum_simc<<endl;
  for(int i=0;i<ENum_simc;i++){
	tree_simc->GetEntry(i);
    TRandom3 ranL_simc;
	//double ran = ranL_simc.Gaus((mm_simc-ML),0.001);
	double ran = mm_simc-ML; 
	//hmm_simc->Fill(ran*1000.+0.052);//-0.50);
	h_mom_dif_simc->Fill(brems_scat);
	h_mom_dif2_simc->Fill(brems_beam);
	}

//--------------------------------//
//             Draw               //
//--------------------------------//

//TH1F* h_eloss_simc->(TH1F*)file_simc->Get("eloss");
double Nsimc_scat = h_mom_dif_simc->Integral(h_mom_dif_simc->FindBin(0.),h_mom_dif_simc->FindBin(200.));
double Nsimc_beam = h_mom_dif2_simc->Integral(h_mom_dif2_simc->FindBin(0.),h_mom_dif2_simc->FindBin(200.));
double NG4_scat = h_mom_dif_G4->Integral(h_mom_dif_G4->FindBin(0.),h_mom_dif_G4->FindBin(200.));
double NG4_beam = h_mom_dif2_G4->Integral(h_mom_dif2_G4->FindBin(0.),h_mom_dif2_G4->FindBin(200.));
double NG4_beam3 = h_mom_dif3_G4->Integral(h_mom_dif3_G4->FindBin(0.),h_mom_dif3_G4->FindBin(200.));
cout<<"Nsimc_scat="<<Nsimc_scat<<endl;
cout<<"Nsimc_beam="<<Nsimc_beam<<endl;
cout<<"NG4_scat="<<NG4_scat<<endl;
cout<<"NG4_beam="<<NG4_beam<<endl;
cout<<"NG4_beam3="<<NG4_beam3<<endl;
//h_mom_dif_simc->Scale(NG4_scat/Nsimc_scat);
//h_mom_dif2_simc->Scale(NG4_beam/Nsimc_beam);
//h_mom_dif3_G4->Scale(NG4_beam/NG4_beam3);
h_mom_dif_G4->Scale(6.);//120000/200000
h_mom_dif2_G4->Scale(6.);
h_mom_dif3_G4->Scale(6.);

	TCanvas* c1 = new TCanvas("c1","if scat",800,800);
c1->cd()->DrawFrame(0.,1.,200.,1000.,"Comparison;Momentum [GeV/c];Counts");
h_mom_dif_simc->Draw("same");
h_mom_dif_G4->Draw("same");
	TCanvas* c2 = new TCanvas("c2","if beam",800,800);
c2->cd()->DrawFrame(0.,1.,200.,1000.,"Comparison;Momentum [GeV/c];Counts");
h_mom_dif2_simc->Draw("same");
h_mom_dif2_G4->Draw("same");
h_mom_dif3_G4->Draw("same");

 file_out->Write();

cout << "Well done!" << endl;
}
