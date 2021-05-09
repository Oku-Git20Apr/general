//-- template for loading files--//
//
//	K. Okuyama, May 9, 2021
//
//taken over from compare.C


void loading(){
  
  TFile *file_G4 = new TFile("/data/41a/ELS/okuyama/g4-gene/root/kout13.root","read");
  TFile *file_simc = new TFile("/data/41a/ELS/okuyama/SIMC_jlab/SIMC/rootfiles/rad_Alupdate.root","read");

  TTree *tree_G4 = (TTree*)file_G4->Get("tree");
  TTree *tree_simc = (TTree*)file_simc->Get("SNT");



	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);



//---Physics Constant---//
 
 const double Mpi = 0.13957018;         // charged pion mass (GeV/c2)
 const double Mp = 0.938272046;         // proton       mass (GeV/c2)
 const double MK = 0.493677;            // charged Kaon mass (GeV/c2)
 const double Me = 0.510998928e-3;      // electron     mass (GeV/c2)
 const double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
 const double PI=3.14159265359;
 const double ML = 1.115683;            // Lambda       mass (GeV/c2)
 const double MS0 = 1.192642;           // Sigma Zero   mass (GeV/c2)
 const double tdc_time=0.056;//ns



//---------------------------------------//
//               Branch                  //
//---------------------------------------//

//GEANT4
  int evid;
  double Mom_mag_i;

	tree_G4->SetBranchStatus("*",0);
	tree_G4->SetBranchStatus("eventid",1);tree_G4->SetBranchAddress("eventid",&evid);
	tree_G4->SetBranchStatus("mom_i",1);tree_G4->SetBranchAddress("mom_i",&Mom_mag_i);

	
//SIMC
	float mm_simc;
	tree_simc->SetBranchStatus("*",0);
	tree_simc->SetBranchStatus("missmass",1);tree_simc->SetBranchAddress("missmass",&mm_simc);




//---------------------------------------//
//             Histrograms               //
//---------------------------------------//

//GEANT4
  TH1D* h_mom_i  = new TH1D("h_mom_i","Initial Momentum;Momentum [MeV/c];Counts/200keV",1000,1800.,2200.);
  h_mom_i->SetLineColor(kRed);

//SIMC
  double xmin = -100., xmax = 200.; int xbin = 300; // 1 MeV / bin
  TH1F* hmm_simc= new TH1F("hmm_simc","hmm_simc",xbin,xmin,xmax);
  hmm_simc->SetLineColor(kAzure);


//--------------------------------//
//             Fill               //
//--------------------------------//

//GEANT4
  int ENum_G4 = tree_G4->GetEntries(); 
cout<<"Entries(G4): "<<ENum_G4<<endl;
  for(int i=0;i<ENum_G4;i++){
	tree_G4->GetEntry(i);
	h_mom_i->Fill(Mom_mag_i);
	}

//SIMC
  int ENum_simc = tree_simc->GetEntries(); 
cout<<"Entries(SIMC): "<<ENum_simc<<endl;
  for(int i=0;i<ENum_simc;i++){
	tree_simc->GetEntry(i);
    TRandom3 ranL_simc;
	double ran = mm_simc-ML; 
	hmm_simc->Fill(ran*1000.+0.052);//-0.50);
	}

//--------------------------------//
//             Draw               //
//--------------------------------//

	TCanvas* c1 = new TCanvas("c1","c1",800,800);
hmm_simc->Draw("");
	TCanvas* c2 = new TCanvas("c2","c2",800,800);
h_mom_i->Draw("");

cout << "Well done!" << endl;
}
