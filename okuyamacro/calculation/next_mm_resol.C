double mm_rand(double mn);
//---Physics Constant---//
 
 const double Mpi = 0.13957018*1000.;         // charged pion mass (MeV/c2)
 const double MP = 0.938272046*1000.;         // proton       mass (MeV/c2)
 const double MK = 0.493677*1000.;            // charged Kaon mass (MeV/c2)
 const double Me = 0.510998928;      // electron     mass (MeV/c2)
 const double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
 const double ML = 1.115683*1000.;            // Lambda       mass (MeV/c2)
 const double MS0 = 1.192642*1000.;           // Sigma Zero   mass (MeV/c2)

 double alpha = 1./137.;

void next_mm_resol()
{
  gROOT->Reset();
  char line[300];



    TCanvas *c1 = new TCanvas("simulation","c1");
    //c1->SetGrid();
    c1->SetFillColor(10);
  
   TH1F *h1 = new TH1F("h1","",200,-3.,3.);
   //TH1F *h1 = new TH1F("Missing Mass","",1000,ML-400.,ML+400.);
   h1->SetXTitle("MM-Mtrue[MeV/c^{2}]");
   h1->SetYTitle("Counts");
   TH1F *hmm_adep = new TH1F("hmm_adep","",200,0.,20.);
   hmm_adep->SetXTitle("A (mass number of the target)");
   hmm_adep->SetYTitle("MM resolution (FWHM) [MeV/c^{2}]");
	TH1F *hist[200];
	TF1 *func[200];
  for(int i=0;i<200;i++){
		hist[i] = new TH1F(Form("hist[%d]",i),"",200,-4.,4.);
		hist[i]->SetXTitle("MM-ML[MeV/c^2]");
  		hist[i]->SetYTitle("Counts");
		func[i] = new TF1(Form("func[%d]",i),"gausn",-4.,4.);
  }
   TH2F *h2 = new TH2F("h2","",200,0,20,1000.,-4.,4.);
  h2->SetXTitle("A (mass number of the target)");
  h2->SetYTitle("MM deviation [MeV/c^{2}]");
   TH2F *h2_pepk = new TH2F("h2_pepk","",200,2740.-124.,2740.+124.,200.,1200.-150.,1200.+150.);
  
  double value;


  for (int i=0; i<1E+4; i++) {
	for(int massnum = 1; massnum<=200;massnum++){
		value = mm_rand(massnum*0.1);
		if(massnum==10){//elementary
			h1->Fill(value);
	double pep0  = 2740.*(1.+0.045*(1.-2.*gRandom->Uniform()));//MeV
	double pk0   = 1200.*(1.+0.125*(1.-2.*gRandom->Uniform()));//MeV
			h2_pepk->Fill(pep0,pk0);
		}
		h2->Fill(massnum*0.1,value);
		hist[massnum-1]->Fill(value);
	}
  }

  

  h1->SetStats( 0 );
  h1->SetLineStyle(2);
  h1->Draw();
  h1->Fit("gausn");

  TCanvas* c2 = new TCanvas("c2","c2",800,800);
  h2->SetStats( 0 );
  h2->Draw("colz");
  TCanvas* c3 = new TCanvas("c3","c3",800,800);
    float pi = 3.141592;
	double Ee   = 4240;//MeV
	double pep  = 2740;//MeV
	double pk   = 1200;//MeV
	double te   = 6.5*pi/180.;//rad
	double tk   = 11.5*pi/180.;//rad
	double tek  = te + tk;
	double pe   = sqrt(Ee*Ee-Me*Me);
	double Eep  = sqrt(pep*pep+Me*Me);
	double Ek   = sqrt(pk*pk+MK*MK);
	//double MM = sqrt(pow((Ee-Eep+MT-Ek),2.)-(pe*pe+pep*pep+pk*pk-2.*pe*pep*cos(te)-2.*pe*pk*cos(tk)+2.*pep*pk*cos(tek)));
  TF1 *f_res_mm = new TF1("f_res_mm","sqrt(pow([0]+x*1000.,2.)-[1])-x*1000.",1.,100.);
	f_res_mm->SetParameter(0,Ee-Eep-Ek);
	f_res_mm->SetParameter(1,pe*pe+pep*pep+pk*pk-2.*pe*pep*cos(te)-2.*pe*pk*cos(tk)+2.*pep*pk*cos(tek));
  //TF1 *f_res_mm = new TF1("f_res_mm","sqrt(pow([0]+x,2.)-[1])-x",1.,100.);
  TProfile *ph2 = h2->ProfileX();
  ph2->Draw();
  //ph2->Fit("f_res_mm");
  f_res_mm->Draw("same");


   int bin1 = h1->FindFirstBinAbove(h1->GetMaximum()/2);
   int bin2 = h1->FindLastBinAbove(h1->GetMaximum()/2);
   double fwhm = h1->GetBinCenter(bin2) - h1->GetBinCenter(bin1);
	cout<<"FWHM = "<<fwhm<< " [MeV/c^2]"<<endl;

	for(int i=0;i<200;i++){
	   int bin1 = hist[i]->FindFirstBinAbove(hist[i]->GetMaximum()/2);
	   int bin2 = hist[i]->FindLastBinAbove(hist[i]->GetMaximum()/2);
	   double fwhm = hist[i]->GetBinCenter(bin2) - hist[i]->GetBinCenter(bin1);
	   hist[i]->Fit(func[i],"Rq0");
		//cout<<"A="<<(i+1.)*0.5<<": FWHM = "<<fwhm<< " [MeV/c^2]"<<endl;
	   //hmm_adep->SetBinContent(i+1,fwhm);
	   hmm_adep->SetBinContent(i+1,2.*sqrt(2.*log(2))*func[i]->GetParameter(2));
	   hmm_adep->SetBinError(i+1,2.*sqrt(2.*log(2))*func[i]->GetParError(2));
	}
  TCanvas *c4 = new TCanvas("c4","c4",800,800);
  hmm_adep->SetStats( 0 );
  hmm_adep->SetLineColor( kAzure );
  hmm_adep->Draw("");
  //c4->Print("nnL_mm_resol.pdf");
  TCanvas* c5 = new TCanvas("c5","c5",1300,600);
  c5->Divide(2,1);
  c5->cd(1);
  h2->Draw("colz");
  c5->cd(2);
  hmm_adep->Draw("");
  //c5->Print("next_mm_resol.pdf");
  //c1->Print("next_mm_resol_A1.pdf");
  TCanvas* c6 = new TCanvas("c6","c6",1300,600);
  h2_pepk->Draw("colz");
}


double mm_rand(double mn)
{
  //first position(x,y) t=theta,p=phi,last position(x_l,y_l)
   float pi = 3.141592;
	double massunit = 931.49410242;//MeV

	double Ee0   = 4240.;//MeV
	double pep0  = 2740.*(1.+0.045*(1.-2.*gRandom->Uniform()));//MeV
	double pk0   = 1200.*(1.+0.125*(1.-2.*gRandom->Uniform()));//MeV
	//double te0   = 13.2*pi/180.+2.*pi/180.*(1.-0.5*gRandom->Uniform());//rad
	//double tk0   = 13.2*pi/180.+2.*pi/180.*(1.-0.5*gRandom->Uniform());//rad
	//
//isomeric (theta_ee'=6.5\pm2.2 deg, theta_ek=11.5\pm2.6 deg)
	double te0 = 0.;
	double tk0 = 0.;
	while(te0*180./pi<4.3||tk0*180./pi<8.9){
	te0   = acos(1.-gRandom->Uniform()*(1.-cos(8.7*pi/180.)));//rad
	tk0   = acos(1.-gRandom->Uniform()*(1.-cos(14.1*pi/180.)));//rad
	}
	double dEe  = Ee0*(0.0001/2.35)*sqrt(-2.*log(gRandom->Uniform()))*cos(2.*pi*gRandom->Uniform());//MeV/c
	double dpep = pep0*(0.0001/2.35)*sqrt(-2.*log(gRandom->Uniform()))*cos(2.*pi*gRandom->Uniform());//MeV/c
	double dpk	= pk0*(0.0002/2.35)*sqrt(-2.*log(gRandom->Uniform()))*cos(2.*pi*gRandom->Uniform());//MeV/c
	double dte  = (0.0006/2.35)*sqrt(-2.*log(gRandom->Uniform()))*cos(2.*pi*gRandom->Uniform());//rad
	double dtk  = (0.0002/2.35)*sqrt(-2.*log(gRandom->Uniform()))*cos(2.*pi*gRandom->Uniform());//rad
	double tek0  = te0 + tk0;//rad
 //cout<<"Ee0="<<Ee0<<", pep0="<<pep0<<", pk0="<<pk0<<", te0="<<te0<<", tk0="<<tk0<<", tek0="<<tek0<<endl;
	//dEe = 0.;
	//dpep= 0.;
	//dpk = 0.;
	//dte = 0.;
	//dtk = 0.;

	double Ee = Ee0 + dEe;
	double pep= pep0 + dpep;
	double pk = pk0 + dpk;
	double te = te0 + dte;
	double tk = tk0 + dtk;
	double tek  = te + tk;//rad
	double dtek = dte + dtk;//rad
 //cout<<"dEe="<<dEe<<", dpep="<<dpep<<", dpk="<<dpk<<", dte="<<dte<<", dtk="<<dtk<<", dtek="<<dtek<<endl;
	

	double pe0  = sqrt(Ee0*Ee0-Me*Me);
	double Eep0 = sqrt(pep0*pep0+Me*Me);
	double Ek0  = sqrt(pk0*pk0+MK*MK);
	double pe   = sqrt(Ee*Ee-Me*Me);
	double Eep  = sqrt(pep*pep+Me*Me);
	double Ek   = sqrt(pk*pk+MK*MK);
	//double MT = MP;
	//double MH = ML-30.;
	//double MT = 37.22*1000.;// 40Ca GeV/c^2
	//double MH = 37.41*1000.+55;// 39K+L
	//double MT = 25.133*1000.;// 27Al
	//double MH = 25.319*1000.+14.;// 26Mg+L
	double MT = mn * massunit;//MeV/c^2
	double MH = MT;

	double MM0= sqrt(pow((Ee0-Eep0+MT-Ek0),2.)-(pe0*pe0+pep0*pep0+pk0*pk0-2.*pe0*pep0*cos(te0)-2.*pe0*pk0*cos(tk0)+2.*pep0*pk0*cos(tek0)));
	double MM = sqrt(pow((Ee-Eep+MT-Ek),2.)-(pe*pe+pep*pep+pk*pk-2.*pe*pep*cos(te)-2.*pe*pk*cos(tk)+2.*pep*pk*cos(tek)));
   return MM - MM0;
	//return 0.5/2.35*sqrt(-2.*log(gRandom->Uniform()))*cos(2.*pi*gRandom->Uniform());
 
}
