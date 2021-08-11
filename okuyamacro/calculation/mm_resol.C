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

void mm_resol()
{
  gROOT->Reset();
  char line[300];



    TCanvas *c1 = new TCanvas("simulation","c1");
    //c1->SetGrid();
    c1->SetFillColor(10);
  
   TH1F *h1 = new TH1F("Missing Mass","",1000,-4.,4.);
  h1->SetXTitle("MM-ML[MeV/c^2]");
  h1->SetYTitle("Counts");
   TH2F *h2 = new TH2F("MM vs Mass number","",100,1,100,1000.,-4.,4.);
  h2->SetXTitle("A (mass number)");
  h2->SetYTitle("MM [MeV/c^2]");
  
  double value;


  for (int i=0; i<1E+4; i++) {
	for(int massnum = 1; massnum<=200;massnum++){
		value = mm_rand(massnum*0.5);
		if(massnum==1)h1->Fill(value);
		h2->Fill(massnum,value);
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
  ph2->Fit("f_res_mm");
  f_res_mm->Draw("same");


   int bin1 = h1->FindFirstBinAbove(h1->GetMaximum()/2);
   int bin2 = h1->FindLastBinAbove(h1->GetMaximum()/2);
   double fwhm = h1->GetBinCenter(bin2) - h1->GetBinCenter(bin1);
	cout<<"FWHM = "<<fwhm<< " [MeV/c^2]"<<endl;

}


double mm_rand(double mn)
{
  //first position(x,y) t=theta,p=phi,last position(x_l,y_l)
   float pi = 3.141592;

	double Ee0   = 4240;//MeV
	//double pep0  = 2740;//MeV
	//double pk0   = 1200;//MeV
	double pep0  = 2740+100.*sqrt(-2.*log(gRandom->Uniform()))*cos(2.*pi*gRandom->Uniform());//MeV
	double pk0   = 1200+100.*sqrt(-2.*log(gRandom->Uniform()))*cos(2.*pi*gRandom->Uniform());//MeV
	double te0   = 6.5*pi/180.;//rad
	double tk0   = 11.5*pi/180.;//rad
	double dEe  = Ee0*0.0001/2.35*sqrt(-2.*log(gRandom->Uniform()))*cos(2.*pi*gRandom->Uniform());//MeV/c
	double dpep = pep0*0.0001/2.35*sqrt(-2.*log(gRandom->Uniform()))*cos(2.*pi*gRandom->Uniform());//MeV/c
	double dpk	= pk0*0.0002/2.35*sqrt(-2.*log(gRandom->Uniform()))*cos(2.*pi*gRandom->Uniform());//MeV/c
	double dte  = 0.0002*sqrt(-2.*log(gRandom->Uniform()))*cos(2.*pi*gRandom->Uniform());//rad
	double dtk  = 0.0006*sqrt(-2.*log(gRandom->Uniform()))*cos(2.*pi*gRandom->Uniform());//rad
	double tek0  = te0 + tk0;//rad

	double Ee = Ee0 + dEe;
	double pep= pep0 + dpep;
	double pk = pk0 + dpk;
	double te = te0 + dte;
	double tk = tk0 + dtk;
	double tek  = te + tk;//rad
	double dtek = dte + dtk;//rad
	

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
	double MT = mn * 1000.;//MeV/c^2
	double MH = MT;

	double MM0= sqrt(pow((Ee0-Eep0+MT-Ek0),2.)-(pe0*pe0+pep0*pep0+pk0*pk0-2.*pe0*pep0*cos(te0)-2.*pe0*pk0*cos(tk0)+2.*pep0*pk0*cos(tek0)));
	double MM = sqrt(pow((Ee-Eep+MT-Ek),2.)-(pe*pe+pep*pep+pk*pk-2.*pe*pep*cos(te)-2.*pe*pk*cos(tk)+2.*pep*pk*cos(tek)));
   return MM-MM0;
	//return 0.5/2./sqrt(2*log(2))*sqrt(-2.*log(gRandom->Uniform()))*cos(2.*pi*gRandom->Uniform());
 
}
