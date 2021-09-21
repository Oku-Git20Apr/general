double fmm_total(double *x, double *par){
  double val = par[0] * TMath::Gaus(x[0],par[1],par[2],kTRUE);//s-orbit
  val += par[3] * TMath::Gaus(x[0],par[4],par[5], kTRUE);//p-orbit
  val += par[6] * TMath::Gaus(x[0],par[7],par[8], kTRUE);//p-orbit
  val += par[9] * TMath::Gaus(x[0],par[10],par[11],kTRUE);//p-orbit
  val += par[12];//Accidentals
  if(x[0]>0)val += par[13] + par[14]*x[0] + par[15]*x[0];//Accidentals

  return val;
}
double fmm_core_total(double *x, double *par){//with core excited states
  double val = 0.;
	for(int i=0;i<19;i++){
		val += par[3*i] * TMath::Gaus(x[0],par[3*i+1],par[3*i+2],kTRUE);
	}
  val += par[60];//Accidentals
  if(x[0]>0)val += par[61] + par[62]*x[0] + par[63]*x[0];//Accidentals

  return val;
}


void fit_mm(){
	string pdfname = "fitting.pdf";
cout << "Output pdf file name is " << pdfname << endl;
  
  //TFile *file = new TFile("250keV_noscale.root","read");
  //TFile *file = new TFile("250keV_10x.root","read");
  //TFile *file = new TFile("250keV_50x.root","read");
  //TFile *file = new TFile("250keV_p532_noscale.root","read");
  TFile *file = new TFile("250keV_single_p532_noscale.root","read");

    
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

  double CSp_tot = 125.*5./6.;//Lambda in p-orbit
  double CSpval1 = CSp_tot/3.;
  double CSpval2 = CSp_tot/3.;
  double CSpval3 = CSp_tot/3.;
  //CSpval1 = CSp_tot;
  //CSpval2 = 0.;
  //CSpval3 = 0.;
  double Expval1 = -6.3;
  double Expval2 = -3.3;
  double Expval3 = -0.8;
  //double Ex[1]={-16.3};
  //double Cs[1]={ 79.0};

//without core excited states
//  double Ex[4] = { -16.3,	   Expval1, Expval2, Expval3};
//  double CS[4] = {  79.0*5./6.,CSpval1, CSpval2, CSpval3,};
//with core excited states
  double Ex[19] = { -16.3, -15.2, -13.5, -12.0, Expval1, Expval1*7.5/8.3, Expval1*5.5/8.3, Expval1*4.0/8.3, Expval1*2.8/8.3, Expval2, (Expval2-Expval1)+Expval1*7.5/8.3, (Expval2-Expval1)+Expval1*5.5/8.3, (Expval2-Expval1)+Expval1*4.0/8.3, (Expval2-Expval1)+Expval1*2.8/8.3, Expval3, (Expval3-Expval1)+Expval1*7.5/8.3, (Expval3-Expval1)+Expval1*5.5/8.3, (Expval3-Expval1)+Expval1*4.0/8.3, (Expval3-Expval1)+Expval1*2.8/8.3};
  double CS[19] = {  79.0*5./6.,  26.3,  13.2,  35.1, CSpval1, CSpval1*100/320, CSpval1*60./320, CSpval1*100/320, CSpval1*60./320, CSpval2, CSpval2*100/320, CSpval2*60./320, CSpval2*100/320, CSpval2*60./320, CSpval3, CSpval3*100/320, CSpval3*60./320, CSpval3*100/320, CSpval3*60./320};


//---Physics Constant---//
 
 const double Mpi = 0.13957018;         // charged pion mass (GeV/c2)
 const double Mp = 0.938272046;         // proton       mass (GeV/c2)
 const double MK = 0.493677;            // charged Kaon mass (GeV/c2)
 const double Me = 0.510998928e-3;      // electron     mass (GeV/c2)
 const double LightVelocity = 0.299792458;          // speed of light in vacuum (m/ns)
 const double PI=3.14159265359;
 const double ML = 1.115683;            // Lambda       mass (GeV/c2)
 const double MS0 = 1.192642;           // Sigma Zero   mass (GeV/c2)

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->Divide(1,1,1E-4,1E-4);
  c1->cd(1);
  //c1->cd(1)->DrawFrame(-20.,0.,20.,1E+4);
    //h1_mm->Draw("E");
	TH1F* h1_mm=(TH1F*)file->Get("h1_mm");
	TH1F* h1_mm2=(TH1F*)file->Get("h1_mm2");
	TH1F* h1_mm_a=(TH1F*)file->Get("h1_mm_a");
	TH1F* h1_mm_a2=(TH1F*)file->Get("h1_mm_a2");
	TH1F* h1_mm_q=(TH1F*)file->Get("h1_mm_q");
    h1_mm->Draw("E");
    //h1_mm->GetXaxis()->SetRangeUser(-25,10);
    //h1_mm->GetXaxis()->SetRangeUser(-30,10);
    //h1_mm_p->Draw("same");
    //h1_mm_q->Draw("same");
    //h1_mm_a->Scale(1./100.);
    h1_mm_a->Draw("same");
  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  c2->Divide(1,2,1E-4,1E-4);
  c2->cd(1);
	TF1* f_acc  = new TF1("f_acc" ,"pol0",-20.,10.);
	TF1* f_qf  = new TF1("f_qf" ,"pol2",-20,10.);
    h1_mm_a2->Fit(f_acc);
    h1_mm_q->Fit(f_qf,"","",0.,10.);
  c2->cd(2);
    h1_mm2->Draw("E");
    h1_mm_a->Draw("same");
	TF1* fs  = new TF1("fs" ,"gausn+pol0(3)",-20.,10.);
	TF1* fp1 = new TF1("fp1","gausn+pol0(3)",-20.,10.);
	TF1* fp2 = new TF1("fp2","gausn+pol0(3)",-20.,10.);
	TF1* fp3 = new TF1("fp3","gausn+pol0(3)",-20.,10.);
	fs ->SetParameter(0,100);
	fs ->SetParameter(1,-16.3);
	fs ->SetParameter(2,0.81);
	fs ->FixParameter(3,f_acc->GetParameter(0));
	fp1->SetParameter(0,50);
	fp1->SetParameter(1,-6.1);
	fp1->SetParameter(2,1.0);
	fp1->FixParameter(3,f_acc->GetParameter(0));
	fp2->SetParameter(0,50);
	fp2->SetParameter(1,-3.1);
	fp2->SetParameter(2,1.0);
	fp2->FixParameter(3,f_acc->GetParameter(0));
	fp3->SetParameter(0,50);
	fp3->SetParameter(1,-0.8);
	fp3->SetParameter(2,1.0);
	fp3->FixParameter(3,f_acc->GetParameter(0));
	h1_mm2->Fit(fs,"","",-20.,-15.);//16.3
	h1_mm2->Fit(fp1,"","",-8.0,-5.5);//6.3
	h1_mm2->Fit(fp2,"","",-5.0,-2.5);//3.3
	h1_mm2->Fit(fp3,"","",-1.5,0.);//0.8
    h1_mm2->Draw("E");
    h1_mm_a->Draw("same");
	fs->Draw("same");
	fp1->Draw("same");
	fp2->Draw("same");
	fp3->Draw("same");
  TCanvas *c4 = new TCanvas("c4","c4",800,600);
  c4->Divide(1,1,1E-4,1E-4);
	TF1* fs_result  = new TF1( "fs_result","gausn",-20.,10.);
	TF1* fp1_result = new TF1("fp1_result","gausn",-20.,10.);
	TF1* fp2_result = new TF1("fp2_result","gausn",-20.,10.);
	TF1* fp3_result = new TF1("fp3_result","gausn",-20.,10.);
	TF1* ftot_result= new TF1("ftot_result",fmm_total,-20.,10.,16);
	TF1* ftot_core_result= new TF1("ftot_core_result",fmm_core_total,-20.,10.,64);
	fs_result->SetNpx(2000);
	fp1_result->SetNpx(2000);
	fp2_result->SetNpx(2000);
	fp3_result->SetNpx(2000);
	ftot_result->SetNpx(2000);
	ftot_core_result->SetNpx(2000);
	fs_result->SetLineColor(kAzure);
	fp1_result->SetLineColor(kAzure);
	fp2_result->SetLineColor(kAzure);
	fp3_result->SetLineColor(kAzure);
	fs_result->SetFillColor(kAzure);
	fp1_result->SetFillColor(kAzure);
	fp2_result->SetFillColor(kAzure);
	fp3_result->SetFillColor(kAzure);
	fs_result->SetFillStyle(3004);
	fp1_result->SetFillStyle(3004);
	fp2_result->SetFillStyle(3004);
	fp3_result->SetFillStyle(3004);
	ftot_result->SetLineColor(kRed);
	ftot_core_result->SetLineColor(kRed);
	ftot_result ->SetParameter(0,fs->GetParameter(0));
	ftot_result ->SetParameter(1,fs->GetParameter(1));
	ftot_result ->SetParameter(2,fs->GetParameter(2));
	ftot_result ->SetParameter(3,fp1->GetParameter(0));
	ftot_result ->SetParameter(4,fp1->GetParameter(1));
	ftot_result ->SetParameter(5,fp1->GetParameter(2));
	ftot_result ->SetParameter(6,fp2->GetParameter(0));
	ftot_result ->SetParameter(7,fp2->GetParameter(1));
	ftot_result ->SetParameter(8,fp2->GetParameter(2));
	ftot_result ->SetParameter(9,fp3->GetParameter(0));
	ftot_result ->SetParameter(10,fp3->GetParameter(1));
	ftot_result ->SetParameter(11,fp3->GetParameter(2));
	ftot_result ->SetParameter(12,f_acc->GetParameter(0));
	ftot_result ->SetParameter(13,f_qf->GetParameter(0));
	ftot_result ->SetParameter(14,f_qf->GetParameter(1));
	ftot_result ->SetParameter(15,f_qf->GetParameter(2));

	//double factor_s = fs->GetParameter(0)/CS[0];
	//double factor_p1 = fp1->GetParameter(0)/CS[4];
	//double factor_p2 = fp2->GetParameter(0)/CS[9];
	//double factor_p3 = fp3->GetParameter(0)/CS[14];
	double factor_s = 32.*7./0.25/CS[0];
	double factor_p1 = 51.*7.*5./10./0.25/CS[4];
	double factor_p2 = 51.*7.*3./10./0.25/CS[9];
	double factor_p3 = 51.*7.*2./10./0.25/CS[14];
	for(int i=0;i<19;i++){
		if(i<4){ftot_core_result ->SetParameter(3*i,CS[i]*factor_s);
		}else if(i<9){ftot_core_result ->SetParameter(3*i,CS[i]*factor_p1);
		}else if(i<14){ftot_core_result ->SetParameter(3*i,CS[i]*factor_p1);
		}else{ftot_core_result ->SetParameter(3*i,CS[i]*factor_p1);}
		ftot_core_result ->SetParLimits(3*i,0.,fs->GetParameter(0));
		ftot_core_result ->FixParameter(3*i+1,Ex[i]);
		//ftot_core_result ->FixParameter(3*i+2,fs->GetParameter(2));
		ftot_core_result ->FixParameter(3*i+2,0.8/2.35);//resolution
	}
		//ftot_core_result ->FixParameter(0,32.*7./0.25/2.);
		//ftot_core_result ->FixParameter(12,51.*7.*5./10./0.25/2.);
		//ftot_core_result ->FixParameter(27,51.*7.*3./10./0.25/2.);
		//ftot_core_result ->FixParameter(42,51.*7.*2./10./0.25/2.);
		ftot_core_result ->SetParameter(60,f_acc->GetParameter(0));
		//ftot_core_result ->SetParameter(61,f_qf->GetParameter(0));
		ftot_core_result ->FixParameter(61,0.);
		ftot_core_result ->SetParameter(62,f_qf->GetParameter(1));
		ftot_core_result ->SetParameter(63,f_qf->GetParameter(2));
	//h1_mm3->Fit(ftot_result,"","",-20.,10.);
	//cout<<"reduced chi-square = "<<ftot_result->GetChisquare()<<"/"<<ftot_result->GetNDF()<<" = "<<ftot_result->GetChisquare()/ftot_result->GetNDF()<<endl;
	//fs_result ->SetParameter(0,ftot_result->GetParameter(0));
	//fs_result ->SetParameter(1,ftot_result->GetParameter(1));
	//fs_result ->SetParameter(2,ftot_result->GetParameter(2));
	//fp1_result->SetParameter(0,ftot_result->GetParameter(3));
	//fp1_result->SetParameter(1,ftot_result->GetParameter(4));
	//fp1_result->SetParameter(2,ftot_result->GetParameter(5));
	//fp2_result->SetParameter(0,ftot_result->GetParameter(6));
	//fp2_result->SetParameter(1,ftot_result->GetParameter(7));
	//fp2_result->SetParameter(2,ftot_result->GetParameter(8));
	//fp3_result->SetParameter(0,ftot_result->GetParameter(9));
	//fp3_result->SetParameter(1,ftot_result->GetParameter(10));
	//fp3_result->SetParameter(2,ftot_result->GetParameter(11));
	TH1F* h1_mm3=(TH1F*)file->Get("h1_mm3");
	h1_mm3->Fit(ftot_core_result,"","",-20.,10.);//with core excited
	cout<<"reduced chi-square = "<<ftot_core_result->GetChisquare()<<"/"<<ftot_core_result->GetNDF()<<" = "<<ftot_core_result->GetChisquare()/ftot_core_result->GetNDF()<<endl;
	fs_result ->SetParameter(0,ftot_core_result->GetParameter(0));
	fs_result ->SetParameter(1,ftot_core_result->GetParameter(1));
	fs_result ->SetParameter(2,ftot_core_result->GetParameter(2));
	fp1_result->SetParameter(0,ftot_core_result->GetParameter(12));
	fp1_result->SetParameter(1,ftot_core_result->GetParameter(13));
	fp1_result->SetParameter(2,ftot_core_result->GetParameter(14));
	fp2_result->SetParameter(0,ftot_core_result->GetParameter(27));
	fp2_result->SetParameter(1,ftot_core_result->GetParameter(28));
	fp2_result->SetParameter(2,ftot_core_result->GetParameter(29));
	fp3_result->SetParameter(0,ftot_core_result->GetParameter(42));
	fp3_result->SetParameter(1,ftot_core_result->GetParameter(43));
	fp3_result->SetParameter(2,ftot_core_result->GetParameter(44));
  h1_mm3->Draw("E");
    h1_mm_a->Draw("same");
	//fs_result->Draw("same");
	//fp1_result->Draw("same");
	//fp2_result->Draw("same");
	//fp3_result->Draw("same");

cout << "Well done!" << endl;
}//fit
