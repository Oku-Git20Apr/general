	double CalcF1TDC(double tdc, double toff) {return (tdc-toff)*0.05623;}
	void temp_ana(){

	const double c = 0.299792458; // speed of light (m/ns)
	const double MK = 0.493677;   // Kaon Mass (GeV)
	const double Me = 0.000511;   // Electron Mass (GeV)
	const double Mp = 0.938272;
	const double ML = 1.115683;
	const double TDCtoT = 0.05623;  //  (ns/ch)
	const double PI = atan(1)*4.;

	double L_tr_n;
	double R_tr_n;
	double L_s2_trpad[10000];
	double R_s2_trpad[10000];
	double L_tr_p[10000];
	double R_tr_p[10000];
	double L_tr_pathl[10000];
	double R_tr_pathl[10000];
	double L_s2_trpath[10000];
	double R_s2_trpath[10000];
	double L_F1Fhit[10000];
	double R_F1Fhit[10000];

	TChain *tr = new TChain("T");
	tr->SetBranchStatus("*",0);
	tr->Add("tritium_111157_1.root");
	tr->SetBranchAddress("L.tr.n" ,&L_tr_n);
	tr->SetBranchAddress("R.tr.n" ,&R_tr_n);
	tr->SetBranchAddress("L.s2.trpad" ,L_s2_trpad);
	tr->SetBranchAddress("L.tr.p" ,L_tr_p);
	tr->SetBranchAddress("R.tr.p" ,R_tr_p);
	tr->SetBranchAddress("L.tr.pathl" ,L_tr_pathl);
	tr->SetBranchAddress("R.tr.pathl" ,R_tr_pathl);
	tr->SetBranchAddress("L.s2.trpath" ,L_s2_trpath);
	tr->SetBranchAddress("R.s2.trpath" ,R_s2_trpath);
	tr->SetBranchAddress("LTDC.F1FirstHit" ,L_F1Fhit);
	tr->SetBranchAddress("RTDC.F1FirstHit" ,R_F1Fhit);
    double ENum = tr->GetEntries();
	cout <<" ENum = " << ENum << endl;
  TH1D * h_ct;
  h_ct       = new TH1D("h_ct" ,"h_ct" ,400, -100., 100.); 
  SetTH1(h_ct,"Coincidence Time" ,"Cointime (ns)" ,"Counts");
	for(int n=0;n<ENum;n++){
		tr->GetEntry(n);
      int NLtr = (int)L_tr_n; 
	  int NRtr = (int)R_tr_n;



////////////F1TDC////////////////
	double L_s2l_t[16];
	double L_s2l_toff[16];
	double L_s2r_toff[16];
	double L_s2r_t[16];
	double L_s2_t[16];	
	double R_s2l_t[16];
	double R_s2l_toff[16];
	double R_s2r_toff[16];
	double R_s2r_t[16];
	double R_s2_t[16];	
	for(int i=0;i<16;i++){
		L_s2l_toff[i] = 0.;
		L_s2r_toff[i] = 0.;
		R_s2l_toff[i] = 0.;
		R_s2r_toff[i] = 0.;
	}
	  for(int i=0;i<16;i++){
		L_s2l_t[i] = CalcF1TDC( L_F1Fhit[i] - L_F1Fhit[30], L_s2l_toff[i] );
		L_s2r_t[i] = CalcF1TDC( L_F1Fhit[i+48] - L_F1Fhit[37]     , L_s2r_toff[i] );
		if( L_F1Fhit[i]>0 && L_F1Fhit[i+48]>0 ){
		L_s2_t[i] = (L_s2l_t[i] + L_s2r_t[i]) / 2.;
		} else{ L_s2_t[i] = -999999; }
	}
	  for(int i=0;i<16;i++){
		R_s2l_t[i] = CalcF1TDC( R_F1Fhit[i] - R_F1Fhit[9], R_s2l_toff[i] );
		R_s2r_t[i] = CalcF1TDC( R_F1Fhit[i+48] - R_F1Fhit[46]     , R_s2r_toff[i] );
		if( R_F1Fhit[i]>0 && R_F1Fhit[i+48]>0 ){
		R_s2_t[i] = (R_s2l_t[i] + R_s2r_t[i]) / 2.;
		} else{ R_s2_t[i] = -999999; }
	}


		for(int lt=0;lt<NLtr;lt++){
//        if( tr->L_tr_chi2[lt]<0.01 ) L_Tr = true;
//        if( tr->L_tr_th[lt]<0.17*tr->L_tr_x[lt]+0.025
//         && tr->L_tr_th[lt]>0.17*tr->L_tr_x[lt]-0.035
//         && tr->L_tr_th[lt]<0.40*tr->L_tr_x[lt]+0.130 ) L_FP = true;

        for(int rt=0;rt<NRtr;rt++){
//          R_Tr = R_FP = false;
//          if( tr->R_tr_chi2[rt]<0.01 ) R_Tr = true;
//          if( tr->R_tr_th[rt]<0.17*tr->R_tr_x[rt]+0.025
//           && tr->R_tr_th[rt]>0.17*tr->R_tr_x[rt]-0.035
//           && tr->R_tr_th[rt]<0.40*tr->R_tr_x[rt]+0.130 ) R_FP = true;
        
//          if(  L_FP && R_FP ){

            int L_s2pad    = (int)L_s2_trpad[lt];
            double L_p     = L_tr_p[lt];//left momentum
            double L_E     = sqrt( Me*Me + L_p*L_p );//left energy
            double L_betae = L_p / sqrt(Me*Me + L_p*L_p);

            int R_s2pad    = (int)R_s2_trpad[rt];
            double R_p     = R_tr_p[rt];//right momentum
            double R_E     = sqrt( MK*MK + R_p*R_p );//right energy
            double R_betaK = R_p / sqrt(MK*MK + R_p*R_p);

            double L_rftime = (L_F1Fhit[47] - L_F1Fhit[40]) * TDCtoT;//internal clock canceling
            double R_rftime = (R_F1Fhit[15] - R_F1Fhit[9] ) * TDCtoT;//internal clock canceling
            double L_tgt = L_s2_t[L_s2pad] -  (L_tr_pathl[lt] + L_s2_trpath[lt])/L_betae/c;
            double R_tgt = R_s2_t[R_s2pad] -  (R_tr_pathl[rt] + R_s2_trpath[rt])/R_betaK/c;
//            double L_tgt = (L_s2_t[L_s2pad] - L_rftime) -  (L_tr_pathl[lt] + L_s2_trpath[lt])/L_betae/c;
//            double R_tgt = (R_s2_t[R_s2pad] - R_rftime) -  (R_tr_pathl[rt] + R_s2_trpath[rt])/R_betaK/c;
            double ct = L_tgt - R_tgt;
            h_ct->Fill( ct );
//            h_ct->Fill( L_tr_n );
		  //}
		}
	  }
	}
  
	h_ct->Draw();

    return 0;
}


