void cs2(){

	string pname = "./cs2.dat"; 
	ifstream ifp(pname.c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << pname.c_str() << endl;

	string buf;
	int nofdata = 24;
	int npoint = 0;
	int npoint2 = 0;
	int npoint3 = 0;
	int flag;
	double cs_val,cs_err_val,mom_i,mom_f;
	double cs[nofdata], cs_err[nofdata], pL[nofdata], mom_width[nofdata], Qsq[nofdata], Qsq_width[nofdata];
	double cs2[nofdata], cs_err2[nofdata], pL2[nofdata], mom_width2[nofdata], Qsq2[nofdata], Qsq_width2[nofdata];
	double cs3[nofdata], cs_err3[nofdata], pL3[nofdata], mom_width3[nofdata], Qsq3[nofdata], Qsq_width3[nofdata];



	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		sbuf >> flag >> mom_i >> mom_f >> cs_val >> cs_err_val;
		cout << mom_i << "/" << mom_f << "/" << cs_val << "/" <<cs_err_val <<endl;

		if(flag==1){
		cs[npoint] = cs_val;
		cs_err[npoint]  = cs_err_val;
		pL[npoint] = (mom_f+mom_i)/2;
		mom_width[npoint] = (mom_f-mom_i)/2;
		Qsq[npoint]=2*4.3*pL[npoint]*(1-cos(0.225));
		Qsq_width[npoint]=mom_width[npoint]*2*4.3*(1-cos(0.225));
		npoint++;
		}
		if(flag==2){
		cs2[npoint2] = cs_val;
		cs_err2[npoint2]  = cs_err_val;
		pL2[npoint2] = (mom_f+mom_i)/2;
		mom_width2[npoint2] = (mom_f-mom_i)/2;
		Qsq2[npoint2]=2*4.3*pL2[npoint2]*(1-cos(0.225));
		Qsq_width2[npoint2]=mom_width2[npoint2]*2*4.3*(1-cos(0.225));
		npoint2++;
		}
		if(flag==3){
		cs3[npoint3] = cs_val;
		cs_err3[npoint3]  = cs_err_val;
		pL3[npoint3] = (mom_f+mom_i)/2;
		mom_width3[npoint3] = (mom_f-mom_i)/2;
		Qsq3[npoint3]=2*4.3*pL3[npoint3]*(1-cos(0.225));
		Qsq_width3[npoint3]=mom_width3[npoint3]*2*4.3*(1-cos(0.225));
		npoint3++;
		}
	}

		TGraphErrors *g = new TGraphErrors( npoint, pL, cs, mom_width, cs_err);
		TGraphErrors *g21 = new TGraphErrors( npoint, Qsq, cs, Qsq_width, cs_err);
		TGraphErrors *g22 = new TGraphErrors( npoint2, Qsq2, cs2, Qsq_width2, cs_err2);
		TGraphErrors *g23 = new TGraphErrors( npoint3, Qsq3, cs3, Qsq_width3, cs_err3);
		g->SetMarkerStyle(21);
		g->SetMarkerColor(kAzure);
		g->SetMarkerSize(1.0);
		g21->SetMarkerStyle(21);
		g21->SetMarkerColor(kAzure);
		g21->SetMarkerSize(1.0);
		g22->SetMarkerStyle(21);
		g22->SetMarkerColor(kGreen);
		g22->SetMarkerSize(1.0);
		g23->SetMarkerStyle(21);
		g23->SetMarkerColor(kRed);
		g23->SetMarkerSize(1.0);
	TCanvas *c = new TCanvas("c","Differential Cross Section",600.,600.);
	c->cd()->SetGrid();
	TH1 *frame = c->DrawFrame(2.08,0.,2.22,1000.);
		g->Draw("psame");
	TCanvas *c2 = new TCanvas("c2","Differential Cross Section",600.,600.);
	c2->cd()->SetGrid();
	TH1 *frame2 = c2->DrawFrame(0.45,0.,0.48,1000.);
		g21->Draw("psame");
		g22->Draw("psame");
		g23->Draw("psame");

}
