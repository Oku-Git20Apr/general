void cs(){

	//string pname = "./cs.dat";//Geant vs Uniformity
	string pname = "./cssimc.dat";//SIMC vs Uniformity 
	ifstream ifp(pname.c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << pname.c_str() << endl;

	string buf;
	int nofdata = 24;
	int npoint = 0;
	int flag=0;
	double cs_val,cs_err_val,mom_i,mom_f,NL;
	double cs[nofdata], cs_err[nofdata], pL[nofdata], mom_width[nofdata], Qsq[nofdata], Qsq_width[nofdata], cs_flat[nofdata], cs_flat_err[nofdata];



	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		sbuf >> flag >> mom_i >> mom_f >> cs_val >> cs_err_val >> NL;
		cout << mom_i << "/" << mom_f << "/" << cs_val << "/" <<cs_err_val <<endl;

		cs[npoint] = cs_val;
		cs_flat_err[npoint]=0.;
		if(flag==1){
		cs_flat[npoint] = 442.536*(NL/663);
		}
		if(flag==2){
		cs_flat[npoint] = 442.536*(NL/663)*2.;
		}
		if(flag==3){
		cs_flat[npoint] = 442.536*(NL/663)*(5/2);
		}
		cout<<cs[npoint]<<endl;
		cs_err[npoint]  = cs_err_val;
		pL[npoint] = (mom_f+mom_i)/2;
		cout<<pL[npoint]<<endl;
		mom_width[npoint] = (mom_f-mom_i)/2;
		Qsq[npoint]=2*4.3*pL[npoint]*(1-cos(0.225));
		cout<<Qsq[npoint]<<endl;
		Qsq_width[npoint]=mom_width[npoint]*2*4.3*(1-cos(0.225));
		npoint++;
	}

		TGraphErrors *g = new TGraphErrors( npoint, pL, cs, mom_width, cs_err);
		TGraphErrors *gg = new TGraphErrors( npoint, pL, cs_flat, mom_width, cs_flat_err);
		TGraphErrors *g2 = new TGraphErrors( npoint, Qsq, cs, Qsq_width, cs_err);
		TGraphErrors *g3 = new TGraphErrors( npoint, Qsq, cs_flat, Qsq_width, cs_flat_err);
		g->SetMarkerStyle(21);
		g->SetMarkerColor(kAzure);
		g->SetMarkerSize(1.0);
		gg->SetMarkerStyle(21);
		gg->SetMarkerColor(kRed);
		gg->SetMarkerSize(1.0);
		g2->SetMarkerStyle(21);
		g2->SetMarkerColor(kAzure);
		g2->SetMarkerSize(1.0);
		g3->SetMarkerStyle(21);
		g3->SetMarkerColor(kRed);
		g3->SetMarkerSize(1.0);
	TCanvas *c = new TCanvas("c","Differential Cross Section",600.,600.);
	c->cd()->SetGrid();
	TH1 *frame = c->DrawFrame(2.08,0.,2.22,1000.);
		g->Draw("psame");
		gg->Draw("psame");
	TCanvas *c2 = new TCanvas("c2","Differential Cross Section",600.,600.);
	c2->cd()->SetGrid();
	TH1 *frame2 = c2->DrawFrame(0.45,0.,0.48,1000.);
		g2->Draw("psame");
		g3->Draw("psame");

}
