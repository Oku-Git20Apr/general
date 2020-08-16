void converge_vperr(){

	string pname = "./vperr.dat";//SIMC vs Uniformity 
	ifstream ifp(pname.c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << pname.c_str() << endl;

	string buf;
	int nofdata = 24;
	int npoint = 0;
	int npoint2 = 0;
	int flag=0;
	double vp_val, vp_err_val, nbin_val;
	double vpflux[nofdata], vperr[nofdata], nbin[nofdata], null[nofdata];
	double vpflux2[nofdata], vperr2[nofdata], nbin2[nofdata], null2[nofdata];



	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		//sbuf >> flag >> mom_i >> mom_f >> cs_val >> cs_err_val >> NL;
		sbuf >>flag >> nbin_val >> vp_val >> vp_err_val;
		cout <<flag << "/" << nbin_val << "/" << vp_val << "/" << vp_err_val <<endl;

		if(flag==1){//Max
		vpflux[npoint] = vp_val;
		vperr[npoint]=vp_err_val;
		nbin[npoint]=nbin_val;
		null[npoint]=0.;
		npoint++;
		}//Max
		if(flag==2){//Min
		vpflux2[npoint2] = vp_val;
		vperr2[npoint2]=vp_err_val;
		nbin2[npoint2]=nbin_val;
		null2[npoint2]=0.;
		npoint2++;
		}//Min
	}

		TGraphErrors *g = new TGraphErrors( npoint, nbin, vpflux, null, vperr);
		TGraphErrors *g2 = new TGraphErrors( npoint2, nbin2, vpflux2, null2, vperr2);
		g->SetMarkerStyle(21);
		g->SetMarkerColor(kRed);
		g->SetMarkerSize(1.0);
		g2->SetMarkerStyle(21);
		g2->SetMarkerColor(kAzure);
		g2->SetMarkerSize(1.0);
	TCanvas *c = new TCanvas("c","Integrated VP Flux",600.,600.);
	c->cd()->SetGrid();
	TH1 *frame = c->DrawFrame(0.,2.2*pow(10,-6),160000,3.*pow(10.,-6));
		g->Draw("psame");
		g2->Draw("psame");

}
