void graph(){

	string pname = "./data.dat"; 
	ifstream ifp(pname.c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << pname.c_str() << endl;

	string buf;
	int ndata = 24;
	int npoint = 0;
	double th, ss, nn, err_val, fom_val;
	double fom[ndata], error[ndata], threshold[ndata], null[ndata];


	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		sbuf >> th >> ss >> nn;

		ss = ss - nn;
		fom_val = sqrt(ss * ss / nn);
		//err_val  = sqrt(ss/nn + ss*nn/4/nn/nn/nn);
		err_val = 0.;
		error[npoint] = err_val;
		fom[npoint]  = fom_val;
		cout << "S = " << ss << "/ N = " << nn << "/E = " << err_val << endl;
		threshold[npoint] = th;
		null[npoint] = 0.;
		npoint++;
	}

		TGraphErrors *g = new TGraphErrors( npoint, threshold, fom, null, error );
		g->SetMarkerStyle(21);
		g->SetMarkerColor(kAzure);
		g->SetMarkerSize(1.0);
	TCanvas *c = new TCanvas("c","a2 threshold min",600.,600.);
	c->cd()->SetGrid();
	TH1 *frame = c->DrawFrame(0.,0.,3500.,10.);
		g->Draw("plsame");

}
