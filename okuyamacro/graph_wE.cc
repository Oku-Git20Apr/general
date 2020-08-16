void graph_wE(){

	string pname = "./gain.dat"; 
	ifstream ifp(pname.c_str(),ios::in);
	if (ifp.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file : " << pname.c_str() << endl;

	string buf;
	int nofsegment = 24;//a1
	int npoint = 0;
	int ac, seg;
	double ped, ope, ped_er, ope_er, gain_val, err_val;
	double gain[nofsegment], err[nofsegment], segment[nofsegment], null[nofsegment];



	while(1){
		getline(ifp,buf);
		if(buf[0]=='#'){continue;}
		if(ifp.eof())break;
		stringstream sbuf(buf);
		sbuf >> ac >> seg >> ped >> ped_er >> ope >> ope_er;
//		cout << ac << "/" << seg << "/" << ped << "/" << ped_er << "/" << ope << "/" << ope_er << endl;

		gain_val = ope - ped;
		err_val  = sqrt(ped_er*ped_er+ope_er*ope_er);
		gain[npoint] = gain_val;
		err[npoint]  = err_val;
		segment[npoint] = (double)seg;
		null[npoint] = 0.;
		npoint++;
	}

		TGraphErrors *g = new TGraphErrors( npoint, segment, gain, null, err );
		g->SetMarkerStyle(21);
		g->SetMarkerColor(kAzure);
		g->SetMarkerSize(1.0);
	TCanvas *c = new TCanvas("c","a1 gain",600.,600.);
	c->cd()->SetGrid();
	TH1 *frame = c->DrawFrame(0.,0.,18.,500.);
		g->Draw("plsame");

}
