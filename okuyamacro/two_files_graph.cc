void two_files_graph(){

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
	TH1 *frame = c->DrawFrame(0.,0.,24.,500.);
		g->Draw("psame");

		ifp.close();

	string pname2= "./offset_ac.dat"; 
	ifstream ifp2(pname2.c_str(),ios::in);
	if (ifp2.fail()){ cout << "Failed" << endl; exit(1);}
cout << "Param file 2 : " << pname2.c_str() << endl;
	
	string buf2;
	int nofsegment2 = 24;//a1
	int npoint2 = 0;
	double gain2[nofsegment], err2[nofsegment], segment2[nofsegment], null2[nofsegment];
	
	
	while(1){
		getline(ifp2,buf2);
//		if(buf2[0]=='#'){continue;}
//		if(buf2[0]=='2'){continue;}
		if(ifp2.eof())break;
		if(buf2[0]=='1'){
		stringstream sbuf2(buf2);
		sbuf2 >> ac >> seg >> ped >> ope;

		gain_val = ope - ped;
		cout << ac << "/" << seg << "/" << ped << "/" <<  ope << "/" << gain_val <<  endl;
		gain2[npoint2] = gain_val;
		err2[npoint2]  = 0.;
		segment2[npoint2] = (double)seg;
		null2[npoint2] = 0.;
		npoint2++;
		}
		}

		TGraphErrors *g2 = new TGraphErrors( npoint2, segment2, gain2, null2, err2 );
		g2->SetMarkerStyle(21);
		g2->SetMarkerColor(kRed);
		g2->SetMarkerSize(1.0);
		g2->Draw("psame");


}
