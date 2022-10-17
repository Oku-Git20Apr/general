void hypchart_log(){
		gStyle->SetTickLength(0.,"X");
		gStyle->SetTickLength(0.,"Y");

		const int NumOfElement = 100;
		const int Range = 131;
		const double aa = 2.;
		int Nnum[NumOfElement];
		int Pnum[NumOfElement];

		int e = -1;
		e++; Nnum[e] = 1; Pnum[e] = 1; //H3L
		e++; Nnum[e] = 2; Pnum[e] = 1; //H4L
		e++; Nnum[e] = 1; Pnum[e] = 2; //He4L
		e++; Nnum[e] = 2; Pnum[e] = 2; //He5L

		TGraph *g = new TGraph(e+1, Nnum, Pnum);
		g->SetMarkerStyle(25);
		g->SetMarkerColor(kAzure);
		g->SetMarkerSize(100/Range);

							// x0, y0, x1, y1
		TLine *h0 = new TLine(aa-0.5,aa-0.5,aa+0.5,aa-0.5);//y=-0.5
		TLine *h1 = new TLine(aa-0.5,aa+0.5,aa+2.5,aa+0.5);//y=0.5
		TLine *h2 = new TLine(aa+0.5,aa+1.5,aa+5.5,aa+1.5);//y=1.5
		TLine *h3 = new TLine(aa+0.5,aa+2.5,aa+6.5,aa+2.5);//y=2.5
		TLine *h4 = new TLine(aa+1.5,aa+3.5,aa+6.5,aa+3.5);//y=3.5
		TLine *h5 = new TLine(aa+1.5,aa+4.5,aa+6.5,aa+4.5);//y=4.5
		TLine *h6 = new TLine(aa+2.5,aa+5.5,aa+7.5,aa+5.5);//y=5.5
		TLine *h7 = new TLine(aa+4.5,aa+6.5,aa+8.5,aa+6.5);//y=6.5
		TLine *h8 = new TLine(aa+5.5,aa+7.5,aa+8.5,aa+7.5);//y=7.5
		TLine *h9 = new TLine(aa+6.5,aa+8.5,aa+7.5,aa+8.5);//y=8.5
		TLine *h10= new TLine(aa+12.5,aa+12.5,aa+13.5,aa+12.5);//y=12.5
		TLine *h11= new TLine(aa+12.5,aa+13.5,aa+13.5,aa+13.5);//y=13.5
		TLine *h12= new TLine(aa+12.5,aa+14.5,aa+13.5,aa+14.5);//y=14.5
		TLine *h13= new TLine(aa+14.5,aa+15.5,aa+15.5,aa+15.5);//y=15.5
		TLine *h14= new TLine(aa+14.5,aa+16.5,aa+15.5,aa+16.5);//y=16.5
		TLine *h15= new TLine(aa+18.5,aa+19.5,aa+19.5,aa+19.5);//y=19.5
		TLine *h16= new TLine(aa+18.5,aa+20.5,aa+19.5,aa+20.5);//y=20.5
		TLine *h17= new TLine(aa+26.5,aa+22.5,aa+28.5,aa+22.5);//y=22.5
		TLine *h18= new TLine(aa+26.5,aa+23.5,aa+28.5,aa+23.5);//y=23.5
		TLine *h19= new TLine(aa+48.5,aa+38.5,aa+49.5,aa+38.5);//y=38.5
		TLine *h20= new TLine(aa+48.5,aa+39.5,aa+49.5,aa+39.5);//y=39.5
		TLine *h21= new TLine(aa+80.5,aa+56.5,aa+81.5,aa+56.5);//y=56.5
		TLine *h22= new TLine(aa+80.5,aa+57.5,aa+81.5,aa+57.5);//y=57.5
		TLine *h23= new TLine(aa+124.5,aa+81.5,aa+125.5,aa+81.5);//y=81.5
		TLine *h24= new TLine(aa+124.5,aa+82.5,aa+125.5,aa+82.5);//y=82.5
		TLine *hh0= new TLine(aa-0.5,aa-0.5,aa+2.5,aa-0.5);//y=-0.5
		TLine *hh1= new TLine(aa+3.5,aa+0.5,aa+4.5,aa+0.5);//y=0.5
		TLine *hh2= new TLine(aa+8.5,aa+7.5,aa+9.5,aa+7.5);//y=7.5
		TLine *hh3= new TLine(aa+8.5,aa+8.5,aa+9.5,aa+8.5);//y=8.5

		TLine *v0 = new TLine(aa-0.5,aa-0.5,aa-0.5,aa+0.5);//x=-0.5
		TLine *v1 = new TLine(aa+0.5,aa-0.5,aa+0.5,aa+2.5);//x=0.5
		TLine *v2 = new TLine(aa+1.5,aa+0.5,aa+1.5,aa+2.5);//x=1.5 part1
		TLine *v3 = new TLine(aa+1.5,aa+3.5,aa+1.5,aa+4.5);//x=1.5 part2
		TLine *v4 = new TLine(aa+2.5,aa+0.5,aa+2.5,aa+5.5);//x=2.5
		TLine *v5 = new TLine(aa+3.5,aa+1.5,aa+3.5,aa+5.5);//x=3.5
		TLine *v6 = new TLine(aa+4.5,aa+1.5,aa+4.5,aa+6.5);//x=4.5
		TLine *v7 = new TLine(aa+5.5,aa+1.5,aa+5.5,aa+7.5);//x=5.5
		TLine *v8 = new TLine(aa+6.5,aa+2.5,aa+6.5,aa+3.5);//x=6.5 part1
		TLine *v9 = new TLine(aa+6.5,aa+4.5,aa+6.5,aa+8.5);//x=6.5 part2
		TLine *v10= new TLine(aa+7.5,aa+5.5,aa+7.5,aa+8.5);//x=7.5
		TLine *v11= new TLine(aa+8.5,aa+6.5,aa+8.5,aa+7.5);//x=8.5
		TLine *v12= new TLine(aa+12.5,aa+12.5,aa+12.5,aa+14.5);//x=12.5
		TLine *v13= new TLine(aa+13.5,aa+12.5,aa+13.5,aa+14.5);//x=13.5
		TLine *v14= new TLine(aa+14.5,aa+15.5,aa+14.5,aa+16.5);//x=14.5
		TLine *v15= new TLine(aa+15.5,aa+15.5,aa+15.5,aa+16.5);//x=15.5
		TLine *v16= new TLine(aa+18.5,aa+19.5,aa+18.5,aa+20.5);//x=18.5
		TLine *v17= new TLine(aa+19.5,aa+19.5,aa+19.5,aa+20.5);//x=19.5
		TLine *v18= new TLine(aa+26.5,aa+22.5,aa+26.5,aa+23.5);//x=26.5
		TLine *v19= new TLine(aa+27.5,aa+22.5,aa+27.5,aa+23.5);//x=27.5
		TLine *v20= new TLine(aa+28.5,aa+22.5,aa+28.5,aa+23.5);//x=28.5
		TLine *v21= new TLine(aa+48.5,aa+38.5,aa+48.5,aa+39.5);//x=48.5
		TLine *v22= new TLine(aa+49.5,aa+38.5,aa+49.5,aa+39.5);//x=49.5
		TLine *v23= new TLine(aa+80.5,aa+56.5,aa+80.5,aa+57.5);//x=80.5
		TLine *v24= new TLine(aa+81.5,aa+56.5,aa+81.5,aa+57.5);//x=81.5
		TLine *v25= new TLine(aa+124.5,aa+81.5,aa+124.5,aa+82.5);//x=124.5
		TLine *v26= new TLine(aa+125.5,aa+81.5,aa+125.5,aa+82.5);//x=125.5
		TLine *vv0= new TLine(aa+1.5,aa-0.5,aa+1.5,aa+0.5);//x=1.5
		TLine *vv1= new TLine(aa+2.5,aa-0.5,aa+2.5,aa+0.5);//x=2.5
		TLine *vv2= new TLine(aa+3.5,aa+0.5,aa+3.5,aa+1.5);//x=3.5
		TLine *vv3= new TLine(aa+4.5,aa+0.5,aa+4.5,aa+1.5);//x=4.5
		TLine *vv4= new TLine(aa+8.5,aa+7.5,aa+8.5,aa+8.5);//x=8.5
		TLine *vv5= new TLine(aa+9.5,aa+7.5,aa+9.5,aa+8.5);//x=9.5
		h0->SetLineColor(kAzure);
		h1->SetLineColor(kAzure);
		h2->SetLineColor(kAzure);
		h3->SetLineColor(kAzure);
		h4->SetLineColor(kAzure);
		h5->SetLineColor(kAzure);
		h6->SetLineColor(kAzure);
		h7->SetLineColor(kAzure);
		h8->SetLineColor(kAzure);
		h9->SetLineColor(kAzure);
		h10->SetLineColor(kAzure);
		h11->SetLineColor(kAzure);
		h12->SetLineColor(kAzure);
		h13->SetLineColor(kAzure);
		h14->SetLineColor(kAzure);
		h15->SetLineColor(kAzure);
		h16->SetLineColor(kAzure);
		h17->SetLineColor(kAzure);
		h18->SetLineColor(kAzure);
		h19->SetLineColor(kAzure);
		h20->SetLineColor(kAzure);
		h21->SetLineColor(kAzure);
		h22->SetLineColor(kAzure);
		h23->SetLineColor(kAzure);
		h24->SetLineColor(kAzure);
		v0->SetLineColor(kAzure);
		v1->SetLineColor(kAzure);
		v2->SetLineColor(kAzure);
		v3->SetLineColor(kAzure);
		v4->SetLineColor(kAzure);
		v5->SetLineColor(kAzure);
		v6->SetLineColor(kAzure);
		v7->SetLineColor(kAzure);
		v8->SetLineColor(kAzure);
		v9->SetLineColor(kAzure);
		v10->SetLineColor(kAzure);
		v11->SetLineColor(kAzure);
		v12->SetLineColor(kAzure);
		v13->SetLineColor(kAzure);
		v14->SetLineColor(kAzure);
		v15->SetLineColor(kAzure);
		v16->SetLineColor(kAzure);
		v17->SetLineColor(kAzure);
		v18->SetLineColor(kAzure);
		v19->SetLineColor(kAzure);
		v20->SetLineColor(kAzure);
		v21->SetLineColor(kAzure);
		v22->SetLineColor(kAzure);
		v23->SetLineColor(kAzure);
		v24->SetLineColor(kAzure);
		v25->SetLineColor(kAzure);
		v26->SetLineColor(kAzure);
		hh0->SetLineColor(kAzure);
		hh1->SetLineColor(kAzure);
		hh2->SetLineColor(kAzure);
		hh3->SetLineColor(kAzure);
		vv0->SetLineColor(kAzure);
		vv1->SetLineColor(kAzure);
		vv2->SetLineColor(kAzure);
		vv3->SetLineColor(kAzure);
		vv4->SetLineColor(kAzure);
		vv5->SetLineColor(kAzure);
		h0->SetLineWidth(3);
		h1->SetLineWidth(3);
		h2->SetLineWidth(3);
		h3->SetLineWidth(3);
		h4->SetLineWidth(3);
		h5->SetLineWidth(3);
		h6->SetLineWidth(3);
		h7->SetLineWidth(3);
		h8->SetLineWidth(3);
		h9->SetLineWidth(3);
		h10->SetLineWidth(3);
		h11->SetLineWidth(3);
		h12->SetLineWidth(3);
		h13->SetLineWidth(3);
		h14->SetLineWidth(3);
		h15->SetLineWidth(3);
		h16->SetLineWidth(3);
		h17->SetLineWidth(3);
		h18->SetLineWidth(3);
		h19->SetLineWidth(3);
		h20->SetLineWidth(3);
		h21->SetLineWidth(3);
		h22->SetLineWidth(3);
		h23->SetLineWidth(3);
		h24->SetLineWidth(3);
		v0->SetLineWidth(3);
		v1->SetLineWidth(3);
		v2->SetLineWidth(3);
		v3->SetLineWidth(3);
		v4->SetLineWidth(3);
		v5->SetLineWidth(3);
		v6->SetLineWidth(3);
		v7->SetLineWidth(3);
		v8->SetLineWidth(3);
		v9->SetLineWidth(3);
		v10->SetLineWidth(3);
		v11->SetLineWidth(3);
		v12->SetLineWidth(3);
		v13->SetLineWidth(3);
		v14->SetLineWidth(3);
		v15->SetLineWidth(3);
		v16->SetLineWidth(3);
		v17->SetLineWidth(3);
		v18->SetLineWidth(3);
		v19->SetLineWidth(3);
		v20->SetLineWidth(3);
		v21->SetLineWidth(3);
		v22->SetLineWidth(3);
		v23->SetLineWidth(3);
		v24->SetLineWidth(3);
		v25->SetLineWidth(3);
		v26->SetLineWidth(3);
		hh0->SetLineWidth(3);
		hh1->SetLineWidth(3);
		hh2->SetLineWidth(3);
		hh3->SetLineWidth(3);
		vv0->SetLineWidth(3);
		vv1->SetLineWidth(3);
		vv2->SetLineWidth(3);
		vv3->SetLineWidth(3);
		vv4->SetLineWidth(3);
		vv5->SetLineWidth(3);
		hh0->SetLineStyle(7);
		hh1->SetLineStyle(7);
		hh2->SetLineStyle(7);
		hh3->SetLineStyle(7);
		vv0->SetLineStyle(7);
		vv1->SetLineStyle(7);
		vv2->SetLineStyle(7);
		vv3->SetLineStyle(7);
		vv4->SetLineStyle(7);
		vv5->SetLineStyle(7);

	TCanvas *c3 = new TCanvas("c3","All",1000,1000);
	TH1 *frame3 = c3->DrawFrame(1,1,130,100);
	frame3->GetXaxis()->SetTitle("N: Neutron Number");
	frame3->GetYaxis()->SetTitle("Z: Proton Number");
	frame3->GetXaxis()->SetTitleOffset(0.5);
	frame3->GetYaxis()->SetTitleOffset(0.5);
	frame3->GetXaxis()->SetLabelOffset(2.);
	frame3->GetYaxis()->SetLabelOffset(2.);
		c3->SetLogx(1);
		c3->SetLogy(1);
		h0->Draw("same");
		h1->Draw("same");
		h2->Draw("same");
		h3->Draw("same");
		h4->Draw("same");
		h5->Draw("same");
		h6->Draw("same");
		h7->Draw("same");
		h8->Draw("same");
		h9->Draw("same");
		v0->Draw("same");
		v1->Draw("same");
		v2->Draw("same");
		v3->Draw("same");
		v4->Draw("same");
		v5->Draw("same");
		v6->Draw("same");
		v7->Draw("same");
		v8->Draw("same");
		v9->Draw("same");
		v10->Draw("same");
		v11->Draw("same");
		hh0->Draw("same");
		hh1->Draw("same");
		hh2->Draw("same");
		hh3->Draw("same");
		vv0->Draw("same");
		vv1->Draw("same");
		vv2->Draw("same");
		vv3->Draw("same");
		vv4->Draw("same");
		vv5->Draw("same");
//Medium-Heavy hypernuclei
		h10->Draw("same");
		h11->Draw("same");
		h12->Draw("same");
		h13->Draw("same");
		h14->Draw("same");
		h15->Draw("same");
		h16->Draw("same");
		h17->Draw("same");
		h18->Draw("same");
		h19->Draw("same");
		h20->Draw("same");
		h21->Draw("same");
		h22->Draw("same");
		h23->Draw("same");
		h24->Draw("same");
		v12->Draw("same");
		v13->Draw("same");
		v14->Draw("same");
		v15->Draw("same");
		v16->Draw("same");
		v17->Draw("same");
		v18->Draw("same");
		v19->Draw("same");
		v20->Draw("same");
		v21->Draw("same");
		v22->Draw("same");
		v23->Draw("same");
		v24->Draw("same");
		v25->Draw("same");
		v26->Draw("same");

}
