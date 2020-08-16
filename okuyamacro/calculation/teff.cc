void fillHistos(TH1 *hcut, TH1 *horg);


void teff() {
	  TEfficiency *heff;
	  //TGraphAsymmErrors *heff;
	  //heff->SetStatisticOption(TEfficiency::kBJeffrey);
	  heff->SetConfidenceLevel(0.95);
	  TH1I *hcut = new TH1I("hcut", "after cut", 100, 0, 100);
	  TH1I *horg = new TH1I("horg", "before cut", 100, 0, 100);
	    
      fillHistos(hcut, horg);

//	  if (TEfficiency::CheckConsistency(*hcut, *horg, "w")){
	      heff = new TEfficiency(*hcut, *horg);
	  //    heff = new TGraphAsymmErrors(hcut, horg, "v,cl=0.95");
	      heff->SetTitle("TEfficiency Test;x-axis;y-axis");
//	    }


	 // Drawing
   
     TH1I *haxis = new TH1I("haxis", "black: before cut, red: after cut;Age;Number", 10, 0, 10);
     haxis->SetMaximum(12);
     haxis->SetMinimum(0);

	 //double xerr, yerr;
	 //for(int n=0;n<10;n++){
	 //    xerr = heff->GetErrorX(n);
	 //    yerr = heff->GetErrorY(n);
	 //    cout<<n<<"(x): "<<xerr<<"   (y): "<<yerr<<endl;
	 //}
           
     TCanvas *c1 = new TCanvas("c1", "c1");
     c1->Divide(1,2);
     c1->cd(1);
     haxis->Draw();
     horg->Draw("same");
     hcut->Draw("same");
     c1->cd(2);
     heff->Draw();
     }


  void fillHistos(TH1 *hcut, TH1 *horg){
     hcut->Sumw2();
     horg->Sumw2();
     hcut->SetLineColor(kRed);
     horg->SetLineColor(kBlack);
                                       
                                             
     int nCut = 100;
	 int nTotal = 100;
	 int val;

//	 for(int i=0;i<95;i++){horg->Fill(5);horg->Fill(7);hcut->Fill(5);hcut->Fill(7);}
//	 horg->Fill(5);
//	 horg->Fill(5);
//	 horg->Fill(5);
//	 horg->Fill(5);
//	 horg->Fill(5);
//	 horg->Fill(7);
//	 horg->Fill(7);
//	 horg->Fill(7);
//	 horg->Fill(7);
//	 horg->Fill(7);
//	 hcut->Fill(7);
//	 hcut->Fill(7);
     for (int i = 0; i< nTotal; i++) {val=i/10;horg->Fill(val);}
     for (int i = 0; i< 10; i++) {for(int n=0;n<i;n++)hcut->Fill(n);}
     }
   
