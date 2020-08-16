void hist(){

 double val[100],small_val[100];
 val[0]=0.;
 for(int i=1;i<100;i++) val[i]=val[i-1]+2.;
 small_val[0]=0.;
 for(int i=1;i<100;i++) small_val[i]=small_val[i-1]+0.002;

 TH1D *h1 = new TH1D("h1", "h1", 100, 0., 200.);
 TH1D *h2 = new TH1D("h2", "h2", 100, 0., 0.2);
 
 int nth = 1000.;

 for(int i=0;i<100;i++){
	for(int k=0;k<nth;k++){
		h1->Fill(val[i]);
	}
 }
 for(int i=0;i<100;i++){
	for(int k=0;k<nth;k++){
		h2->Fill(small_val[i]+0.001);
	}
 }

 TCanvas *c1 = new TCanvas("c1","c1",800,600);
 c1->Divide(1,2);
 c1->cd(1);h1->Draw("");
 c1->cd(2);h2->Draw("");
 
}
