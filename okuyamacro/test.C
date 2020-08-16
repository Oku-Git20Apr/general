void hist(){

const Double_t a = 1000.; //scale
const Double_t b = 500.; //mean
const Double_t c = 100.; //sigma
 

 TH1D *h1 = new TH1D("h1", "h1", 100, -5., 5.);
 h1->FillRandom("gaus",1000.);

 h1->Draw("");
 
}
