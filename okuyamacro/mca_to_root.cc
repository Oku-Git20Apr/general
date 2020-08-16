Double_t F_lambda( Double_t *x, Double_t *par )
{
 Double_t fitpoly;
fitpoly = par[0] * pow(x[0], 3.0)  + par[1] * pow(x[0], 2.0) + par[2] * x[0] + par[3];
 return fitpoly;
}

void draw(){

const Double_t a = 1.; //x^3
const Double_t b = -1.; //x^2
const Double_t c = 1.; //x^1
const Double_t d = 1.; //x^0
 
 TGraph *g = new TGraph("out1.dat");
 g->Draw("A*");

 TF1 *f1 = new TF1("f1", F_lambda, 10., 300., 4);
 f1->SetParameter(0, a);
 f1->SetParameter(1, b);
 f1->SetParameter(2, c);
 f1->SetParameter(3, d);

 g->Fit(f1,"0","",10.,300.); //fitting

 f1->Draw("same");
 
}
