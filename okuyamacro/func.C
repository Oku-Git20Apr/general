double F_Gauss(double *t, double *par){

	double x = t[0];
	double scale = par[0];
	double mean = par[1];
	double sigma = par[2];

	double y = 1/sqrt(2.*TMath::Pi()*sigma*sigma)*TMath::Exp(-pow((x-mean),2.)/2/sigma/sigma);
	return y;
}

double F_Voigt( double *x, double *par )
{
  // par[0] : area
  // par[1] : location
  // par[2] : gaussian sigma
  // par[3] : lorentz fwhm
  double val = par[0] * TMath::Voigt(x[0]-par[1],par[2],par[3],4);
  return val;
}

void func(){

const Double_t a = -5.; //scale
const Double_t b = 0.; //mean
const Double_t c = 0.05; //sigma
const Double_t d = -0.01; //sigma
 

// TF1 *f1 = new TF1("f1", F_Gauss, 0., 1000., 3);
 TF1 *f1 = new TF1("f1", F_Voigt, -0.1, 0.15, 4);
 f1->SetParameter(0, a);
 f1->SetParameter(1, b);
 f1->SetParameter(2, c);
 f1->SetParameter(3, d);

 f1->Draw("");
 
}
