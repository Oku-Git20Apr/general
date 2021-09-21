void param_calc(){
  const double Na = 6.02214129*1e+23;    // Avogadro constant
  const double hbar = 6.58211928*1e-22;  // Planc constant (reduced) (MeV x s)
  const double PI = 4.*atan(1.);
  double Z = 13;
  double A = 27;
  double ro = 2.7;              // density (g/cm3)
  double delta;                          // density correction
  double C = 0;                              // shell correction
//  double W;                              // maximum energy transfer in a single collision
  double I;                              // mean excitation potential
  double C0, a, m, X0, X1;               // parameters of delta
  double Ne = Na * ro * Z / A;           // density of electrons
  double mup = sqrt(80.617 * 1e+6 * Ne); // plasma frequency (cm^3/2 x Hz) = sqrt(Ne*e^2/pi/me)

//--  I  --//
  if(Z<13) I = (12.*Z + 7.) * 1e-6;
  else I = (9.76 + 58.8*pow(Z,-1.19)) * Z * 1e-6; // MeV

//--  C0  --//
  C0 = - (2. * log(I/(2.*PI) / hbar / mup) + 1.); // MeV

//--  delta  --//
//  if(X<X0) delta = 0;
//  else if(X>=X0 && X<X1) delta = 4.6052*X + C0 + a*pow(X1-X,m);
//  else delta = 4.6052*X + C0;




cout<<"I [eV] = "<<I*1E+6<<endl;
cout<<"-C0    = "<<C0<<endl;
//cout<<"a      = "<<a<<endl;
//cout<<"m      = "<<m<<endl;
//cout<<"X1     = "<<X1<<endl;
//cout<<"X0     = "<<X0<<endl;
}
