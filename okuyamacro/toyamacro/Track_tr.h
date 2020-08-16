#ifndef Track_tr_h
#define Track_tr_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

class Track_tr
{
public:
  Track_tr();
  ~Track_tr();

public:
  int runnum, evnum;             // run info
  int spill, spillend, radflag;  // spill info
  int evsize[3];

  int    nt, np, nhits, nshits;

  double chi2v, chi2, trmom;

  int charge;
  double mom , momh , momv;
  double dmom, dmomh, dmomv;
  double ft, fl, beta, mass2;
  double poshIH, poshOH, posvIH, posvOH;
  double zIH, zOH;
  double passIH, passOH;

  int    ihlr, ihseg;
  double ihwid, ihthick;
  int    ihtdc1, ihtdc2, ihadc1, ihadc2;
  double iht1  , iht2  , ihde1 , ihde2 , ihct1 , ihct2;
  double iht, ihct, ihde;
  double ihdedx1, ihdedx2, ihdedx;
  double ihposvt, ihposva;

  int    ohlr, ohseg;
  double ohwid, ohthick;
  int    ohtdc1, ohtdc2, ohadc1, ohadc2;
  double oht1  , oht2  , ohde1 , ohde2 , ohct1 , ohct2;
  double oht, ohct, ohde;
  double ohdedx1, ohdedx2, ohdedx;
  double ohposvt, ohposva;


  int tbseg;
  double tbt,tbct,tbw;//time, ctime, time width
  double tbeg,tbdeg;//E_gamma, dE_gamma

public:
  TChain *tree;

  void readtree();
  void add(string ifname);
  int GetEntries()    const { return tree->GetEntries(); }

private:

};

#endif

