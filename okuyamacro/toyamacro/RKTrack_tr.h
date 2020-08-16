#ifndef RKTrack_tr_h
#define RKTrack_tr_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

const int NTDL = 15;
const int NIH = 10;
const int NOHV = 12;
const int NOHH = 9;
const int NB   = 40;
const int NF   = 160;
const int NEV   = 4;
const int TDLtrig = 0;
const int Btrig   = 38;

class RKTrack_tr
{
public:
  RKTrack_tr();
  ~RKTrack_tr();

public:
  int runnum, evnum;             // run info
  int spill, spillend, radflag;  // spill info
  int evsize[3];

  double ihx,ihy,ihz;
  double ohhx,ohhy,ohhz;
  double ohvx,ohvy,ohvz;
  double l13x,l13y,l13z;

  int    nt, np, nhits, nshits;

  double chi2v, chi2, trmom;

  int charge;
  double mom;
  double fl_ih, fl_oh;
  double dl_ih, dl_oh;

  int    ass_nih;
  double ihwid, ihthick;
  int    ihtdc1, ihtdc2, ihadc1, ihadc2;
  double iht1  , iht2  , ihde1 , ihde2 , ihct1 , ihct2;
  double iht, ihct, ihde;
  double ihdedx1, ihdedx2, ihdedx;
  double ihposvt, ihposva;

  int    ass_noh, ohlr, ohseg;
  double ohwid, ohthick;
  int    ohtdc1, ohtdc2, ohadc1, ohadc2;
  double oht1  , oht2  , ohde1 , ohde2 , ohct1 , ohct2;
  double oht, ohct, ohde;
  double ohdedx1, ohdedx2, ohdedx;
  double ohposvt, ohposva;


  int tbseg;
  double tbt,tbct,tbw;//time, ctime, time width
  double tbeg,tbdeg;//E_gamma, dE_gamma

// TDL 
  int tdllseg, tdlrseg;                // Hit segment
  int ntdl, ntdll, ntdlr;    // Multiplicity
  double tdllutime_w[NTDL], tdlrutime_w[NTDL], tdlldtime_w[NTDL],   tdlrdtime_w[NTDL];   // time (width)
  double tdlluctime[NTDL], tdlructime[NTDL], tdlldctime[NTDL],   tdlrdctime[NTDL];   // time (width)
  double tdllmtime[NTDL], tdlrmtime[NTDL];   // time (mean)
  double tdllmctime[NTDL], tdlrmctime[NTDL];   // time (mean)

// OHV
  int noh;//Multiplicity
  int ohvlseg, ohvrseg;                // Hit segment
  int nohvl, nohvr;    // Multiplicity
  double ohvlutime[NOHV],  ohvrutime[NOHV],  ohvldtime[NOHV],  ohvrdtime[NOHV];  // Time
  double ohvluctime[NOHV], ohvructime[NOHV], ohvldctime[NOHV], ohvrdctime[NOHV]; // Time w/PHC
  double ohvlmctime[NOHV], ohvrmctime[NOHV];                                     // mean Time w/PHC
  double ohvlude[NOHV],    ohvrude[NOHV],    ohvldde[NOHV],    ohvrdde[NOHV];    // dE
  double ohvlde[NOHV],    ohvrde[NOHV];                                          // mean dE

// OHH
  int ohhlseg, ohhrseg;                 // Hit segment
  int nohhl, nohhr;     // Multiplicity
  double ohhlutime[NOHH],  ohhrutime[NOHH],  ohhldtime[NOHH],  ohhrdtime[NOHH];  // Time
  double ohhluctime[NOHH], ohhructime[NOHH], ohhldctime[NOHH], ohhrdctime[NOHH]; // Time w/PHC
  double ohhlmctime[NOHH], ohhrmctime[NOHH];                                     // mean Time w/PHC
  double ohhlude[NOHH],    ohhrude[NOHH],    ohhldde[NOHH],    ohhrdde[NOHH];    // dE
  double ohhlde[NOHH],    ohhrde[NOHH];                                         // mean dE
// Tagger
  int ntagb;
  double tagbtime[NB],tagbctime[NB];
  double tagbtime_w[NB];//V1290

///2018.01.10 by Gaio
  double xpos_ih[10], ypos_ih[10];
  double xpos_oh[10], ypos_oh[10];
  int    ihlr[10], ihseg[10];
  int    best_ihseg, best_ihlr;

public:
  TChain *tree;

  void readtree();
  void add(string ifname);
  int GetEntries()    const { return tree->GetEntries(); }

private:

};

#endif

