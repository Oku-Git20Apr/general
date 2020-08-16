#ifndef RKVertex_tr_h
#define RKVertex_tr_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

const int NTDL=10;

class RKVertex_tr
{
public:
  RKVertex_tr();
  ~RKVertex_tr();

public:
  int runnum, evnum;             // run info
  int spill, spillend, radflag;  // spill info
  int evsize;

  int    nt, np, nhits, nshits;

  double chi2v[3], chi2[3], trmom[3];

  int charge[3];
  double mom[3];
  double momx[3], momy[3], momz[3];
  double fl_ih[3], fl_oh[3];

  int    ass_nih[3], ihlr[3], ihseg[3];
  double ih_dydx[3], ih_dzdx[3];
  double ihwid[3], ihthick[3];
  int    ihtdc1[3], ihtdc2[3], ihadc1[3], ihadc2[3];
  double iht1  [3], iht2  [3], ihde1 [3], ihde2 [3], ihct1 [3], ihct2[3];
  double iht[3], ihct[3], ihde[3];
  double ihdedx1[3], ihdedx2[3], ihdedx[3];
  double ihposvt[3], ihposva[3];
  double ih_x[3] , ih_y[3] ;
  double ih_nx[3], ih_ny[3];
  double ih_dx[3], ih_dy[3];

  int    ass_noh[3], ohlr[3], ohseg[3];
  double oh_dydx[3], oh_dzdx[3];
  double ohwid[3], ohthick[3];
  int    ohtdc1[3], ohtdc2[3], ohadc1[3], ohadc2[3];
  double oht1  [3], oht2  [3], ohde1 [3], ohde2 [3], ohct1 [3], ohct2[3];
  double oht[3], ohct[3], ohde[3];
  double ohdedx1[3], ohdedx2[3], ohdedx[3];
  double ohposvt[3], ohposva[3];
  double oh_x[3]  , oh_y[3] ;
  double oh_nx[3] , oh_ny[3];
  double oh_dx[3] , oh_dy[3];


  int NumOfTagBHit;
  int tagbtdc[40];
  double tagbctime[40], tagbtime_w[40], tagbtime[40];
  double tagbtime_l[40][10];

  int NumOfTagFHit;
  int tagftdc[160];
  double tagftime[160];

  double vertex[3], pvertex[3], vdist, dlength, dca;
  int vtype;
  int ntr;
   /*
    |vertex position    : vertex[3]
    |momentum @ vertex  : pvertex[3]
    |opning angle       : oa
    |vertex distance    : vdist
    |decay length       : dlength
    |number of track    : ntr
    */

// TDL 
  int ntdl, ntdll, ntdlr;    // Multiplicity
  double tdllutime_w[NTDL], tdlrutime_w[NTDL], tdlldtime_w[NTDL],   tdlrdtime_w[NTDL];   // time (width)
  double tdlluctime[NTDL], tdlructime[NTDL], tdlldctime[NTDL],   tdlrdctime[NTDL];   // time (width)
  double tdllmtime[NTDL], tdlrmtime[NTDL];   // time (mean)
  double tdllmctime[NTDL], tdlrmctime[NTDL];   // time (mean)

public:
  TChain *tree;

  void readtree();
  void add(string ifname);
  int GetEntries()    const { return tree->GetEntries(); }

private:

};

#endif

