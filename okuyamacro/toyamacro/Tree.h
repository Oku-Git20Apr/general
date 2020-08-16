//Class to read out TTree make by UserAnalysis.cpp 

#ifndef Tree_h
#define Tree_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

const int NTDL = 24;
const int NIH  = 10;
const int NOHV = 12;
const int NOHH =  9;
const int NB   = 40;
const int NF   =160;
const int NEV  =  4;
const int NMRPC =2;

class Tree
{
public:
  Tree();
  ~Tree();

public:
  int runnum, evnum;             // run info
  int tdllutdc_l[NTDL][10], tdlrutdc_l[NTDL][10], tdlldtdc_l[NTDL][10],   tdlrdtdc_l[NTDL][10];   // TDC (leading)
  int tdllutdc_t[NTDL][10], tdlrutdc_t[NTDL][10], tdlldtdc_t[NTDL][10],   tdlrdtdc_t[NTDL][10];   // TDC (trailing)
  double tdllmtime[NTDL], tdlrmtime[NTDL];   // time (mean)
  double tdllmctime[NTDL], tdlrmctime[NTDL];   // time (mean)
  int ihlutdc[NIH], ihldtdc[NIH], ihrutdc[NIH], ihrdtdc[NIH];
  int ohvlutdc[NOHV], ohvldtdc[NOHV], ohvrutdc[NOHV], ohvrdtdc[NOHV];
  int ohhlutdc[NOHV], ohhldtdc[NOHV], ohhrutdc[NOHV], ohhrdtdc[NOHV];
  int ohvluadc[NOHH], ohvldadc[NOHH], ohvruadc[NOHH], ohvrdadc[NOHH];
  int ohhluadc[NOHH], ohhldadc[NOHH], ohhruadc[NOHH], ohhrdadc[NOHH];
  int evltdc[NEV], evrtdc[NEV];
  int evladc[NEV], evradc[NEV];
  int evlgltdc, evlgrtdc;
  int evlgladc, evlgradc;
  int tagbtdc_lm[NB][10], tagftdc_lm[NF][10];
  int tagbtdc_tm[NB][10], tagftdc_tm[NF][10];
  int tagbtdc_l[NB], tagftdc_l[NF];
  int tagf_lsize[NF];
  double tdllutime_l[NTDL][10], tdlrutime_l[NTDL][10], tdlldtime_l[NTDL][10],   tdlrdtime_l[NTDL][10];   // time (leading)
  double tdllutime_t[NTDL][10], tdlrutime_t[NTDL][10], tdlldtime_t[NTDL][10],   tdlrdtime_t[NTDL][10];   // time (leading)
  double tdllutime_w[NTDL], tdlrutime_w[NTDL], tdlldtime_w[NTDL],   tdlrdtime_w[NTDL];   // time (width)
  double tdllutime_wn[NTDL], tdlrutime_wn[NTDL], tdlldtime_wn[NTDL],   tdlrdtime_wn[NTDL];   // time (width)
  double tdlluctime[NTDL], tdlructime[NTDL], tdlldctime[NTDL],   tdlrdctime[NTDL];   // time (width)
  double tagbtime_l[NB][10], tagbtime_t[NB][10],tagftime_l[NF][10],tagftime_t[NF];
  double tagbtime_w[NB], tagftime_w[NF];
  double tagbtime[NB], tagbctime[NB];
  double tagb_trig_time, tagf_trig_time1, tagf_trig_time2;//trigger time
  double trig_time1, trig_time2;//trigger time
  double ihlutime[NIH], ihldtime[NIH], ihrutime[NIH], ihrdtime[NIH];
  double ihluctime[NIH], ihldctime[NIH], ihructime[NIH], ihrdctime[NIH];
  double ihlmctime[NIH], ihrmctime[NIH];
  double ihlude[NIH], ihldde[NIH], ihrude[NIH], ihrdde[NIH];
  double ihlde[NIH], ihrde[NIH];
  double ohvlutime[NOHV],ohvldtime[NOHV],ohvrutime[NOHV], ohvrdtime[NOHV];
  double ohvluctime[NOHV], ohvldctime[NOHV], ohvructime[NOHV], ohvrdctime[NOHV];
  double ohvlmctime[NIH], ohvrmctime[NIH];
  double ohvlude[NOHV], ohvldde[NOHV], ohvrude[NOHV], ohvrdde[NOHV];
  double ohvlde[NOHV], ohvrde[NOHV];
  double ohhlutime[NOHH],ohhldtime[NOHH],ohhrutime[NOHH], ohhrdtime[NOHH];
  double ohhluctime[NOHH], ohhldctime[NOHH], ohhructime[NOHH], ohhrdctime[NOHH];
  double ohhlude[NOHH], ohhldde[NOHH], ohhrude[NOHH], ohhrdde[NOHH];
  double ohhlde[NOHH], ohhrde[NOHH];

// MRPC 
  int nmrpc, mrpcseg;
  int mrpcutdc_l[NMRPC][10],mrpcdtdc_l[NMRPC][10];
  int mrpcutdc_t[NMRPC][10],mrpcdtdc_t[NMRPC][10];
  double mrpcutime_l[NMRPC][10],mrpcdtime_l[NMRPC][10];
  double mrpcutime_t[NMRPC][10],mrpcdtime_t[NMRPC][10];
  double mrpcutime_w[NMRPC],mrpcdtime_w[NMRPC];
  double mrpcuctime[NMRPC], mrpcdctime[NMRPC];   // time (width)
  double mrpcmctime[NMRPC];
  double mrpcmtime[NMRPC];

public:
  TChain *tree;

  void readtree();
  void add(string ifname);
  int GetEntries()    const { return tree->GetEntries(); }

private:

};

#endif

