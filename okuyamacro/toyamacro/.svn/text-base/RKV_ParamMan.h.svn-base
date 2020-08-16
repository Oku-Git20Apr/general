#ifndef RKV_ParamMan_h
#define RKV_ParamMan_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <sstream>
using namespace std;

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

class RKV_ParamMan
{
public:
  RKV_ParamMan();
  ~RKV_ParamMan();
  void read_param(const char* InputFileName);
  void write_param(const char* OutputFileName);
  void SetT0( int lr, int cid, int seg, double t0 );
  double GetT0( int lr, int cid, int seg );

private:
  int CID_IH, CID_OH, CID_TagB, CID_TDL;
  double IHL_t0[10], IHR_t0[10];
  double TDLL_t0[10], TDLR_t0[10];
  double OHVL_t0[12], OHVR_t0[12], OHHL_t0[9], OHHR_t0[9];
  double TagB_t0[40];
};

#endif
