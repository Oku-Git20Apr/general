#ifndef ToF_ParamMan_h
#define ToF_ParamMan_h

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

class ToF_ParamMan
{
public:
  ToF_ParamMan();
  ~ToF_ParamMan();
  void read_param(const char* InputFileName);
  void write_param(const char* OutputFileName);
  void SetToF( int lr, int cid, int seg, double ToF );
  double GetToF( int lr, int cid, int seg );
  double GetWidthLimit( int lr, int cid, int tw, int ud, int seg );

private:
  int CID_TagB, CID_TDL;
  double TDLL_ToF[10], TDLR_ToF[10];
  double TDLLU_Wmin[10], TDLLU_Wmax[10];
  double TDLLD_Wmin[10], TDLLD_Wmax[10];
  double TDLRU_Wmin[10], TDLRU_Wmax[10];
  double TDLRD_Wmin[10], TDLRD_Wmax[10];
  double TagB_ToF[40];
  double TagB_Wmin[10], TagB_Wmax[10];
};

#endif
