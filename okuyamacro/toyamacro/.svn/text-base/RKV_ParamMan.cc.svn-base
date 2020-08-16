#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;
#include "RKV_ParamMan.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
RKV_ParamMan::RKV_ParamMan()
{
  CID_IH  =1;
  CID_OH  =2;
  CID_TDL =15;
  CID_TagB=18;
  for(int i=0; i<10; i++){
    IHL_t0[i] = 0.;
    IHR_t0[i] = 0.;
    TDLL_t0[i] = 0.;
    TDLR_t0[i] = 0.;
  }
  for(int i=0; i<12; i++){
    OHVL_t0[i] = 0.;
    OHVR_t0[i] = 0.;
  }
  for(int i=0; i< 9; i++){
    OHHL_t0[i] = 0.;
    OHHR_t0[i] = 0.;
  }
  for(int i=0; i<40; i++){
    TagB_t0[i] = 0.;
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
RKV_ParamMan::~RKV_ParamMan()
{
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void RKV_ParamMan::read_param(const char* InputFileName){
  static const std::string funcname = "RKV_ParamMan::Initialize";
  FILE *fp;
  char str[200];
  
  if((fp=fopen(InputFileName,"r"))==0){
    std::cerr << "[" << funcname << "]: file open fail" << std::endl;
    exit(-1);
  }
  
  int lr,cid,seg;
  double p0;
  
  while(fgets(str,200,fp)!=0){
    if(str[0]=='#') continue;
    else if(sscanf(str,"%d %d %d %lf",&lr,&cid,&seg,&p0)==4){
        if(cid==CID_IH && lr==-1)        IHL_t0[seg-1] = p0;
        else if(cid==CID_IH && lr== 1)   IHR_t0[seg-1] = p0;
        else if(cid==CID_TDL&& lr==-1)  TDLL_t0[seg-1] = p0;
        else if(cid==CID_TDL&& lr== 1)  TDLR_t0[seg-1] = p0;
        else if(cid==CID_OH && lr==-1)   OHVL_t0[seg-1]= p0;
        else if(cid==CID_OH && lr== 1)   OHVR_t0[seg-1]= p0;
        else if(cid==CID_OH && lr==-2)   OHHL_t0[seg-1]= p0;
        else if(cid==CID_OH && lr== 2)   OHHR_t0[seg-1]= p0;
        else if(cid==CID_TagB        )   TagB_t0[seg-1]= p0;
        else   cerr << "[" << funcname << "]: fail" << endl;
      }
  }

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void RKV_ParamMan::write_param(const char* OutputFileName){

  ofstream fout;
  if( fout.is_open() ) fout.close();
  fout.open(OutputFileName, ios::out|ios::trunc);
  fout.setf(ios_base::fixed);
  fout << "#" << endl
       << "#  "  << OutputFileName << endl
       << "#" << endl;
  fout << "# LR CID SEG t0[ns]" << endl;
  fout << "#IHL" << endl;
    for(int i=0; i<10; i++){
      fout << std::setw(3) << -1
	   << std::setw(4) << CID_IH
	   << std::setw(4) << i+1
	   << std::setw(13) << std::setprecision(6)
	   << IHL_t0[i]<<std::endl;
    }

  fout << "#IHR" << endl;
    for(int i=0; i<10; i++){
      fout << std::setw(3) << 1
	   << std::setw(4) << CID_IH
	   << std::setw(4) << i+1
	   << std::setw(13) << std::setprecision(6)
	   << IHR_t0[i]<<std::endl;
    }

  fout << "#OHVL" << endl;
    for(int i=0; i<12; i++){
      fout << std::setw(3) << -1
	   << std::setw(4) << CID_OH
	   << std::setw(4) << i+1
	   << std::setw(13) << std::setprecision(6)
	   << OHVL_t0[i]<<std::endl;
    }

  fout << "#OHVR" << endl;
    for(int i=0; i<12; i++){
      fout << std::setw(3) << 1
	   << std::setw(4) << CID_OH
	   << std::setw(4) << i+1
	   << std::setw(13) << std::setprecision(6)
	   << OHVR_t0[i]<<std::endl;
    }

  fout << "#OHHL" << endl;
    for(int i=0; i<9; i++){
      fout << std::setw(3) << -2
	   << std::setw(4) << CID_OH
	   << std::setw(4) << i+1
	   << std::setw(13) << std::setprecision(6)
	   << OHHL_t0[i]<<std::endl;
    }

  fout << "#OHHR" << endl;
    for(int i=0; i<9; i++){
      fout << std::setw(3) << 2
	   << std::setw(4) << CID_OH
	   << std::setw(4) << i+1
	   << std::setw(13) << std::setprecision(6)
	   << OHHR_t0[i]<<std::endl;
    }

  fout << "#TDLL" << endl;
    for(int i=0; i<10; i++){
      fout << std::setw(3) << -1
	   << std::setw(4) << CID_TDL
	   << std::setw(4) << i+1
	   << std::setw(13) << std::setprecision(6)
	   << TDLL_t0[i]<<std::endl;
    }

  fout << "#TDLR" << endl;
    for(int i=0; i<10; i++){
      fout << std::setw(3) << 1
	   << std::setw(4) << CID_TDL
	   << std::setw(4) << i+1
	   << std::setw(13) << std::setprecision(6)
	   << TDLR_t0[i]<<std::endl;
    }
  fout << "#TagB" << endl;
    for(int i=0; i<40; i++){
      fout << std::setw(3) << 2
	   << std::setw(4) << CID_TagB
	   << std::setw(4) << i+1
	   << std::setw(13) << std::setprecision(6)
	   << TagB_t0[i]<<std::endl;
    }

std::cout<<"write:"<<OutputFileName<<std::endl;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void RKV_ParamMan::SetT0( int lr, int cid, int seg, double t0 )
{
  static const std::string funcname = "RKV_ParamMan::SetT0";

  if(cid==CID_IH && lr==0)         IHL_t0[seg-1] = t0;
  else if(cid==CID_IH && lr== 1)   IHR_t0[seg-1] = t0;
  else if(cid==CID_TDL&& lr==-1)  TDLL_t0[seg-1] = t0;
  else if(cid==CID_TDL&& lr== 1)  TDLR_t0[seg-1] = t0;
  else if(cid==CID_OH && lr==-1)   OHVL_t0[seg-1]= t0;
  else if(cid==CID_OH && lr== 1)   OHVR_t0[seg-1]= t0;
  else if(cid==CID_OH && lr==-2)   OHHL_t0[seg-1]= t0;
  else if(cid==CID_OH && lr== 2)   OHHR_t0[seg-1]= t0;
  else if(cid==CID_TagB        )   TagB_t0[seg-1]= t0;
  else   cerr << "[" << funcname << "]: fail" << endl;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
double RKV_ParamMan::GetT0( int lr, int cid, int seg )
{
  static const std::string funcname = "RKV_ParamMan::GetT0";

  if(cid==CID_IH && lr==-1)        return IHL_t0[seg-1] ;
  else if(cid==CID_IH && lr== 1)   return IHR_t0[seg-1] ;
  else if(cid==CID_TDL&& lr==-1)   return TDLL_t0[seg-1];
  else if(cid==CID_TDL&& lr== 1)   return TDLR_t0[seg-1];
  else if(cid==CID_OH && lr==-1)   return OHVL_t0[seg-1];
  else if(cid==CID_OH && lr== 1)   return OHVR_t0[seg-1];
  else if(cid==CID_OH && lr==-2)   return OHHL_t0[seg-1];
  else if(cid==CID_OH && lr== 2)   return OHHR_t0[seg-1];
  else if(cid==CID_TagB        )   return TagB_t0[seg-1];
  else   cerr << "[" << funcname << "]: fail CID:" <<cid<<"  seg:" <<seg<<"  lr:"<<lr<< endl;

  return 0.;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
