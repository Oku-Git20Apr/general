#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;
#include "ToF_ParamMan.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
ToF_ParamMan::ToF_ParamMan()
{
  CID_TDL =15;
  CID_TagB=18;
  for(int i=0; i<10; i++){
    TDLL_ToF[i] = 0.;
    TDLR_ToF[i] = 0.;
    TDLLU_Wmin[i] = 0.;
    TDLLU_Wmax[i] =50.;
    TDLLD_Wmin[i] = 0.;
    TDLLD_Wmax[i] =50.;
    TDLRU_Wmin[i] = 0.;
    TDLRU_Wmax[i] =50.;
    TDLRD_Wmin[i] = 0.;
    TDLRD_Wmax[i] =50.;
  }
  for(int i=0; i<40; i++){
    TagB_ToF[i] = 0.;
    TagB_Wmin[i] = 0.;
    TagB_Wmax[i] = 0.;
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
ToF_ParamMan::~ToF_ParamMan()
{
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void ToF_ParamMan::read_param(const char* InputFileName){
  static const std::string funcname = "ToF_ParamMan::Initialize";
  FILE *fp;
  char str[200];
  
  if((fp=fopen(InputFileName,"r"))==0){
    std::cerr << "[" << funcname << "]: file open fail" << std::endl;
    exit(-1);
  }
  
  int lr,cid,seg,ud,tw;
  //tw : 0 ToF, 1 widthmin, 2 widthmax
  double p0;
  
  while(fgets(str,200,fp)!=0){
    if(str[0]=='#') continue;
    else if(sscanf(str,"%d %d %d %d %d %lf",&lr,&cid,&seg,&tw,&ud,&p0)==6){
             if(cid==CID_TDL && tw==0&& lr==-1)  TDLL_ToF[seg-1] = p0;
        else if(cid==CID_TDL && tw==0&& lr== 1)  TDLR_ToF[seg-1] = p0;
        else if(cid==CID_TagB&& tw==0        )   TagB_ToF[seg-1]= p0;
        else if(cid==CID_TDL && tw==1&& ud==0 && lr==-1)  TDLLU_Wmin[seg-1] = p0;
        else if(cid==CID_TDL && tw==1&& ud==1 && lr==-1)  TDLLD_Wmin[seg-1] = p0;
        else if(cid==CID_TDL && tw==2&& ud==0 && lr==-1)  TDLLU_Wmax[seg-1] = p0;
        else if(cid==CID_TDL && tw==2&& ud==1 && lr==-1)  TDLLD_Wmax[seg-1] = p0;
        else if(cid==CID_TDL && tw==1&& ud==0 && lr== 1)  TDLRU_Wmin[seg-1] = p0;
        else if(cid==CID_TDL && tw==1&& ud==1 && lr== 1)  TDLRD_Wmin[seg-1] = p0;
        else if(cid==CID_TDL && tw==2&& ud==0 && lr== 1)  TDLRU_Wmax[seg-1] = p0;
        else if(cid==CID_TDL && tw==2&& ud==1 && lr== 1)  TDLRD_Wmax[seg-1] = p0;
        else   cerr << "[" << funcname << "]: fail" << endl;
      }
  }

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void ToF_ParamMan::write_param(const char* OutputFileName){

  ofstream fout;
  if( fout.is_open() ) fout.close();
  fout.open(OutputFileName, ios::out|ios::trunc);
  fout.setf(ios_base::fixed);
  fout << "#" << endl
       << "#  "  << OutputFileName << endl
       << "#" << endl;
  fout << "# LR CID SEG tw ToF[ns]" << endl;
  fout << "#TDLL ToF" << endl;
    for(int i=0; i<10; i++){
      fout << std::setw(3) << -1
	   << std::setw(4) << CID_TDL
	   << std::setw(4) << i+1
	   << std::setw(4) << 0
	   << std::setw(4) << 0
	   << std::setw(13) << std::setprecision(6)
	   << TDLL_ToF[i]<<std::endl;
    }

  fout << "#TDLR ToF" << endl;
    for(int i=0; i<10; i++){
      fout << std::setw(3) << 1
	   << std::setw(4) << CID_TDL
	   << std::setw(4) << i+1
	   << std::setw(4) << 0
	   << std::setw(4) << 0
	   << std::setw(13) << std::setprecision(6)
	   << TDLR_ToF[i]<<std::endl;
    }
  fout << "#TagB ToF" << endl;
    for(int i=0; i<40; i++){
      fout << std::setw(3) << 2
	   << std::setw(4) << CID_TagB
	   << std::setw(4) << i+1
	   << std::setw(4) << 0
	   << std::setw(4) << 0
	   << std::setw(13) << std::setprecision(6)
	   << TagB_ToF[i]<<std::endl;
    }
  fout << "# LR CID SEG tw width[ns]" << endl;
  fout << "#TDLL Width min" << endl;
    for(int i=0; i<10; i++){
      fout << std::setw(3) << -1
	   << std::setw(4) << CID_TDL
	   << std::setw(4) << i+1
	   << std::setw(4) << 1
	   << std::setw(4) << 0
	   << std::setw(13) << std::setprecision(6)
	   << TDLLU_Wmin[i]<<std::endl;
    }
    for(int i=0; i<10; i++){
      fout << std::setw(3) << -1
	   << std::setw(4) << CID_TDL
	   << std::setw(4) << i+1
	   << std::setw(4) << 1
	   << std::setw(4) << 1//ud
	   << std::setw(13) << std::setprecision(6)
	   << TDLLD_Wmin[i]<<std::endl;
    }
  fout << "#TDLL Width max" << endl;
    for(int i=0; i<10; i++){
      fout << std::setw(3) << -1
	   << std::setw(4) << CID_TDL
	   << std::setw(4) << i+1
	   << std::setw(4) << 2
	   << std::setw(4) << 0
	   << std::setw(13) << std::setprecision(6)
	   << TDLLU_Wmax[i]<<std::endl;
    }
    for(int i=0; i<10; i++){
      fout << std::setw(3) << -1
	   << std::setw(4) << CID_TDL
	   << std::setw(4) << i+1
	   << std::setw(4) << 2
	   << std::setw(4) << 1
	   << std::setw(13) << std::setprecision(6)
	   << TDLLD_Wmax[i]<<std::endl;
    }
  fout << "# LR CID SEG tw width[ns]" << endl;
  fout << "#TDLR Width min" << endl;
    for(int i=0; i<10; i++){
      fout << std::setw(3) <<  1
	   << std::setw(4) << CID_TDL
	   << std::setw(4) << i+1
	   << std::setw(4) << 1
	   << std::setw(4) << 0
	   << std::setw(13) << std::setprecision(6)
	   << TDLRU_Wmin[i]<<std::endl;
    }
    for(int i=0; i<10; i++){
      fout << std::setw(3) <<  1
	   << std::setw(4) << CID_TDL
	   << std::setw(4) << i+1
	   << std::setw(4) << 1
	   << std::setw(4) << 1
	   << std::setw(13) << std::setprecision(6)
	   << TDLRD_Wmin[i]<<std::endl;
    }
  fout << "#TDLR Width max" << endl;
    for(int i=0; i<10; i++){
      fout << std::setw(3) <<  1
	   << std::setw(4) << CID_TDL
	   << std::setw(4) << i+1
	   << std::setw(4) << 2
	   << std::setw(4) << 0
	   << std::setw(13) << std::setprecision(6)
	   << TDLRU_Wmax[i]<<std::endl;
    }
    for(int i=0; i<10; i++){
      fout << std::setw(3) <<  1
	   << std::setw(4) << CID_TDL
	   << std::setw(4) << i+1
	   << std::setw(4) << 2
	   << std::setw(4) << 1
	   << std::setw(13) << std::setprecision(6)
	   << TDLRD_Wmax[i]<<std::endl;
    }

std::cout<<"write:"<<OutputFileName<<std::endl;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void ToF_ParamMan::SetToF( int lr, int cid, int seg, double ToF )
{
  static const std::string funcname = "ToF_ParamMan::SetToF";

       if(cid==CID_TDL&& lr==-1)  TDLL_ToF[seg-1] = ToF;
  else if(cid==CID_TDL&& lr== 1)  TDLR_ToF[seg-1] = ToF;
  else if(cid==CID_TagB        )   TagB_ToF[seg-1]= ToF;
  else   cerr << "[" << funcname << "]: fail" << endl;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
double ToF_ParamMan::GetToF( int lr, int cid, int seg )
{
  static const std::string funcname = "ToF_ParamMan::GetToF";

       if(cid==CID_TDL&& lr==-1)   return TDLL_ToF[seg-1];
  else if(cid==CID_TDL&& lr== 1)   return TDLR_ToF[seg-1];
  else if(cid==CID_TagB        )   return TagB_ToF[seg-1];
  else   cerr << "[" << funcname << "]: fail CID:" <<cid<<"  seg:" <<seg<<"  lr:"<<lr<< endl;

  return 0.;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
double ToF_ParamMan::GetWidthLimit( int lr, int cid, int tw, int ud, int seg )
{
  static const std::string funcname = "ToF_ParamMan::GetToF";

       if(cid==CID_TDL && tw==1&& ud==0 && lr==-1)  return TDLLU_Wmin[seg-1];
  else if(cid==CID_TDL && tw==1&& ud==1 && lr==-1)  return TDLLD_Wmin[seg-1];
  else if(cid==CID_TDL && tw==2&& ud==0 && lr==-1)  return TDLLU_Wmax[seg-1];
  else if(cid==CID_TDL && tw==2&& ud==1 && lr==-1)  return TDLLD_Wmax[seg-1];
  else if(cid==CID_TDL && tw==1&& ud==0 && lr== 1)  return TDLRU_Wmin[seg-1];
  else if(cid==CID_TDL && tw==1&& ud==1 && lr== 1)  return TDLRD_Wmin[seg-1];
  else if(cid==CID_TDL && tw==2&& ud==0 && lr== 1)  return TDLRU_Wmax[seg-1];
  else if(cid==CID_TDL && tw==2&& ud==1 && lr== 1)  return TDLRD_Wmax[seg-1];
  else   cerr << "[" << funcname << "]: fail CID:" <<cid<<"  seg:" <<seg<<"  lr:"<<lr<< endl;

  return 0.;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
