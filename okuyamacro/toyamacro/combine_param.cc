//combine_param.cc  2017.8.3 Y.Toyama
//combine parameters


#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
using namespace std;


#include "RKV_ParamMan.h"
#include "HodoParamMan.hh"
#include "TagParamMan.h"

static const double PI = 4.0*atan(1.);
static const double mrad_to_deg = 1./1000*180./PI;
const double Mp = 938.272046;          // proton       mass (MeV/c2)
const double Mpi = 139.57018;          // charged pion mass (MeV/c2)
const double MK = 493.677;             // charged Kaon mass (MeV/c2)
const double c = 0.299792458;          // speed of light in vacuum (m/ns)
const int NIH   = 20;  // No. of IH
const int NTB   = 40;  // No. of TagB
const int NTF   =160;  // No. of TagB
const int NMRPC =  6;  // No. of MRPC 

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class combine_param
{
 public:
         combine_param();
        ~combine_param();
  void SetParamHodo(string ifname); 
  void SetParamTag(string ifname); 
  void SetParamRK(string ifname); 
  void SetValue(); 
  void WriteParam(string hodoname, string tagname); 
  RKV_ParamMan *ParamManRKV;
  HodoParamMan *ParamManHodo;
  TagParamMan  *ParamManTag;

  private:
    int run_num;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
combine_param::combine_param()
{
  ParamManRKV  = new RKV_ParamMan();
  ParamManTag  = new TagParamMan();
}
////////////////////////////////////////////////////////////////////////////
combine_param::~combine_param(){
}
////////////////////////////////////////////////////////////////////////////
void combine_param::SetParamHodo(string ifname){
  ParamManHodo = new HodoParamMan(ifname.c_str() );
  ParamManHodo -> Initialize();
} 
////////////////////////////////////////////////////////////////////////////
void combine_param::SetParamTag(string ifname){
  ParamManTag->ReadParamFile(ifname);
} 
////////////////////////////////////////////////////////////////////////////
void combine_param::SetParamRK(string ifname){
  ParamManRKV->read_param(ifname.c_str());
} 
////////////////////////////////////////////////////////////////////////////
void combine_param::WriteParam(string hodoname, string tagname){
  ParamManHodo -> WriteToFile(hodoname.c_str());
  ParamManTag  -> RewriteParam(tagname);
} 
////////////////////////////////////////////////////////////////////////////
void combine_param::SetValue(){

  int lr,cid;
  for(int i=0;i<12;i++){//OHV
    cid=2;
    lr = 2;ParamManHodo ->   SetTimeTune(lr,cid,i+1,0, -1.*ParamManRKV->GetT0(-1,cid,i+1)); //OHVLU
    lr = 2;ParamManHodo ->   SetTimeTune(lr,cid,i+1,1, -1.*ParamManRKV->GetT0(-1,cid,i+1)); //OHVLD
    lr = 3;ParamManHodo ->   SetTimeTune(lr,cid,i+1,0, -1.*ParamManRKV->GetT0( 1,cid,i+1)); //OHVRU
    lr = 3;ParamManHodo ->   SetTimeTune(lr,cid,i+1,1, -1.*ParamManRKV->GetT0( 1,cid,i+1)); //OHVRD
  }

  for(int i=0;i<9 ;i++){//OHH
    cid=2;
    lr = 0;ParamManHodo ->   SetTimeTune(lr,cid,i+1,0, -1.*ParamManRKV->GetT0(-2,cid,i+1)); //OHHLU
    lr = 0;ParamManHodo ->   SetTimeTune(lr,cid,i+1,1, -1.*ParamManRKV->GetT0(-2,cid,i+1)); //OHHLD
    lr = 1;ParamManHodo ->   SetTimeTune(lr,cid,i+1,0, -1.*ParamManRKV->GetT0( 2,cid,i+1)); //OHHRU
    lr = 1;ParamManHodo ->   SetTimeTune(lr,cid,i+1,1, -1.*ParamManRKV->GetT0( 2,cid,i+1)); //OHHRD
  }

  for(int i=0;i<40;i++){//TagB
    cid = 18;
    ParamManTag -> SetTimeTune("TagB",i+1,  -1*ParamManRKV->GetT0(0,cid,i+1));
  }


  for(int i=0;i<10;i++){//TDL
    cid = 15;
    lr=0;  ParamManHodo ->   SetTimeTune(lr,cid,i+1,0, -1.*ParamManRKV->GetT0(-1,cid,i+1));//TDLLU
    lr=0;  ParamManHodo ->   SetTimeTune(lr,cid,i+1,1, -1.*ParamManRKV->GetT0(-1,cid,i+1));//TDLLD
    lr=1;  ParamManHodo ->   SetTimeTune(lr,cid,i+1,0, -1.*ParamManRKV->GetT0( 1,cid,i+1));//TDLRU
    lr=1;  ParamManHodo ->   SetTimeTune(lr,cid,i+1,1, -1.*ParamManRKV->GetT0( 1,cid,i+1));//TDLRD
  }


}
////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname_hodo = "input_hodo.dat";
  string ifname_tag  = "input_tag.dat";
  string ifname_rkv  = "input_rkv.dat";
  string ofname_hodo = "output_hodo.dat";
  string ofname_tag  = "output_tag.dat";
  int ch;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:t:r:w:g:"))!=-1){
    switch(ch){
    case 'f':
      ifname_hodo = optarg;
      cout<<"input hodo param : "<<ifname_hodo<<endl;
      break;
    case 't':
      ifname_tag = optarg;
      cout<<"input tag param : "<<ifname_tag<<endl;
      break;
    case 'r':
      ifname_rkv = optarg;
      cout<<"input RKV param : "<<ifname_rkv<<endl;
      break;
    case 'w':
      ofname_hodo = optarg;
      cout<<"output hodo param : "<<ofname_hodo<<endl;
      break;
    case 'g':
      ofname_tag = optarg;
      cout<<"output tag param : "<<ofname_tag<<endl;
      break;
    case 'h':
      cout<<"-f : input hodo param filename"<<endl;
      cout<<"-t : input tag param filename"<<endl;
      cout<<"-r : input RKV param filename"<<endl;
      cout<<"-w : output hodo param filename"<<endl;
      cout<<"-g : output tag  param filename"<<endl;
      //cout<<"-p : print pdf file"<<endl;
      return 0;
      break;
    case '?':
      cout<<"unknown option...."<<endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;
    }
  }

  combine_param *combine = new combine_param();

  combine->SetParamHodo(ifname_hodo);
  combine->SetParamRK(ifname_rkv);
  combine->SetParamTag(ifname_tag);
  combine->SetValue();
  combine->WriteParam(ofname_hodo,ofname_tag);
  delete combine;

  return 0;
}

