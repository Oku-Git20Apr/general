#include "RKV_ParamMan.cc"
using namespace std;

void copy_param(string input_param){

  RKV_ParamMan *ParamMan=new RKV_ParamMan();
  RKV_ParamMan *ParamBW =new RKV_ParamMan();
  //ParamBW->read_param("param/OHBW.param");
  //ParamBW->read_param("param/10226_10238bw.param");
  //ParamBW->read_param("param/10296_10350bw.param");
  //ParamBW->read_param("param/group1_10226_10237BW.param");
  //ParamBW->read_param("param/group2_10238_10245BW.param");
  //ParamBW->read_param("param/group3_10246_10256BW.param");
  //ParamBW->read_param("param/group4_10257_10269BW.param");
  //ParamBW->read_param("param/group5_10270_10280BW.param");
  //ParamBW->read_param("param/group6_10282_10293BW.param");
  //ParamBW->read_param("param/group7_10294_10298BW.param");
  //ParamBW->read_param("param/group8_10304_10320BW.param");
  ParamBW->read_param("param/group9_10321_10350BW.param");
  ParamMan->read_param(input_param.c_str());
  int lr,cid;
  //cid=2;
  //for(int i=8; i<12; i++){//OHV8-12
  //  for(int j=0; j<2; j++){//lr
  //    lr=2*j-1;
  //    ParamMan->SetT0( lr,cid,i+1,ParamBW->GetT0(lr,cid,i+1));
  //  }
  //}

  //for(int i=0; i<9; i++){//OHHL1-3,7-9
  //  for(int j=0; j<2; j++){//lr
  //    lr=4*j-2;
  //    if(i>2&&i<6)continue;
  //    ParamMan->SetT0( lr,cid,i+1,ParamBW->GetT0(lr,cid,i+1));
  //  }
  //}


  cid=15;
  for(int i=6; i<9; i++){//TDL7-9
    for(int j=0; j<2; j++){//lr
      lr=2*j-1;
      ParamMan->SetT0( lr,cid,i+1,ParamBW->GetT0(lr,cid,i+1));
    }
  }


  ParamMan->write_param(input_param.c_str());

}
