#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

#include "Tree.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
Tree::Tree()
{
  tree = new TChain("tree");
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
Tree::~Tree()
{
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void Tree::add(string ifname)
{
  tree->Add(ifname.c_str());
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void Tree::readtree()
{
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("runnum"         ,1);  tree->SetBranchAddress("runnum"         ,& runnum      );
  tree->SetBranchStatus("evnum"          ,1);  tree->SetBranchAddress("evnum"          ,& evnum       );
  tree->SetBranchStatus("tdllutdc_l"     ,1);  tree->SetBranchAddress("tdllutdc_l"     , tdllutdc_l   );
  tree->SetBranchStatus("tdlldtdc_l"     ,1);  tree->SetBranchAddress("tdlldtdc_l"     , tdlldtdc_l   );
  tree->SetBranchStatus("tdlrutdc_l"     ,1);  tree->SetBranchAddress("tdlrutdc_l"     , tdlrutdc_l   );
  tree->SetBranchStatus("tdlrdtdc_l"     ,1);  tree->SetBranchAddress("tdlrdtdc_l"     , tdlrdtdc_l   );
  tree->SetBranchStatus("tdllutdc_t"     ,1);  tree->SetBranchAddress("tdllutdc_t"     , tdllutdc_t   );
  tree->SetBranchStatus("tdlldtdc_t"     ,1);  tree->SetBranchAddress("tdlldtdc_t"     , tdlldtdc_t   );
  tree->SetBranchStatus("tdlrutdc_t"     ,1);  tree->SetBranchAddress("tdlrutdc_t"     , tdlrutdc_t   );
  tree->SetBranchStatus("tdlrdtdc_t"     ,1);  tree->SetBranchAddress("tdlrdtdc_t"     , tdlrdtdc_t   );
  //tree->SetBranchStatus("ihlutdc"        ,1);  tree->SetBranchAddress("ihlutdc"        , ihlutdc      );
  //tree->SetBranchStatus("ihldtdc"        ,1);  tree->SetBranchAddress("ihldtdc"        , ihldtdc      );
  //tree->SetBranchStatus("ihrutdc"        ,1);  tree->SetBranchAddress("ihrutdc"        , ihrutdc      );
  //tree->SetBranchStatus("ihrdtdc"        ,1);  tree->SetBranchAddress("ihrdtdc"        , ihrdtdc      );
  tree->SetBranchStatus("ohvlutdc"       ,1);  tree->SetBranchAddress("ohvlutdc"       , ohvlutdc     );
  tree->SetBranchStatus("ohvldtdc"       ,1);  tree->SetBranchAddress("ohvldtdc"       , ohvldtdc     );
  tree->SetBranchStatus("ohvrutdc"       ,1);  tree->SetBranchAddress("ohvrutdc"       , ohvrutdc     );
  tree->SetBranchStatus("ohvrdtdc"       ,1);  tree->SetBranchAddress("ohvrdtdc"       , ohvrdtdc     );
  tree->SetBranchStatus("ohvluadc"       ,1);  tree->SetBranchAddress("ohvluadc"       , ohvluadc     );
  tree->SetBranchStatus("ohvldadc"       ,1);  tree->SetBranchAddress("ohvldadc"       , ohvldadc     );
  tree->SetBranchStatus("ohvruadc"       ,1);  tree->SetBranchAddress("ohvruadc"       , ohvruadc     );
  tree->SetBranchStatus("ohvrdadc"       ,1);  tree->SetBranchAddress("ohvrdadc"       , ohvrdadc     );
//  tree->SetBranchStatus("ohhlutdc"       ,1);  tree->SetBranchAddress("ohhlutdc"       , ohhlutdc     );
//  tree->SetBranchStatus("ohhldtdc"       ,1);  tree->SetBranchAddress("ohhldtdc"       , ohhldtdc     );
//  tree->SetBranchStatus("ohhrutdc"       ,1);  tree->SetBranchAddress("ohhrutdc"       , ohhrutdc     );
//  tree->SetBranchStatus("ohhrdtdc"       ,1);  tree->SetBranchAddress("ohhrdtdc"       , ohhrdtdc     );
//  tree->SetBranchStatus("ohhluadc"       ,1);  tree->SetBranchAddress("ohhluadc"       , ohhluadc     );
//  tree->SetBranchStatus("ohhldadc"       ,1);  tree->SetBranchAddress("ohhldadc"       , ohhldadc     );
//  tree->SetBranchStatus("ohhruadc"       ,1);  tree->SetBranchAddress("ohhruadc"       , ohhruadc     );
//  tree->SetBranchStatus("ohhrdadc"       ,1);  tree->SetBranchAddress("ohhrdadc"       , ohhrdadc     );
//  tree->SetBranchStatus("evltdc"         ,1);  tree->SetBranchAddress("evltdc"         , evltdc     );
//  tree->SetBranchStatus("evrtdc"         ,1);  tree->SetBranchAddress("evrtdc"         , evrtdc     );
//  tree->SetBranchStatus("evladc"         ,1);  tree->SetBranchAddress("evladc"         , evladc     );
//  tree->SetBranchStatus("evradc"         ,1);  tree->SetBranchAddress("evradc"         , evradc     );
//  tree->SetBranchStatus("evlgltdc"       ,1);  tree->SetBranchAddress("evlgltdc"       ,&evlgltdc     );
//  tree->SetBranchStatus("evlgrtdc"       ,1);  tree->SetBranchAddress("evlgrtdc"       ,&evlgrtdc     );
//  tree->SetBranchStatus("evlgladc"       ,1);  tree->SetBranchAddress("evlgladc"       ,&evlgladc     );
//  tree->SetBranchStatus("evlgradc"       ,1);  tree->SetBranchAddress("evlgradc"       ,&evlgradc     );
  tree->SetBranchStatus("tagbtdc_lm"     ,1);  tree->SetBranchAddress("tagbtdc_lm"     , tagbtdc_lm     );
  tree->SetBranchStatus("tagbtdc_tm"     ,1);  tree->SetBranchAddress("tagbtdc_tm"     , tagbtdc_tm     );
  tree->SetBranchStatus("tagftdc_lm"     ,1);  tree->SetBranchAddress("tagftdc_lm"     , tagftdc_lm     );
  tree->SetBranchStatus("tagftdc_tm"     ,1);  tree->SetBranchAddress("tagftdc_tm"     , tagftdc_tm     );
  tree->SetBranchStatus("tagf_lsize"     ,1);  tree->SetBranchAddress("tagf_lsize"     , tagf_lsize     );
  tree->SetBranchStatus("tagftime_l"     ,1);  tree->SetBranchAddress("tagftime_l"     , tagftime_l     );
  tree->SetBranchStatus("tagbtime_l"     ,1);  tree->SetBranchAddress("tagbtime_l"     , tagbtime_l );
  tree->SetBranchStatus("tagbtime_w"     ,1);  tree->SetBranchAddress("tagbtime_w"     , tagbtime_w );
  tree->SetBranchStatus("tagbtime"       ,1);  tree->SetBranchAddress("tagbtime"       , tagbtime  );
  tree->SetBranchStatus("tagbctime"      ,1);  tree->SetBranchAddress("tagbctime"      , tagbctime );
  tree->SetBranchStatus("tagbtdc_l"      ,1);  tree->SetBranchAddress("tagbtdc_l"      , tagbtdc_l );
  tree->SetBranchStatus("tagbtdc_lm"     ,1);  tree->SetBranchAddress("tagbtdc_lm"     , tagbtdc_lm     );

  tree->SetBranchStatus("tdllutime_l"    ,1); tree->SetBranchAddress("tdllutime_l"  , tdllutime_l );
  tree->SetBranchStatus("tdlrutime_l"    ,1); tree->SetBranchAddress("tdlrutime_l"  , tdlrutime_l );
  tree->SetBranchStatus("tdlldtime_l"    ,1); tree->SetBranchAddress("tdlldtime_l"  , tdlldtime_l );
  tree->SetBranchStatus("tdlrdtime_l"    ,1); tree->SetBranchAddress("tdlrdtime_l"  , tdlrdtime_l );
  tree->SetBranchStatus("tdllutime_w"    ,1); tree->SetBranchAddress("tdllutime_w"  , tdllutime_w );
  tree->SetBranchStatus("tdlrutime_w"    ,1); tree->SetBranchAddress("tdlrutime_w"  , tdlrutime_w );
  tree->SetBranchStatus("tdlldtime_w"    ,1); tree->SetBranchAddress("tdlldtime_w"  , tdlldtime_w );
  tree->SetBranchStatus("tdlrdtime_w"    ,1); tree->SetBranchAddress("tdlrdtime_w"  , tdlrdtime_w );
  tree->SetBranchStatus("tdllutime_wn"    ,1); tree->SetBranchAddress("tdllutime_wn"  , tdllutime_wn );
  tree->SetBranchStatus("tdlrutime_wn"    ,1); tree->SetBranchAddress("tdlrutime_wn"  , tdlrutime_wn );
  tree->SetBranchStatus("tdlldtime_wn"    ,1); tree->SetBranchAddress("tdlldtime_wn"  , tdlldtime_wn );
  tree->SetBranchStatus("tdlrdtime_wn"    ,1); tree->SetBranchAddress("tdlrdtime_wn"  , tdlrdtime_wn );
  tree->SetBranchStatus("tdlluctime"     ,1); tree->SetBranchAddress("tdlluctime"   , tdlluctime  );
  tree->SetBranchStatus("tdlructime"     ,1); tree->SetBranchAddress("tdlructime"   , tdlructime  );
  tree->SetBranchStatus("tdlldctime"     ,1); tree->SetBranchAddress("tdlldctime"   , tdlldctime  );
  tree->SetBranchStatus("tdlrdctime"     ,1); tree->SetBranchAddress("tdlrdctime"   , tdlrdctime  );
  tree->SetBranchStatus("tdllmtime"      ,1); tree->SetBranchAddress("tdllmtime"    , tdllmtime );
  tree->SetBranchStatus("tdllmctime"     ,1); tree->SetBranchAddress("tdllmctime"   , tdllmctime);
  tree->SetBranchStatus("tdlrmtime"      ,1); tree->SetBranchAddress("tdlrmtime"    , tdlrmtime );
  tree->SetBranchStatus("tdlrmctime"     ,1); tree->SetBranchAddress("tdlrmctime"   , tdlrmctime);

  tree->SetBranchStatus("ohvlutime" , 1);  tree->SetBranchAddress("ohvlutime" , ohvlutime );
  tree->SetBranchStatus("ohvrutime" , 1);  tree->SetBranchAddress("ohvrutime" , ohvrutime );
  tree->SetBranchStatus("ohvldtime" , 1);  tree->SetBranchAddress("ohvldtime" , ohvldtime );
  tree->SetBranchStatus("ohvrdtime" , 1);  tree->SetBranchAddress("ohvrdtime" , ohvrdtime );
  tree->SetBranchStatus("ohvluctime", 1);  tree->SetBranchAddress("ohvluctime", ohvluctime);
  tree->SetBranchStatus("ohvructime", 1);  tree->SetBranchAddress("ohvructime", ohvructime);
  tree->SetBranchStatus("ohvldctime", 1);  tree->SetBranchAddress("ohvldctime", ohvldctime);
  tree->SetBranchStatus("ohvrdctime", 1);  tree->SetBranchAddress("ohvrdctime", ohvrdctime);
  tree->SetBranchStatus("ohvlmctime", 1);  tree->SetBranchAddress("ohvlmctime", ohvlmctime);
  tree->SetBranchStatus("ohvrmctime", 1);  tree->SetBranchAddress("ohvrmctime", ohvrmctime);
  tree->SetBranchStatus("ohvlude"   , 1);  tree->SetBranchAddress("ohvlude"   , ohvlude   );
  tree->SetBranchStatus("ohvrude"   , 1);  tree->SetBranchAddress("ohvrude"   , ohvrude   );
  tree->SetBranchStatus("ohvldde"   , 1);  tree->SetBranchAddress("ohvldde"   , ohvldde   );
  tree->SetBranchStatus("ohvrdde"   , 1);  tree->SetBranchAddress("ohvrdde"   , ohvrdde   );
  tree->SetBranchStatus("ohvlde"    , 1);  tree->SetBranchAddress("ohvlde"    , ohvlde    );
  tree->SetBranchStatus("ohvrde"    , 1);  tree->SetBranchAddress("ohvrde"    , ohvrde    );

//  tree->SetBranchStatus("ohhlutime" , 1);  tree->SetBranchAddress("ohhlutime" , ohhlutime );
//  tree->SetBranchStatus("ohhrutime" , 1);  tree->SetBranchAddress("ohhrutime" , ohhrutime );
//  tree->SetBranchStatus("ohhldtime" , 1);  tree->SetBranchAddress("ohhldtime" , ohhldtime );
//  tree->SetBranchStatus("ohhrdtime" , 1);  tree->SetBranchAddress("ohhrdtime" , ohhrdtime );
//  tree->SetBranchStatus("ohhluctime", 1);  tree->SetBranchAddress("ohhluctime", ohhluctime);
//  tree->SetBranchStatus("ohhructime", 1);  tree->SetBranchAddress("ohhructime", ohhructime);
//  tree->SetBranchStatus("ohhldctime", 1);  tree->SetBranchAddress("ohhldctime", ohhldctime);
//  tree->SetBranchStatus("ohhrdctime", 1);  tree->SetBranchAddress("ohhrdctime", ohhrdctime);
//  tree->SetBranchStatus("ohhlude"   , 1);  tree->SetBranchAddress("ohhlude"   , ohhlude   );
//  tree->SetBranchStatus("ohhrude"   , 1);  tree->SetBranchAddress("ohhrude"   , ohhrude   );
//  tree->SetBranchStatus("ohhldde"   , 1);  tree->SetBranchAddress("ohhldde"   , ohhldde   );
//  tree->SetBranchStatus("ohhrdde"   , 1);  tree->SetBranchAddress("ohhrdde"   , ohhrdde   );
//  tree->SetBranchStatus("ohhlde"    , 1);  tree->SetBranchAddress("ohhlde"    , ohhlde    );
//  tree->SetBranchStatus("ohhrde"    , 1);  tree->SetBranchAddress("ohhrde"    , ohhrde    );

//  tree->SetBranchStatus("ihlutime"  , 1);  tree->SetBranchAddress("ihlutime"  , ihlutime  );
//  tree->SetBranchStatus("ihrutime"  , 1);  tree->SetBranchAddress("ihrutime"  , ihrutime  );
//  tree->SetBranchStatus("ihldtime"  , 1);  tree->SetBranchAddress("ihldtime"  , ihldtime  );
//  tree->SetBranchStatus("ihrdtime"  , 1);  tree->SetBranchAddress("ihrdtime"  , ihrdtime  );
//  tree->SetBranchStatus("ihluctime" , 1);  tree->SetBranchAddress("ihluctime" , ihluctime );
//  tree->SetBranchStatus("ihructime" , 1);  tree->SetBranchAddress("ihructime" , ihructime );
//  tree->SetBranchStatus("ihldctime" , 1);  tree->SetBranchAddress("ihldctime" , ihldctime );
//  tree->SetBranchStatus("ihrdctime" , 1);  tree->SetBranchAddress("ihrdctime" , ihrdctime );
//  tree->SetBranchStatus("ihlmctime" , 1);  tree->SetBranchAddress("ihlmctime" , ihlmctime );
//  tree->SetBranchStatus("ihrmctime" , 1);  tree->SetBranchAddress("ihrmctime" , ihrmctime );
//  tree->SetBranchStatus("ihlude"    , 1);  tree->SetBranchAddress("ihlude"    , ihlude    );
//  tree->SetBranchStatus("ihrude"    , 1);  tree->SetBranchAddress("ihrude"    , ihrude    );
//  tree->SetBranchStatus("ihldde"    , 1);  tree->SetBranchAddress("ihldde"    , ihldde    );
//  tree->SetBranchStatus("ihrdde"    , 1);  tree->SetBranchAddress("ihrdde"    , ihrdde    );
//  tree->SetBranchStatus("ihlde"     , 1);  tree->SetBranchAddress("ihlde"     , ihlde     );
//  tree->SetBranchStatus("ihrde"     , 1);  tree->SetBranchAddress("ihrde"     , ihrde     );
  tree->SetBranchStatus("mrpcseg"    ,1);  tree->SetBranchAddress("mrpcseg"    ,&mrpcseg      );
  tree->SetBranchStatus("nmrpc"      ,1);  tree->SetBranchAddress("nmrpc"      ,&nmrpc        );
  tree->SetBranchStatus("mrpcutdc_l" ,1);  tree->SetBranchAddress("mrpcutdc_l" , mrpcutdc_l   );
  tree->SetBranchStatus("mrpcdtdc_l" ,1);  tree->SetBranchAddress("mrpcdtdc_l" , mrpcdtdc_l   );
  tree->SetBranchStatus("mrpcutdc_t" ,1);  tree->SetBranchAddress("mrpcutdc_t" , mrpcutdc_t   );
  tree->SetBranchStatus("mrpcdtdc_t" ,1);  tree->SetBranchAddress("mrpcdtdc_t" , mrpcdtdc_t   );
  tree->SetBranchStatus("mrpcutime_l",1);  tree->SetBranchAddress("mrpcutime_l", mrpcutime_l  );
  tree->SetBranchStatus("mrpcdtime_l",1);  tree->SetBranchAddress("mrpcdtime_l", mrpcdtime_l  );
  tree->SetBranchStatus("mrpcutime_t",1);  tree->SetBranchAddress("mrpcutime_t", mrpcutime_t  );
  tree->SetBranchStatus("mrpcdtime_t",1);  tree->SetBranchAddress("mrpcdtime_t", mrpcdtime_t  );
  tree->SetBranchStatus("mrpcutime_w",1);  tree->SetBranchAddress("mrpcutime_w", mrpcutime_w  );
  tree->SetBranchStatus("mrpcdtime_w",1);  tree->SetBranchAddress("mrpcdtime_w", mrpcdtime_w  );
  tree->SetBranchStatus("mrpcmtime"  ,1);  tree->SetBranchAddress("mrpcmtime"  , mrpcmtime    );
  tree->SetBranchStatus("mrpcuctime" ,1);  tree->SetBranchAddress("mrpcuctime" , mrpcuctime   );
  tree->SetBranchStatus("mrpcdctime" ,1);  tree->SetBranchAddress("mrpcdctime" , mrpcdctime   );
  tree->SetBranchStatus("mrpcmctime" ,1);  tree->SetBranchAddress("mrpcmctime" , mrpcmctime   );
}

