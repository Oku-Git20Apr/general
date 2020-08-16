#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

#include "RKVertex_tr.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
RKVertex_tr::RKVertex_tr()
{
  tree = new TChain("tree");
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
RKVertex_tr::~RKVertex_tr()
{
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void RKVertex_tr::add(string ifname)
{
  tree->Add(ifname.c_str());
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void RKVertex_tr::readtree()
{
 tree->SetBranchStatus("*",0);
 tree->SetBranchStatus("runnum"     ,1); tree->SetBranchAddress( "runnum"    , &runnum   );
 tree->SetBranchStatus("evnum"      ,1); tree->SetBranchAddress( "evnum"     , &evnum    );


//vertex
  tree->SetBranchStatus("vtype"     ,1); tree->SetBranchAddress( "vtype"     , &vtype    );
  tree->SetBranchStatus("ntr"       ,1); tree->SetBranchAddress( "ntr"       , &ntr      );
  tree->SetBranchStatus("vertex"    ,1); tree->SetBranchAddress( "vertex"    ,  vertex   );
  tree->SetBranchStatus("pvertex"   ,1); tree->SetBranchAddress( "pvertex"   ,  pvertex  );
  tree->SetBranchStatus("vdist"     ,1); tree->SetBranchAddress( "vdist"     , &vdist    );
  tree->SetBranchStatus("dlength"   ,1); tree->SetBranchAddress( "dlength"   , &dlength  );
  tree->SetBranchStatus("dca"       ,1); tree->SetBranchAddress( "dca"       , &dca      );

//track
  tree->SetBranchStatus("nt"        ,1); tree->SetBranchAddress( "nt"        , &nt       );
  tree->SetBranchStatus("np"        ,1); tree->SetBranchAddress( "np"        , &np       );
  tree->SetBranchStatus("nhits"     ,1); tree->SetBranchAddress( "nhits"     , &nhits    );
  tree->SetBranchStatus("nshits"    ,1); tree->SetBranchAddress( "nshits"    , &nshits   );
                                                                              
  tree->SetBranchStatus("chi2v"     ,1); tree->SetBranchAddress( "chi2v"     , chi2v     );
  tree->SetBranchStatus("chi2"      ,1); tree->SetBranchAddress( "chi2"      , chi2      );
  tree->SetBranchStatus("trmom"     ,1); tree->SetBranchAddress( "trmom"     , trmom     );
                                                                              
  tree->SetBranchStatus("charge"    ,1); tree->SetBranchAddress( "charge"    , charge    );
  tree->SetBranchStatus("mom"       ,1); tree->SetBranchAddress( "mom"       , mom       );
  tree->SetBranchStatus("momx"      ,1); tree->SetBranchAddress( "momx"      , momx       );
  tree->SetBranchStatus("momy"      ,1); tree->SetBranchAddress( "momy"      , momy       );
  tree->SetBranchStatus("momz"      ,1); tree->SetBranchAddress( "momz"      , momz       );
                                                                              
  //IH
  tree->SetBranchStatus("ass_nih"   ,1); tree->SetBranchAddress( "ass_nih"   , ass_nih   );
  tree->SetBranchStatus("fl_ih"     ,1); tree->SetBranchAddress( "fl_ih"     , fl_ih     );
  tree->SetBranchStatus("ihlr"      ,1); tree->SetBranchAddress( "ihlr"      , ihlr      );
  tree->SetBranchStatus("ihseg"     ,1); tree->SetBranchAddress( "ihseg"     , ihseg     );
  tree->SetBranchStatus("ih_dydx"   ,1); tree->SetBranchAddress( "ih_dydx"   , ih_dydx   );
  tree->SetBranchStatus("ih_dzdx"   ,1); tree->SetBranchAddress( "ih_dzdx"   , ih_dzdx   );
  tree->SetBranchStatus("ihwid"     ,1); tree->SetBranchAddress( "ihwid"     , ihwid     );
  tree->SetBranchStatus("ihthick"   ,1); tree->SetBranchAddress( "ihthick"   , ihthick   );
  tree->SetBranchStatus("ihtdc1"    ,1); tree->SetBranchAddress( "ihtdc1"    , ihtdc1    );
  tree->SetBranchStatus("ihtdc2"    ,1); tree->SetBranchAddress( "ihtdc2"    , ihtdc2    );
  tree->SetBranchStatus("ihadc1"    ,1); tree->SetBranchAddress( "ihadc1"    , ihadc1    );
  tree->SetBranchStatus("ihadc2"    ,1); tree->SetBranchAddress( "ihadc2"    , ihadc2    );
  tree->SetBranchStatus("iht1"      ,1); tree->SetBranchAddress( "iht1"      , iht1      );
  tree->SetBranchStatus("iht2"      ,1); tree->SetBranchAddress( "iht2"      , iht2      );
  tree->SetBranchStatus("ihde1"     ,1); tree->SetBranchAddress( "ihde1"     , ihde1     );
  tree->SetBranchStatus("ihde2"     ,1); tree->SetBranchAddress( "ihde2"     , ihde2     );
  tree->SetBranchStatus("ihdedx1"   ,1); tree->SetBranchAddress( "ihdedx1"   , ihdedx1   );
  tree->SetBranchStatus("ihdedx2"   ,1); tree->SetBranchAddress( "ihdedx2"   , ihdedx2   );
  tree->SetBranchStatus("ihdedx"    ,1); tree->SetBranchAddress( "ihdedx"    , ihdedx    );
  tree->SetBranchStatus("ihct1"     ,1); tree->SetBranchAddress( "ihct1"     , ihct1     );
  tree->SetBranchStatus("ihct2"     ,1); tree->SetBranchAddress( "ihct2"     , ihct2     );
  tree->SetBranchStatus("iht"       ,1); tree->SetBranchAddress( "iht"       , iht       );
  tree->SetBranchStatus("ihct"      ,1); tree->SetBranchAddress( "ihct"      , ihct      );
  tree->SetBranchStatus("ihde"      ,1); tree->SetBranchAddress( "ihde"      , ihde      );
  tree->SetBranchStatus("ihposvt"   ,1); tree->SetBranchAddress( "ihposvt"   , ihposvt   );
  //tree->SetBranchStatus("ihposva"   ,1); tree->SetBranchAddress( "ihposva"   , ihposva   );
  tree->SetBranchStatus("ih_x"      ,1); tree->SetBranchAddress( "ih_x"      , ih_x      );
  tree->SetBranchStatus("ih_y"      ,1); tree->SetBranchAddress( "ih_y"      , ih_y      );
  tree->SetBranchStatus("ih_nx"     ,1); tree->SetBranchAddress( "ih_nx"     , ih_nx     );
  tree->SetBranchStatus("ih_ny"     ,1); tree->SetBranchAddress( "ih_ny"     , ih_ny     );
  tree->SetBranchStatus("ih_dx"     ,1); tree->SetBranchAddress( "ih_dx"     , ih_dx     );
  tree->SetBranchStatus("ih_dy"     ,1); tree->SetBranchAddress( "ih_dy"     , ih_dy     );
                                                                              
  //OH
  tree->SetBranchStatus("ass_noh"   ,1); tree->SetBranchAddress( "ass_noh"   , ass_noh   );
  tree->SetBranchStatus("fl_oh"     ,1); tree->SetBranchAddress( "fl_oh"     , fl_oh     );
  tree->SetBranchStatus("ohlr"      ,1); tree->SetBranchAddress( "ohlr"      , ohlr      );
  tree->SetBranchStatus("ohseg"     ,1); tree->SetBranchAddress( "ohseg"     , ohseg     );
  tree->SetBranchStatus("oh_dydx"   ,1); tree->SetBranchAddress( "oh_dydx"   , oh_dydx   );
  tree->SetBranchStatus("oh_dzdx"   ,1); tree->SetBranchAddress( "oh_dzdx"   , oh_dzdx   );
  tree->SetBranchStatus("ohwid"     ,1); tree->SetBranchAddress( "ohwid"     , ohwid     );
  tree->SetBranchStatus("ohthick"   ,1); tree->SetBranchAddress( "ohthick"   , ohthick   );
  tree->SetBranchStatus("ohtdc1"    ,1); tree->SetBranchAddress( "ohtdc1"    , ohtdc1    );
  tree->SetBranchStatus("ohtdc2"    ,1); tree->SetBranchAddress( "ohtdc2"    , ohtdc2    );
  tree->SetBranchStatus("ohadc1"    ,1); tree->SetBranchAddress( "ohadc1"    , ohadc1    );
  tree->SetBranchStatus("ohadc2"    ,1); tree->SetBranchAddress( "ohadc2"    , ohadc2    );
  tree->SetBranchStatus("oht1"      ,1); tree->SetBranchAddress( "oht1"      , oht1      );
  tree->SetBranchStatus("oht2"      ,1); tree->SetBranchAddress( "oht2"      , oht2      );
  tree->SetBranchStatus("ohde1"     ,1); tree->SetBranchAddress( "ohde1"     , ohde1     );
  tree->SetBranchStatus("ohde2"     ,1); tree->SetBranchAddress( "ohde2"     , ohde2     );
  tree->SetBranchStatus("ohdedx1"   ,1); tree->SetBranchAddress( "ohdedx1"   , ohdedx1   );
  tree->SetBranchStatus("ohdedx2"   ,1); tree->SetBranchAddress( "ohdedx2"   , ohdedx2   );
  tree->SetBranchStatus("ohdedx"    ,1); tree->SetBranchAddress( "ohdedx"    , ohdedx    );
  tree->SetBranchStatus("ohct1"     ,1); tree->SetBranchAddress( "ohct1"     , ohct1     );
  tree->SetBranchStatus("ohct2"     ,1); tree->SetBranchAddress( "ohct2"     , ohct2     );
  tree->SetBranchStatus("oht"       ,1); tree->SetBranchAddress( "oht"       , oht       );
  tree->SetBranchStatus("ohct"      ,1); tree->SetBranchAddress( "ohct"      , ohct      );
  tree->SetBranchStatus("ohde"      ,1); tree->SetBranchAddress( "ohde"      , ohde      );
  tree->SetBranchStatus("ohposvt"   ,1); tree->SetBranchAddress( "ohposvt"   , ohposvt   );
  tree->SetBranchStatus("ohposva"   ,1); tree->SetBranchAddress( "ohposva"   , ohposva   );
  tree->SetBranchStatus("oh_x"      ,1); tree->SetBranchAddress( "oh_x"      , oh_x      );
  tree->SetBranchStatus("oh_y"      ,1); tree->SetBranchAddress( "oh_y"      , oh_y      );
  tree->SetBranchStatus("oh_nx"     ,1); tree->SetBranchAddress( "oh_nx"     , oh_nx     );
  tree->SetBranchStatus("oh_ny"     ,1); tree->SetBranchAddress( "oh_ny"     , oh_ny     );
  tree->SetBranchStatus("oh_dx"     ,1); tree->SetBranchAddress( "oh_dx"     , oh_dx     );
  tree->SetBranchStatus("oh_dy"     ,1); tree->SetBranchAddress( "oh_dy"     , oh_dy     );

  //Tagger
  tree->SetBranchStatus("NumOfTagBHit",1); tree->SetBranchAddress("NumOfTagBHit", &NumOfTagBHit);
  tree->SetBranchStatus("tagbtdc"     ,1); tree->SetBranchAddress("tagbtdc"     ,  tagbtdc     );
  tree->SetBranchStatus("tagbtime"    ,1); tree->SetBranchAddress("tagbtime"    ,  tagbtime    );
  tree->SetBranchStatus("tagbtime_l"  ,1); tree->SetBranchAddress("tagbtime_l"  ,  tagbtime_l  );
  tree->SetBranchStatus("tagbctime"   ,1); tree->SetBranchAddress("tagbctime"   ,  tagbctime   );
  tree->SetBranchStatus("tagbtime_w"  ,1); tree->SetBranchAddress("tagbtime_w"  ,  tagbtime_w  );
  tree->SetBranchStatus("NumOfTagFHit",1); tree->SetBranchAddress("NumOfTagFHit", &NumOfTagFHit);
  tree->SetBranchStatus("tagftdc"     ,1); tree->SetBranchAddress("tagftdc"     ,  tagftdc     );
  tree->SetBranchStatus("tagftime"    ,1); tree->SetBranchAddress("tagftime"    ,  tagftime    );
//TDL
  tree->SetBranchStatus("ntdll"        ,1);tree->SetBranchAddress("ntdll"       ,&ntdll      );
  tree->SetBranchStatus("ntdlr"        ,1);tree->SetBranchAddress("ntdlr"       ,&ntdlr      );
  tree->SetBranchStatus("tdllutime_w"  ,1);tree->SetBranchAddress("tdllutime_w" , tdllutime_w);
  tree->SetBranchStatus("tdlrutime_w"  ,1);tree->SetBranchAddress("tdlrutime_w" , tdlrutime_w);
  tree->SetBranchStatus("tdlldtime_w"  ,1);tree->SetBranchAddress("tdlldtime_w" , tdlldtime_w);
  tree->SetBranchStatus("tdlrdtime_w"  ,1);tree->SetBranchAddress("tdlrdtime_w" , tdlrdtime_w);
  tree->SetBranchStatus("tdlluctime"   ,1);tree->SetBranchAddress("tdlluctime"  , tdlluctime );
  tree->SetBranchStatus("tdlructime"   ,1);tree->SetBranchAddress("tdlructime"  , tdlructime );
  tree->SetBranchStatus("tdlldctime"   ,1);tree->SetBranchAddress("tdlldctime"  , tdlldctime );
  tree->SetBranchStatus("tdlrdctime"   ,1);tree->SetBranchAddress("tdlrdctime"  , tdlrdctime );
  tree->SetBranchStatus("tdllmtime"    ,1);tree->SetBranchAddress("tdllmtime"   , tdllmtime  );
  tree->SetBranchStatus("tdlrmtime"    ,1);tree->SetBranchAddress("tdlrmtime"   , tdlrmtime  );
  tree->SetBranchStatus("tdllmctime"   ,1);tree->SetBranchAddress("tdllmctime"  , tdllmctime );
  tree->SetBranchStatus("tdlrmctime"   ,1);tree->SetBranchAddress("tdlrmctime"  , tdlrmctime );
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
