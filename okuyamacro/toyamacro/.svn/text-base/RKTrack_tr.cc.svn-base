#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

#include "RKTrack_tr.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
RKTrack_tr::RKTrack_tr()
{
  tree = new TChain("tree");
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
RKTrack_tr::~RKTrack_tr()
{
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void RKTrack_tr::add(string ifname)
{
  tree->Add(ifname.c_str());
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void RKTrack_tr::readtree()
{
  tree->SetBranchStatus("*",0);

  tree->SetBranchStatus("runnum"       , 1);  tree->SetBranchAddress("runnum"       , &runnum       );
  tree->SetBranchStatus("evnum"        , 1);  tree->SetBranchAddress("evnum"        , &evnum        );
  tree->SetBranchStatus("spill"        , 1);  tree->SetBranchAddress("spill"        , &spill        );
  tree->SetBranchStatus("spillend"     , 1);  tree->SetBranchAddress("spillend"     , &spillend     );
  tree->SetBranchStatus("radflag"      , 1);  tree->SetBranchAddress("radflag"      , &radflag      );
  tree->SetBranchStatus("evsize"       , 1);  tree->SetBranchAddress("evsize"       ,  evsize       );
                                                                                                    
  tree->SetBranchStatus("ihx"          , 1);  tree->SetBranchAddress("ihx"          , &ihx          );
  tree->SetBranchStatus("ihy"          , 1);  tree->SetBranchAddress("ihy"          , &ihy          );
  tree->SetBranchStatus("ihz"          , 1);  tree->SetBranchAddress("ihz"          , &ihz          );
  tree->SetBranchStatus("ohhx"         , 1);  tree->SetBranchAddress("ohhx"         , &ohhx         );
  tree->SetBranchStatus("ohhy"         , 1);  tree->SetBranchAddress("ohhy"         , &ohhy         );
  tree->SetBranchStatus("ohhz"         , 1);  tree->SetBranchAddress("ohhz"         , &ohhz         );
  tree->SetBranchStatus("ohvx"         , 1);  tree->SetBranchAddress("ohvx"         , &ohvx         );
  tree->SetBranchStatus("ohvy"         , 1);  tree->SetBranchAddress("ohvy"         , &ohvy         );
  tree->SetBranchStatus("ohvz"         , 1);  tree->SetBranchAddress("ohvz"         , &ohvz         );
  tree->SetBranchStatus("l13x"         , 1);  tree->SetBranchAddress("l13x"         , &l13x         );
  tree->SetBranchStatus("l13y"         , 1);  tree->SetBranchAddress("l13y"         , &l13y         );
  tree->SetBranchStatus("l13z"         , 1);  tree->SetBranchAddress("l13z"         , &l13z         );
                                                                                                    
  tree->SetBranchStatus("nt"           , 1);  tree->SetBranchAddress("nt"           , &nt           );
  tree->SetBranchStatus("np"           , 1);  tree->SetBranchAddress("np"           , &np           );
  tree->SetBranchStatus("nhits"        , 1);  tree->SetBranchAddress("nhits"        , &nhits        );
  tree->SetBranchStatus("nshits"       , 1);  tree->SetBranchAddress("nshits"       , &nshits       );
                                                                                                    
  tree->SetBranchStatus("chi2v"        , 1);  tree->SetBranchAddress("chi2v"        , &chi2v        );
  tree->SetBranchStatus("chi2"         , 1);  tree->SetBranchAddress("chi2"         , &chi2         );
  tree->SetBranchStatus("trmom"        , 1);  tree->SetBranchAddress("trmom"        , &trmom        );
                                                                                                    
  tree->SetBranchStatus("charge"       , 1);  tree->SetBranchAddress("charge"       , &charge       );
  tree->SetBranchStatus("mom"          , 1);  tree->SetBranchAddress("mom"          , &mom          );
                                                                                                   
  tree->SetBranchStatus("ass_nih"      , 1);  tree->SetBranchAddress("ass_nih"      , &ass_nih      );
  tree->SetBranchStatus("fl_ih"        , 1);  tree->SetBranchAddress("fl_ih"        , &fl_ih        );
  tree->SetBranchStatus("dl_ih"        , 1);  tree->SetBranchAddress("dl_ih"        , &dl_ih        );
  tree->SetBranchStatus("ihwid"        , 1);  tree->SetBranchAddress("ihwid"        , &ihwid        );
  tree->SetBranchStatus("ihthick"      , 1);  tree->SetBranchAddress("ihthick"      , &ihthick      );
  tree->SetBranchStatus("ihtdc1"       , 1);  tree->SetBranchAddress("ihtdc1"       , &ihtdc1       );
  tree->SetBranchStatus("ihtdc2"       , 1);  tree->SetBranchAddress("ihtdc2"       , &ihtdc2       );
  tree->SetBranchStatus("ihadc1"       , 1);  tree->SetBranchAddress("ihadc1"       , &ihadc1       );
  tree->SetBranchStatus("ihadc2"       , 1);  tree->SetBranchAddress("ihadc2"       , &ihadc2       );
  tree->SetBranchStatus("iht1"         , 1);  tree->SetBranchAddress("iht1"         , &iht1         );
  tree->SetBranchStatus("iht2"         , 1);  tree->SetBranchAddress("iht2"         , &iht2         );
  tree->SetBranchStatus("ihde1"        , 1);  tree->SetBranchAddress("ihde1"        , &ihde1        );
  tree->SetBranchStatus("ihde2"        , 1);  tree->SetBranchAddress("ihde2"        , &ihde2        );
  tree->SetBranchStatus("ihct1"        , 1);  tree->SetBranchAddress("ihct1"        , &ihct1        );
  tree->SetBranchStatus("ihct2"        , 1);  tree->SetBranchAddress("ihct2"        , &ihct2        );
  tree->SetBranchStatus("iht"          , 1);  tree->SetBranchAddress("iht"          , &iht          );
  tree->SetBranchStatus("ihct"         , 1);  tree->SetBranchAddress("ihct"         , &ihct         );
  tree->SetBranchStatus("ihde"         , 1);  tree->SetBranchAddress("ihde"         , &ihde         );
  tree->SetBranchStatus("ihdedx1"      , 1);  tree->SetBranchAddress("ihdedx1"      , &ihdedx1      );
  tree->SetBranchStatus("ihdedx2"      , 1);  tree->SetBranchAddress("ihdedx2"      , &ihdedx2      );
  tree->SetBranchStatus("ihdedx"       , 1);  tree->SetBranchAddress("ihdedx"       , &ihdedx       );
  tree->SetBranchStatus("ihposvt"      , 1);  tree->SetBranchAddress("ihposvt"      , &ihposvt      );
  tree->SetBranchStatus("ihposva"      , 1);  tree->SetBranchAddress("ihposva"      , &ihposva      );
                                                                                                    
  tree->SetBranchStatus("ass_noh"      , 1);  tree->SetBranchAddress("ass_noh"      , &ass_noh      );
  tree->SetBranchStatus("fl_oh"        , 1);  tree->SetBranchAddress("fl_oh"        , &fl_oh        );
  tree->SetBranchStatus("dl_oh"        , 1);  tree->SetBranchAddress("dl_oh"        , &dl_oh        );
  tree->SetBranchStatus("ohlr"         , 1);  tree->SetBranchAddress("ohlr"         , &ohlr         );
  tree->SetBranchStatus("ohseg"        , 1);  tree->SetBranchAddress("ohseg"        , &ohseg        );
  tree->SetBranchStatus("ohwid"        , 1);  tree->SetBranchAddress("ohwid"        , &ohwid        );
  tree->SetBranchStatus("ohthick"      , 1);  tree->SetBranchAddress("ohthick"      , &ohthick      );
  tree->SetBranchStatus("ohtdc1"       , 1);  tree->SetBranchAddress("ohtdc1"       , &ohtdc1       );
  tree->SetBranchStatus("ohtdc2"       , 1);  tree->SetBranchAddress("ohtdc2"       , &ohtdc2       );
  tree->SetBranchStatus("ohadc1"       , 1);  tree->SetBranchAddress("ohadc1"       , &ohadc1       );
  tree->SetBranchStatus("ohadc2"       , 1);  tree->SetBranchAddress("ohadc2"       , &ohadc2       );
  tree->SetBranchStatus("oht1"         , 1);  tree->SetBranchAddress("oht1"         , &oht1         );
  tree->SetBranchStatus("oht2"         , 1);  tree->SetBranchAddress("oht2"         , &oht2         );
  tree->SetBranchStatus("ohde1"        , 1);  tree->SetBranchAddress("ohde1"        , &ohde1        );
  tree->SetBranchStatus("ohde2"        , 1);  tree->SetBranchAddress("ohde2"        , &ohde2        );
  tree->SetBranchStatus("ohct1"        , 1);  tree->SetBranchAddress("ohct1"        , &ohct1        );
  tree->SetBranchStatus("ohct2"        , 1);  tree->SetBranchAddress("ohct2"        , &ohct2        );
  tree->SetBranchStatus("oht"          , 1);  tree->SetBranchAddress("oht"          , &oht          );
  tree->SetBranchStatus("ohct"         , 1);  tree->SetBranchAddress("ohct"         , &ohct         );
  tree->SetBranchStatus("ohde"         , 1);  tree->SetBranchAddress("ohde"         , &ohde         );
  tree->SetBranchStatus("ohdedx1"      , 1);  tree->SetBranchAddress("ohdedx1"      , &ohdedx1      );
  tree->SetBranchStatus("ohdedx2"      , 1);  tree->SetBranchAddress("ohdedx2"      , &ohdedx2      );
  tree->SetBranchStatus("ohdedx"       , 1);  tree->SetBranchAddress("ohdedx"       , &ohdedx       );
  tree->SetBranchStatus("ohposvt"      , 1);  tree->SetBranchAddress("ohposvt"      , &ohposvt      );
  tree->SetBranchStatus("ohposva"      , 1);  tree->SetBranchAddress("ohposva"      , &ohposva      );
                                                                                                   
  tree->SetBranchStatus("tbseg"        , 1);  tree->SetBranchAddress("tbseg"        , &tbseg        );
  tree->SetBranchStatus("tbt"          , 1);  tree->SetBranchAddress("tbt"          , &tbt          );
  tree->SetBranchStatus("tbct"         , 1);  tree->SetBranchAddress("tbct"         , &tbct         );
  tree->SetBranchStatus("tbw"          , 1);  tree->SetBranchAddress("tbw"          , &tbw          );
  tree->SetBranchStatus("tbeg"         , 1);  tree->SetBranchAddress("tbeg"         , &tbeg         );
  tree->SetBranchStatus("tbdeg"        , 1);  tree->SetBranchAddress("tbdeg"        , &tbdeg        );
                                                                                                   
  tree->SetBranchStatus("tdllseg"      , 1);  tree->SetBranchAddress("tdllseg"      , &tdllseg      ); 
  tree->SetBranchStatus("tdlrseg"      , 1);  tree->SetBranchAddress("tdlrseg"      , &tdlrseg      ); 
  tree->SetBranchStatus("ntdll"        , 1);  tree->SetBranchAddress("ntdll"        , &ntdll        ); 
  tree->SetBranchStatus("ntdlr"        , 1);  tree->SetBranchAddress("ntdlr"        , &ntdlr        ); 
  tree->SetBranchStatus("tdllutime_w"  , 1);  tree->SetBranchAddress("tdllutime_w"  ,  tdllutime_w  );
  tree->SetBranchStatus("tdlrutime_w"  , 1);  tree->SetBranchAddress("tdlrutime_w"  ,  tdlrutime_w  );
  tree->SetBranchStatus("tdlldtime_w"  , 1);  tree->SetBranchAddress("tdlldtime_w"  ,  tdlldtime_w  );
  tree->SetBranchStatus("tdlrdtime_w"  , 1);  tree->SetBranchAddress("tdlrdtime_w"  ,  tdlrdtime_w  );
  tree->SetBranchStatus("tdlluctime"   , 1);  tree->SetBranchAddress("tdlluctime"   ,  tdlluctime   );
  tree->SetBranchStatus("tdlructime"   , 1);  tree->SetBranchAddress("tdlructime"   ,  tdlructime   );
  tree->SetBranchStatus("tdlldctime"   , 1);  tree->SetBranchAddress("tdlldctime"   ,  tdlldctime   );
  tree->SetBranchStatus("tdlrdctime"   , 1);  tree->SetBranchAddress("tdlrdctime"   ,  tdlrdctime   );
  tree->SetBranchStatus("tdllmtime"    , 1);  tree->SetBranchAddress("tdllmtime"    ,  tdllmtime    );
  tree->SetBranchStatus("tdlrmtime"    , 1);  tree->SetBranchAddress("tdlrmtime"    ,  tdlrmtime    );
  tree->SetBranchStatus("tdllmctime"   , 1);  tree->SetBranchAddress("tdllmctime"   ,  tdllmctime   );
  tree->SetBranchStatus("tdlrmctime"   , 1);  tree->SetBranchAddress("tdlrmctime"   ,  tdlrmctime   );
                                                                                                  
  tree->SetBranchStatus("ohvlseg"      , 1);  tree->SetBranchAddress("ohvlseg"      , &ohvlseg      );
  tree->SetBranchStatus("ohvrseg"      , 1);  tree->SetBranchAddress("ohvrseg"      , &ohvrseg      );
  tree->SetBranchStatus("noh"          , 1);  tree->SetBranchAddress("noh"          , &noh          );
  tree->SetBranchStatus("nohvl"        , 1);  tree->SetBranchAddress("nohvl"        , &nohvl        );
  tree->SetBranchStatus("nohvr"        , 1);  tree->SetBranchAddress("nohvr"        , &nohvr        );
  tree->SetBranchStatus("ohvlutime"    , 1);  tree->SetBranchAddress("ohvlutime"    ,  ohvlutime    );
  tree->SetBranchStatus("ohvrutime"    , 1);  tree->SetBranchAddress("ohvrutime"    ,  ohvrutime    );
  tree->SetBranchStatus("ohvldtime"    , 1);  tree->SetBranchAddress("ohvldtime"    ,  ohvldtime    );
  tree->SetBranchStatus("ohvrdtime"    , 1);  tree->SetBranchAddress("ohvrdtime"    ,  ohvrdtime    );
  tree->SetBranchStatus("ohvluctime"   , 1);  tree->SetBranchAddress("ohvluctime"   ,  ohvluctime   );
  tree->SetBranchStatus("ohvructime"   , 1);  tree->SetBranchAddress("ohvructime"   ,  ohvructime   );
  tree->SetBranchStatus("ohvldctime"   , 1);  tree->SetBranchAddress("ohvldctime"   ,  ohvldctime   );
  tree->SetBranchStatus("ohvrdctime"   , 1);  tree->SetBranchAddress("ohvrdctime"   ,  ohvrdctime   );
  tree->SetBranchStatus("ohvlmctime"   , 1);  tree->SetBranchAddress("ohvlmctime"   ,  ohvlmctime   );
  tree->SetBranchStatus("ohvrmctime"   , 1);  tree->SetBranchAddress("ohvrmctime"   ,  ohvrmctime   );
  tree->SetBranchStatus("ohvlude"      , 1);  tree->SetBranchAddress("ohvlude"      ,  ohvlude      );
  tree->SetBranchStatus("ohvrude"      , 1);  tree->SetBranchAddress("ohvrude"      ,  ohvrude      );
  tree->SetBranchStatus("ohvldde"      , 1);  tree->SetBranchAddress("ohvldde"      ,  ohvldde      );
  tree->SetBranchStatus("ohvrdde"      , 1);  tree->SetBranchAddress("ohvrdde"      ,  ohvrdde      );
  tree->SetBranchStatus("ohvlde"       , 1);  tree->SetBranchAddress("ohvlde"       ,  ohvlde       );
  tree->SetBranchStatus("ohvrde"       , 1);  tree->SetBranchAddress("ohvrde"       ,  ohvrde       );
  
  tree->SetBranchStatus("ohhlseg"      , 1);  tree->SetBranchAddress("ohhlseg"      , &ohhlseg      );
  tree->SetBranchStatus("ohhrseg"      , 1);  tree->SetBranchAddress("ohhrseg"      , &ohhrseg      );
  tree->SetBranchStatus("nohhl"        , 1);  tree->SetBranchAddress("nohhl"        , &nohhl        );
  tree->SetBranchStatus("nohhr"        , 1);  tree->SetBranchAddress("nohhr"        , &nohhr        );
  tree->SetBranchStatus("ohhlutime"    , 1);  tree->SetBranchAddress("ohhlutime"    ,  ohhlutime    );
  tree->SetBranchStatus("ohhrutime"    , 1);  tree->SetBranchAddress("ohhrutime"    ,  ohhrutime    );
  tree->SetBranchStatus("ohhldtime"    , 1);  tree->SetBranchAddress("ohhldtime"    ,  ohhldtime    );
  tree->SetBranchStatus("ohhrdtime"    , 1);  tree->SetBranchAddress("ohhrdtime"    ,  ohhrdtime    );
  tree->SetBranchStatus("ohhluctime"   , 1);  tree->SetBranchAddress("ohhluctime"   ,  ohhluctime   );
  tree->SetBranchStatus("ohhructime"   , 1);  tree->SetBranchAddress("ohhructime"   ,  ohhructime   );
  tree->SetBranchStatus("ohhldctime"   , 1);  tree->SetBranchAddress("ohhldctime"   ,  ohhldctime   );
  tree->SetBranchStatus("ohhrdctime"   , 1);  tree->SetBranchAddress("ohhrdctime"   ,  ohhrdctime   );
  tree->SetBranchStatus("ohhlmctime"   , 1);  tree->SetBranchAddress("ohhlmctime"   ,  ohhlmctime   );
  tree->SetBranchStatus("ohhrmctime"   , 1);  tree->SetBranchAddress("ohhrmctime"   ,  ohhrmctime   );
  tree->SetBranchStatus("ohhlude"      , 1);  tree->SetBranchAddress("ohhlude"      ,  ohhlude      );
  tree->SetBranchStatus("ohhrude"      , 1);  tree->SetBranchAddress("ohhrude"      ,  ohhrude      );
  tree->SetBranchStatus("ohhldde"      , 1);  tree->SetBranchAddress("ohhldde"      ,  ohhldde      );
  tree->SetBranchStatus("ohhrdde"      , 1);  tree->SetBranchAddress("ohhrdde"      ,  ohhrdde      );
  tree->SetBranchStatus("ohhlde"       , 1);  tree->SetBranchAddress("ohhlde"       ,  ohhlde       );
  tree->SetBranchStatus("ohhrde"       , 1);  tree->SetBranchAddress("ohhrde"       ,  ohhrde       );
 
  tree->SetBranchStatus("ntagb"        , 1);  tree->SetBranchAddress("ntagb"        , &ntagb        );
  tree->SetBranchStatus("tagbtime"     , 1);  tree->SetBranchAddress("tagbtime"     ,  tagbtime     );
  tree->SetBranchStatus("tagbctime"    , 1);  tree->SetBranchAddress("tagbctime"    ,  tagbctime    );
  tree->SetBranchStatus("tagbtime_w"   , 1);  tree->SetBranchAddress("tagbtime_w"   ,  tagbtime_w   );

///2018.01.10 by Gaio
  tree->SetBranchStatus("xpos_ih"      , 1);  tree->SetBranchAddress("xpos_ih"      ,  xpos_ih      );
  tree->SetBranchStatus("ypos_ih"      , 1);  tree->SetBranchAddress("ypos_ih"      ,  ypos_ih      );
  tree->SetBranchStatus("xpos_oh"      , 1);  tree->SetBranchAddress("xpos_oh"      ,  xpos_oh      );
  tree->SetBranchStatus("ypos_oh"      , 1);  tree->SetBranchAddress("ypos_oh"      ,  ypos_oh      );
  tree->SetBranchStatus("best_ihseg"   , 1);  tree->SetBranchAddress("best_ihseg"   , &best_ihseg   );
  tree->SetBranchStatus("best_ihlr"    , 1);  tree->SetBranchAddress("best_ihlr"    , &best_ihlr    );
  tree->SetBranchStatus("ihlr"         , 1);  tree->SetBranchAddress("ihlr"         ,  ihlr         );
  tree->SetBranchStatus("ihseg"        , 1);  tree->SetBranchAddress("ihseg"        ,  ihseg        );
}

