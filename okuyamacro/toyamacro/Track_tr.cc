#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

#include "Track_tr.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
Track_tr::Track_tr()
{
  tree = new TChain("tree");
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
Track_tr::~Track_tr()
{
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void Track_tr::add(string ifname)
{
  tree->Add(ifname.c_str());
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void Track_tr::readtree()
{
 tree->SetBranchStatus("*",0);
 tree->SetBranchStatus("runnum"    ,1);  tree->SetBranchAddress("runnum"    , &runnum   );
 tree->SetBranchStatus("evnum"     ,1);  tree->SetBranchAddress("evnum"     , &evnum    );
 tree->SetBranchStatus("spill"     ,1);  tree->SetBranchAddress("spill"     , &spill    );
 tree->SetBranchStatus("spillend"  ,1);  tree->SetBranchAddress("spillend"  , &spillend );
 tree->SetBranchStatus("radflag"   ,1);  tree->SetBranchAddress("radflag"   , &radflag  );
 tree->SetBranchStatus("evsize"    ,1);  tree->SetBranchAddress("evsize"    ,  evsize   );
 tree->SetBranchStatus("nt"        ,1);  tree->SetBranchAddress("nt"        , &nt       );
 tree->SetBranchStatus("np"        ,1);  tree->SetBranchAddress("np"        , &np       );
 tree->SetBranchStatus("nhits"     ,1);  tree->SetBranchAddress("nhits"     , &nhits    );
 tree->SetBranchStatus("nshits"    ,1);  tree->SetBranchAddress("nshits"    , &nshits   );
 tree->SetBranchStatus("chi2v"     ,1);  tree->SetBranchAddress("chi2v"     , &chi2v    );
 tree->SetBranchStatus("chi2"      ,1);  tree->SetBranchAddress("chi2"      , &chi2     );
 tree->SetBranchStatus("trmom"     ,1);  tree->SetBranchAddress("trmom"     , &trmom       );
 tree->SetBranchStatus("charge"    ,1);  tree->SetBranchAddress("charge"    , &charge   );
 tree->SetBranchStatus("mom"       ,1);  tree->SetBranchAddress("mom"       , &mom      );
 tree->SetBranchStatus("momh"      ,1);  tree->SetBranchAddress("momh"      , &momh     );
 tree->SetBranchStatus("momv"      ,1);  tree->SetBranchAddress("momv"      , &momv     );
 tree->SetBranchStatus("dmom"      ,1);  tree->SetBranchAddress("dmom"      , &dmom     );
 tree->SetBranchStatus("dmomh"     ,1);  tree->SetBranchAddress("dmomh"     , &dmomh    );
 tree->SetBranchStatus("dmomv"     ,1);  tree->SetBranchAddress("dmomv"     , &dmomv    );
 tree->SetBranchStatus("ft"        ,1);  tree->SetBranchAddress("ft"        , &ft       );
 tree->SetBranchStatus("fl"        ,1);  tree->SetBranchAddress("fl"        , &fl       );
 tree->SetBranchStatus("beta"      ,1);  tree->SetBranchAddress("beta"      , &beta     );
 tree->SetBranchStatus("mass2"     ,1);  tree->SetBranchAddress("mass2"     , &mass2    );
 tree->SetBranchStatus("poshIH"    ,1);  tree->SetBranchAddress("poshIH"    , &poshIH   );
 tree->SetBranchStatus("poshOH"    ,1);  tree->SetBranchAddress("poshOH"    , &poshOH   );
 tree->SetBranchStatus("posvIH"    ,1);  tree->SetBranchAddress("posvIH"    , &posvIH   );
 tree->SetBranchStatus("posvOH"    ,1);  tree->SetBranchAddress("posvOH"    , &posvOH   );
 tree->SetBranchStatus("zIH"       ,1);  tree->SetBranchAddress("zIH"       , &zIH      );
 tree->SetBranchStatus("zOH"       ,1);  tree->SetBranchAddress("zOH"       , &zOH      );
 tree->SetBranchStatus("passIH"    ,1);  tree->SetBranchAddress("passIH"    , &passIH   );
 tree->SetBranchStatus("passOH"    ,1);  tree->SetBranchAddress("passOH"    , &passOH   );
 tree->SetBranchStatus("ihlr"      ,1);  tree->SetBranchAddress("ihlr"      , &ihlr     );
 tree->SetBranchStatus("ihseg"     ,1);  tree->SetBranchAddress("ihseg"     , &ihseg    );
 tree->SetBranchStatus("ihwid"     ,1);  tree->SetBranchAddress("ihwid"     , &ihwid    );
 tree->SetBranchStatus("ihthick"   ,1);  tree->SetBranchAddress("ihthick"   , &ihthick  );
 tree->SetBranchStatus("ihtdc1"    ,1);  tree->SetBranchAddress("ihtdc1"    , &ihtdc1   );
 tree->SetBranchStatus("ihtdc2"    ,1);  tree->SetBranchAddress("ihtdc2"    , &ihtdc2   );
 tree->SetBranchStatus("ihadc1"    ,1);  tree->SetBranchAddress("ihadc1"    , &ihadc1   );
 tree->SetBranchStatus("ihadc2"    ,1);  tree->SetBranchAddress("ihadc2"    , &ihadc2   );
 tree->SetBranchStatus("iht1"      ,1);  tree->SetBranchAddress("iht1"      , &iht1     );
 tree->SetBranchStatus("iht2"      ,1);  tree->SetBranchAddress("iht2"      , &iht2     );
 tree->SetBranchStatus("ihde1"     ,1);  tree->SetBranchAddress("ihde1"     , &ihde1    );
 tree->SetBranchStatus("ihde2"     ,1);  tree->SetBranchAddress("ihde2"     , &ihde2    );
 tree->SetBranchStatus("ihct1"     ,1);  tree->SetBranchAddress("ihct1"     , &ihct1    );
 tree->SetBranchStatus("ihct2"     ,1);  tree->SetBranchAddress("ihct2"     , &ihct2    );
 tree->SetBranchStatus("iht"       ,1);  tree->SetBranchAddress("iht"       , &iht      );
 tree->SetBranchStatus("ihct"      ,1);  tree->SetBranchAddress("ihct"      , &ihct     );
 tree->SetBranchStatus("ihde"      ,1);  tree->SetBranchAddress("ihde"      , &ihde     );
 tree->SetBranchStatus("ihdedx1"   ,1);  tree->SetBranchAddress("ihdedx1"   , &ihdedx1  );
 tree->SetBranchStatus("ihdedx2"   ,1);  tree->SetBranchAddress("ihdedx2"   , &ihdedx2  );
 tree->SetBranchStatus("ihdedx"    ,1);  tree->SetBranchAddress("ihdedx"    , &ihdedx   );
 tree->SetBranchStatus("ihposvt"   ,1);  tree->SetBranchAddress("ihposvt"   , &ihposvt  );
 tree->SetBranchStatus("ihposva"   ,1);  tree->SetBranchAddress("ihposva"   , &ihposva  );

 tree->SetBranchStatus("ohlr"      ,1);  tree->SetBranchAddress("ohlr"      , &ohlr     );
 tree->SetBranchStatus("ohseg"     ,1);  tree->SetBranchAddress("ohseg"     , &ohseg    );
 tree->SetBranchStatus("ohwid"     ,1);  tree->SetBranchAddress("ohwid"     , &ohwid    );
 tree->SetBranchStatus("ohthick"   ,1);  tree->SetBranchAddress("ohthick"   , &ohthick  );
 tree->SetBranchStatus("ohtdc1"    ,1);  tree->SetBranchAddress("ohtdc1"    , &ohtdc1   );
 tree->SetBranchStatus("ohtdc2"    ,1);  tree->SetBranchAddress("ohtdc2"    , &ohtdc2   );
 tree->SetBranchStatus("ohadc1"    ,1);  tree->SetBranchAddress("ohadc1"    , &ohadc1   );
 tree->SetBranchStatus("ohadc2"    ,1);  tree->SetBranchAddress("ohadc2"    , &ohadc2   );
 tree->SetBranchStatus("oht1"      ,1);  tree->SetBranchAddress("oht1"      , &oht1     );
 tree->SetBranchStatus("oht2"      ,1);  tree->SetBranchAddress("oht2"      , &oht2     );
 tree->SetBranchStatus("ohde1"     ,1);  tree->SetBranchAddress("ohde1"     , &ohde1    );
 tree->SetBranchStatus("ohde2"     ,1);  tree->SetBranchAddress("ohde2"     , &ohde2    );
 tree->SetBranchStatus("ohct1"     ,1);  tree->SetBranchAddress("ohct1"     , &ohct1    );
 tree->SetBranchStatus("ohct2"     ,1);  tree->SetBranchAddress("ohct2"     , &ohct2    );
 tree->SetBranchStatus("oht"       ,1);  tree->SetBranchAddress("oht"       , &oht      );
 tree->SetBranchStatus("ohct"      ,1);  tree->SetBranchAddress("ohct"      , &ohct     );
 tree->SetBranchStatus("ohde"      ,1);  tree->SetBranchAddress("ohde"      , &ohde     );
 tree->SetBranchStatus("ohdedx1"   ,1);  tree->SetBranchAddress("ohdedx1"   , &ohdedx1  );
 tree->SetBranchStatus("ohdedx2"   ,1);  tree->SetBranchAddress("ohdedx2"   , &ohdedx2  );
 tree->SetBranchStatus("ohdedx"    ,1);  tree->SetBranchAddress("ohdedx"    , &ohdedx   );
 tree->SetBranchStatus("ohposvt"   ,1);  tree->SetBranchAddress("ohposvt"   , &ohposvt  );
 tree->SetBranchStatus("ohposva"   ,1);  tree->SetBranchAddress("ohposva"   , &ohposva  );

 tree->SetBranchStatus("tbseg"     ,1);  tree->SetBranchAddress("tbseg"     , &tbseg    );
 tree->SetBranchStatus("tbt"       ,1);  tree->SetBranchAddress("tbt"       , &tbt      );
 tree->SetBranchStatus("tbct"      ,1);  tree->SetBranchAddress("tbct"      , &tbct     );
 tree->SetBranchStatus("tbw"       ,1);  tree->SetBranchAddress("tbw"       , &tbw      );
 tree->SetBranchStatus("tbeg"      ,1);  tree->SetBranchAddress("tbeg"      , &tbeg     );
 tree->SetBranchStatus("tbdeg"     ,1);  tree->SetBranchAddress("tbdeg"     , &tbdeg    );

}

