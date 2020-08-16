//RKVertex_selectevent.cc  2017.8.3 Y.Toyama
//remake same tree with UserRKVertexAnalysis.cpp after apply additional parameters.
//just for time consuming...

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
using namespace std;

#include "TApplication.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGaxis.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"

#include "RKVertex_tr.h"

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

const double MaxCoinTime =  20.;//IH-Tagger Coincidence Timing[ns]
const double MinCoinTime = -20.;//IH-Tagger Coincidence Timing[ns]
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct TreeBranch{
  int runnum, evnum;             // run info
  int spill, spillend, radflag;  // spill info
  int evsize;

  int    nt, np, nhits, nshits;

  double chi2v[3], chi2[3], trmom[3];

  int charge[3];
  double mom[3];
  double momx[3], momy[3], momz[3];
  double fl_ih[3], fl_oh[3];

  int    ass_nih[3], ihlr[3], ihseg[3];
  double ih_dydx[3], ih_dzdx[3];
  double ihwid[3], ihthick[3];
  int    ihtdc1[3], ihtdc2[3], ihadc1[3], ihadc2[3];
  double iht1  [3], iht2  [3], ihde1 [3], ihde2 [3], ihct1 [3], ihct2[3];
  double ihdedx1[3], ihdedx2[3], ihdedx[3];
  double iht[3], ihct[3], ihde[3];
  double ihposvt[3], ihposva[3];
  double ih_x[3] , ih_y[3] ;
  double ih_nx[3], ih_ny[3];
  double ih_dx[3], ih_dy[3];

  int    ass_noh[3], ohlr[3], ohseg[3];
  double oh_dydx[3], oh_dzdx[3];
  double ohwid[3], ohthick[3];
  int    ohtdc1[3], ohtdc2[3], ohadc1[3], ohadc2[3];
  double oht1  [3], oht2  [3], ohde1 [3], ohde2 [3], ohct1 [3], ohct2[3];
  double ohdedx1[3], ohdedx2[3], ohdedx[3];
  double oht[3], ohct[3], ohde[3];
  double ohposvt[3], ohposva[3];
  double oh_x[3]  , oh_y[3] ;
  double oh_nx[3] , oh_ny[3];
  double oh_dx[3] , oh_dy[3];


  int NumOfTagBHit;
  int tagbtdc[40];
  double tagbctime[40], tagbtime_w[40], tagbtime[40];
  double tagbtime_l[40][10];

  int NumOfTagFHit;
  int tagftdc[160];
  double tagftime[160];

  double vertex[3], pvertex[3], vdist, dlength;
  int vtype;
  int ntr;

// TDL 
  int ntdl, ntdll, ntdlr;    // Multiplicity
  double tdllutime_w[NTDL], tdlrutime_w[NTDL], tdlldtime_w[NTDL],   tdlrdtime_w[NTDL];   // time (width)
  double tdlluctime[NTDL], tdlructime[NTDL], tdlldctime[NTDL],   tdlrdctime[NTDL];   // time (width)
  double tdllmtime[NTDL], tdlrmtime[NTDL];   // time (mean)
  double tdllmctime[NTDL], tdlrmctime[NTDL];   // time (mean)
};

static TreeBranch tr;
//*****************************************************************//
void Initialization(){
  tr.runnum =  tr.evnum = -99;             // run info
  tr.spill =  tr.spillend =  tr.radflag = -9;  // spill info
  tr.evsize = -99;

  tr.nt =  tr.np =  tr.nhits =  tr.nshits = -99;
  tr.vdist =  tr.dlength = -999.;
  tr.vtype = -99;
  tr.ntr   = -99;

  for(int i=0;i<3;i++){
    tr.chi2v[i] =  tr.chi2[i] =  tr.trmom[i];

    tr.charge[i]= -99;
    tr.mom[i]   = -999.;
    tr.momx[i] = tr.momy[i]   = tr.momz[i] = -999.;
    tr.fl_ih[i] =  tr.fl_oh[i] = -999.;

    tr.ass_nih[i] =  tr.ihlr[i] =  tr.ihseg[i] = -99;
    tr.ih_dydx[i] =  tr.ih_dzdx[i] = -999.;
    tr.ihwid[i] =  tr.ihthick[i] = -999.;
    tr.ihtdc1[i] =  tr.ihtdc2[i] =  tr.ihadc1[i] =  tr.ihadc2[i]= -999;
    tr.iht1  [i] =  tr.iht2  [i] =  tr.ihde1 [i] =  tr.ihde2 [i] =  tr.ihct1 [i] =  tr.ihct2[i]=-999.;
    tr.iht[i] =  tr.ihct[i] =  tr.ihde[i] = -999.;
    tr.ihposvt[i] =  tr.ihposva[i] = -999.;
    tr.ih_x[i]  =  tr.ih_y[i]  = -999.;
    tr.ih_nx[i] =  tr.ih_ny[i] = -999.;
    tr.ih_dx[i] =  tr.ih_dy[i] = -999.;

    tr.ass_noh[i] =  tr.ohlr[i] =  tr.ohseg[i]=-99;
    tr.oh_dydx[i] =  tr.oh_dzdx[i] = -999.;
    tr.ohwid[i] =  tr.ohthick[i]=-999.;
    tr.ohtdc1[i] =  tr.ohtdc2[i] =  tr.ohadc1[i] =  tr.ohadc2[i]=-999;
    tr.oht1  [i] =  tr.oht2  [i] =  tr.ohde1 [i] =  tr.ohde2 [i] =  tr.ohct1 [i] =  tr.ohct2[i]=-999.;
    tr.oht[i] =  tr.ohct[i] =  tr.ohde[i]=-999.;
    tr.ohposvt[i] =  tr.ohposva[i]=-999.;
    tr.oh_x[i]  =  tr.oh_y[i]  = -999.;
    tr.oh_nx[i] =  tr.oh_ny[i] = -999.;
    tr.oh_dx[i] =  tr.oh_dy[i] = -999.;

    tr.vertex[i] =  tr.pvertex[i]  =-999.;
  }


  tr.NumOfTagBHit = -100;
  for(int i=0; i<40; i++){
    tr.tagbtdc[i] = -100; tr.tagbtime[i] = tr.tagbtime_w[i] = tr.tagbctime[i] = -1000.;
    for(int j=0;j<10;j++){
      tr.tagbtime_l[i][j]=-1000.;
    }
  }

  tr.NumOfTagFHit = -100;
  for(int i=0; i<160; i++){tr.tagftdc[i] = -100; tr.tagftime[i] = -1000.;}


  tr.vtype = -99;
  tr.ntr   = -99;

//////////////
//TDL(V1290)//
//////////////
  tr.ntdll   = tr.ntdlr   = -9;
  for(int i=0; i<NTDL; i++){
    tr.tdllutime_w[i]= tr.tdlrutime_w[i] = tr.tdlldtime_w[i] = tr.tdlrdtime_w[i] = -9999.;
    tr.tdlluctime[i] = tr.tdlructime[i]  = tr.tdlldctime[i]  = tr.tdlrdctime[i]  = -9999.;
    tr.tdllmtime[i]  = tr.tdlrmtime[i]   = -9999;
    tr.tdllmctime[i] = tr.tdlrmctime[i]  = -9999;
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class RKVertex_selectevent : public RKVertex_tr
{
 public:
         RKVertex_selectevent();
        ~RKVertex_selectevent();
  void loop();
  void SetRoot(string ifname); 
  void SetList(string ifname); 
  void SetBranch(string ofname); 
  void SetMaxEvent( int N )  { ENumMax = N; }

  private:
    TFile *ofp ;
    TTree *tree_out;
    string pdf_name;
    string input_param;
    int GetMaxEvent() { return ENumMax; }
    int ENumMax;
    int ENum;
    int runmin,runmax;

    int run_num;
    double beta[3];

    vector<int> run, event;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
RKVertex_selectevent::RKVertex_selectevent()
{
  run.clear();
  event.clear();
}
////////////////////////////////////////////////////////////////////////////
RKVertex_selectevent::~RKVertex_selectevent(){
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_selectevent::SetList(string ifname){
  input_param = ifname;

  ifstream ifs( ifname.c_str() );
  if(ifs.fail()){ cout<<"file open fail!"<<endl; exit(1); }

  string line;
  int val1,val2;
  while(1){
    getline(ifs,line);
    if(line[0]=='#') continue;
    if( ifs.eof() ) break;
    istringstream sline(line);
    sline >> val1 >> ws;
    sline >> val2 >> ws;

    run.push_back(val1);
    event.push_back(val2);
  }
  
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_selectevent::SetBranch(string ofname){
cout<<"make branch"<<endl;
  ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  tree_out = new TTree("tree","tree");
 
//vertex
  tree_out->Branch("runnum"    , &tr.runnum    , "runnum/I"       );
  tree_out->Branch("evnum"     , &tr.evnum     , "evnum/I"        );
  tree_out->Branch("vtype"     , &tr.vtype     , "VertexType/I"   );
  tree_out->Branch("ntr"       , &tr.ntr       , "NTrack/I"       );
  tree_out->Branch("vertex"    ,  tr.vertex    , "VertexPos[3]/D" );
  tree_out->Branch("pvertex"   ,  tr.pvertex   , "VertexMom[3]/D" );
  tree_out->Branch("vdist"     , &tr.vdist     , "VertexDist/D"   );
  tree_out->Branch("dlength"   , &tr.dlength   , "DecayLength/D"  );

//track
  tree_out->Branch("nt"        , tr.nt        , "nt[3]/I"          );
  tree_out->Branch("np"        , tr.np        , "np[3]/I"          );
  tree_out->Branch("nhits"     , tr.nhits     , "nhits[3]/I"       );
  tree_out->Branch("nshits"    , tr.nshits    , "nshits[3]/I"      );

  tree_out->Branch("chi2v"     , tr.chi2v     , "chi2v[3]/D"       );
  tree_out->Branch("chi2"      , tr.chi2      , "chi2[3]/D"        );
  tree_out->Branch("trmom"     , tr.trmom     , "trmom[3]/D"       );

  tree_out->Branch("charge"    , tr.charge    , "charge[3]/I"      );
  tree_out->Branch("mom"       , tr.mom       , "mom[3]/D"         );
  tree_out->Branch("momx"      , tr.momx      , "momx[3]/D"        );
  tree_out->Branch("momy"      , tr.momy      , "momy[3]/D"        );
  tree_out->Branch("momz"      , tr.momz      , "momz[3]/D"        );

//IH
  tree_out->Branch("ass_nih"   , tr.ass_nih   , "ass_nih[3]/I"     );
  tree_out->Branch("fl_ih"     , tr.fl_ih     , "fl_ih[3]/D"       );
  tree_out->Branch("ihlr"      , tr.ihlr      , "ihlr[3]/I"        );
  tree_out->Branch("ihseg"     , tr.ihseg     , "ihseg[3]/I"       );
  tree_out->Branch("ih_dydx"   , tr.ih_dydx   , "ih_dydx[3]/D"     );
  tree_out->Branch("ih_dzdx"   , tr.ih_dzdx   , "ih_dzdx[3]/D"     );
  tree_out->Branch("ihwid"     , tr.ihwid     , "ihwid[3]/D"       );
  tree_out->Branch("ihthick"   , tr.ihthick   , "ihthick[3]/D"     );
  tree_out->Branch("ihtdc1"    , tr.ihtdc1    , "ihtdc1[3]/I"      );
  tree_out->Branch("ihtdc2"    , tr.ihtdc2    , "ihtdc2[3]/I"      );
  tree_out->Branch("ihadc1"    , tr.ihadc1    , "ihadc1[3]/I"      );
  tree_out->Branch("ihadc2"    , tr.ihadc2    , "ihadc2[3]/I"      );
  tree_out->Branch("iht1"      , tr.iht1      , "iht1[3]/D"        );
  tree_out->Branch("iht2"      , tr.iht2      , "iht2[3]/D"        );
  tree_out->Branch("ihde1"     , tr.ihde1     , "ihde1[3]/D"       );
  tree_out->Branch("ihde2"     , tr.ihde2     , "ihde2[3]/D"       );
  tree_out->Branch("ihdedx1"   , tr.ihdedx1   , "ihdedx1[3]/D"     );
  tree_out->Branch("ihdedx2"   , tr.ihdedx2   , "ihdedx2[3]/D"     );
  tree_out->Branch("ihdedx"    , tr.ihdedx    , "ihdedx[3]/D"      );
  tree_out->Branch("ihct1"     , tr.ihct1     , "ihct1[3]/D"       );
  tree_out->Branch("ihct2"     , tr.ihct2     , "ihct2[3]/D"       );
  tree_out->Branch("iht"       , tr.iht       , "iht[3]/D"         );
  tree_out->Branch("ihct"      , tr.ihct      , "ihct[3]/D"        );
  tree_out->Branch("ihde"      , tr.ihde      , "ihde[3]/D"        );
  tree_out->Branch("ihposvt"   , tr.ihposvt   , "ihposvt[3]/D"     );
  tree_out->Branch("ihposva"   , tr.ihposva   , "ihposva[3]/D"     );
  tree_out->Branch("ih_x"      , tr.ih_x      , "ih_x[3]/D"        );
  tree_out->Branch("ih_y"      , tr.ih_y      , "ih_y[3]/D"        );
  tree_out->Branch("ih_nx"     , tr.ih_nx     , "ih_nx[3]/D"       );
  tree_out->Branch("ih_ny"     , tr.ih_ny     , "ih_ny[3]/D"       );
  tree_out->Branch("ih_dy"     , tr.ih_dy     , "ih_dy[3]/D"       );

//OH
  tree_out->Branch("ass_noh"   , tr.ass_noh   , "ass_noh[3]/I"     );
  tree_out->Branch("fl_oh"     , tr.fl_oh     , "fl_oh[3]/D"       );
  tree_out->Branch("ohlr"      , tr.ohlr      , "ohlr[3]/I"        );
  tree_out->Branch("ohseg"     , tr.ohseg     , "ohseg[3]/I"       );
  tree_out->Branch("oh_dydx"   , tr.oh_dydx   , "oh_dydx[3]/D"     );
  tree_out->Branch("oh_dzdx"   , tr.oh_dzdx   , "oh_dzdx[3]/D"     );
  tree_out->Branch("ohwid"     , tr.ohwid     , "ohwid[3]/D"       );
  tree_out->Branch("ohthick"   , tr.ohthick   , "ohthick[3]/D"     );
  tree_out->Branch("ohtdc1"    , tr.ohtdc1    , "ohtdc1[3]/I"      );
  tree_out->Branch("ohtdc2"    , tr.ohtdc2    , "ohtdc2[3]/I"      );
  tree_out->Branch("ohadc1"    , tr.ohadc1    , "ohadc1[3]/I"      );
  tree_out->Branch("ohadc2"    , tr.ohadc2    , "ohadc2[3]/I"      );
  tree_out->Branch("oht1"      , tr.oht1      , "oht1[3]/D"        );
  tree_out->Branch("oht2"      , tr.oht2      , "oht2[3]/D"        );
  tree_out->Branch("ohde1"     , tr.ohde1     , "ohde1[3]/D"       );
  tree_out->Branch("ohde2"     , tr.ohde2     , "ohde2[3]/D"       );
  tree_out->Branch("ohdedx1"   , tr.ohdedx1   , "ohdedx1[3]/D"     );
  tree_out->Branch("ohdedx2"   , tr.ohdedx2   , "ohdedx2[3]/D"     );
  tree_out->Branch("ohdedx"    , tr.ohdedx    , "ohdedx[3]/D"      );
  tree_out->Branch("ohct1"     , tr.ohct1     , "ohct1[3]/D"       );
  tree_out->Branch("ohct2"     , tr.ohct2     , "ohct2[3]/D"       );
  tree_out->Branch("oht"       , tr.oht       , "oht[3]/D"         );
  tree_out->Branch("ohct"      , tr.ohct      , "ohct[3]/D"        );
  tree_out->Branch("ohde"      , tr.ohde      , "ohde[3]/D"        );
  tree_out->Branch("ohposvt"   , tr.ohposvt   , "ohposvt[3]/D"     );
  tree_out->Branch("ohposva"   , tr.ohposva   , "ohposva[3]/D"     );
  tree_out->Branch("oh_x"      , tr.oh_x      , "oh_x[3]/D"        );
  tree_out->Branch("oh_y"      , tr.oh_y      , "oh_y[3]/D"        );
  tree_out->Branch("oh_nx"     , tr.oh_nx     , "oh_nx[3]/D"       );
  tree_out->Branch("oh_ny"     , tr.oh_ny     , "oh_ny[3]/D"       );
  tree_out->Branch("oh_dx"     , tr.oh_dx     , "oh_dx[3]/D"       );
  tree_out->Branch("oh_dy"     , tr.oh_dy     , "oh_dy[3]/D"       );


//TDL
  tree_out->Branch("ntdll"      ,&tr.ntdll      ,"ntdll/I");
  tree_out->Branch("ntdlr"      ,&tr.ntdlr      ,"ntdlr/I");
  tree_out->Branch("tdllutime_w"  , tr.tdllutime_w  ,Form("tdllutime_w[%d]/D"    ,NTDL));
  tree_out->Branch("tdlrutime_w"  , tr.tdlrutime_w  ,Form("tdlrutime_w[%d]/D"    ,NTDL));
  tree_out->Branch("tdlldtime_w"  , tr.tdlldtime_w  ,Form("tdlldtime_w[%d]/D"    ,NTDL));
  tree_out->Branch("tdlrdtime_w"  , tr.tdlrdtime_w  ,Form("tdlrdtime_w[%d]/D"    ,NTDL));
  tree_out->Branch("tdlluctime"   , tr.tdlluctime   ,Form("tdlluctime[%d]/D"     ,NTDL));
  tree_out->Branch("tdlructime"   , tr.tdlructime   ,Form("tdlructime[%d]/D"     ,NTDL));
  tree_out->Branch("tdlldctime"   , tr.tdlldctime   ,Form("tdlldctime[%d]/D"     ,NTDL));
  tree_out->Branch("tdlrdctime"   , tr.tdlrdctime   ,Form("tdlrdctime[%d]/D"     ,NTDL));
  tree_out->Branch("tdllmtime"    , tr.tdllmtime    ,Form("tdllmtime[%d]/D"      ,NTDL));
  tree_out->Branch("tdlrmtime"    , tr.tdlrmtime    ,Form("tdlrmtime[%d]/D"      ,NTDL));
  tree_out->Branch("tdllmctime"   , tr.tdllmctime   ,Form("tdllmctime[%d]/D"     ,NTDL));
  tree_out->Branch("tdlrmctime"   , tr.tdlrmctime   ,Form("tdlrmctime[%d]/D"     ,NTDL));

//Tagger
  tree_out->Branch("NumOfTagBHit",  &tr.NumOfTagBHit  ,      "NumOfTagBHit/I"       );
  tree_out->Branch("tagbtdc"     ,   tr.tagbtdc       ,      "TagBTDC[40]/I"        );
  tree_out->Branch("tagbtime"    ,   tr.tagbtime      ,      "TagBTime[40]/D"       );
  tree_out->Branch("tagbtime_l"  ,   tr.tagbtime_l    ,      "TagBTime_l[40][10]/D" );
  tree_out->Branch("tagbtime_w"  ,   tr.tagbtime_w    ,      "TagBTime_w[40]/D"     );
  tree_out->Branch("tagbctime"   ,   tr.tagbctime     ,      "TagBCTime[40]/D"      );
  tree_out->Branch("NumOfTagFHit",  &tr.NumOfTagFHit  ,      "NumOfTagFHit/I"       );
  tree_out->Branch("tagftdc"     ,   tr.tagftdc       ,      "TagFTDC[160]/I"       );
  tree_out->Branch("tagftime"    ,   tr.tagftime      ,       "TagFTime[160]/D"     );
}
////////////////////////////////////////////////////////////////////////////
void RKVertex_selectevent::loop(){
cout<<"start loop"<<endl;
runmin=99999;runmax=-99;
int j=0;

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){

    Initialization();

    tree->GetEntry(n);
    if(runnum>runmax)runmax=runnum;
    if(runnum<runmin)runmin=runnum;

    if(j>run.size()){
      break;
    }
    else if((runnum==run[j]&&evnum<event[j])||runnum>run[j]){
      continue;
    }
    else if(runnum<run[j]){
      continue;
    }
    else if((runnum==run[j]&&evnum>event[j])||runnum>run[j]){
      j++;
      continue;
    }
    else if(runnum>run[j]){
      j++;
      continue;
    }

    if(runnum==run[j]&&evnum==event[j]){
      cout<<"event match! run"<<runnum<<"  event"<<evnum<<endl;
    }

    if(n%100000==0) cout<<n<<" / "<<ENum<<endl;

    tr.runnum   =  runnum  ; 
    tr.evnum    =  evnum   ; 
    tr.spill    =  spill   ; 
    tr.spillend =  spillend; 
    tr.radflag  =  radflag ; 
    tr.evsize   =  evsize  ; 

    tr.nt        =  nt      ; 
    tr.np        =  np      ; 
    tr.nhits     =  nhits   ; 
    tr.nshits    =  nshits  ; 
    tr.vdist     =  vdist   ; 
    tr.dlength   =  dlength ; 

    tr.vtype     = vtype;
    tr.ntr       = ntr;

    for(int i=0;i<3;i++){
      tr.vertex[i] =  vertex[i] ;
      tr.pvertex[i]=  pvertex[i];
    }

      //////////////
      //each track//
      //////////////
    for(int i=0;i<3;i++){
      tr.chi2v[i]   = chi2v[i]  ; 
      tr.chi2[i]    = chi2[i]   ; 
      tr.trmom[i]   = trmom[i]  ; 

      tr.charge[i]  = charge[i] ; 
      tr.mom[i]     = mom[i]    ; 
      tr.momx[i]    = momx[i]   ; 
      tr.momy[i]    = momy[i]   ; 
      tr.momz[i]    = momz[i]   ; 
      tr.fl_ih[i]   = fl_ih[i]  ; 
      tr.fl_oh[i]   = fl_oh[i]  ;  
      tr.ass_nih[i] = ass_nih[i]; 
      tr.ass_noh[i] = ass_noh[i]; 

      //* IH *//
      tr.ihlr[i]    = ihlr[i]   ; 
      tr.ihseg[i]   = ihseg[i]  ; 
      tr.ih_dydx[i] = ih_dydx[i]; 
      tr.ih_dzdx[i] = ih_dzdx[i];  
      tr.ihwid[i]   = ihwid[i]  ; 
      tr.ihthick[i] = ihthick[i]; 
      tr.ihtdc1[i]  = ihtdc1[i] ; 
      tr.ihtdc2[i]  = ihtdc2[i] ; 
      tr.ihadc1[i]  = ihadc1[i] ; 
      tr.ihadc2[i]  = ihadc2[i] ;  
      tr.iht1[i]    = iht1[i]   ; 
      tr.iht2[i]    = iht2[i]   ; 
      tr.ihde1[i]   = ihde1[i]  ; 
      tr.ihde2[i]   = ihde2[i]  ; 
      tr.ihdedx1[i] = ihdedx1[i]; 
      tr.ihdedx2[i] = ihdedx2[i]; 
      tr.ihdedx[i]  = ihdedx[i] ; 
      if(ihseg[i]>0){
        tr.ihct1[i]   = ihct1[i]; 
        tr.ihct2[i]   = ihct2[i]; 
        tr.ihct[i]    = ihct[i]; 
      }
      else{
        tr.ihct1[i]   = ihct1[i]; 
        tr.ihct2[i]   = ihct2[i]; 
        tr.ihct[i]    = ihct[i] ; 
      }
      tr.iht[i]     = iht[i]    ; 
      tr.ihde[i]    = ihde[i]   ; 
      tr.ihposvt[i] = ihposvt[i]; 
      tr.ihposva[i] = ihposva[i];  
      tr.ih_x[i]    = ih_x[i]   ; 
      tr.ih_y[i]    = ih_y[i]   ; 
      tr.ih_nx[i]   = ih_nx[i]  ; 
      tr.ih_ny[i]   = ih_ny[i]  ;  
      tr.ih_dx[i]   = ih_dx[i]  ; 
      tr.ih_dy[i]   = ih_dy[i]  ; 

      //* OH *//
      tr.ohlr[i]    = ohlr[i]   ; 
      tr.ohseg[i]   = ohseg[i]  ; 
      tr.oh_dydx[i] = oh_dydx[i]; 
      tr.oh_dzdx[i] = oh_dzdx[i]; 
      tr.ohwid[i]   = ohwid[i]  ; 
      tr.ohthick[i] = ohthick[i]; 
      tr.ohtdc1[i]  = ohtdc1[i] ; 
      tr.ohtdc2[i]  = ohtdc2[i] ; 
      tr.ohadc1[i]  = ohadc1[i] ; 
      tr.ohadc2[i]  = ohadc2[i] ; 
      tr.oht1[i]    = oht1[i]   ; 
      tr.oht2[i]    = oht2[i]   ; 
      tr.ohde1[i]   = ohde1[i]  ; 
      tr.ohde2[i]   = ohde2[i]  ; 
      tr.ohdedx1[i] = ohdedx1[i]; 
      tr.ohdedx2[i] = ohdedx2[i]; 
      tr.ohdedx[i]  = ohdedx[i] ; 
      if(ohseg[i]>0){
        tr.ohct1[i]   = ohct1[i]; 
        tr.ohct2[i]   = ohct2[i]; 
        tr.ohct[i]    = ohct[i]; 
      }
      else{
        tr.ohct1[i]   = ohct1[i];
        tr.ohct2[i]   = ohct2[i];
        tr.ohct[i]    = ohct[i] ;
      }
      tr.oht[i]     = oht[i]    ; 
      tr.ohde[i]    = ohde[i]   ; 
      tr.ohposvt[i] = ohposvt[i]; 
      tr.ohposva[i] = ohposva[i]; 
      tr.oh_x[i]    = oh_x[i]   ; 
      tr.oh_y[i]    = oh_y[i]   ; 
      tr.oh_nx[i]   = oh_nx[i]  ; 
      tr.oh_ny[i]   = oh_ny[i]  ; 
      tr.oh_dx[i]   = oh_dx[i]  ; 
      tr.oh_dy[i]   = oh_dy[i]  ; 

    }//for each track


    //////////////
    //Tag info.//
    //////////////
    tr.NumOfTagBHit  = NumOfTagBHit;
    tr.NumOfTagFHit  = NumOfTagFHit;
    for(int i=0;i<40;i++){
      tr.tagbtdc[i]     = tagbtdc[i]   ;
      tr.tagbtime[i]    = tagbtime[i] ;
      tr.tagbctime[i]   = tagbctime[i];
      tr.tagbtime_w[i]  = tagbtime_w[i];
      for(int j=0;j<10;j++){
        tr.tagbtime_l[i][j] = tagbtime_l[i][j];
      }
    }
    for(int i=0;i<160;i++){
      tr.tagftdc[i]       = tagftdc[i]  ;
      tr.tagftime[i]      = tagftime[i] ;
    }

    ///////
    //TDL//
    ///////
    tr.ntdll     = ntdll;
    tr.ntdlr     = ntdlr;
    for(int i=0;i<NTDL;i++){
      tr.tdllutime_w[i] = tdllutime_w[i];
      tr.tdlrutime_w[i] = tdlrutime_w[i];
      tr.tdlldtime_w[i] = tdlldtime_w[i];
      tr.tdlrdtime_w[i] = tdlrdtime_w[i];
      tr.tdlluctime[i]  = tdlluctime[i];
      tr.tdlructime[i]  = tdlructime[i];
      tr.tdlldctime[i]  = tdlldctime[i];
      tr.tdlrdctime[i]  = tdlrdctime[i];
      tr.tdllmtime[i]   = tdllmtime[i] ;
      tr.tdlrmtime[i]   = tdlrmtime[i] ;
      tr.tdllmctime[i]  = tdllmctime[i];
      tr.tdlrmctime[i]  = tdlrmctime[i];
    }



    tree_out ->Fill();

  }//event loop

  ofp->Write();
}
////////////////////////////////////////////////
void RKVertex_selectevent::SetRoot(string ifname)
{
  cout<<"SetRoot"<<endl;
  std::ifstream fp(ifname.c_str());
  if(fp.fail()){ cout<<"file open fail!"<<endl; exit(1); }
  vector<string> runname;
  string line, aa;
  while(1){
    getline(fp,line);
    if(line[0]=='#') continue;
    if( fp.eof() ) break;
    istringstream sline(line);
    sline >> aa >> ws;
    runname.push_back(line);
  }

  cout<<"Run list"<<endl;
  for(int i=0;i<(int)runname.size();i++){
    cout<<runname[i]<<endl;
    add(Form("%s",runname[i].c_str()));
  }

   readtree();
  ENum = GetEntries();
}
////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "output0000.dat";
  string ofname = "root/hoge.root";
  string pdfname = "pdf/track/hoge.pdf";
  string input_param  = "param/default.param";
  int ch;
  int MaxNum = 0;
  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = true;
  bool coin_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:w:n:bcop:i:"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;
    case 'w':
      output_flag = true;
      draw_flag = false;
      ofname = optarg;
      cout<<"output filename : "<<ofname<<endl;
      break;
    case 'n':
      MaxNum = atoi(optarg);
      break;
    case 'c':
      coin_flag = true;
      break;
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
    case 'p':
      pdfname = optarg;
      break;
    case 'i':
      input_param= optarg;
      break;
    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output root filename"<<endl;
      cout<<"-n : maximum number of analysed events"<<endl;
      cout<<"-i : input event list file"<<endl;
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

  TApplication *theApp = new TApplication("App", &argc, argv);
  RKVertex_selectevent *remake = new RKVertex_selectevent();

  remake->SetMaxEvent(MaxNum);
  remake->SetList(input_param);
  remake->SetRoot(ifname);
  remake->SetBranch(ofname);
  remake->loop();
  delete remake;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}

