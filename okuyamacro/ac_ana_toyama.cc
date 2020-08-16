#include <iostream>
#include <iomanip>
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

#include "TRandom.h"
#include "Tree.h"
#include "Setting.h"

const int NCanvas = 4;

class ana_example : public Tree
{
 public:
         ana_example();
        ~ana_example();
  void makehist();
  void savehist(string ifname);
  void loop();
  void draw(); 
  void SetRoot(string ifname); 
  void savecanvas(string PdfFileName); 
  void SetMaxEvent( int N )  { ENumMax = N; }

  private:
    Setting *set;
    string rfile;
    int GetMaxEvent() { return ENumMax; }
    int ENumMax;
    //int ENum;
    TH2F *h2_qq , *h2_qq_wT1 , *h2_qq_wT2 , *h2_qq_wTT;
    TH2F *h2_npe, *h2_npe_wT1, *h2_npe_wT2, *h2_npe_wTT;
    TH1F *h_q[2],*h_t[2];
    TH1F *h_q_wT[2],*h_q_wTT[2];
    TH1F *h_npe[2],*h_npe_wT[2],*h_npe_wTT[2];
    TH2F *h_frame;
    TF1  *f1;


    TFile *file_out;
    TLatex *latex;
    TLine *line;
    TGraph *gr;

    TCanvas *c[NCanvas];

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ana_example::ana_example()
{

  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");
  gROOT->SetBatch(1);

  gStyle->SetOptDate(0);
  gStyle->SetOptFit(1);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetLineWidth(0);
  gStyle->SetOptDate(0);
//  gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleTextColor(1);
  gStyle->SetGridWidth(1);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetNdivisions(510); // tertiary*10000 + secondary*100 + first
  gStyle->SetOptStat("iMen");
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);

  //const Int_t NRGBs = 5;
  //const Int_t NCont = 255;
  //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  //Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  //Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  //Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  //TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  //gStyle->SetNumberContours(NCont);
  //  cout<<"&ACL_TDC:"<<&ACL_TDC  <<endl;
  //  cout<<"&ACL_QDC:"<<&ACL_QDC  <<",  &ACL_TDC + 1:"<<&ACL_TDC+1<<endl;
  //  cout<<"&ACL_Time:"<<&ACL_Time<<",  &ACL_TDC + 2:"<<&ACL_TDC+2<<endl;
  //  cout<<"&ACL_Edep:"<<&ACL_Edep<<",  &ACL_TDC + 4:"<<&ACL_TDC+4<<endl;
      

  for(int i=0;i<NCanvas;i++)c[i]= new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1200,800 );

  latex = new TLatex();
  line  = new TLine();
}
////////////////////////////////////////////////////////////////////////////
ana_example::~ana_example(){
}
////////////////////////////////////////////////////////////////////////////
void ana_example::makehist(){
  cout<<"makehist"<<endl;

  h2_qq      = new TH2F("h2_qq"     ,"h2_qq"     ,1000,0,4000,1000,0,4000);
  h2_qq_wT1  = new TH2F("h2_qq_wT1" ,"h2_qq_wT1" ,1000,0,4000,1000,0,4000);
  h2_qq_wT2  = new TH2F("h2_qq_wT2" ,"h2_qq_wT2" ,1000,0,4000,1000,0,4000);
  h2_qq_wTT  = new TH2F("h2_qq_wTT" ,"h2_qq_wTT" ,1000,0,4000,1000,0,4000);
  set->SetTH2(h2_qq      , "ACL Q vs ACR Q"              , "ACL QDC [ch]", "ACR QDC [ch]");
  set->SetTH2(h2_qq_wT1  , "ACL Q vs ACR Q (w/ ACL TDC)" , "ACL QDC [ch]", "ACR QDC [ch]");
  set->SetTH2(h2_qq_wT2  , "ACL Q vs ACR Q (w/ ACR TDC)" , "ACL QDC [ch]", "ACR QDC [ch]");
  set->SetTH2(h2_qq_wTT  , "ACL Q vs ACR Q (w/ ACLR TDC)", "ACL QDC [ch]", "ACR QDC [ch]");

  h2_npe     = new TH2F("h2_npe"    ,"h2_npe"    ,1000,0, 200,1000,0, 200);
  h2_npe_wT1 = new TH2F("h2_npe_wT1","h2_npe_wT1",1000,0, 200,1000,0, 200);
  h2_npe_wT2 = new TH2F("h2_npe_wT2","h2_npe_wT2",1000,0, 200,1000,0, 200);
  h2_npe_wTT = new TH2F("h2_npe_wTT","h2_npe_wTT",1000,0, 200,1000,0, 200);
  set->SetTH2(h2_npe     , "ACL NPE vs ACR NPE"              , "ACL PE ", "ACR PE");
  set->SetTH2(h2_npe_wT1 , "ACL NPE vs ACR NPE (w/ ACL TDC)" , "ACL PE ", "ACR PE");
  set->SetTH2(h2_npe_wT2 , "ACL NPE vs ACR NPE (w/ ACR TDC)" , "ACL PE ", "ACR PE");
  set->SetTH2(h2_npe_wTT , "ACL NPE vs ACR NPE (w/ ACLR TDC)", "ACL PE ", "ACR PE");

  for(int i=0;i<2;i++){
    h_t[i]       = new TH1F(Form("h_t%d"       ,i+1),Form("h_t%d"      ,i+1), 1000,1,4001);
    h_q[i]       = new TH1F(Form("h_q%d"       ,i+1),Form("h_q%d"      ,i+1), 1000,0,4000);
    h_q_wT[i]    = new TH1F(Form("h_q_wT%d"    ,i+1),Form("h_q_wT%d"   ,i+1), 1000,0,4000);
    h_q_wTT[i]   = new TH1F(Form("h_q_wTT%d"   ,i+1),Form("h_q_wTT%d"  ,i+1), 1000,0,4000);
    h_npe[i]     = new TH1F(Form("h_npe%d"     ,i+1),Form("h_npe%d"    ,i+1), 1000,0, 200);
    h_npe_wT[i]  = new TH1F(Form("h_npe_wT%d"  ,i+1),Form("h_npe_wT%d" ,i+1), 1000,0, 200);
    h_npe_wTT[i] = new TH1F(Form("h_npe_wTT%d" ,i+1),Form("h_npe_wTT%d",i+1), 1000,0, 200);
    set->SetTH1(h_t[i]      ,Form("AC%d TDC",i+1)                  ,"TDC [ch]","Counts/4ch",1,3001,4);
    set->SetTH1(h_q[i]      ,Form("AC%d QDC",i+1)                  ,"QDC [ch]","Counts/4ch",1,3001,4);
    set->SetTH1(h_q_wT[i]   ,Form("AC%d QDC (w/ AC%d TDC)",i+1,i+1),"QDC [ch]","Counts/4ch",1,3001,2);
    set->SetTH1(h_q_wTT[i]  ,Form("AC%d QDC (w/ ACLR TDC)",i+1)    ,"QDC [ch]","Counts/4ch",1,3001,6);
    set->SetTH1(h_npe[i]      ,Form("AC%d NPE",i+1)                  ,"PE","Counts",1,3001,4);
    set->SetTH1(h_npe_wT[i]   ,Form("AC%d NPE (w/ AC%d TDC)",i+1,i+1),"PE","Counts",1,3001,2);
    set->SetTH1(h_npe_wTT[i]  ,Form("AC%d NPE (w/ ACLR TDC)",i+1)    ,"PE","Counts",1,3001,6);
  }
}
////////////////////////////////////////////////////////////////////////////
void ana_example::loop(){

  for(int n=0;n<ENum;n++){
    if(n%100000==0) cout<<n<<" / "<<ENum<<endl;
    tree->GetEntry(n);
    h_q[0]       ->Fill(ACL_QDC);
    h_q[1]       ->Fill(ACR_QDC);
    h_t[0]       ->Fill(ACL_TDC);
    h_t[1]       ->Fill(ACR_TDC);
    h_npe[0]     ->Fill(ACL_Edep);
    h_npe[1]     ->Fill(ACR_Edep);
    h2_qq        ->Fill(ACL_QDC , ACR_QDC);
    h2_npe       ->Fill(ACL_Edep, ACR_Edep);

    if(ACL_TDC>0){
      h_q_wT[0]    ->Fill(ACL_QDC);
      h_npe_wT[0]  ->Fill(ACL_Edep);
      h2_qq_wT1    ->Fill(ACL_QDC , ACR_QDC);
      h2_npe_wT1   ->Fill(ACL_Edep, ACR_Edep);
    }

    if(ACR_TDC>0){
      h_q_wT[1]    ->Fill(ACR_QDC);
      h_npe_wT[1]  ->Fill(ACR_Edep);
      h2_qq_wT2    ->Fill(ACL_QDC , ACR_QDC);
      h2_npe_wT2   ->Fill(ACL_Edep, ACR_Edep);
    }

    if(ACL_TDC>0 && ACR_TDC>0){
      h_q_wTT[0]    ->Fill(ACL_QDC);
      h_npe_wTT[0]  ->Fill(ACL_Edep);
      h_q_wTT[1]    ->Fill(ACR_QDC);
      h_npe_wTT[1]  ->Fill(ACR_Edep);
      h2_qq_wTT     ->Fill(ACL_QDC , ACR_QDC);
      h2_npe_wTT    ->Fill(ACL_Edep, ACR_Edep);
    }

  }

}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void ana_example::draw(){

  c[0]->Clear();c[0]->Divide(2,4);
  c[0]->cd(1);gPad->SetLogy(1); h_q[0]      ->Draw("");
  c[0]->cd(2);gPad->SetLogy(1); h_q[1]      ->Draw("");
  c[0]->cd(3);gPad->SetLogy(1); h_t[0]      ->Draw("");
  c[0]->cd(4);gPad->SetLogy(1); h_t[1]      ->Draw("");
  c[0]->cd(5);gPad->SetLogy(1); h_npe[0]    ->Draw("");
  c[0]->cd(6);gPad->SetLogy(1); h_npe[1]    ->Draw("");
  c[0]->cd(7);gPad->SetLogz(1); h2_qq       ->Draw("colz");
  c[0]->cd(8);gPad->SetLogz(1); h2_npe      ->Draw("colz");

  c[1]->Clear();c[1]->Divide(2,2);
  c[1]->cd(1);gPad->SetLogy(1);h_q_wT[0]    ->Draw("");
  c[1]->cd(2);gPad->SetLogy(1);h_npe_wT[0]  ->Draw("");
  c[1]->cd(3);gPad->SetLogz(1);h2_qq_wT1    ->Draw("colz");
  c[1]->cd(4);gPad->SetLogz(1);h2_npe_wT1   ->Draw("colz");

  c[2]->Clear();c[2]->Divide(2,2);
  c[2]->cd(1);gPad->SetLogy(1);h_q_wT[1]    ->Draw("");
  c[2]->cd(2);gPad->SetLogy(1);h_npe_wT[1]  ->Draw("");
  c[2]->cd(3);gPad->SetLogz(1);h2_qq_wT2    ->Draw("colz");
  c[2]->cd(4);gPad->SetLogz(1);h2_npe_wT2   ->Draw("colz");

  c[3]->Clear();c[3]->Divide(2,2);
  c[3]->cd(1);gPad->SetLogy(1);h_q_wTT[0]   ->Draw("");        
  c[3]->cd(2);gPad->SetLogy(1);h_npe_wTT[0] ->Draw("");     
  c[3]->cd(3);gPad->SetLogz(1);h_q_wTT[1]   ->Draw("");     
  c[3]->cd(4);gPad->SetLogz(1);h_npe_wTT[1] ->Draw("");     
  c[3]->cd(3);gPad->SetLogz(1);h2_qq_wTT    ->Draw("colz"); 
  c[3]->cd(4);gPad->SetLogz(1);h2_npe_wTT   ->Draw("colz"); 
}
////////////////////////////////////////////////////////////////////////////
void ana_example::savecanvas(string PdfFileName){
  c[0]->Print(Form("%s[",PdfFileName.c_str() ));
  for(int i=0;i<NCanvas;i++)c[i]->Print(Form("%s" ,PdfFileName.c_str() ));
  c[NCanvas-1]->Print(Form("%s]",PdfFileName.c_str() ));
  cout<<"saved : "<<Form("%s",PdfFileName.c_str() )<<endl;
}
////////////////////////////////////////////////
void ana_example::SetRoot(string ifname)
{
  cout<<"SetRoot"<<endl;
  rfile = ifname;
  readtree(rfile);
  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  //ENum = tr->GetEntries();
}
////////////////////////////////////////////////
void ana_example::savehist(string ifname)
{
  cout<<"save histogram to root file"<<endl;
  file_out = new TFile(ifname.c_str() ,"RECREATE");
  file_out->cd();
  //h      ->Write();
  file_out->Close();
}
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "tmp.root";
  string ofname = "tmp.param";
  string ofname_root = "tmp.root";
  string ofname_pdf = "tmp.pdf";
  TString ParamFileName= "default.param";
  int ch;
  int MaxNum = 0;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:w:n:bc:i:p:t"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;
    case 'w':
      ofname = optarg;
      cout<<"output filename : "<<ofname<<endl;
      break;
    case 'n':
      MaxNum = atoi(optarg);
      break;
    case 'b':
      cout<<"BACH MODE!"<<endl;
      break;
    case 'p':
      ParamFileName = optarg;
      break;
    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output root filename"<<endl;
      cout<<"-n : maximum number of events to be analysed "<<endl;
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

  ofname_root = ofname;

  ofname_pdf = ofname;
  ofname_pdf.erase(ofname_pdf.size()-5);
  ofname_pdf.append(".pdf");


  TApplication *theApp = new TApplication("App", &argc, argv);
  ana_example *ana = new ana_example();
  ana->SetMaxEvent(MaxNum);
  ana->SetRoot(ifname);
  ana->makehist();
  ana->loop();
  ana->savehist(ofname_root);
  ana->draw();
  ana->savecanvas(ofname_pdf);
  delete ana;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}
