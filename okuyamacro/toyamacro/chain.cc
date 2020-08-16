#include<iostream>
#include<fstream>
#include<math.h>
#include<string>
using namespace std;

#include "TApplication.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"

int main(int argc, char** argv){

  TApplication theApp("App", &argc, argv);
  //TFile f(str);
//  cout<<str<<" will be analyzed "<<endl;
  if(argc!=3){
		cout<<"Usage: ./bin/chain runnum runnum"<<endl;
    exit(1);
  }

  int runnum1 = atoi(argv[1]);
  int runnum2 = atoi(argv[2]);

  if(runnum1>runnum2){
		cout<<"WORNING argv[1] > argv[2]"<<endl;
    exit(1);
  }

  int num = runnum2-runnum1+1;

  TFile *ifp;
  TChain *chain = new TChain("tree");
  for(int i=0;i<num;i++){
    chain->Add(Form("../../root/all/%d.root",runnum1+i));
  }

  TFile *ofp = new TFile(Form("../../root/all/%d_%d.root",runnum1,runnum2),"recreate");


  TTree *tree_out = new TTree();
  tree_out = chain->CloneTree(0);
  int ENum = chain->GetEntries();
  for(int i=0;i<ENum;i++){
    chain->GetEntry(i);
    tree_out->Fill();
    if(i%10000==0){ cout<<i<<" / "<<ENum<<"\r"<<flush; }
  }
  cout<<endl;

//  chain->Merge(Form("rootfiles/output%d_%d.root",runnum1,runnum2),"keep");
  tree_out->AutoSave();
  ofp->Write();
  ofp->Close();

//  theApp.Run();
  return 0;
  exit(1);

}
