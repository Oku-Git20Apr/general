#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"

#include "ReadTxt.cc"
using namespace RooFit ;

void rf_lifetime(){

  std::vector<double> imass,time,vx,vy;
  ReadTxt *readtxt = new ReadTxt();
  readtxt->ReadFile("txt/Lambda_lifetime.txt",imass,time,vx,vy);

  RooRealVar invm("invm","invariant mass",0,2);
  RooRealVar dtime("dtime","decaytime",-2,2);
  RooRealVar xpos("xpos","vertex pos(x)",-20,20);
  RooRealVar ypos("ypos","vertex pos(y)",-20,20);

  //RooDataSet *data = RooDataSet::read("txt/Lambda_lifetime.txt",RooArgList(invm,dtime,xpos,ypos));
  RooDataSet data("data","data",dtime);
  for(int i=0;i<time.size();i++){
    if(i%100==0)std::cout<<i<<" / "<<time.size()<<std::endl;
    dtime = time[i];
    data.add(dtime);
  }

  RooPlot *tframe = dtime.frame(Title("decay time"));
  data.plotOn(tframe);

  TCanvas* c = new TCanvas("c","c",800,400) ;
  tframe->Draw();
}
