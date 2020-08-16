#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
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
  RooPlot *tframe2 = dtime.frame(Title("decay time"));
  RooPlot *tframe3 = dtime.frame(Title("decay time"));
  RooPlot *fframe = dtime.frame(Title("function"));
  data.plotOn(tframe);
  data.plotOn(tframe2);
  data.plotOn(tframe3);

  RooRealVar ag("ag","amp gaus"  ,1000);
  RooRealVar mg("mg","mean gaus" ,0.00);
  RooRealVar sg("sg","sigma gaus",0.185);
  RooRealVar tau("tau","tau", 0.25, 0.0, 10);
  //RooGaussian   gaus("gaus","gaus",dtime,mg,sg);
  //RooGenericPdf exp("exp","exp","TMath::Exp(-1.*dtime/tau)",RooArgSet(dtime,tau));
  //RooExponential exp("exp","exp",dtime,tau);
  //RooFFTConvPdf expgaus("expgaus","expgaus",dtime,gaus,exp) ;
  
  //RooGExpModel expgaus("expgaus","exp(x)gaus",dtime,mg,sg,tau);
  
  RooRealVar bias1("bias1","bias1",0) ;
  RooRealVar sigma1("sigma1","sigma1",0.185) ;
  RooGaussModel gm1("gm1","gauss model 1",dtime,bias1,sigma1) ;
  // Construct decay(t) (x) gauss1(t)
  RooDecay expgaus("decay_gm1","decay",dtime,tau,gm1,RooDecay::SingleSided) ;

  RooGaussian gm2("gm2","gauss model 2",dtime,bias1,sigma1) ;
  double in_c11 = 0.5;
  double in_c12 = 0.1;
  double in_c13 = 0.8;
  RooRealVar c11("c11","c11",in_c11);
  RooRealVar c12("c12","c12",in_c12);
  RooRealVar c13("c13","c13",in_c13);
  //RooRealVar c11("c11","c11",in_c11,lo_c11,hi_c11);
  RooAddPdf total1("total1","total1",RooArgList(gm2,expgaus),c11);
  RooAddPdf total2("total2","total2",RooArgList(gm2,expgaus),c12);
  RooAddPdf total3("total3","total3",RooArgList(gm2,expgaus),c13);

  //RooNumConvPdf expgaus("expgaus","expgaus",dtime,exp,gaus) ;
  //expgaus.setConvolutionWindow(mg,sg,6.);
  

  total1.fitTo(data);
  std::cout<<"*************"<<std::endl;
  std::cout<<"============="<<std::endl;
  std::cout<<"total1 tau: "<<tau<<std::endl;
  std::cout<<"============="<<std::endl;
  std::cout<<"*************"<<std::endl;
  total2.fitTo(data);
  std::cout<<"*************"<<std::endl;
  std::cout<<"============="<<std::endl;
  std::cout<<"total2 tau: "<<tau<<std::endl;
  std::cout<<"============="<<std::endl;
  std::cout<<"*************"<<std::endl;
  total3.fitTo(data);
  std::cout<<"*************"<<std::endl;
  std::cout<<"============="<<std::endl;
  std::cout<<"total3 tau: "<<tau<<std::endl;
  std::cout<<"============="<<std::endl;
  std::cout<<"*************"<<std::endl;

  total1.plotOn(tframe);
  total2.plotOn(tframe2);
  total3.plotOn(tframe3);
  //gm2.plotOn(tframe);
  //expgaus.plotOn(tframe);
  total1.plotOn(fframe);

  TCanvas* c[2];
  for(int i=0;i<2;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),800,800) ;
  }
  c[0]->Divide(2,2);
  c[0]->cd(1);
  tframe->Draw();
  c[0]->cd(2);
  fframe->Draw();
  c[0]->cd(3);
  tframe2->Draw();
  c[0]->cd(4);
  tframe3->Draw();

  c[1]->Divide(2,2);
  c[1]->cd(1);
  gPad->SetLogy(1);
  tframe->Draw();
  c[1]->cd(2);
  gPad->SetLogy(1);
  fframe->Draw();
  c[1]->cd(3);
  gPad->SetLogy(1);
  tframe2->Draw();
  c[1]->cd(4);
  gPad->SetLogy(1);
  tframe3->Draw();
}
