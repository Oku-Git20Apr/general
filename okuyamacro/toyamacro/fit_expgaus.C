//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetTF1(TF1 *f, int LColor=2, int LWidth=0, int LStyle=1, int Npx=1000){
  f->SetLineColor(LColor);
  f->SetLineWidth(LWidth);
  f->SetLineStyle(LStyle);
  f->SetNpx(Npx);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expgaus(double *x, double *par) {
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;

// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
  }
//  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));

  return val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expgaus_gaus(double *x, double *par){
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  return expgaus(x,&par[0]) + par[0]*par[4]*TMath::Gaus(x[0],par[3],par[2],1);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double expexpgaus_gaus(double *x, double *par){
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  return expgaus(x,&par[0]) + expgaus(x,&par[4]) + par[8]*TMath::Gaus(x[0],par[10],par[9],1);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void fit_expgaus(){
gStyle->SetPadTopMargin(0.04);
gStyle->SetPadBottomMargin(0.20);
gStyle->SetPadRightMargin(0.04);
gStyle->SetPadLeftMargin(0.07);
//TString ifname1=Form("./root/ana10226_10350.root");//
//TString ifname1=Form("./root/all.root");//
TString ifname1=Form("./root/RKV_ana/ppi_g1to9.root");//
TFile *ifp1  = new TFile(ifname1  );

  double mean = 0.0;
 // double Ltau = 270.;  // Lambda lifetime
  double Ltau = 263.2;  // Lambda lifetime PDG
  double Tau, Sigma, Mean;
  double eTau, eSigma, eMean;

  double Sigma_ga[2], Mean_ga[2];
  double eSigma_ga[2], eMean_ga[2];

  double Sigma_bg, Mean_bg;
  double eSigma_bg, eMean_bg;
  vector<double> vTau, veTau;

  double sigmaBG ;
  double sigmaH3L;
  double sigmaL  ;

 //  sigma = sigmaH3L;

  double fluc = 0.0;//fluctuation of num of event of Lambda and H3Lambda (+/- hoge[ratio])
  int nn = 0, NN = 0;  // No. of Requested events to get statistics
#if 0
  gROOT->SetBatch(1);
  if(N<400)        NN = 4000;
  else if(N<1000)  NN = 10000;
  else if(N<4000)  NN = 20000;
  else if(N<10000) NN = 40000;
  else NN = N;
#endif

  TH1F *h, *h_BG;
  h = (TH1F*)ifp1->Get("cointime/h_ctime_Lam");
  h ->Rebin(5);
  h->SetStats(kFALSE);
  h->SetMinimum(0.1);
  h->SetLineWidth(3);

  h_BG = (TH1F*)ifp1->Get("cointime/h_ctime_ppipi");
  //h_BG = (TH1F*)ifp1->Get("cointime/h_ctime_BG");
  h_BG ->Rebin(1);
  h_BG->SetStats(kFALSE);
  //h_BG->SetMinimum(0.1);
  h_BG->SetLineWidth(3);


  TH1F *h_cor, *h_BG_cor;//w/ flight time correction
  h_cor = (TH1F*)ifp1->Get("cointime/h_ctime_Lam_cor");
  h_cor ->Rebin(4);
  h_cor->SetStats(kFALSE);
  h_cor->SetMinimum(0.1);
  h_cor->SetLineWidth(3);

  h_BG_cor = (TH1F*)ifp1->Get("cointime/h_ctime_ppipi_cor");
  h_BG_cor ->Rebin(1);
  h_BG_cor->SetStats(kFALSE);
  h_BG_cor->SetLineWidth(3);

  TF1 *f, *f_H3L, *f_Lam, *f_bg;
  TF1 *f_ga[2];

//initial parameters
double area  = 200;
double tau   = 0.25;
double sigma = 0.2;
double mean  = 0.;
double areaBG  = 0.5;

 
// Fitting //
    double par[11],par_bg[11];
    f = new TF1("f",expgaus_gaus,-1,3,5);
    f->SetParNames("area","tau","#sigma","mean","ratio B.G.");
    //f->SetParNames("area","tau","#sigma","mean");
    f->SetParameter(0,area); //20 bin num.
    f->SetParameter(1,tau);
    f->SetParameter(2,sigma);
    f->SetParameter(3,mean);
    f->SetParameter(4,areaBG);
    f->FixParameter(2,0.18);
    f->FixParameter(3,0.);
    f->FixParameter(4,0.4);
    //f->SetParLimits(1,0.20,0.39);
    //f->SetParLimits(2,0.18,0.19);
    //f->SetParLimits(3,0.1,0.18);
    //f->FixParameter(3,-0.1);

    for(int i =0;i<2;i++){
      f_ga[i] = new TF1(Form("f_ga%d",i+1),"gaus",-1,3);
      f_ga[i]->SetParameter(0,area); //20 bin num.
      f_ga[i]->SetParameter(1,mean);  f_ga[i]->SetParameter(2,sigma);
      SetTF1(f_ga[i] ,2,3,1,2000);
    }
  
    f_bg = new TF1("f_bg","gaus",-1,3);
    f_bg->SetParNames("area","tau","#sigma","mean");
    f_bg->SetParameter(0,area); //20 bin num.
    f_bg->SetParameter(1,0.1);  f_bg->SetParameter(2,sigma);  f_bg->SetParameter(3,mean);
#if 1
//    f->FixParameter(0,N*20.);
//    f->SetParLimits(0,N*20.-N*fluc*20.,N*20.+N*fluc*20.);
//    f->FixParameter(2,sigmaH3L);  f->FixParameter(3,mean);
//    f->FixParameter(4,N*20.*lratio);
//    f->SetParLimits(4,N*20.-N*lratio*fluc*20.,N*20.+N*lratio*fluc*20.);
//    f->FixParameter(5,Ltau);  f->FixParameter(6,sigmaL);  f->FixParameter(7,mean);
//    f->FixParameter(8,N*20.*bratio);
//    f->FixParameter(9,sigmaBG);  f->SetParameter(10,mean);
//    h->Fit(f,"LL","",-1000,3000);
    h_cor->Fit(f,"0LL","",-1,2);
    h_cor->Fit(f_ga[0],"0QLL","",-1,2);

    h_BG_cor->Fit(f_bg,"0QLL","",-1,2);
    h_BG_cor->Fit(f_ga[1],"0QLL","",-1,2);
  
// GetParameters //
    f->GetParameters(&par[0]);
    Tau  = f->GetParameter(1); Sigma  = f->GetParameter(2); Mean  = f->GetParameter(3);
    eTau = f->GetParError(1);  eSigma = f->GetParError(2);  eMean = f->GetParError(3);
  
    for(int i =0;i<2;i++){
      Sigma_ga[i]  = f_ga[i]->GetParameter(2); Mean_ga[i]  = f_ga[i]->GetParameter(1);
      eSigma_ga[i] = f_ga[i]->GetParError(2);  eMean_ga[i] = f_ga[i]->GetParError(1);
    }
 #endif 

// Drawing //
//           ,tau,N,Ltau,int(N*lratio),sigma,int(N*bratio)),"#it{c#tau} (ps)","counts / 20ps");
  SetTF1(f    ,1,3,1,2000);

  TLatex *t_tau   = new TLatex(0.9,h_cor->GetMaximum()*exp(0.)  ,Form("#tau = %.3lf #pm %.3lf ps",Tau,eTau));
  TLatex *t_sig   = new TLatex(0.9,h_cor->GetMaximum()*exp(-0.7),Form("#sigma = %.3lf #pm %.3lf ps",Sigma,eSigma));
  TLatex *t_chi   = new TLatex(0.9,h_cor->GetMaximum()*exp(-1.5),Form("#chi^{2} = %.1lf/%d = %.3lf",f->GetChisquare(),f->GetNDF(),f->GetChisquare()/f->GetNDF()));
 // t_tau->SetTextSize(0.035);
 // t_chi->SetTextSize(0.035);
  t_tau->SetTextSize(0.05);
  t_sig->SetTextSize(0.05);
  t_chi->SetTextSize(0.05);
  t_tau->SetTextAlign(13);
  t_sig->SetTextAlign(13);
  t_chi->SetTextAlign(13);


  TCanvas *c[6];
  for(int i=0;i<6;i++){
    c[i] =new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1400,700);
    c[i]->Divide(1,1,1E-5,1E-5);
  }
 // TCanvas *c1 =new TCanvas("c1","c1",900,700);
  c[0]->cd(1)->SetMargin(0.10,0.05,0.18,0.07);
  gPad->SetLogy(1);
  //h->SetTitle("");
  h->SetTitleSize(0.07,"x");
  h->SetTitleSize(0.07,"y");
  h->GetYaxis()->SetTitleOffset(0.47);
  h->Draw();
//  f->Draw("same");
  //f_H3L->Draw("same");
  //f_Lam->Draw("same");


  c[1]->cd(1)->SetMargin(0.10,0.05,0.18,0.07);
  gPad->SetLogy(1);
  h->Draw();
//  f_ga[0]->Draw("same");
  //f_H3L->Draw("same");
  //f_Lam->Draw("same");
  //t_tau->Draw();
                     
  c[2]->cd(1)->SetMargin(0.10,0.05,0.18,0.07);
  gPad->SetLogy(1);
  //h_BG->SetTitle("");
  h_BG->SetTitleSize(0.07,"x");
  h_BG->SetTitleSize(0.07,"y");
  h_BG->GetYaxis()->SetTitleOffset(0.47);
  h_BG ->Draw();
//  f_bg->Draw("same");
  //t_tau->DrawLatex(0.9,h_BG->GetMaximum()*exp(0.)  ,Form("#tau = %.3lf #pm %.3lf ps",Tau,eTau));
  //t_sig->DrawLatex(0.9,h_BG->GetMaximum()*exp(-0.7),Form("#sigma = %.3lf #pm %.3lf ps",Sigma,eSigma));
  //t_chi->DrawLatex(0.9,h_BG->GetMaximum()*exp(-1.5),Form("#chi^{2} = %.1lf/%d = %.3lf",f_bg->GetChisquare(),f_bg->GetNDF(),f_bg->GetChisquare()/f_bg->GetNDF()));

  c[3]->cd(1)->SetMargin(0.10,0.05,0.18,0.07);
  gPad->SetLogy(1);
  //h_BG->SetTitle("");
  h_BG->SetTitleSize(0.07,"x");
  h_BG->SetTitleSize(0.07,"y");
  h_BG->GetYaxis()->SetTitleOffset(0.47);
  h_BG ->Draw();
//  f_ga[1]->Draw("same");

  c[4]->cd(1)->SetMargin(0.10,0.05,0.18,0.07);
  gPad->SetLogy(1);
  h_cor->SetTitleSize(0.07,"x");
  h_cor->SetTitleSize(0.07,"y");
  h_cor->GetYaxis()->SetTitleOffset(0.47);
  h_cor->Draw();
  t_tau->Draw();
  t_sig->Draw();
  t_chi->Draw();
  f->Draw("same");
  //t_sig->DrawLatex(0.9,h_cor->GetMaximum()*exp(-0.7),Form("#sigma = %.3lf #pm %.3lf ps",Sigma_ga[0],eSigma_ga[0]));
  //t_chi->DrawLatex(0.9,h_cor->GetMaximum()*exp(-1.5),Form("#chi^{2} = %.1lf/%d = %.3lf",f_ga[0]->GetChisquare(),f_ga[0]->GetNDF(),f_ga[0]->GetChisquare()/f_ga[0]->GetNDF()));

  c[5]->cd(1)->SetMargin(0.10,0.05,0.18,0.07);
    Tau  = f_bg->GetParameter(1); Sigma  = f_bg->GetParameter(2); Mean  = f_bg->GetParameter(3);
    eTau = f_bg->GetParError(1);  eSigma = f_bg->GetParError(2);  eMean = f_bg->GetParError(3);
  gPad->SetLogy(1);
  h_BG_cor->SetTitleSize(0.07,"x");
  h_BG_cor->SetTitleSize(0.07,"y");
  h_BG_cor->GetYaxis()->SetTitleOffset(0.47);
  h_BG_cor ->Draw();

  f_bg ->Draw("same");
  t_sig->DrawLatex(0.9,h_BG_cor->GetMaximum()*exp(-0.7),Form("#sigma = %.3lf #pm %.3lf ps",Sigma_ga[1],eSigma_ga[1]));
  t_chi->DrawLatex(0.9,h_BG_cor->GetMaximum()*exp(-1.5),Form("#chi^{2} = %.1lf/%d = %.3lf",f_ga[1]->GetChisquare(),f_ga[1]->GetNDF(),f_ga[1]->GetChisquare()/f_ga[1]->GetNDF()));

  c[0]->Print("pdf/draw_expgaus.pdf[");
  c[0]->Print("pdf/draw_expgaus.pdf" );
 // c[1]->Print("pdf/draw_expgaus.pdf" );
  c[2]->Print("pdf/draw_expgaus.pdf" );
 // c[3]->Print("pdf/draw_expgaus.pdf" );
  c[4]->Print("pdf/draw_expgaus.pdf" );
  c[5]->Print("pdf/draw_expgaus.pdf" );
  c[5]->Print("pdf/draw_expgaus.pdf]");
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
