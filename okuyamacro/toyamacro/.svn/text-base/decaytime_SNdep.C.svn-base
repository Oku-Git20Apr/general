using namespace std;
#include "Settings.cc"
#define pol1fit

void SetTH1org(TH1F *h1, TString hname, TString xname, TString yname, int LColor, int FStyle, int FColor){
  h1->SetTitle(hname);
  h1->GetXaxis()->SetTitle(xname);
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->SetTitle(yname);
  h1->GetYaxis()->CenterTitle();
  h1->SetMinimum(0.0);
  h1->SetLineColor(LColor);
  h1->SetFillStyle(FStyle);
  h1->SetFillColor(FColor);
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->SetLineWidth(0);
  h1->SetTitleSize(0.05,"");
  h1->SetTitleSize(0.05,"x");
  h1->SetTitleSize(0.05,"y");
  h1->SetLabelSize(0.05,"x");
  h1->SetLabelSize(0.05,"y");
  ((TGaxis*)h1->GetYaxis())->SetMaxDigits(3);
  ((TGaxis*)h1->GetXaxis())->SetMaxDigits(3);
}
////////////////
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gaus_pol1bg(double *x, double *par) {
  //par[0]=area (gaus)
  //par[1]=mean (gaus)
  //par[2]=sigma (gaus)
  //par[3]=pol const
  //par[4]=pol tilt
  double val;
  double bg;
  double ga;

  ga=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
  bg =  par[3]+ x[0]*par[4] ;
  val = ga + bg;
  return val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gaus_pol0bg(double *x, double *par) {
  //par[0]=area (gaus)
  //par[1]=mean (gaus)
  //par[2]=sigma (gaus)
  //par[3]=pol const
  //par[4]=pol tilt
  double val;
  double bg;
  double ga;

  ga=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
  bg =  par[3];
  val = ga + bg;
  return val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gaus_pol3bg(double *x, double *par) {
  //par[0]=area (gaus)
  //par[1]=mean (gaus)
  //par[2]=sigma (gaus)
  //par[3]=pol tilt
  //par[4]=pol tilt
  //par[5]=pol tilt
  //par[6]=pol const
  double val;
  double bg;
  double ga;

  ga=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
  bg =  par[3]*pow(x[0]-par[4],3.) + par[5]*(x[0]-par[4])+par[6];
  val = ga + bg;
  return val;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double get_1st_mom(TH1 *h,double low = 0.,double high = 0.){
  double rel_freq;//soutai dosuu
  double xval,mom=0;
  //int first_bin =h->FindFirstBinAbove(low_th); 
  //int last_bin  =h->FindLastBinAbove(high_th);
  int first_bin =h->FindBin(low);
  int last_bin  =h->FindBin(high);
  int total_event = h->Integral(first_bin,last_bin);
  for(int i=first_bin;i<last_bin;i++){
    rel_freq = (double)h->GetBinContent(i)/total_event;
    xval = h->GetBinCenter(i);
    mom += xval*rel_freq;
  }
  return mom;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//____________________________________________________________________________________________

void decaytime_SNdep(double fmin =1.08, double fmax =1.16, int pol3flag=1 ){
  gStyle->SetOptStat("ie");
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TString ifname1=Form("root/dtime/true.root");
  TFile *ifp1  = new TFile(ifname1  );
  TTree *tree = (TTree*)ifp1->Get("tree");
  Settings* set = new Settings();
  cout<<"input filename : "<<ifname1<<endl;

  TLine *line = new TLine(0,0,0,0);
  line->SetLineColor(1);
  line->SetLineWidth(3);
  line->SetLineStyle(5);

  TArrow *arrow = new TArrow(0,0,0,0,0.01,"");//x1,y1,x2,y2,arrow size,option
  arrow->SetAngle(45);
  arrow->SetLineColor(1);
  arrow->SetLineWidth(3);

  ofstream fout;
  if( fout.is_open() ) fout.close();
      fout.open(Form("dtime_SNdep_result.dat"), ios::out|ios::app);
      fout.setf(ios::left, ios_base::adjustfield);


  int pol =3;
  double MassRange_Lam_min[4]   =  { 1.105, 1.106, 1.110, 1.110};
  double MassRange_Lam_max[4]   =  { 1.125, 1.121, 1.120, 1.116};

  double param_bg[4];
  double param[7], para_error[7];
  double GeVperBin = 0.002;//Gev/bin
  TH1F *h_invm = new TH1F("h_invm","h_invm",150,1.0,1.3);//
  TH1F *h_dtime[4];
  TH1F *h_invm_cut[4];
  SetTH1org(h_invm,"","Invariant Mass[GeV/#it{c}^{2}]",Form("Counts/%.1lfMeV/#it{c}^{2}",1000.*GeVperBin),1,3000,0);
  for(int i=0;i<4;i++){
    h_dtime[i] = new TH1F(Form("h_dtime%d",i+1),Form("h_dtime%d",i+1),160/4,-2.,2.);//
    SetTH1org(h_dtime[i],"Decay time","dtime[ns]","Counts/100ps",1,3001,i+1);
    h_dtime[i] ->SetMinimum(0.8);//

    h_invm_cut[i] = new TH1F(Form("h_invm_cut%d",i+1),Form("h_invm_cut%d",i+1),150,1.,1.3);//
    SetTH1org(h_invm_cut[i],"Invariant Mass","Invariant Mass[GeV/#it{c}^{2}]",Form("Counts/%.1lfMeV/#it{c}^{2}",1000.*GeVperBin)    ,i+1,3000,0);
    h_invm_cut[i] -> SetLineWidth(1);
  }
  h_invm -> SetLineWidth(1);
  h_invm -> SetStats(0);
  h_invm->GetXaxis()->SetRangeUser(1.06,1.26);

  double invm,dtime;
  tree->SetBranchAddress("inv_mass" , &invm);
  tree->SetBranchAddress("dtime"    , &dtime);

  double mom[4]={0.,0.,0.,0.};
  int counter[4]={0,0,0,0};
  for(int n=0;n<tree->GetEntries();n++){
    tree->GetEntry(n);
    h_invm->Fill(invm);
    for(int i=0;i<4;i++){
      if(invm>MassRange_Lam_min[i] && invm<MassRange_Lam_max[i]){
        h_dtime[i]->Fill(dtime);
        h_invm_cut[i]->Fill(invm);
        if(dtime>-0.6){mom[i] += dtime;counter[i]++;}
      }
    }
  }

  cout<<"********"<<endl;
  for(int i=0;i<4;i++){
    //cout<<get_1st_mom(h_dtime[i],-0.7,1.5)<<endl;
    cout<<mom[i]/counter[i]<<" +/- "<<mom[i]/(counter[i]*sqrt(counter[i]))<<endl;
  }
  cout<<"********"<<endl;

  TF1 *f_bg = new TF1("f_bg","[0]*(x-[1])*(x-[1])*(x-[1])+[2]*(x-[1])+[3]",fmin,fmax);
  f_bg -> SetParameter(0,100);
  f_bg -> SetParameter(1,1.3);
  //f_bg -> FixParameter(1,1.3);
  f_bg -> SetParameter(2,-10);
  f_bg -> SetParameter(3,1);
  h_invm ->Fit(f_bg,"0QR","",1.06,1.26);
  f_bg -> GetParameters(&param_bg[0]);


  TF1 *fadd = new TF1("fadd",gaus_pol3bg,fmin,fmax,7);
  fadd -> SetParameter(0,150.);
  fadd -> SetParameter(1,1.115);
  fadd -> SetParameter(2,0.004);
  fadd -> SetParameter(3,param_bg[0]);
  fadd -> SetParameter(4,param_bg[1]);
  fadd -> SetParameter(5,param_bg[2]);
  fadd -> SetParameter(6,param_bg[3]);
  if( pol3flag!=1){fadd ->FixParameter(3,0.);pol=1;}
  fadd ->SetLineWidth(2);  fadd ->SetLineColor(4);fadd ->SetLineStyle(1);fadd ->SetNpx(1000);
  h_invm->Fit(fadd,"0QR","",fmin,fmax);
  fadd ->GetParameters(&param[0]);
  double chi2 = fadd->GetChisquare();
  double ndf  = fadd->GetNDF();
  double g_area  = param[0];
  double g_mean  = param[1];
  double g_sigma = param[2];
  double bg_p0   = param[3];
  double bg_p1   = param[4];
  double bg_p2   = param[5];
  double bg_p3   = param[6];
  //cout<<g_area/GeVperBin<<endl;
  cout<<"BG param(2nd step) ="<<bg_p0<<" "<<bg_p1<<" "<<bg_p2<<" "<<bg_p3<<" "<<endl;//bg_p1<<endl;

  double erg_area =fadd ->GetParError(0);
  double erg_mean =fadd ->GetParError(1);
  double erg_sigma=fadd ->GetParError(2);
  double erbg_p0  =fadd ->GetParError(3);
  //double erbg_p1  =fadd ->GetParError(4);

  TF1 *ga = new TF1("ga","gausn",fmin,fmax);
  ga ->SetLineWidth(2);  ga ->SetLineColor(2);ga ->SetLineStyle(2);
  ga -> SetParameter(0, g_area);
  ga -> SetParameter(1, g_mean);
  ga -> SetParameter(2, g_sigma);


  f_bg->SetParameter(0,bg_p0);
  f_bg->SetParameter(1,bg_p1);
  f_bg->SetParameter(2,bg_p2);
  f_bg->SetParameter(3,bg_p3);

  double x_min=g_mean-3.*g_sigma;
  double x_max=g_mean+3.*g_sigma;

    cout<<setw(16)<<"MassRange_Lam"<<flush;
    cout<<setw(20)<<"NLambda"<<flush;
    cout<<setw(20)<<"mu_Lambda"<<endl;
  for(int i=0;i<4;i++){
    double NLambda = (ga  ->Integral(MassRange_Lam_min[i], MassRange_Lam_max[i]))/GeVperBin;
    double er_NLambda = erg_area*(ga  ->Integral(MassRange_Lam_min[i], MassRange_Lam_max[i]))/g_area   /GeVperBin;
    double bg_3s   = (f_bg->Integral(MassRange_Lam_min[i], MassRange_Lam_max[i]))/GeVperBin;
    cout<<setw(8)<<MassRange_Lam_min[i]<<flush;
    cout<<setw(8)<<MassRange_Lam_max[i]<<flush;
    cout<<setw(8)<<NLambda<<flush;
    cout<<" +/- "<<flush;
    cout<<setw(8)<<er_NLambda<<flush;
    cout<<" "<<flush;
    cout<<setw(8)<<mom[i]/NLambda<<" +/- "<<sqrt(pow(er_NLambda*mom[i]/(NLambda*NLambda),2.) + pow( mom[i]/(counter[i]*sqrt(coun    ter[i])),2.) )<<endl;
    //cout<<setw(8)<<bg_3s<<endl;

    fout<<setw(4)<<pol<<flush;
    fout<<setw(10)<<setprecision(6)<<fmin<<flush;
    fout<<setw(10)<<setprecision(6)<<fmax<<flush;
    fout<<setw(10)<<setprecision(6)<<MassRange_Lam_min[i]<<flush;
    fout<<setw(10)<<setprecision(6)<<MassRange_Lam_max[i]<<flush;
    fout<<setw(10)<<setprecision(6)<<NLambda<<flush;
    fout<<setw(10)<<setprecision(6)<<er_NLambda<<flush;
    fout<<setw(10)<<setprecision(6)<<bg_3s<<flush;
    fout<<setw(10)<<setprecision(6)<<fadd->GetChisquare()/fadd->GetNDF()<<endl;
  }
  double bg_in3s=f_bg->Integral(x_min,x_max)/GeVperBin;

  int p_flag=1;
  //cout<<"Do you want to save canvas?"<<endl;
  //cout<<"yes :1, no :0"<<endl;
  //cin>>p_flag;

////////////

  TH2F *h_frame = new TH2F("h_frame","h_frame",10,0,1,10,0,1);
  SetTH2(h_frame,"","","");
  h_frame->GetXaxis()->SetNdivisions(000);
  h_frame->GetYaxis()->SetNdivisions(000);
  h_frame->SetStats(0);

TCanvas *c[5];
for(int i=0;i<3;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1800/2,1500/2);}
c[3] = new TCanvas(Form("c%d",4),Form("c%d",4),1500,1500/2);
c[4] = new TCanvas(Form("c%d",5),Form("c%d",5),1500,1500/2);
c[0]->Clear();
c[0]->Divide(2,2);
c[0]->cd(1);
gPad->SetLogy(0);h_invm ->Draw("");
c[0]->cd(2);
gPad->SetLogy(0);h_invm ->Draw("PE");f_bg->Draw("same");
c[0]->cd(3);
gPad->SetLogy(0);h_invm ->Draw("");  fadd->Draw("same");ga->Draw("same");f_bg->Draw("same");
c[0]->cd(4);
gPad->SetLogy(0);h_invm ->Draw("PE");fadd->Draw("same");ga->Draw("same");f_bg->Draw("same");

c[1]->Clear();
c[1]->Divide(2,2);
  c[1] -> cd(1);
TH1F *h_zoom = (TH1F*)h_invm->Clone("h_zoom");
h_zoom->GetXaxis()->SetNdivisions(505);
h_zoom->GetXaxis()->SetRangeUser(x_min,x_max);
gPad->SetLogy(0);h_zoom->Draw("PE");fadd->Draw("same");ga->Draw("same");f_bg->Draw("same");
  c[1] -> cd(2);
h_frame->Draw("");
TLatex *tex  = new TLatex(0.5,0.7,"");
tex  -> SetTextSize(0.060);
tex  -> SetTextAlign(22);
                tex->DrawLatex(0.5,0.9,Form("mean: %.01lf#pm%.01lf[MeV/#it{c}^{2}]"          ,g_mean*1000.,erg_mean*1000.                ) );
                tex->DrawLatex(0.5,0.7,Form("Num. of #Lambda : %.01lf#pm%.01lf"              ,g_area/GeVperBin, erg_area/GeVperB    in   ) );
                tex->DrawLatex(0.5,0.5,Form("width(#sigma): %.01lf#pm%.01lf[MeV/#it{c}^{2}]" ,g_sigma*1000.   , erg_sigma*1000.          ) );
                tex->DrawLatex(0.5,0.3,Form("B.G. (#pm 3#sigma): %.01lf"                     ,bg_in3s                                      ) );
                //tex->DrawLatex(0.5,0.1,Form("S/#sqrt{S+N}: %.01lf"                           ,g_area/0.0025/sqrt(g_area/0.0025    +bg_3s)) );
                tex->DrawLatex(0.5,0.1,Form("#chi^{2}/ndf= %.1lf / %d = %.3lf"                           ,fadd->GetChisquare(),f    add->GetNDF(), fadd->GetChisquare()/fadd->GetNDF()));

  c[2]->Clear();
  c[2]->Divide(1,2);
  c[2]->cd(1);
  h_invm->Draw();
    for(int i=1;i<4;i++){
      h_invm_cut[i]->Draw("same");
    }
  c[2]->cd(2);
    gPad->SetLogy(1);
    h_dtime[0]->Draw("");
    for(int i=1;i<4;i++){
      h_dtime[i]->Draw("same");
    }

  c[3]->Clear();
  c[3]->cd(1);
  gPad->SetLogy(0);h_invm ->Draw("PE");fadd->Draw("same");ga->Draw("same");f_bg->Draw("same");

  c[4]->Clear();
  c[4]->cd(1);
  gPad->SetLogy(0);h_invm ->Draw("PE");fadd->Draw("same");ga->Draw("same");f_bg->Draw("same");
  line->DrawLine(MassRange_Lam_min[1],0,MassRange_Lam_min[1],h_invm->GetMaximum());
  line->DrawLine(MassRange_Lam_max[1],0,MassRange_Lam_max[1],h_invm->GetMaximum());
  arrow->DrawArrow(MassRange_Lam_min[1]-0.01, 0.8*h_invm->GetMaximum(),MassRange_Lam_min[1] ,0.8*h_invm->GetMaximum(),0.01,">");
  arrow->DrawArrow(MassRange_Lam_max[1]+0.01, 0.8*h_invm->GetMaximum(),MassRange_Lam_max[1] ,0.8*h_invm->GetMaximum(),0.01,">");

if(p_flag==1){
string ofname_pdf="dtime";
    c[0]->Print(Form("./pdf/dtime/fit_%s_%.0lf_%.0lf_pol%d.pdf[",ofname_pdf.c_str(), fmin*1000, fmax*1000,pol)  );
  for(int i=0;i<5;i++){
    c[i]->Print(Form("./pdf/dtime/fit_%s_%.0lf_%.0lf_pol%d.pdf" ,ofname_pdf.c_str(), fmin*1000, fmax*1000,pol)  );
  }
    c[4]->Print(Form("./pdf/dtime/fit_%s_%.0lf_%.0lf_pol%d.pdf]",ofname_pdf.c_str(), fmin*1000, fmax*1000,pol)  );
 }
}
