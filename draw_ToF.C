#include "Setting.cc"
const int nCanvas = 5;
//////////////////////ooooooooOOOOOOOooooooo/////////////
void draw_ToF(string rname = "../root/r0000024.root"){
  TFile *ifp  = new TFile(rname.c_str(),"READONLY");
  TTree *tree = (TTree*)ifp->Get("tree");
  Setting *set = new Setting();

  TLine *line = new TLine(); 
  line ->SetLineColor(6);
  line ->SetLineWidth(1);
  line ->SetLineStyle(1);

  int tdc[4],qdc[4];
  tree->SetBranchAddress("tdc", tdc);
  tree->SetBranchAddress("qdc", qdc);

  int ENum = tree->GetEntries();

  TH1D *h_tof,*h_tof_cut;
  h_tof      = new TH1D("h_tof"    ,"h_tof"    ,4000,-100,100);
  h_tof_cut  = new TH1D("h_tof_cut","h_tof_cut",4000,-100,100);
  set->SetTH1(h_tof    ,"ToF w/o Q cut","ToF [ns]","Counts/50ps",1,1,0);
  set->SetTH1(h_tof_cut,"ToF w/  Q cut","ToF [ns]","Counts/50ps",1,1,0);
  int q_min[4],q_max[4];
  q_min[0] = 400; q_max[0] = 700;
  q_min[1] = 400; q_max[1] = 700;
  q_min[2] = 550; q_max[2] =1000;
  q_min[3] = 400; q_max[3] = 700;
  

  TH1D *h_qdc[4], *h_tdc[4];
  TH2D  *h2_q_tof[4], *h2_q_tof_cut[4];
  for(int i=0;i<4;i++){
    h_qdc[i] = new TH1D(Form("h_qdc%d",i),Form("h_qdc%d",i),3000,0,3000);
    h_tdc[i] = new TH1D(Form("h_tdc%d",i),Form("h_tdc%d",i),2000,0,4000);
    h2_q_tof[i] = new TH2D(Form("h2_q%d_tof",i),Form("h2_q%d_tof",i), 3000,0,3000,4000,-140,140);
    h2_q_tof_cut[i] = new TH2D(Form("h2_q%d_tof_cut",i),Form("h2_q%d_tof_cut",i), 3000,0,3000,4000,-140,140);
    set->SetTH1(h_qdc[i],Form("QDC KT%d",i+1),"QDC [ch]","Counts/ch",2,1,0);
    set->SetTH1(h_tdc[i],Form("TDC KT%d",i+1),"TDC [ch]","Counts/2ch",4,1,0);
    set->SetTH2(h2_q_tof[i]    ,Form("TDC KT%d",i+1),"QDC [ch]","ToF [ns]");
    set->SetTH2(h2_q_tof_cut[i],Form("TDC KT%d",i+1),"QDC [ch]","ToF [ns]");
  }

  for(int n=0;n<ENum;n++){
    bool q_flag = false;
    if(n%1000==0)cout<<n<<" / "<<ENum<<endl;
    tree->GetEntry(n);

    if(tdc[0]>0&&tdc[1]>0&&tdc[2]>0&&tdc[3]>0){ // w/ all TDC ch has a hit
      double tof = 0.035*( 0.5*(tdc[0] + tdc[1]) - 0.5*(tdc[2] + tdc[3]) );//[ns]
      h_tof->Fill(tof);

      if(qdc[0]>q_min[0] && qdc[0]<q_max[0] && 
         qdc[1]>q_min[1] && qdc[1]<q_max[1] && 
         qdc[2]>q_min[2] && qdc[2]<q_max[2] && 
         qdc[3]>q_min[3] && qdc[3]<q_max[3] )q_flag=true;
      if(q_flag)h_tof_cut->Fill(tof);

      for(int i=0;i<4;i++){
        h_qdc[i] ->Fill(qdc[i]);
        h_tdc[i] ->Fill(tdc[i]);
        h2_q_tof[i] ->Fill(qdc[i],tof);
        if(q_flag)h2_q_tof_cut[i] ->Fill(qdc[i],tof);
      }
    }
    
  }

  double tof_min = h_tof->GetBinCenter(h_tof->GetMaximumBin())-3;
  double tof_max = h_tof->GetBinCenter(h_tof->GetMaximumBin())+3;
  h_tof    ->GetXaxis()->SetRangeUser(tof_min,tof_max);
  h_tof_cut->GetXaxis()->SetRangeUser(tof_min,tof_max);
  for(int i=0;i<4;i++){
    h2_q_tof[i]->GetYaxis()->SetRangeUser(tof_min,tof_max);
    h2_q_tof_cut[i]->GetYaxis()->SetRangeUser(tof_min,tof_max);
  }

  TF1 *ga = new TF1("ga","gausn",tof_min,tof_max);
  ga->SetLineColor(2);
  ga->SetNpx(2000);
  for(int i=0;i<5;i++)h_tof_cut->Fit(ga);

  TCanvas *c[nCanvas];
  for(int i=0;i<nCanvas;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),800,800);
  }
  c[0]->Clear();
  c[0]->Divide(2,2);
  for(int i=0;i<4;i++){
    c[0]->cd(i+1);h_qdc[i]->Draw();
  }
  c[1]->Clear();
  c[1]->Divide(2,2);
  for(int i=0;i<4;i++){
    c[1]->cd(i+1);h_tdc[i]->Draw();
  }
  c[2]->Clear();
  c[2]->Divide(2,2);
  for(int i=0;i<4;i++){
    c[2]->cd(i+1);h2_q_tof[i]->Draw("colz");line->DrawLine(q_max[i],tof_min,q_max[i],tof_max);line->DrawLine(q_min[i],tof_min,q_min[i],tof_max);
  }

  c[3]->Clear();
  c[3]->Divide(2,2);
  for(int i=0;i<4;i++){
    c[3]->cd(i+1);h2_q_tof_cut[i]->Draw("colz");line->DrawLine(q_max[i],tof_min,q_max[i],tof_max);line->DrawLine(q_min[i],tof_min,q_min[i],tof_max);
  }

  c[4]->Clear();
  c[4]->Divide(1,2);
  c[4]->cd(1);h_tof    ->Draw("");
  c[4]->cd(2);h_tof_cut->Draw("");

  string pdf_name=rname;
  pdf_name.erase(pdf_name.size()-5,5);
  pdf_name.erase(0,8);
  pdf_name.append(".pdf");
  //cout<<pdf_name<<endl;
  //cout<<Form("%s[",pdf_name.c_str())<<endl;
  c[0]->Print(Form("%s[",pdf_name.c_str()));
  for(int i=0;i<nCanvas;i++)c[i]->Print(Form("%s",pdf_name.c_str()));
  c[nCanvas-1]->Print(Form("%s]",pdf_name.c_str()));
}
