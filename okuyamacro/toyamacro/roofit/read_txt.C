#include "ReadTxt.cc"

void read_txt(){
  std::vector<double> imass,time,xpos,ypos;

  ReadTxt *readtxt = new ReadTxt();
  readtxt->ReadFile("txt/Lambda_lifetime.txt",imass,time,xpos,ypos);

  
  TH1F* h_invmass = new TH1F("h_invmass","h_invmass",100,1.0,1.3);

  for(int i=0;i<imass.size();i++){
    //std::cout<<Form("%.01lf, %.01lf, %.01lf, %.01lf", imass[i], time[i], xpos[i],ypos[i])<<std::endl;
    h_invmass ->Fill(imass[i]);

  }

    h_invmass ->Draw();
}
