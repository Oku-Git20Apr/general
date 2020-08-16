///////////////
//
// TDL geometry file create
//
// each scintillator size 60h x 20w x 0.5t [cm3]
//
//
//# LR CID SEG        X          Y        Z    Width  Thick  ZLeng   Angle   Tofs  sigma
//
//
///////////////
#include <math.h>
#include <fstream>
#include <iomanip>
using namespace std;
void tdl_geo(){
ofstream ofs("tdl_geo.param");
ofs<<" # LR CID SEG        X          Y        Z    Width  Thick  ZLeng   Angle   Tofs  sigma "<<endl;
  for(int j=0;j<2;j++){ //L R
   for(int i=0;i<10;i++){
     double r=8.25;//radius [cm]
     double x;//x axis is beam direction[cm]
     double y;//y axis LR direction (L + side, R - side)[cm]
     double z = 0.;//z axis vertical direction[cm]
     double w = 2.02;//width
     double mergine = 0.3;//mergine for width
     double t = 0.5;//thickness
     double zl= 6.;//length
     double angle=2.*(0.5-(double)j)*w/r*((double)i +0.5 );
     double tof=0.;
     double sigma=0.;//[cm] I don't know what this value means.
     w += mergine;

     x=r*cos(angle);
     y=r*sin(angle);
     //if(j==1)y*=-1.;
  
     //ofs<<"   "<<j<<"   "<<1<<"   "<<i+1<<"    "<<Form("%.03lf",x)<<"    "<<Form("%.04lf",y) <<"    "
     //<<Form("%.04lf",z) <<"    "<<Form("%.03lf",w)    <<"    "<<Form("%.03lf",t)  <<"    "
     //<<Form("%.03lf",zl)<<"    "<<Form("%.03lf",angle*rad_to_deg)<<"    "<<Form("%.02lf",tof)<<"    "<<Form("%.02lf",sigma)<<endl;
  
     ofs<<std::setw(4)<<j
        <<std::setw(4)<<15
        <<std::setw(4)<<i+1
        <<std::setw(10)<<std::setprecision(6)<<x
        <<std::setw(10)<<std::setprecision(7)<<y
        <<std::setw(9) <<std::setprecision(7)<<z
        <<std::setw(9) <<std::setprecision(7)<<w
        <<std::setw(7) <<std::setprecision(5)<<t
        <<std::setw(7) <<std::setprecision(6)<<zl
        <<std::setw(7) <<std::setprecision(4)<<angle*rad_to_deg
        <<std::setw(7) <<std::setprecision(4)<<tof
        <<std::setw(7) <<std::setprecision(4)<<sigma
        <<endl;
   }
  }



}
