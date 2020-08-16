#include "ReadTxt.h"
//using namespace std;

ReadTxt::ReadTxt()
{
}
ReadTxt::~ReadTxt()
{
}
////////////////
void ReadTxt::ReadFile(std::string ifname,  std::vector<double> &value){ 
  std::ifstream ifs( ifname.c_str() );
  if ( ifs.fail() ) {
    std::cout << "file open fail : " << ifname << std::endl;
    return;
  }
  std::string line;
  std::cout<<"file name :"<<ifname<<std::endl;
  int count=0;
  double tmp;
    while(!ifs.eof()){
      std::getline(ifs,line);
      if(line[0]=='#') continue;
      std::istringstream sline(line);
      sline >>tmp ;
      value.push_back(tmp);
      count++;
    }
  std::cout<<"end of file"<<std::endl;
} 
////////////////
void ReadTxt::ReadFile(std::string ifname,  std::vector<double> &value1,  std::vector<double> &value2,  std::vector<double> &value3,  std::vector<double> &value4){ 
  std::ifstream ifs( ifname.c_str() );
  if ( ifs.fail() ) {
    std::cout << "file open fail : " << ifname << std::endl;
    return;
  }
  std::string line;
  std::cout<<"file name :"<<ifname<<std::endl;
  int count=0;
  double tmp[4];
    while(!ifs.eof()){
      std::getline(ifs,line);
      if(line[0]=='#') continue;
      std::istringstream sline(line);
      sline >>tmp[0]>>tmp[1]>>tmp[2]>>tmp[3] ;
      value1.push_back(tmp[0]);
      value2.push_back(tmp[1]);
      value3.push_back(tmp[2]);
      value4.push_back(tmp[3]);
      count++;
    }
  std::cout<<"end of file"<<std::endl;
} 
////////////////
