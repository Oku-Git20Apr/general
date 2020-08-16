#include "ReadCSV.h"
using namespace std;

ReadCSV::ReadCSV()
{
}
ReadCSV::~ReadCSV()
{
}
vector<string> ReadCSV::splitline(string &line){
  char delimiter = ',';
  istringstream stream(line);
  string field;
  vector<string> result;
  while(getline(stream,field,delimiter)){
    result.push_back(field);
  }
  return result;
}
////////////////
void ReadCSV::ReadFile(string ifname,  vector<double> &value){ 
  ifstream ifs( ifname.c_str() );
  if ( ifs.fail() ) {
    cout << "file open fail : " << ifname << endl;
    return;
  }
  string line;
  vector<string> sline;
  cout<<"file name :"<<ifname<<endl;
  int count=0;
  double tmp;
    while(!ifs.eof()){
      getline(ifs,line);
      if(line[0]=='#') continue;
      sline = splitline(line);
      for(int i=0;i<sline.size();i++){
        tmp = atof(sline[i].c_str()) ;
        value.push_back(tmp);
      //  cout<<sline[i]<<endl;
      }
      count++;
    }
  cout<<"end of file"<<endl;
} 
////////////////
