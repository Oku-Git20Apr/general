#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXLEN 255

int main( int argc, char **argv ){

  char inFileName[MAXLEN]="countermap.param";
  char outFileName[MAXLEN]="tmp_vmethre.normal";

  int i,j,k;
  int At[16];
  int Thre[16][32];
  FILE *fpin, *fpout;

  char str[MAXLEN];
  int c,n,a,lr,cid,seg,at,ud,th;

  if(argc!=2){
    printf("Input the countermap file name\n");
    return 0;
  }

  strcpy(inFileName, argv[1]);
  printf("Input  File Name : %s \n",inFileName);
  printf("Output File Name : %s \n",outFileName);
  
  for(i=0; i<16; i++){
    for(j=0; j<32; j++){
      Thre[i][j]=4000;
    }
    At[i]=0;
  }
  
  if((fpin=fopen( inFileName, "r" ))==0){
    printf(" File Open Fail : [%s]\n", inFileName );
    return 0;
  }
  if((fpout=fopen( outFileName, "w" ))==0){
    printf(" File Open Fail : [%s]\n", outFileName);
    fclose(fpin);
    return 0;
  }
  
  while(fgets(str,MAXLEN,fpin)!=0){
    if(str[0]!='#'){
      if(sscanf(str, "%d %d %d %d %d %d %d %d %d",
		&c,&n,&a,&lr,&cid,&seg,&at,&ud,&th)==9){
	if(c==2) Thre[n][a]=th;
	if(c==2 && a==0) At[n]=at;
      }
    }
  }
  
  fprintf(fpout, "#CR SLOT SA AT THRESHOLD\n");
  for(i=0; i<16; i++){
    for(j=0; j<32; j++){
      fprintf(fpout, "%d %3d %3d %2d %4d\n",
	      2, i, j, At[i], Thre[i][j]);
      /* cr slot(geoaddr) sa(ch) at threshold */
    }
  }
  
  fclose(fpin);
  fclose(fpout);
  return 0;
}
