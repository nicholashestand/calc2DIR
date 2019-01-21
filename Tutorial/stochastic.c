#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "randomlib.h"
#include "stochastic.h"

int main(int argc,char *argv[]){
  int length;
  float w1,wa,wao;
  float time,sigma;
  float deltat,w0;
  float a,b;
  float mux,muy,muz;
  int i;
  FILE *E_FH,*mu_FH;

  length=atoi(argv[1]);
  deltat=atof(argv[2]);
  sigma=atof(argv[3]);
  time=atof(argv[4]);
  w0=atof(argv[5]);

  a=exp(-deltat/time);
  b=sqrt(1-a*a);

  RandomInitialise(2511,1974);
  wao=RandomGaussian(0,sigma);

  E_FH=fopen("Energy.bin","wb");
  mu_FH=fopen("Dipole.bin","wb");

  for(i=0;i<length;i++){
    wa=wao*a+b*RandomGaussian(0,sigma);
    w1=wa+w0;
    fwrite(&i,sizeof(i),1,E_FH);
    fwrite(&w1,sizeof(w1),1,E_FH);
    fwrite(&i,sizeof(i),1,mu_FH);
    mux=1.; fwrite(&mux,sizeof(mux),1,mu_FH);
    muy=1.; fwrite(&muy,sizeof(muy),1,mu_FH);
    muz=1.; fwrite(&muz,sizeof(muz),1,mu_FH);
    wao=wa;
  } 

  return 0;
}
