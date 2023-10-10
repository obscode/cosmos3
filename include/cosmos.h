#ifndef  INCLUDE_COSMOS_H
#define  INCLUDE_COSMOS_H

#include "cosdat.h"

int   OpenParm(char*);
int   ReadParm_s(char*,char*);
int   ReadParm_i(char*,int*);
int   ReadParm_b(char*,int*);
int   ReadParm_r(char*,float*);
int   Readcdf(char*);
int   ReadDistor(char*);
int   fp2ccd(float,float,float*,float*,float);
int   ccd8(int,float,float,float*,float*);
int   Readcof(char*);
void  GetCof(float*, float*, float*, float*);
void  SetCof(float, float, float, float);
int   GetChipdat_f(char*, float*);
int   GetChipdat_i(char*, int*);
int   ReadGain(char*,char*,float*);
void  order(int*, float*, int*, float*);
void  dorder(int*, double*, int*, double*);
float polyvalue(float,float[],int);
void  polyval(float,float[],int,float*);
void  Getchipdat(dewdat*);
int   ReadObsDef(char*, obsdef*);
int   GetParam_s(char*, char*, char*);
int   GetParam_i(char*, char*, int*);
int   GetParam_f(char*, char*, float*);
void  plyfit(float*,float*,int,int,float*,double*,double*,double**);
void  plyfit_w(float*,float*,float*,int,int,float*,double*,double*,double**);
int   OpenCosParm(char*);
void  die(char*);
void  die2(char*,int);
void  fits_die(char*,int);
void  addbar(char*);
void  subbars(char*);
int   mspos(int,float,float,float*,float*);
void subbias(float**,long[],int,int,int,int,int,int);
void subbias_i(int**,long[],int,int,int,int,int,int);
float f_interpol(float**,long[],float,float);
float i_interpol(int**,long[],float,float);
float e_interpol(float**,float**,long[],float,float);
int e2v1( double*, double*);
int e2v2(double*, double*);
float   zernike(float,float,int,int*, float*);

#ifdef INCLUDE_OPTUTILS_H
int   SetupInstr(obsdef*,element**);
int   ReadSMFfile(char*,Obs**);
int   SetupCamera(obsdef*);
#endif

#ifdef _FITSIO_H
int   OpenFitsFile(char*,fitsfile**,fitsdef*);
int   ReadFitsFile(fitsfile**,fitsdef*,int,int*);
#endif

#endif
