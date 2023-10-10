 /*
 *  ReadObsDef   reads obs definition files
 *
 * GetParam[ifs] returns a single parameter value if type =[int float string]from obs definition files
 *
 *  VERSION      10 Sept 03
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cosdat.h"


int ReadObsDef(char filename[], obsdef* obsd){

  FILE   *parmfile;
  char   file[80],line[133],parm[80],value[80];
  int    npar,i,n;
  float  f;
  obsdef obsdata;

  strcpy(file,filename);
  if(strstr(file,".obsdef")==NULL) strcat(file,".obsdef");
  if((parmfile=fopen(file,"r")) == NULL) return 99;
  npar=0;
  obsdata.nshuffle=obsdata.ranod=obsdata.decnod=0;
  obsdata.alignrot=0.0;
  strcpy(obsdata.instrument,"IMACS");
  while(fgets(line,133,parmfile)){
    sscanf(line,"%s",parm);
    if(!strcmp(parm,"INSTRUMENT")){
      sscanf(line,"%s %s",parm,value);
      strcpy(obsdata.instrument,value);
      npar++;
      continue;}
    if(!strcmp(parm,"GR_ORDER")){
      sscanf(line,"%s %d",parm,&i);
      obsdata.gr_order=i;
      npar++;
      continue;}
    if(!strcmp(parm,"GR_ANGLE")){
      sscanf(line,"%s %f",parm,&f);
      obsdata.gr_angle=f;
       npar++;
     continue;}
    if(!strcmp(parm,"DEWOFF")){
      sscanf(line,"%s %s",parm,value);
      strcpy(obsdata.dewoff,value);
      npar++;
      continue;}
    if(!strcmp(parm,"CAMERA")){
      n=sscanf(line,"%s %s",parm,value);
      if(n==2){
	strcpy(obsdata.camera,value);}
      else{
	strcpy(obsdata.camera,"");}
      npar++;
      continue;}
    if(!strcmp(parm,"MODE")){
      sscanf(line,"%s %s",parm,value);
      strcpy(obsdata.mode,value);
      npar++;
      continue;}
    if(!strcmp(parm,"MASK")){
      sscanf(line,"%s %s",parm,value);
      strcpy(obsdata.mask,value);
      npar++;
      continue;}
    if(!strcmp(parm,"DEWAR")){
      sscanf(line,"%s %s",parm,value);
      strcpy(obsdata.dewar,value);
      npar++;
      continue;}
      if(!strcmp(parm,"DISTORTION")){
        sscanf(line,"%s %s",parm,value);
        strcpy(obsdata.distor,value);
        npar++;
        continue;}
    if(!strcmp(parm,"GRATING")){
      sscanf(line,"%s %s",parm,value);
      strcpy(obsdata.grating,value);
      npar++;
      continue;}
    if(!strcmp(parm,"D_ALIGNROT")){
      sscanf(line,"%s %f",parm,&f);
      obsdata.alignrot=f;
      npar++;
      continue;}
    }
  fclose(parmfile);
  if(npar<8){
    printf("Error: missing obsdef parameter\n");
    return 1;}
  *obsd=obsdata;
  return 0;}
