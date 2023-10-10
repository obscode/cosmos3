/*****************************************************************************\
*                                                                             *
*  COSMOSUTILS - utility functions for COSMOS programs                        *
*                                                                             *
*  VERSION       28 Feb 2004                                                  *
*                                                                             *
\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cosdat.h"
#include "cosmos.h"


/*
 *   die- print message and quit
 *
 */

void die(char message[]){

  printf("%s\n",message);
  exit(1);}



/*
 *   addbar -add underline to file names containing _ that don't end in _
 *
 */

void addbar(char file[]){

  unsigned long len;
  int i,index;
  char *p;

  len=strlen(file);
  if(*(file+len-1) == '_') return;

  index=0;
  for(i=0;i<len;i++){
    if(file[i]=='/') index=i;}
  if(strstr(file+index,"_")==NULL) return;
  strcat(file,"_");
  return;}

/*
 *  subbars- deletes _b_ _f_ _s_ appendages from file names
 *
 */

void subbars(char file[]){

  char *p;
  unsigned long len;

  p=strstr(file,"_b_");
  if(p==NULL){
    p=strstr(file,"_f_");
    if(p==NULL){
      p=strstr(file,"_s_");
      if(p==NULL) return;}
    }
  len=p-file;
  *(file+len)='\0';
  return;}

/*
 *
 * interpolate- interpolate within array
 *
 */

/* interpolate floating point numbers */

float f_interpol(float **array, long naxes[], float xpos, float ypos){

  int   ixpos,iypos;
  float dix,diy,value,d11,d21,d22,d12;

  ixpos=(int)xpos;
  iypos=(int)ypos;
  dix=xpos-(int)xpos;
  diy=ypos-(int)ypos;
  if(ixpos<0 || ixpos>=naxes[0]-2 || iypos<0 || iypos>=naxes[1]-2){
     //value=nanf(&empty);
     value=0.0;
     }
  else{
    if(dix==0){
      value= *(*(array+iypos)+ixpos)*(1-diy) + *(*(array+iypos+1)+ixpos)*diy;}
    else{
      if(diy==0){
	value= *(*(array+iypos)+ixpos)*(1-dix) + *(*(array+iypos)+ixpos+1)*dix;}
      else{
	d22=dix*diy;
	d12=(1-dix)*diy;
	d21=dix*(1-diy);
	d11=(1-dix)*(1-diy);
	value= *(*(array+iypos)+ixpos)*d11 +
	  *(*(array+iypos)+ixpos+1)*d21 +
	  *(*(array+iypos+1)+ixpos)*d12 +
	  *(*(array+iypos+1)+ixpos+1)*d22;
	value/=(d11+d12+d21+d22);}
    }
  }
  return value;}

/* interpolate integer numbers */

float i_interpol(int **array, long naxes[], float xpos, float ypos){

  int   ixpos,iypos;
  float dix,diy,value,d22,d21,d11,d12;
  char empty;

  ixpos=(int)xpos;
  iypos=(int)ypos;
  dix=xpos-ixpos;
  diy=ypos-iypos;
  if(ixpos<0 || ixpos>=naxes[0]-2 || iypos<0 || iypos>=naxes[1]-2){
     //value=nanf(&empty);
     value=0;
     }
  else{
    if(dix==0){
      value= *(*(array+iypos)+ixpos)*(1-diy) + *(*(array+iypos+1)+ixpos)*diy;}
    else{
      if(diy==0){
	value= *(*(array+iypos)+ixpos)*(1-dix) + *(*(array+iypos)+ixpos+1)*dix;}
      else{
	d22=dix*diy;
	d12=(1-dix)*diy;
	d21=dix*(1-diy);
	d11=(1-dix)*(1-diy);
	value= *(*(array+iypos)+ixpos)*d11 +
	  *(*(array+iypos)+ixpos+1)*d21 +
	  *(*(array+iypos+1)+ixpos)*d12 +
	  *(*(array+iypos+1)+ixpos+1)*d22;
	value/=(d11+d12+d21+d22);}
    }
  }
  return value;}

/* interpolate floating point numbers, with zero weight for bad values */

float e_interpol(float **array, float **earray, long naxes[], float xpos,
                 float ypos){

  int   ixpos,iypos;
  float dix,diy,value,d11,d21,d22,d12,wt,e00,e01,e10,e11;

  ixpos=(int)xpos;
  iypos=(int)ypos;
  dix=xpos-(int)xpos;
  diy=ypos-(int)ypos;
  if(ixpos<0 || ixpos>=naxes[0]-2 || iypos<0 || iypos>=naxes[1]-2){
    //value=nanf(&empty);
    value=0.0;
    }
  else{
    //any bad values?
    e00=e01=e10=e11=1;
    if(*(*(earray+iypos)+ixpos)<=0) e00=0;
    if(*(*(earray+iypos+1)+ixpos)<=0) e01=0;
    if(*(*(earray+iypos)+ixpos+1)<=0) e10=0;
    if(*(*(earray+iypos+1)+ixpos+1)<=0) e11=0;
       if(dix==0){
         value=*(*(array+iypos)+ixpos)*(1-diy)*e00 +
                *(*(array+iypos+1)+ixpos)*diy*e01;
         wt=(1-diy)*e00+diy*e01;
         if(wt>0) value/=wt;}
    else{
      if(diy==0){
        value= *(*(array+iypos)+ixpos)*(1-dix)*e00 +
               *(*(array+iypos)+ixpos+1)*dix*e10;
        wt=(1-dix)*e00+dix*e10;
        if(wt>0)value/=wt;}
      else{
        d22=dix*diy*e11;
        d12=(1-dix)*diy*e01;
        d21=dix*(1-diy)*e10;
        d11=(1-dix)*(1-diy)*e00;
        value= *(*(array+iypos)+ixpos)*d11 +
          *(*(array+iypos)+ixpos+1)*d21 +
          *(*(array+iypos+1)+ixpos)*d12 +
          *(*(array+iypos+1)+ixpos+1)*d22;
        wt=(d11+d12+d21+d22);
        if(wt>0) value/=wt;
        else value=0;}
    }
  }
  return value;}




/*
 *  subbias  bias subtraction routine
 */


void subbias(float **array,long naxes[],int biasx0,int biasx1,int biasy0,
             int biasy1,int xsize,int ysize){

  int j,k,n,nbiasx2,nbiasy2,irnk[1024];
  float f,rank[1024],bmean,bscat,bsum;

  nbiasx2=(biasx1-biasx0+1)/2;
  nbiasy2=(biasy1-biasy0+1)/2;

  //bias row correction
  for(j=0;j<naxes[1];j++){
    n=0;
    for(k=biasx0;k<biasx1;k++){
      f=(float)(*(*(array+j)+k));
      order(&n,&f,irnk,rank);
      n++;}
    bmean=rank[nbiasx2];
    n=0;
    for(k=biasx0;k<biasx1;k++){
      f=fabs((float)(*(*(array+j)+k)) - bmean);
      order(&n,&f,irnk,rank);
      n++;}
    bscat=1.49*rank[nbiasx2];
    n=0;
    bsum=0.0;
    for(k=biasx0;k<biasx1;k++){
      f=(float)(*(*(array+j)+k));
      if (fabs(f-bmean) < 5*bscat) {
        bsum += f;
        n++;}
      }
    if (n) {
      bmean = bsum/n;
        for(k=0;k<naxes[0];k++) *(*(array+j)+k)-=bmean;
      } else {
        for(k=0;k<naxes[0];k++) *(*(array+j)+k)=0.0;
      }
   }

  //bias column correction
  bmean=0.;
  if(nbiasy2){
    for(k=0;k<xsize;k++){
      n=0;
      for(j=biasy0;j<biasy1;j++){
	f=(float)(*(*(array+j)+k));
	order(&n,&f,irnk,rank);
	n++;}
      bmean=rank[nbiasy2];
      n=0;
      for(j=biasy0;j<biasy1;j++){
	f=fabs((float)(*(*(array+j)+k))-bmean);
	order(&n,&f,irnk,rank);
	n++;}
      bscat=1.49*rank[nbiasy2];
      n=0;
      bsum=0.0;
      for(j=biasy0;j<biasy1;j++){
	f=(float)(*(*(array+j)+k));
        if (fabs(f-bmean) < 5*bscat) {
	  bsum += f;
	  n++;}
        }
      if (n) {
        bmean = bsum/n;
        for(j=0;j<ysize;j++) *(*(array+j)+k)-=bmean;
      } else {
        bmean = bsum/n;
        for(j=0;j<ysize;j++) *(*(array+j)+k)=0.0;
        }
    }
    }

  return;}

int ReadGain(char dewar[], char speed[], float gain[]){

  FILE     *gainfile;
  char     CHAR,dwr[10],spd[4],line[133],fpath[130];
  int      i,nopar;
  char     *COS_HOME;

  COS_HOME=malloc(sizeof(CHAR)*100);
  COS_HOME=getenv("COSMOS_HOME");
  if(COS_HOME==NULL){
    printf("COSMOS_HOME undefined!\n");
    return 1;}
  strcpy(fpath,COS_HOME);
  strcat(fpath,"/sdata/gains.dat");
  if((gainfile=fopen(fpath,"r")) == NULL){
    printf("%s does not exist!\n",fpath);
    return 1;}
//  printf("fpath = %s\n",fpath);
  while(fgets(line, 133, gainfile)){
    if(*line=='#') continue;
    sscanf(line,"%s %s",dwr,spd);
    if( !strcmp(dwr,dewar) && !strcmp(spd,speed) ){
      sscanf(line,"%s %s %f %f %f %f %f %f %f %f",dwr,spd,&gain[1],&gain[2],&gain[3],&gain[4],&gain[5],&gain[6],&gain[7],&gain[8]);
      break;}
//    printf("%s\n",line);
    }
  fclose(gainfile);
  return 0;}
