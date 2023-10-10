/*****************************************************************************\
*
*   GetSlitXY8   calculates the postions of mask holes and slits on each ccd
*               from mask definition file for direct or dispersed images, and combines, for IMACS,
*               8-chip coord system used by mosaic
*
*   VERSION 6 Nov 2017
*
*   USAGE:   GetSlitXY8(obsdef_file,nlines,wavelength list, plotpos,which, X8,Y8)
                     return value: number of slit images
\*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ioutils.h"
#include "kdcutil.h"
#include "mgutils.h"
#include "optutils.h"
#include "cosmos.h"
#include "clardy.h"

static    GotCof=0;

int  GetSlitXY8(char* obsdefile, int nlines, float* lines,int ibin, int plotpos, int which, 
				int* X8, int* Y8){

  Obs     *obset = NULL;
  slit    slitdat;
  element *instrument,*optics;
  vect3   mskpos,ccdpos;
  vect2   mpend;
  double  cur_wavl, cur_temp;
  int     cen,xsize,ysize,msize,ixcd,iycd,narg,n,i,chip,chip1,ends,ybin,nim,imacs,xcd,ycd;
  float   FLOAT,xfp,yfp,x8cd,y8cd,xccd,yccd,xrccd,yrccd;
  char    *env;
  char    CHAR,flnm[80],file[80],maskfile[80],linfile[80],
          line[133],SMF_FILE[80],imacsdir[133];
  objq    *oq, *oq1;
  obsdef  obsdata;
  objdat  objdata;
  FILE    *linefile,*outfile;

  /*------------------------get observ data----------------------------------*/
  imacs=0;
  cen=ends=0;
  if(plotpos>0) ends=1;
  else if(plotpos==0) cen=1;

  if(ReadObsDef(obsdefile,&obsdata)!=0){
      printf("Error reading observation definition file! %s\n",obsdefile);
      return 0;}
  //get mask data
  strcpy(flnm,obsdata.mask);
  strcpy(SMF_FILE,flnm);
  strcat(SMF_FILE,".SMF");
  if(ReadSMFfile(SMF_FILE,&obset)) return 1;
  //set up instrument
  if(!strcmp(obsdata.mode,"SPEC") ){
    if(nlines==0){
      printf("dispersed images require a line list\n");
      return 0;}
    }
  SetupInstr(&obsdata,&instrument);
  if(!GotCof){
      GotCof=1;
      if(SetupCamera(&obsdata)==1) return 0;}
  if(!strcmp(obsdata.instrument,"IMACS") ) imacs=1;

  /*_______________________use mask file parameters--------------------------*/

  cur_temp=obset->temp;
  //calculate slit positions
  oq=obset->ob;
  oq1=oq;
  n=1;
  nim=0;
  while(1){
    while(2){
      slitdat=oq->slit;
      //     if(slitdat.shape !=2) break;
      objdata=oq->dat;
      mskpos = get3de2v(oq->smpos,0.0);
      for(i=0;i<nlines;i++){
        cur_wavl=lines[i];
        //"left" end if not plotting object position
        if(cen || ends){
          mpend=sum2vect(oq->smpos,lslit(oq->slit));
          mskpos = get3de2v(mpend,0.0);}
        ccdpos = Op_transform(mskpos,instrument,cur_wavl,cur_temp);
        xfp=ccdpos.x;
        yfp=ccdpos.y;
        chip=fp2ccd(xfp,yfp,&xccd,&yccd,cur_wavl);
        if(!chip) continue;
        //"right" end if not plotting object position
        if(cen || ends){
          mpend=sum2vect(oq->smpos,rslit(oq->slit));
          mskpos = get3de2v(mpend,0.0);
          ccdpos = Op_transform(mskpos,instrument,cur_wavl,cur_temp);
          xfp=ccdpos.x;
          yfp=ccdpos.y;
          chip1=fp2ccd(xfp,yfp,&xrccd,&yrccd,cur_wavl);
          if(!chip1 || abs(chip)!=abs(chip1)) continue;
          if(!ends){
            xccd=(xccd+xrccd)/2;
            yccd=(yccd+yrccd)/2;}
          }
        mspos(chip,xccd,yccd,&x8cd,&y8cd);
        xccd/=ibin;
        yccd/=ibin;
        x8cd/=ibin;
        y8cd/=ibin;
        xcd=floor(xccd+1.0);
        ycd=floor(yccd+1.0);
        ixcd=floor(x8cd+1.0);
        iycd=floor(y8cd+1.0);
        X8[nim]= imacs ? ixcd : xcd;
        Y8[nim]= imacs ? iycd: ycd;
        nim++;
        if(ends){
          mspos(chip,xrccd,yrccd,&x8cd,&y8cd);
          xrccd/=ibin;
          yrccd/=ibin;
          x8cd/=ibin;
          y8cd/=ibin;
          xcd=floor(xrccd+1.0);
          ycd=floor(yrccd+1.0);
          ixcd=floor(x8cd+1.0);
          iycd=floor(y8cd+1.0);
          X8[nim]= imacs? ixcd: xcd;
          Y8[nim]= imacs ? iycd:ycd;
          nim++;
                }
      }
      n++;
      break;}
    oq=oq->next;
    if(oq==oq1) break;}
  return nim;}
