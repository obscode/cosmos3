/*****************************************************************************\
*                                                                             *
*   BADORDERS   calculates the postions of spectral lines on each ccd         *
*               from mask definition file for dispersed images, in undesirec  *
*               orders and adds positions to a bad pixel file                 *
*                                                                             *
*   VERSION 17 Jan 06                                                         *
*                                                                             *
*   USAGE:   badorders -o observeset -b badpixelfile                          *
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
#include "mgfeats.h"

void      writeit(void);
int       xsize,ysize,border,xmin,xmax,ymin,ymax,chippy;
FILE      *outfile;

int  main(int argc,char *argv[]){

  Obs     *obset = NULL;
  slit    slitdat;
  element *instrument,*optics;
  vect3   mskpos,ccdpos;
  vect2   mpend;
  double  cur_wavl, cur_temp;
  int     ocen,msize,ixcd,iycd,narg,ibin,n,i,chip[5],nlines,chip1,zerord,tword,
          ends,xsize2,ysize2,chgrp[5],nchips,nchp,j,order,nl,nll,n0,n1,nord,
          ordernum[4],no,shuffle,shf,xbin,ybin,k,l,m;
  float   xcd[4],ycd[4],xccd,yccd,xfp,yfp,lines[500],x8cd,y8cd,lambda[2],
          hwidth;
  char    *env,*COS_HOME,*point;
  char    CHAR,flnm[80],file[80],alphabet[53],linfile[80],maskfile[80],
          line[133],SMF_FILE[80],camdef[80],camoff[80],imacsdir[133],
          obadfile[80],ibadfile[80],orders[10];
  objq    *oq, *oq1;
  obsdef  obsdata;
  objdat  objdata;
  dewdat  dewinfo;
  FILE    *linefile,*badfile;

  point=NULL;
  COS_HOME=malloc(sizeof(CHAR)*80);
  COS_HOME=getenv("COSMOS_HOME");
  if(COS_HOME==NULL){
    printf("COSMOS_HOME undefined!\n");
    return 1;}

  /*------------------------get observ data----------------------------------*/
  i=1;
  ocen=ends=0;
  nlines=0;
  strcpy(ibadfile,"");
  if(argc<5){
    for(;;){
      printf("Enter observation name: ");
      scanf("%s",maskfile);
      if(ReadObsDef(maskfile,&obsdata)==0) break;
      printf("Cannot find obsdef file %s\n",maskfile);}
    strcpy(line,COS_HOME);
    strcat(line,"/sdata/");
    strcat(line,obsdata.dewar);
    printf("Enter input bad pixel file [%s]: ",line);
    scanf("%s",ibadfile);
    printf("Enter output bad pixel file: ");
    scanf("%s",obadfile);
    }

  else{
    if(!(strcmp(argv[1],"-o"))){
      narg=2;}
    else{
      if(!(strcmp(argv[3],"-o"))){
        narg=4;}
      else{
        if(argc==7 && !(strcmp(argv[5],"-o"))){
          narg=6;}
        else{
          printf("proper invocation:badorders -o obserset [-i input_badpixel_file] -f  output_badpixel_file\n");
          return 1;}
      }
    }
    strcpy(maskfile,argv[narg]);;
    if(ReadObsDef(maskfile,&obsdata)!=0){
      printf("Error reading observation definition file! %s\n",maskfile);
      return 1;}
    if(!(strcmp(argv[1],"-f"))){
      narg=2;}
    else{
      if(!(strcmp(argv[3],"-f"))){
        narg=4;}
      else{
        if(argc==7 && !(strcmp(argv[5],"-f"))){
          narg=6;}
        else{
          printf("proper invocation:badorders -o obserset [-i input_badpixel_file] -f  output_badpixel_file\n");
          return 1;}
      }
    }
    strcpy(obadfile,argv[narg]);
    narg=0;
    if(!(strcmp(argv[1],"-i"))){
      narg=2;}
    else{
      if(!(strcmp(argv[3],"-i"))){
        narg=4;}
      else{
        if(argc==7 && !(strcmp(argv[5],"-i"))){
          narg=6;}
      }
    }
    if(narg) strcpy(ibadfile,argv[narg]);
  }

  //get program parameters

  nlines=0;
  if(OpenCosParm("badorders")) die("Cannot open parameter file");
  if(ReadParm_i("BORDER",&border)) die("parameter file error");
  if(ReadParm_s("ORDERS",orders))die("parameter file error");
  nord=strlen(orders);
  for(i=0;i<nord;i++){
    line[0]=orders[i];
    line[1]='\0';
    sscanf(line,"%d",ordernum+i);
    if(ordernum[i]==0){
      if(ReadParm_r("LAMBDA0",&lambda[0])) die("parameter file error");
      if(ReadParm_r("LAMBDA1",&lambda[1])) die("parameter file error");}
    if(ReadParm_i("SHUFFLE",&shuffle)) shuffle=0;
    if(ordernum[i]>1 && nlines==0){
      if(ReadParm_r("HWIDTH",&hwidth)) die("Parameter file error");
      if(ReadParm_s("LINEFILE",linfile)) die("parameter file error");
      linefile=fopen(linfile,"r");
      if(linefile==NULL) die("cannot open 2nd order line file");
      while(nlines<500){
        if(fgets(line,133,linefile)==NULL) break;
        sscanf(line,"%f",&lines[nlines]);
        nlines++;}
      if(nlines==500) printf("2nd order line list limited to 500 lines\n");
      }
    }

  //copy existing bad pixel file to new file

  strcat(obadfile,".badpix");
  if((outfile=fopen(obadfile,"w"))==NULL)die("Cannot open new bad pixel file");
  //default input bad pixel file
  if(strlen(ibadfile)==0){
    strcpy(ibadfile,COS_HOME);
    strcat(ibadfile,"/sdata/badpix/");
    strcat(ibadfile,obsdata.dewar);}
  ibin=xbin=ybin=1;
  //if using input bad pixel file
  if(strcasecmp(ibadfile,"none")){
    strcat(ibadfile,".badpix");
    if((badfile=fopen(ibadfile,"r"))==NULL){
      printf ("cannot open bad pixel file %s\n",ibadfile);}
    else{
      if(!(fgets(line,133,badfile))){
        printf ("Input bad pixel file is empty\n");}
      else{
        fputs(line,outfile);
        //binning
        if(sscanf(line,"%d %d %d %d %d %d",&i,&j,&k,&l,&m,&ibin)<6){
          ibin=xbin=ybin=1;}
        else{
          if(ibin<10){
            xbin=ybin=ibin;}
          else{
            i=(int) ibin/10;
            ybin=ibin-10*i;
            xbin=i;}
        }
      }
      while(1){
        if(!fgets(line,133,badfile)) break;
        fputs(line,outfile);}
      }
    }
  printf("producing bad pixel file for %d x %d binning\n",xbin,ybin);

  //for each order analyzed

  for(no=0;no<nord;no++){
    order=ordernum[no];
    //get obsdef data
    strcpy(flnm,maskfile);
    strcat(flnm,"_");
    n=strlen(flnm);
    flnm[n]=orders[no];
    flnm[n+1]='\0';
    if(ReadObsDef(flnm,&obsdata)!=0){
      printf("Error reading observation definition file! %s\n",flnm);
      return 1;}
    //get mask data
    strcpy(flnm,obsdata.mask);
    strcpy(SMF_FILE,flnm);
    strcat(SMF_FILE,".SMF");
    if(ReadSMFfile(SMF_FILE,&obset)) return 1;
    //set up instrument
    if(strcmp(obsdata.mode,"SPEC")){
      printf("spectral-lines requires a dispersed image\n");
      return 1;}
    SetupInstr(&obsdata,&instrument);
    if(SetupCamera(&obsdata)==1) return 1;
    Getchipdat(&dewinfo);
    xsize=dewinfo.xchip/xbin;
    ysize=dewinfo.ychip/ybin;
    xsize2=xsize/2;
    ysize2=ysize/2;
    oq=obset->ob;
    oq1=oq;
    n=1;

    //loop through slits
    while(1){
      while(2){
      slitdat=oq->slit;
      objdata=oq->dat;
      mskpos = get3de2v(oq->smpos,0.0);
      nchips=0;
      //end wavelengths
      nll = (order==0) ? 1 : nlines;
      //for each line
      for(nl=0;nl<nll;nl++){
        for(shf=0;shf<=shuffle;shf+=shuffle){
          nchp=0;
          for(i=0;i<2;i++){
            if(order>0){
              cur_wavl=lines[nl]+hwidth*(2*i-1);}
            else{
              cur_wavl=lambda[i];}
            //"left" end
            mpend=sum2vect(oq->smpos,lslit(oq->slit));
            mskpos = get3de2v(mpend,0.0);
            ccdpos = Op_transform(mskpos,instrument,cur_wavl,cur_temp);
            xfp=ccdpos.x;
            yfp=ccdpos.y;
            chip[nchp]=fp2ccd(xfp,yfp,&xccd,&yccd,cur_wavl);
            xcd[nchp]=xccd/xbin;
            ycd[nchp]=(yccd+shf*dewinfo.sy[chip[nchp]])/ybin;
            if(ycd[nchp]<0 || ycd[nchp]>=ysize) chip[nchp]=0;
            nchips=(chip[nchp]>0)?chip[nchp]:0;
            nchp++;
            //"right" end
            mpend=sum2vect(oq->smpos,rslit(oq->slit));
            mskpos = get3de2v(mpend,0.0);
            ccdpos = Op_transform(mskpos,instrument,cur_wavl,cur_temp);
            xfp=ccdpos.x;
            yfp=ccdpos.y;
            chip[nchp]=fp2ccd(xfp,yfp,&xccd,&yccd,cur_wavl);
            xcd[nchp]=xccd/xbin;
            ycd[nchp]=(yccd+shf*dewinfo.sy[chip[nchp]])/ybin;
            if(ycd[nchp]<0 || ycd[nchp]>=ysize) chip[nchp]=0;
            nchips=(chip[nchp]>0)?chip[nchp]:0;
            nchp++;}

          //no corner on a chip?
          if(!nchips) break;

          //all 4 corners on one chip?
          if(chip[0]>0 && chip[0]==chip[1] && chip[0]==chip[2] &&
             chip[0]==chip[3]){
            chippy=chip[0];
            xmin=xsize;
            ymin=ysize;
            xmax=ymax=0;
            for(i=0;i<4;i++){
              xmin = (xmin<xcd[i]) ? xmin : xcd[i];
              ymin = (ymin<ycd[i]) ? ymin : ycd[i];
              xmax = (xmax>xcd[i]) ? xmax : xcd[i];
              ymax = (ymax>ycd[i]) ? ymax : ycd[i];}
            writeit();}

          //other cases
          else{
            for(i=0;i<4;i++){
              if(chip[i]<=0) continue;
              chippy=chip[i];
              chgrp[0]=i;
              nchips=1;
              for(j=i+1;j<4;j++){
                if(chip[j]==chip[i]){
                  chgrp[nchips]=j;
                  chip[j]=0;
                  nchips++;}
              }

              //one corner case
              if(nchips==1){
                if(xcd[chgrp[0]]<xsize2){
                  xmin=1;
                  xmax=xcd[chgrp[0]];}
                else{
                  xmin=xcd[chgrp[0]];
                  xmax=xsize;}
                if(ycd[chgrp[0]]<ysize2){
                  ymin=1;
                  ymax=ycd[chgrp[0]];}
                else{
                  ymin=ycd[chgrp[0]];
                  ymax=ysize;}
              }
              else{
                // two corner case
                if(nchips==2){
                  if(fabs(xcd[chgrp[0]]-xcd[chgrp[1]])<fabs(ycd[chgrp[0]]-
                                                            ycd[chgrp[1]])){
                    //vertical pair
                    if(xcd[chgrp[0]]<xsize2){
                      xmin=1;
                      xmax= (xcd[chgrp[0]]>xcd[chgrp[1]]) ? xcd[chgrp[0]] :
                        xcd[chgrp[1]];}
                    else{
                      xmin= (xcd[chgrp[0]]<xcd[chgrp[1]]) ? xcd[chgrp[0]] :
                      xcd[chgrp[1]];
                      xmax=xsize;}
                    if(ycd[chgrp[0]]<ycd[chgrp[1]]){
                      ymin=ycd[chgrp[0]];
                      ymax=ycd[chgrp[1]];}
                    else{
                      ymin=ycd[chgrp[1]];
                      ymax=ycd[chgrp[0]];}
                  }
                  else{
                    //horizontal pair
                    if(ycd[chgrp[0]]<ysize2){
                      ymin=1;
                      ymax= (ycd[chgrp[0]]>ycd[chgrp[1]]) ? ycd[chgrp[0]] :
                        ycd[chgrp[1]];}
                    else{
                      ymin= (ycd[chgrp[0]]<ycd[chgrp[1]]) ? ycd[chgrp[0]] :
                      ycd[chgrp[1]];
                      ymax=ysize;}
                    if(xcd[chgrp[0]]<xcd[chgrp[1]]){
                      xmin=xcd[chgrp[0]];
                      xmax=xcd[chgrp[1]];}
                    else{
                      xmin=xcd[chgrp[1]];
                      xmax=xcd[chgrp[0]];}
                  }
                }
                //three corner case
                else{
                  xmin= (xcd[chgrp[0]]<xcd[chgrp[1]]) ? xcd[chgrp[0]] : xcd[chgrp[1]];
                  xmin= (xmin < xcd[chgrp[2]]) ? xmin : xcd[chgrp[2]];
                  ymin= (ycd[chgrp[0]]<ycd[chgrp[1]]) ? ycd[chgrp[0]] : ycd[chgrp[1]];
                  ymin= (ymin < ycd[chgrp[2]]) ? ymin : ycd[chgrp[2]];
                  xmax= (xcd[chgrp[0]]>xcd[chgrp[1]]) ? xcd[chgrp[0]] : xcd[chgrp[1]];
                  xmax= (xmax>xcd[chgrp[2]]) ? xmax : xcd[chgrp[2]];
                  ymax= (ycd[chgrp[0]]>ycd[chgrp[1]]) ? ycd[chgrp[0]] : ycd[chgrp[1]];
                  ymax= (ymax>ycd[chgrp[2]]) ? ymax : ycd[chgrp[2]];
                }
              }
              writeit();}
          }
          if(!shuffle) break;}
      }
      break;}
      oq=oq->next;
      if(oq==oq1) break;}
  }
  fclose(outfile);
  return 0;}


void writeit(){

  xmin-=border;
  xmin = (xmin>0) ? xmin : 1;
  ymin-=border;
  ymin = (ymin>0) ? ymin : 1;
  xmax+=border+1;
  xmax = (xmax<xsize) ? xmax : xsize;
  ymax+=border+1;
  ymax = (ymax<ysize) ? ymax : ysize;
  fprintf(outfile,"%1d %4d %4d %4d %4d\n",chippy,xmin,xmax,ymin,ymax);
  return;}
