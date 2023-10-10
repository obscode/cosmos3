/*****************************************************************************\
*                                                                             *
*   SPECTRAL-LINES   calculates the postions of spectral lines on each ccd    *
*               from mask definition file for dispersed images, combines on   *
*               8-chip coord system used by mosaic, and writes out a list     *
*               of apertues with chip and mosaic coordinates and aper type    *
*                                                                             *
*   VERSION 18 Mar 04                                                         *
*                                                                             *
*   USAGE:   spectral-lines -o obserevset -b binning -l linefile              *
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

int  main(int argc,char *argv[]){

  Obs     *obset = NULL;
  slit    slitdat;
  element *instrument,*optics;
  vect3   mskpos,ccdpos;
  vect2   mpend;
  double  cur_wavl, cur_temp;
  int     ocen,xsize,ysize,msize,ixcd,iycd,narg,ibin,n,i,chip,nlines,chip1,
          ends,ybin;
  float   FLOAT,xccd,yccd,xrccd,yrccd,xfp,yfp,*lines,x8cd,y8cd;
  char    *env;
  char    CHAR,flnm[80],file[80],maskfile[80],alphabet[53],linfile[80],
          line[133],SMF_FILE[80],imacsdir[133];
  objq    *oq, *oq1;
  obsdef  obsdata;
  objdat  objdata;
  FILE    *linefile,*outfile;


  /*------------------------get observ data----------------------------------*/
  ibin=1;
  ocen=ends=0;
  nlines=0;
  lines=malloc(sizeof(FLOAT));
  strcpy(alphabet,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
  if(argc<7){
    for(;;){
      printf("Enter observation name: ");
      scanf("%s",maskfile);
      if(ReadObsDef(maskfile,&obsdata)==0) break;}
    printf("Enter binning: ");
    scanf("%d",&ibin);
    while(1){
      printf("Enter spectral line file  ");
      scanf("%s",linfile);
      linefile=fopen(linfile,"r");
      if(linefile!=NULL) break;
      printf("Cannot open line file\n");}
  }

  else{
    if(!(strcmp(argv[1],"-o"))){
      narg=2;}
    else{
      if(!(strcmp(argv[3],"-o"))){
        narg=4;}
      else{
        if(!(strcmp(argv[5],"-o"))){
          narg=6;}
        else{
          printf("proper invocation:spectral-lines -o obserset -b binning -l linefile\n");
          return 1;}
      }
    }
    strcpy(maskfile,argv[narg]);;
    if(ReadObsDef(maskfile,&obsdata)!=0){
      printf("Error reading observation definition file! %s\n",maskfile);
      return 1;}
    if(!(strcmp(argv[1],"-b"))){
      narg=2;}
    else{
      if(!(strcmp(argv[3],"-b"))){
        narg=4;}
      else{
        if(!(strcmp(argv[5],"-b"))){
          narg=6;}
        else{
          printf("usage: spectral-lines -o obserset -b binning -l linefile\n");
          return 1;}
      }
    }
    sscanf(argv[narg],"%d",&ibin);
    if(!(strcmp(argv[1],"-l"))){
      narg=2;}
    else{
      if(!(strcmp(argv[3],"-l"))){
        narg=4;}
      else{
        if(!(strcmp(argv[5],"-l"))){
          narg=6;}
        else{
          printf("usage: spectral-lines -o obserset -b binning -l linefile\n");
          return 1;}
      }
    }

    //binning
    if(ibin<10){
      ybin=ibin;}
    else{
      i=(int) ibin/10;
      ybin=ibin-10*i;
      ibin=i;}


    //is linelist a file or wavelength?
    if(strpbrk(argv[narg],alphabet)==NULL){
      nlines=1;
      sscanf(argv[narg],"%f",&lines[0]);}
    else{
      sscanf(argv[narg],"%s",linfile);
      linefile=fopen(linfile,"r");
      if(linefile==NULL) die("Cannot open line file");}
    if(argc==8){
      if(!strcmp(argv[7],"-c")) ocen=1;
      if(!strcmp(argv[7],"-e")) ends=1;
       }
  }
  //read spectral line file

  //count lines
  if(!nlines){
    while(1){
      if(fgets(line,133,linefile)==NULL) break;
      nlines++;}
    rewind(linefile);
    lines=realloc(lines,nlines*sizeof(FLOAT));
    for(i=0;i<nlines;i++){
      if(fgets(line,133,linefile)==NULL) break;
      sscanf(line,"%f",&lines[i]);}
    close(linefile);
        }

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

  /*_______________________use mask file parameters--------------------------*/

  cur_temp=obset->temp;

  //print object list

  strcpy(file,maskfile);
  strcat(file,".xy");
    outfile=fopen(file,"w");
    if(outfile == NULL){
      printf("Error opening output file\n");
      return 1;}

  oq=obset->ob;
  oq1=oq;
  n=1;
  while(1){
    while(2){
      slitdat=oq->slit;
      //     if(slitdat.shape !=2) break;
      objdata=oq->dat;
      mskpos = get3de2v(oq->smpos,0.0);
      for(i=0;i<nlines;i++){
        cur_wavl=lines[i];
        //"left" end if not plotting object position
        if(!ocen){
          mpend=sum2vect(oq->smpos,lslit(oq->slit));
          mskpos = get3de2v(mpend,0.0);}
        ccdpos = Op_transform(mskpos,instrument,cur_wavl,cur_temp);
        xfp=ccdpos.x;
        yfp=ccdpos.y;
        chip=fp2ccd(xfp,yfp,&xccd,&yccd,cur_wavl);
//        printf("%f %f %f\n",cur_wavl,xfp,yfp);
        if(!chip) continue;
        //"right" end if not plotting object position
        if(!ocen){
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
        yccd/=ybin;
        x8cd/=ibin;
        y8cd/=ybin;
        xccd=floor(xccd+1.5);
        yccd=floor(yccd+1.5);
        ixcd=floor(x8cd+1.5);
        iycd=floor(y8cd+1.5);
        fprintf(outfile,"%4d %4d %5d %s %d %6.1f %6.1f %6.1f\n",ixcd,iycd,n,
                objdata.name,chip,xccd,yccd,cur_wavl);
        if(ends){
          mspos(chip,xrccd,yrccd,&x8cd,&y8cd);
          xrccd/=ibin;
          yrccd/=ybin;
          x8cd/=ibin;
          y8cd/=ybin;
         xrccd=floor(xrccd+1.5);
         yrccd=floor(yrccd+1.5);
          ixcd=floor(x8cd+1.5);
          iycd=floor(y8cd+1.5);
          fprintf(outfile,"%4d %4d %5d %s %d %6.1f %6.1f %6.1f\n",ixcd,iycd,n,
                  objdata.name,chip,xrccd,yrccd,cur_wavl);
                }
      }
      n++;
      break;}
    oq=oq->next;
    if(oq==oq1) break;}
  fclose(outfile);
  return 0;}
