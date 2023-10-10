/*****************************************************************************\
*                                                                             *
*   APERTURES   calculates the postions of mask holes and slits on each ccd   *
*               from mask definition file for direct images, combines on      *
*               8-chip coord system used by mosaic, and writes out a list     *
*               of apertues with chip and mosaic coordinates and aper type    *
*                                                                             *
*   VERSION  15 Mar 2005                                                      *
*                                                                             *
*   USAGE:   apertures -o obserevset -b binning                               *
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
  objdat  objdata;
  slit    slitdat;
  element *instrument;
  vect2   mpend;
  vect3   mskpos,ccdpos;
  double  cur_wavl, cur_temp;
  int     ixcd,iycd,narg,ibin,n,i,chip,chip1,ybin;
  float   xcd,ycd,xccd,yccd,xrccd,yrccd,xfp,yfp;
  char    CHAR,flnm[80],file[80],maskfile[80],SMF_FILE[80],camoff[80],imacsdir[133];
  objq    *oq, *oq1;
  obsdef  obsdata;
  FILE    *outfile;


  /*------------------------get observ data----------------------------------*/
  ibin=1;
  if(argc<5){
    for(;;){
      printf("Enter observation name: ");
      scanf("%s",maskfile);
      if(ReadObsDef(maskfile,&obsdata)==0) break;
      printf("Error reading obsdef file %s\n",maskfile);}
    printf("Enter binning: ");
    scanf("%d",&ibin);}
  else{
    if(!(strcmp(argv[1],"-o"))){
      narg=2;}
    else{
      if(!(strcmp(argv[3],"-o"))){
	narg=4;}
      else{
	printf("proper invocation: apertures -o obserset -b binning\n");
	return 1;}
    }
    strcpy(maskfile,argv[narg]);
    if(ReadObsDef(maskfile,&obsdata)!=0){
      printf("Error reading observation definition file %s!\n",maskfile);
      return 1;}
    if(!(strcmp(argv[1],"-b"))){
      narg=2;}
    else{
      if(!(strcmp(argv[3],"-b"))){
	narg=4;}
      else{
	printf("proper invocation: apertures -o obserset -b binning\n");
	return 1;}
    }
    sscanf(argv[narg],"%d",&ibin);}

  //binning
  if(ibin<10){
    ybin=ibin;}
  else{
    i=(int) ibin/10;
    ybin=ibin-10*i;
    ibin=i;}

  //read obs def file data

  if(strcmp(obsdata.mode,"DIRECT")){
    printf("apertures requires a direct image\n");
    return 1;}
  SetupInstr(&obsdata,&instrument);
  if(SetupCamera(&obsdata)==1) return 1;
   strcpy(flnm,obsdata.mask);

  /*------------------------- get mask data--------------------------------- */

   strcpy(SMF_FILE,flnm);
   strcat(SMF_FILE,".SMF");
  obset=read_smdf(SMF_FILE,NULL,0);
  if(obset == NULL){
    printf ("No data in SMF file %s\n",SMF_FILE);
    return 1;}

  /*_______________________use mask file parameters--------------------------*/

  cur_wavl=  obset->cw;
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
    while(1){
      slitdat=oq->slit;
      objdata=oq->dat;
//     mskpos = get3de2v(oq->smpos,0.0);
      mpend=sum2vect(oq->smpos,lslit(oq->slit));
      mskpos = get3de2v(mpend,0.0);
      ccdpos = Op_transform(mskpos,instrument,cur_wavl,cur_temp);
      xfp=ccdpos.x;
      yfp=ccdpos.y;
      chip=fp2ccd(xfp,yfp,&xccd,&yccd,cur_wavl);
      if(chip<1) break;;
      mpend=sum2vect(oq->smpos,rslit(oq->slit));
      mskpos = get3de2v(mpend,0.0);
      ccdpos = Op_transform(mskpos,instrument,cur_wavl,cur_temp);
      xfp=ccdpos.x;
      yfp=ccdpos.y;
      chip1=fp2ccd(xfp,yfp,&xrccd,&yrccd,cur_wavl);
      if(!chip1 || abs(chip)!=abs(chip1)) break;
      xccd=(xccd+xrccd)/2;
      yccd=(yccd+yrccd)/2;
      mspos(chip,xccd,yccd,&xcd,&ycd);
      xccd/=ibin;
      yccd/=ybin;
      ixcd=xcd/ibin;
      iycd=ycd/ybin;
      fprintf(outfile,"%4d %4d %5d %s %d %6.1f %6.1f %d %f %f\n",ixcd,iycd,n,
              objdata.name,chip,xccd,yccd,slitdat.shape,xfp,yfp);
      break;}
    oq=oq->next;
    n++;
    if(oq==oq1) break;}
  fclose(outfile);
  return 0;}
