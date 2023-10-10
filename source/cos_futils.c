/*
 * OpenFitsFile opens a fits image file and retrieves primary information
 *
 * VERSION 5 Nov '03
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fitsio.h"
#include "cosdat.h"

void fits_die(char *mes,int status){

  char  line[80],out[80];

  strcpy(out,mes);
  strcat(out,".  Error: ");
  fits_get_errstatus(status,line);
  strcat(out,line);
  printf("%s\n",out);
  exit(1);}

int OpenFitsFile(char *filename,fitsfile **fptr,fitsdef *fitsinfo){
  int biaslins,overscan,ibin,i,j,k,l,binning,bitpix,status,naxis,i1,i2,i3,i4,
      nshuffle;
  float ranod,decnod,exptime;
  long naxes[3];
  char param[20],line[80],dwri[3],utdate[12];

  //try to open file
  status=0;
  fits_open_file(fptr,filename,READONLY,&status);
  if(status) return status;

  //get parameters

  fits_get_img_param(*fptr,3,&bitpix,&naxis,naxes,&status);
  if(status) return status;
  fitsinfo->naxis=naxis;
  fitsinfo->bitpix=bitpix;
  for(i=0;i<naxis;i++) *(fitsinfo->naxes+i)=naxes[i];
  fits_read_key(*fptr,TSTRING,"BINNING",param,line,&status);
  //binning
  fitsinfo->ybinning=0;
  if(status){
    ibin=1;
    fitsinfo->binning=ibin;}
  else{
    switch(param[0]){
      case '1':
        ibin=1;
        break;
      case '2':
        ibin=2;
        break;
      case '4':
        ibin=4;
        break;
      default:
        printf("inappropriate binning %s\n",binning);
        return 1;}
    fitsinfo->binning=ibin;}
  //is binning square?
  if(param[1] == 'x'){
    switch(param[2]){
      case '1':
        ibin=1;
        break;
      case '2':
        ibin=2;
        break;
      case '4':
        ibin=4;
        break;
      default:
        printf("inappropriate y binning %s\n",binning);
        return 1;}
    fitsinfo->ybinning=ibin;}
  fits_read_key(*fptr,TFLOAT,"EXPTIME",&exptime,line,&status);
  if(status) exptime=0;
  fitsinfo->exptime=exptime;
  //what is the date of the observations
  fits_read_key(*fptr,TSTRING,"UT-DATE",utdate,line,&status);
  if(status) strcpy(utdate,"");
  strcpy(fitsinfo->date,utdate);
  //is this a subraster?
  fitsinfo->subrstr=0;
  status=0;
  fits_read_key(*fptr,TSTRING,"SUBRASTR",param,line,&status);
  if(!status && strcmp(param,"none")){
    fitsinfo->subrstr=1;
    sscanf(param,"%d:%d,%d:%d",&i,&j,&k,&l);
    *fitsinfo->subx=i;
    *(fitsinfo->subx+1)=j;
    *fitsinfo->suby=k;
    *(fitsinfo->suby+i)=l;}
  status=0;
  fits_read_key(*fptr,TINT,"OVERSCAN",&overscan,line,&status);
  if(status){
    status=0;
    fits_read_key(*fptr,TINT,"NOVERSCN",&overscan,line,&status);}
  if(status) overscan=0;
  fitsinfo->overscan=overscan;
  status=0;
  fits_read_key(*fptr,TSTRING,"DEWARORI",dwri,line,&status);
  if(status)strcpy(fitsinfo->dewarori,"N");
  else{
    if(!strcmp(dwri,"Nod&Shuffle")){
      strcpy(fitsinfo->dewarori,"NS");}
    else strcpy(fitsinfo->dewarori,"N");}
    fitsinfo->overscan=overscan;
    status=0;
  fits_read_key(*fptr,TINT,"BIASLINS",&biaslins,line,&status);
  if(status){
    status=0;
    fits_read_key(*fptr,TINT,"NBIASLNS",&biaslins,line,&status);}
  if(status)biaslins=0;
  fitsinfo->biaslins=biaslins;
  status=0;
  fits_read_key(*fptr,TINT,"NSHUFFLE",&nshuffle,line,&status);
  if(status) nshuffle=0;
  fitsinfo->nshuffle=nshuffle;
  if(!status){
    fits_read_key(*fptr,TFLOAT,"NOD-RA",&ranod,line,&status);
    if(status) fits_die("FITS header error",status);
    fitsinfo->ranod=ranod;
    fits_read_key(*fptr,TFLOAT,"NOD-DEC",&decnod,line,&status);
    if(status)  fits_die("FITS header error",status);
    fitsinfo->decnod=decnod;}
  return 0;}

int ReadFitsFile(fitsfile **fppr, fitsdef *fitsinfo, int type, int *image){
  
  int   point,anynl,status,nulvl,nx0,nx1,nchip,nelem,i,j,fxpxl,fypxl;
  int   **array, *INTP;  
  long  naxes[2],firste[2];
  long  firstel=1;
  float *fimage,**farray,*FLOATP;
  
  //set data type
  if(type == TFLOAT){
    //point=(int) image;
    fimage=(float*) image;}
  nx0=naxes[0]=fitsinfo->naxes[0];
  nx1=naxes[1]=fitsinfo->naxes[1];
  nelem=naxes[0]*naxes[1];
  status=0;
  if(type==TINT){
    fits_read_img(*fppr,TINT,firstel,nelem,&nulvl,image,&anynl,&status);}
  else{
    fits_read_img(*fppr,TFLOAT,firstel,nelem,&nulvl,fimage,&anynl,&status);}
  if(status) return status;
  return 0;}

