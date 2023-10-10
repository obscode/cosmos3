/*****************************************************************************\
*                                                                             *
*  MOSAIC combines  chip images into one large FITS image for viewing only   *
*                                                                             *
*  USAGE: mosaic framename                                                    *
*                                                                             *
*  OUTPUT: framename_m.fits a 8192x8192 (or smaller if input data binned) fits*
*          file                                                               *
*                                                                             *
*  VERSION: 05 Nov 03                                                         *
*                                                                             *
\*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cosmos.h"
#include "fitsio.h"

int main(int argc, char *argv[]){

  int      outelem,status,i,j,k,ii,jj,kk,chip,anynl,ibin,ibiny,bitpix,
           nelem,thisone,subr,xaxis,yaxis;
  unsigned short INT, **array,*image,**outarray,*outimage,*INTP;
  float    FLOAT,*FLOATP,**farray,*fimage,**foutarray,*foutimage;
  long     firstpt[2],npts[2],firstcn[3],lastcn[3],naxes[2],outaxes[2],firstel,
           offset,inc[2],firstelem[3];
  char     CHAR,file[80],ifile[80],dfile[80],line[80],string[133];
  char     *HOME,*DATA_DIR,*COS_HOME;
  fitsfile *fptr[9],*outptr;
  fitsdef  fitsinfo;
  int      nulvl=1;

  inc[0]=inc[1]=1;
  firstel=firstelem[0]=firstelem[1]=firstelem[2]=firstcn[2]=lastcn[2]=1;
  for(i=1;i<=8;i++) fptr[i]=0;
  //directories
  COS_HOME=malloc(sizeof(CHAR)*80);
  HOME=malloc(sizeof(CHAR)*80);
  DATA_DIR=malloc(sizeof(CHAR)*80);
  DATA_DIR=getenv("COSMOS_IMAGE_DIR");
  if(DATA_DIR==NULL){
    printf("COSMOS_IMAGE_DIR undefined!\n");
    return 1;}
  strcat(DATA_DIR,"/");


  //input images
  thisone=0;
  if(argc>1){
    strcpy(dfile,argv[1]);}
  else{
    printf("\nenter file to process: ");
    fgets(string,133,stdin);
    sscanf(string,"%s",dfile);} 
  strcpy(ifile,DATA_DIR);
  strcat(ifile,dfile);
  for(i=1;i<=8;i++){
    strcpy(file,ifile);
    addbar(file);
    strcat(file,"c");
    sprintf(line,"%d.fits",i);
    strcat(file,line);
    status=0;
    status=OpenFitsFile(file,&fptr[i],&fitsinfo);
    if(!status) thisone=i;}
  if(!thisone) die("Cannot find image files ");
      
  //image properties

  subr=fitsinfo.subrstr;
  ibin=ibiny=fitsinfo.binning;
  if(fitsinfo.ybinning) ibiny=fitsinfo.ybinning;
  xaxis=2048/ibin;
  yaxis=4096/ibiny;
  naxes[0]=2048/ibin+fitsinfo.overscan;
  naxes[1]=4096/ibiny+fitsinfo.biaslins;
  bitpix=fitsinfo.bitpix;
  bitpix=abs(bitpix);
  nelem=naxes[0]*naxes[1];
  if(bitpix==16){
    image=malloc(sizeof(INT)*nelem);
    array=malloc(sizeof(INTP)*naxes[1]);
    for(j=0;j<naxes[1];j++) *(array+j)=image+naxes[0]*j;}
  if(bitpix==32){
    fimage=malloc(sizeof(FLOAT)*nelem);
    farray=malloc(sizeof(FLOATP)*naxes[1]);
    for(j=0;j<naxes[1];j++) *(farray+j)=fimage+naxes[0]*j;}
  if(subr){
    firstpt[0]=fitsinfo.subx[0];
    firstpt[1]=fitsinfo.suby[0];
    npts[0]=fitsinfo.subx[1]-fitsinfo.subx[0]+1+fitsinfo.overscan;
    npts[1]=fitsinfo.suby[1]-fitsinfo.suby[0]+1+fitsinfo.biaslins;}


  //output file

  outaxes[0]=8192/ibin;;
  outaxes[1]=8192/ibiny;;
  outelem=outaxes[0]*outaxes[1];
  if(bitpix==16){
    outimage=malloc(sizeof(INT)*outelem);
    for(i=0;i<outelem;i++) *(outimage+i)=0;
    outarray=malloc(sizeof(INTP)*outaxes[1]);
    for(i=0;i<outaxes[1];i++) outarray[i]=outimage+outaxes[0]*i;}
  if(bitpix==32){
    foutimage=malloc(sizeof(FLOAT)*outelem);
    for(i=0;i<outelem;i++) *(foutimage+i)=0;
    foutarray=malloc(sizeof(FLOATP)*outaxes[1]);
    for(i=0;i<outaxes[1];i++) foutarray[i]=foutimage+outaxes[0]*i;}
  status=0;
  strcpy(ifile,"!");
  strcat(ifile,DATA_DIR);
  strcat(ifile,dfile);
  
  strcat(ifile,"_mos.fits");
  fits_create_file(&outptr,ifile,&status);
  if(bitpix==16)fits_create_img(outptr,USHORT_IMG,2,outaxes,&status);
  else fits_create_img(outptr,FLOAT_IMG,2,outaxes,&status);
  if(status)fits_die("Output file error",status);

  //read in image,copy to output image

  for(chip=1;chip<=4;chip++){
    if(!fptr[chip]) continue;
    printf("Reading chip %d\r",chip);
    fflush(stdout);
    status=0;
    if(bitpix==16){
      if(!subr){
	fits_read_pix(fptr[chip],TUSHORT,firstelem,nelem,&nulvl,image,&anynl,
                      &status);}
      else{
	firstcn[0]=1;
	lastcn[0]=npts[0];
	for(i=1;i<=npts[1];i++){
	  firstcn[1]=lastcn[1]=i;
	  offset=(naxes[0]*(i+firstpt[1]-1))+firstpt[0]-1;
	  fits_read_subset(fptr[chip],TUSHORT,firstcn,lastcn,inc,&nulvl,
                           image+offset,&anynl,&status);}
        }
      }
    else{
      if(!subr){
	fits_read_pix(fptr[chip],TFLOAT,firstelem,nelem,&nulvl,fimage,&anynl,
		      &status);}
      else{
	firstcn[0]=1;
	lastcn[0]=npts[0];
	for(i=1;i<=npts[1];i++){
	  firstcn[1]=lastcn[1]=i;
	  offset=(naxes[0]*(i+firstpt[1]-1))+firstpt[0]-1;
	  fits_read_subset(fptr[chip],TFLOAT,firstcn,lastcn,inc,&nulvl,
                           fimage+offset,&anynl,&status);}
        }
      }
    if(status)fits_die("File error",status);

    //copy image, flipping in y and offsetting
    for(i=0;i<yaxis;i++){
      kk=outaxes[1]-1-i;
      jj=(chip-1)*xaxis;
	if(bitpix==16){
	  for(j=0;j<xaxis;j++){
	    ii=jj+j;
	    *(outarray[kk]+ii)= *(array[i]+j);}
	  }
	else{
	  for(j=0;j<xaxis;j++){
	    ii=jj+j;
	    *(foutarray[kk]+ii)= *(farray[i]+j);}
	  }
      }
    }

  for(chip=5;chip<=8;chip++){
    if(!fptr[chip]) continue;
    printf("Reading chip %d\r",chip);
    fflush(stdout);
    status=0;
    if(bitpix==16){
      if(!subr){
	fits_read_pix(fptr[chip],TUSHORT,firstelem,nelem,&nulvl,image,&anynl,
                      &status);}
      else{
	firstcn[0]=1;
	lastcn[0]=npts[0];
	for(i=1;i<=npts[1];i++){
	  firstcn[1]=lastcn[1]=i;
	  offset=(naxes[0]*(i+firstpt[1]-1))+firstpt[0]-1;
	  fits_read_subset(fptr[chip],TUSHORT,firstcn,lastcn,inc,&nulvl,
                           image+offset,&anynl,&status);}
        }
      }
    else{
      if(!subr){
	fits_read_pix(fptr[chip],TFLOAT,firstelem,nelem,&nulvl,fimage,&anynl,
		      &status);}
      else{
	firstcn[0]=1;
	lastcn[0]=npts[0];
	for(i=1;i<=npts[1];i++){
	  firstcn[1]=lastcn[1]=i;
	  offset=(naxes[0]*(i+firstpt[1]-1))+firstpt[0]-1;
	  fits_read_subset(fptr[chip],TFLOAT,firstcn,lastcn,inc,&nulvl,
                           fimage+offset,&anynl,&status);}
        }
      }
    if(status)fits_die("Input file error",status);
    fits_close_file(fptr[chip],&status);
    //copy image, flipping in y and offsetting
    for(i=0;i<yaxis;i++){
	switch(chip){
	case 6: jj=xaxis-1;
	  break;
	case 5: jj=2*xaxis-1;
	  break;
	case 8: jj=3*xaxis-1;
	  break;
	case 7: jj=4*xaxis-1;}
	if(bitpix==16){
	  for(j=0;j<xaxis;j++){
	    ii=jj-j;
	    *(outarray[i]+ii)= *(array[i]+j);}
	  }
	else{
	  for(j=0;j<xaxis;j++){
	    ii=jj-j;
	    *(foutarray[i]+ii)= *(farray[i]+j);}
	  }
      }
    }
  //write output file
  printf("Writing output file\n");
  status=0;
  if(bitpix==16){
    fits_write_pix(outptr,TUSHORT,firstelem,outelem,outimage,&status);}
  else{
    fits_write_pix(outptr,TFLOAT,firstelem,outelem,foutimage,&status);}
  if(status)fits_die("Output file error",status);
  fits_close_file(outptr,&status);
  return  0;
}

 

