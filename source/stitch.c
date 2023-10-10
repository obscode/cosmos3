/*****************************************************************************\
 *                                                                            *
 *  STITCH combines individual amplifier readouts of chip into one ccd frame  *
 *         with overscan and bias lines subtracted, and preliminary gain      *
 *         correction                                                         *
 *  VERSION  08 April 2005                                                    *
 *                                                                            *
\*****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fitsio.h"
#include "cosmos.h"

#define DELIMITERS ",:[]"

int main(int argc, char *argv[]){

  int      xsize,ysize,biasx0,biasx1,biasy0,biasy1,bitpix,nbfile,nffile,*INTP,
           nulvl,totsizx,totsizy,nflatx,nflaty,n,naxis,nelem,i,j,k,l,nosub,
           oldbsfl,ncard,overscan,zero,biaslins,xsz,ysz,status,INT,
           jl,kl,i1,i2,i3,i4,fxpxl,fypxl,anynl,ibin,ibiny,xbin,ybin;
  float    FLOAT,gain[9],avg;
  float    *FLOATP;
  char     CHAR,*datadir,dfile[80],file[1000],bfile[80],ffile[80],
           string[133],cardline[81],subrastr[81],binline[81],line[133],
           indewar[80],outdewar[80],*homedir,datasec[80],ndatasec[80],
           *secstr;
  long     naxes[3],noxes[3],firstel,firstelem[3],val[4];
  float    *bimage,*dimage,**barray,**darray,*fimage,**farray,*outimage,
           **outarray;
  fitsfile *fptr,*bptr,*fptr_o;
  fitsdef  fitsinfo;
  dewdat   indewinfo,outdewinfo;

  avg=1.;
  nulvl=0;
  zero=0;
  firstelem[0]=firstelem[1]=firstelem[2]=1;
  firstel=1;
  datadir=malloc(sizeof(CHAR)*80);
  datadir=getenv("COSMOS_IMAGE_DIR");
  if(datadir==NULL){
    printf("COSMOS_IMAGE_DIR undefined!\n");
    return 1;}
  strcat(datadir,"/");

  //get parameters
  if(OpenCosParm("stitch")) die("Cannot open parameter file");
  if(ReadParm_r("gain1",&gain[1])==1) die("Parameter file error");
  if(ReadParm_r("gain2",&gain[2])==1) die("Parameter file error");
  if(ReadParm_r("gain3",&gain[3])==1) die("Parameter file error");
  if(ReadParm_r("gain4",&gain[4])==1) die("Parameter file error");
  if(ReadParm_s("dewar",indewar)==1) die("Parameter file error");

  //read input data
  nbfile=0;
  if(argc>1){
    strcpy(dfile,argv[1]);}
  else{
    printf("\nenter file to process: ");
      fgets(string,133,stdin);
      sscanf(string,"%s",dfile);}

  //look in parameter file for bias file
  if(ReadParm_s("bias_file",bfile)==1) die("Parameter file error");
  if(strcasecmp(bfile,"none") && strncmp(bfile,"\0",1)) nbfile=1;

  //get dewar information
  if(Readcdf(indewar)) die("error reading dewar definition file");
  Getchipdat(&indewinfo);

  //determine binning
  strcpy(file,datadir);
  strcat(file,dfile);
  //addbar(file);
  strcat(file,"c");
  sprintf(string,"%d.fits",1);
  strcat(file,string);
  status=0;

  status=OpenFitsFile(file,&fptr,&fitsinfo);
  if(status){
      printf("file %s does not exist!",file);
      return 1;}
  ibin=ibiny=fitsinfo.binning;
  if(fitsinfo.ybinning) ibiny=fitsinfo.ybinning;
  xbin=ibin;
  ybin=ibiny;
  fits_close_file(fptr,&status);

  //setup output buffer
  outimage=(float *) malloc(sizeof(float)*indewinfo.xarsize/xbin*indewinfo.yarsize/ybin);
  outarray=(float **) malloc(sizeof(float *)*indewinfo.yarsize/ybin);
  for(j=0;j<indewinfo.yarsize/ybin;j++) *(outarray+j)=outimage+indewinfo.xarsize/xbin*j;


  //Work chip by chip
  for(i=1;i<=indewinfo.nchip;i++){
    //read in bias file if exist
    if(nbfile){
      strcpy(file,datadir);
      strcat(file,bfile);
      strcat(file,"c");
      sprintf(string,"%d.fits",i);
      strcat(file,string);
      status=0;
      status=OpenFitsFile(file,&bptr,&fitsinfo);
      if(status) die("file %s does not exist!");
      bitpix=fitsinfo.bitpix;
      totsizx=naxes[0]=fitsinfo.naxes[0];
      totsizy=naxes[1]=fitsinfo.naxes[1];
      nelem=totsizx*totsizy;
      //bias border
      overscan=fitsinfo.overscan;
      biaslins=fitsinfo.biaslins;



      status=0;
      if(i==1){
        bimage=(float *) malloc(sizeof(float)*naxes[0]*naxes[1]);
        barray=(float **) malloc(sizeof(float *)*naxes[1]);
        for(j=0;j<naxes[1];j++) *(barray+j)=bimage+naxes[0]*j;}
      fits_read_img(bptr,TFLOAT,firstel,nelem,&nulvl,bimage,&anynl,&status);
      if(status)fits_die("Bias file error",status);
      fits_close_file(bptr,&status);}

    //Read in data file, correct

    strcpy(file,datadir);
    strcat(file,dfile);
//    addbar(file);
    strcat(file,"c");
    sprintf(string,"%d.fits",i);
    strcat(file,string);
    status=0;

    fits_open_file(&fptr,file,READONLY,&status);
    if(status)fits_die("data file error",status);
    fits_get_img_param(fptr,2,&bitpix,&naxis,naxes,&status);
    //check for consistent sizes
    if(nbfile){
      if(naxes[0] != totsizx || naxes[1] != totsizy){
        die("Incompatible images sizes!\n");}
    }
      totsizx=naxes[0];
      totsizy=naxes[1];
      //bias border
      oldbsfl=1;
      fits_read_key(fptr,TSTRING,"FILENAME",string,line,&status);
      printf("file: %s\n",string);
      fits_read_key(fptr,TINT,"OVERSCAN",&overscan,line,&status);
      if(status){
        oldbsfl=0;
        status=0;
        fits_read_key(fptr,TINT,"NOVERSCN",&overscan,line,&status);}
      if(status) overscan=0;
      status=0;
      fits_read_key(fptr,TINT,"BIASLINS",&biaslins,line,&status);
      if(status){
        oldbsfl=0;
        status=0;
        fits_read_key(fptr,TINT,"NBIASLNS",&biaslins,line,&status);}
      if(status)biaslins=0;
      xsize=naxes[0]-overscan;
      ysize=naxes[1]-biaslins;
      xsz=xsize;
      ysz=ysize;
      biasx0=xsz+3;
      biasx1=xsize+overscan-3;
      biasy0=ysz+3;
      biasy1=ysize+biaslins-3;
      nelem=naxes[0]*naxes[1];
      status=0;
      if(i==1){
        dimage=(float *) malloc(sizeof(float)*naxes[0]*naxes[1]);
        darray=(float **) malloc(sizeof(float *)*naxes[1]);}
      for(j=0;j<naxes[1];j++) *(darray+j)=dimage+naxes[0]*j;
      status=0;
      fits_read_img(fptr,TFLOAT,firstel,nelem,&nulvl,dimage,&anynl,&status);
      if(status) fits_die("Error opening data file",status);

      //subtract bias
      if(nbfile){
        for(j=0;j<naxes[1];j++){
          for(k=0;k<naxes[0];k++) *(*(darray+j)+k)-= *(*(barray+j)+k);}
      }
      if(biaslins&&overscan){
        subbias(darray,naxes,biasx0,biasx1,biasy0,biasy1,xsz,ysz);}

      //correct for gain
      for(j=0;j<ysz;j++){
        for(k=0;k<xsz;k++){
          *(*(darray+j)+k)=gain[i]*((float)(*(*(darray+j)+k)));}
      }

      //add to output buffer
      i1=indewinfo.xz[i]/xbin+indewinfo.sx[i]-1;
      i2= (i1==0) ? indewinfo.xchip: -1;
      i3=indewinfo.yz[i]/ybin+indewinfo.sy[i]-1;
      i4= (i3==0) ? indewinfo.ychip : -1;
      fxpxl=indewinfo.xloc[i]*indewinfo.xchip;
      fypxl=indewinfo.yloc[i]*indewinfo.ychip;
      jl=indewinfo.yloc[i]*indewinfo.ychip-1;
      for(j=i3;j!=i4;j+=indewinfo.sy[i]){
        jl++;
        kl=indewinfo.xloc[i]*indewinfo.xchip-1;
        for(k=i1;k!=i2;k+=indewinfo.sx[i]){
          kl++;
          *(*(outarray+jl)+kl)=*(*(darray+j)+k);}
      }
  }

  //write output file

  strcpy(file,"!");
  strcat(file,datadir);
  strcat(file,dfile);
  strcat(file,".fits");
  status=0;
  naxes[0]=indewinfo.xarsize/xbin;
  naxes[1]=indewinfo.yarsize/ybin;
  bitpix=-32;
  status=0;
  fits_create_file(&fptr_o,file,&status);
  if(status)fits_die("Output file error",status);
  fits_create_img(fptr_o,bitpix,naxis,naxes,&status);
  nelem=naxes[0]*naxes[1];
  status=0;
  fits_write_img(fptr_o,TFLOAT,firstel,nelem,outimage,&status);

  //fix up header

  //   fits_delete_record(fptr_o,6,&status);
  fits_delete_record(fptr_o,7,&status);
  fits_delete_record(fptr_o,7,&status);

  status=0;
  ncard=0;
  while(1){
    ncard++;
    fits_read_record(fptr,ncard,cardline,&status);
    if(status){
      printf("Error reading input file header\n");
      status=0;
      break;}
    if(strstr(cardline,"END")==cardline) break;
    if(strstr(cardline,"NAXIS")!=NULL) continue;
    if(strstr(cardline,"SIMPLE")!=NULL) continue;
    if(strstr(cardline,"BITPIX")!=NULL) continue;
    if(strstr(cardline,"BSCALE")!=NULL) continue;
    if(strstr(cardline,"BZERO")!=NULL) continue;
    fits_write_record(fptr_o,cardline,&status);
    if(status) die("Error writing output file");}

  if(oldbsfl){
    fits_update_key(fptr_o,TINT,"BIASLINS",&zero,"biaslines",&status);
    if(status) die("Error writing output file");
    fits_update_key(fptr_o,TINT,"OVERSCAN",&zero,"overscan",&status);
    if(status) die("Error writing output file");
  } else {
    fits_update_key(fptr_o,TINT,"NBIASLNS",&zero,"biaslines",&status);
    if(status) die("Error writing output file");
    fits_update_key(fptr_o,TINT,"NOVERSCN",&zero,"overscan",&status);
    if(status) die("Error writing output file");
  }


  /* update datasec */
  status=0;
  fits_read_key(fptr,TSTRING,"DATASEC",&datasec,line,&status);
  if(!status) {
    secstr = strtok(datasec,DELIMITERS);
    val[0] = atol(secstr);
    for(l=1;l<4;l++){
       secstr = strtok(NULL,DELIMITERS);
       val[l] = atol(secstr);
       if(l==1) {
          val[l] = naxes[0];
       } else {
          val[l] = naxes[1] - val[l] + 1;
       }
    }
    sprintf(ndatasec,"[%i:%i,%i:%i]",val[0],val[1],val[3],val[2]);
    fits_update_key(fptr_o,TSTRING,"DATASEC",&ndatasec,"new datasec",&status);
    if(status) die("Error writing output file");
  }
  status=0;

  /* remove biassec */
  fits_delete_key(fptr_o,"BIASSEC",&status);
  if(status) die("Error writing output file");

  fits_write_key(fptr_o,TINT,"NMOSAIC",&indewinfo.nchip,"Number of chips mosaiced",&status);
  fits_close_file(fptr,&status);
  fits_close_file(fptr_o,&status);
  printf("file %s written\r",file);
  fflush(stdout);
  printf("\n");
  return 0;}
