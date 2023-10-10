/*****************************************************************************\
 *                                                                                                                                        *
 *  BIASFLAT does bias and flatfiled corrections to ccd data sets                                             *
 *                                                                                                                                        *
 *  VERSION  25 Oct 2017                                                                                                  *
 *                                                                                                                                        *
\*****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fitsio.h"
#include "cosmos.h"


int main(int argc, char *argv[]){

  int      xsize,ysize,biasx0,biasx1,biasy0,biasy1,bitpix,nbfile,nffile,*INTP,
           nulvl,totsizx,totsizy,nflatx,nflaty,n,naxis,nelem,i,j,k,nosub,nms,
           oldbsfl,ncard,overscan,zero,biaslins,xsz,ysz,status,ndfile,INT,
           nchips,imacs,anynl,nkeys,end;
  float    FLOAT,gain[9],avg;
  float    *FLOATP;
  char     CHAR,*datadir,dfile[50][80],file[1000],bfile[80],ffile[80],dewar[10],
              peed[10],string[133],cardline[81],subrastr[81],binline[81],line[133],
             *homedir,*pardir,speed[10],*headstr,*epoint;
  long     naxes[3],noxes[3],firstel,firstelem[3];
  float    *bimage,*dimage,**barray,**darray,*fimage,**farray,*outimage,
            **outarray;
  fitsfile *fptr,*bptr,*fptr_o;
  fitsdef  fitsinfo;
  dewdat   dewinfo;

  avg=1.;
  imacs=nulvl=0;
  zero=0;
  firstelem[0]=firstelem[1]=firstelem[2]=1;
  firstel=1;
  datadir=malloc(sizeof(CHAR)*80);
  epoint=getenv("COSMOS_IMAGE_DIR");
  if(epoint==NULL){
    printf("COSMOS_IMAGE_DIR undefined!\n");
    return 1;}
  strcpy(datadir,epoint);
  strcat(datadir,"/");
  homedir=malloc(sizeof(CHAR)*80);
  epoint=getenv("COSMOS_HOME");
  if(epoint==NULL){
    printf("COSMOS_HOME undefined!\n");
    return 1;}
  strcpy(homedir,epoint);
  strcat(homedir,"/");

  //get parameters

  if(OpenCosParm("biasflat")) die("Cannot open parameter file");
  if(ReadParm_s("dewar",dewar)==1) die("parameter file error");
  strcpy(string,"unable to read dewar definitionn file ");
  strcat(string,dewar);
  if(Readcdf(dewar)) die(string);
  Getchipdat(&dewinfo);
  nchips=dewinfo.nchip;
  if(!strcmp(dewar,"SITE"))imacs=1;
  if(ReadParm_s("speed",speed)==1) die("parameter file error");
  if(nchips==1){
    gain[1]=1;}
  else{
    if(ReadGain(dewar,speed,gain)==1) die("error reading gain values");}
  ndfile=nffile=nbfile=0;
  //read input data
  if(argc>1){
    i=1;
    while(i<argc){
      if(!strcmp(argv[i],"-b")){
	if(strcasecmp(argv[i+1],"none")){
	  nbfile=1;
	  strcpy(bfile,argv[i+1]);}
	i+=2;
	continue;}
      if(!strcmp(argv[i],"-f")){
	nffile=1;
	strcpy(ffile,argv[i+1]);
	i+=2;
	continue;}
      strcpy(&dfile[ndfile][0],argv[i]);
      ndfile++;
      i+=1;}
    }
  else{
    printf("enter bias file name: ");
    fgets(string,133,stdin);
    if(strncmp(string,"\n",1)){
      nbfile=1;
      sscanf(string,"%s",bfile);}
    printf("enter flatfield file name: ");
    fgets(string,133,stdin);
    if(strncmp(string,"\n",1)){
      nffile=1;
      sscanf(string,"%s",ffile);}
    }
    if(!ndfile){
      printf("\nenter files to process; CR => end\n\n");
      while(1){
	fgets(string,133,stdin);
	if(!strncmp(string,"\n",1)) break;
	sscanf(string,"%s",&dfile[ndfile][0]);
	ndfile++;}
      if(!ndfile){
	printf("No data files specified!\n");
	return 1;}
      }

  //if no bias and/or flatfiled files specified, look in parameter file

  if(!nbfile){
    if(ReadParm_s("bias_file",bfile)==1){
      printf("parameter file error bf!\n");
      return 1;}
    if(strcasecmp(bfile,"none") && strncmp(bfile,"\0",1))nbfile=1;}
  if(!nffile){
    if(ReadParm_s("flat_file",ffile)==1){
      printf("parameter file error ff!\n");
      return 1;}
    if(strcasecmp(ffile,"none") && strncmp(ffile,"\0",1))nffile=1;}

  //Work chip by chip

  for(i=1;i<=nchips;i++){

    //read in bias file if exist

    if(nbfile){
      strcpy(file,datadir);
      strcat(file,bfile);
      if(nchips>1){
        strcat(file,"c");
        sprintf(string,"%d",i);
        strcat(file,string);}
      strcat(file,".fits");
      status=0;
      status=OpenFitsFile(file,&bptr,&fitsinfo);
      if(status){
        printf("file %s does not exist! \n",file);
        continue;}
      if(!imacs && nchips==1){
        fits_read_key(bptr,TINT,"NMOSAIC",&nms,line,&status);
        if(status) die("bias file was not properly stitched");}
      bitpix=fitsinfo.bitpix;
      totsizx=naxes[0]=fitsinfo.naxes[0];
      totsizy=naxes[1]=fitsinfo.naxes[1];
      nelem=totsizx*totsizy;
    //bias border
    overscan=fitsinfo.overscan;
    biaslins=fitsinfo.biaslins;
    status=0;
    if(i==1){
      bimage=malloc(sizeof(FLOAT)*naxes[0]*naxes[1]);
      barray=malloc(sizeof(FLOATP)*naxes[1]);
      for(j=0;j<naxes[1];j++) *(barray+j)=bimage+naxes[0]*j;}
    fits_read_img(bptr,TFLOAT,firstel,nelem,&nulvl,bimage,&anynl,&status);
    if(status)fits_die("Bias file error",status);
    fits_close_file(bptr,&status);}

    //read in flatfield file if exists

    if(nffile){
      strcpy(file,datadir);
      strcat(file,ffile);
      if(nchips>1){
        strcat(file,"_c");
        sprintf(string,"%d",i);
        strcat(file,string);}
      strcat(file,".fits");
      status=0;
      status=OpenFitsFile(file,&fptr,&fitsinfo);
      if(status) fits_die(file,status);
      if(!imacs && nchips==1){
        fits_read_key(fptr,TINT,"NMOSAIC",&nms,line,&status);
        if(status) die("flat file was not properly stitched");}
      bitpix=fitsinfo.bitpix;
      nflatx=noxes[0]=fitsinfo.naxes[0];
      nflaty=noxes[1]=fitsinfo.naxes[1];
      nelem=nflatx*nflaty;
      status=0;
      //bias border
      overscan=fitsinfo.overscan;
      biaslins=fitsinfo.biaslins;
      fimage=malloc(sizeof(FLOAT)*noxes[0]*noxes[1]);
      farray=malloc(sizeof(FLOATP)*noxes[1]);
      for(j=0;j<noxes[1];j++) *(farray+j)=fimage+noxes[0]*j;
      status=0;
      fits_read_img(fptr,TFLOAT,firstel,nelem,&nulvl,fimage,&anynl,&status);
      if(status)fits_die("Flat field file error",status);
      fits_close_file(fptr,&status);
      xsize=noxes[0]-overscan;
      ysize=noxes[1]-biaslins;
      //subtract bias if necessary
      if(overscan && biaslins){
	if(nbfile){
	  for(j=0;j<naxes[1];j++){
	    for(k=0;k<naxes[0];k++) *(*(farray+j)+k)-= *(*(barray+j)+k);}
	  }
	xsz=xsize;
	ysz=ysize;
	biasx0=xsz+3;
	biasx1=(xsize+overscan)-3;
	biasy0=ysz+3;
	biasy1=(ysize+biaslins)-3;
	subbias(farray,naxes,biasx0,biasx1,biasy0,biasy1,xsz,ysz);}
      //normalize to unity
      /*      avg=0.;
      for(j=0;j<ysize;j++){
	for(k=0;k<xsize;k++){
	  avg+= *(*(farray+j)+k);}
      }
      avg/=xsize*ysize;}*/
      }
    //Read in data files, correct

    for(n=0;n<ndfile;n++){
      strcpy(file,datadir);
      strcat(file,dfile[n]);
      //addbar(file);
      if(nchips>1){
        strcat(file,"c");
        sprintf(string,"%d",i);
        strcat(file,string);}
      strcat(file,".fits");
      status=0;
      fits_open_file(&fptr,file,READONLY,&status);
      if(status)fits_die(file,status);
      if(!imacs && nchips==1){
        fits_read_key(fptr,TINT,"NMOSAIC",&nms,line,&status);
        if(status) die("data file was not properly stitched");}
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
      if(i==1 && n==0){
	dimage=malloc(sizeof(FLOAT)*naxes[0]*naxes[1]);
	darray=malloc(sizeof(FLOATP)*naxes[1]);
	for(j=0;j<naxes[1];j++) *(darray+j)=dimage+naxes[0]*j;
    outimage=malloc(sizeof(FLOAT)*xsz*ysz);
    outarray=malloc(sizeof(FLOATP)*ysz);
	for(j=0;j<ysz;j++) *(outarray+j)=outimage+xsz*j;}
      status=0;
      fits_read_img(fptr,TFLOAT,firstel,nelem,&nulvl,dimage,&anynl,&status);
      if(status)fits_die("data file error",status);

      //create new output file

      strcpy(file,"!");
      strcat(file,datadir);
      strcat(file,dfile[n]);
      if(nffile){
        strcat(file,"_f");}
      else{
        strcat(file,"_b");}
      if(nchips>1){
        strcat(file,"_c");
        sprintf(string,"%d",i);
        strcat(file,string);}
      strcat(file,".fits");
      status=0;

      //subtract bias

      if(nbfile){
	for(j=0;j<naxes[1];j++){
	  for(k=0;k<naxes[0];k++) *(*(darray+j)+k)-= *(*(barray+j)+k);}
        }
      if(biaslins&&overscan){
	subbias(darray,naxes,biasx0,biasx1,biasy0,biasy1,xsz,ysz);}

      //if flatfiled, flatten

      if(nffile){
	for(j=0;j<ysz;j++){
	  for(k=0;k<xsz;k++){
         if (*(*(farray+j)+k) > 0.0) {
	        *(*(outarray+j)+k)=((float)(*(*(darray+j)+k)))*avg*gain[i]/
                               (float)(*(*(farray+j)+k)); }
        else
            { *(*(outarray+j)+k)= 0.0; }
	  }
        }
      }else{
	for(j=0;j<ysz;j++){
	  for(k=0;k<xsz;k++){
	    *(*(outarray+j)+k)=gain[i]*((float)(*(*(darray+j)+k)));}
	  }
      }

      //write output file

     status=0;
      naxes[0]=xsz;
      naxes[1]=ysz;
      bitpix=-32;
      status=0;
      fits_create_file(&fptr_o,file,&status);
      if(status)fits_die("Output file error",status);
      fits_create_img(fptr_o,bitpix,naxis,naxes,&status);
      //      fits_write_record(fptr_o,binline,&status);
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
	  printf("Error reading input file header, line %d\n",ncard);
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
	if(status) die("Error writing output file");}
      else{
	fits_update_key(fptr_o,TINT,"NBIASLNS",&zero,"biaslines",&status);
	if(status) die("Error writing output file");
	fits_update_key(fptr_o,TINT,"NOVERSCN",&zero,"overscan",&status);
	    if(status) die("Error writing output file");}
        fits_close_file(fptr,&status);
        fits_close_file(fptr_o,&status);
         printf("file %s written\r",file);
        fflush(stdout);}
    }
  printf("\n");
  return 0;}
