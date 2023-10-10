/*****************************************************************************\
 *                                                                            *
 *   SUMSPEC  sums 1-d and 2-d spectral stacks, with optional CR cleaning     *
 *            for 2-d stacks                                                  *
 *                                                                            *
 *   VERSION  04 Jun 2004                                                     *
 *                                                                            *
\*****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fitsio.h"
#include "cosmos.h"

int main(int argc,char *argv[]){

  char     name[80],infile[20][80],outfile[80],flnm[1000],string[133],line[80],
           file[80],CHAR,lane[80],cardline[80];
  char     *DATA_DIR,*HOME,*COS_HOME,*append;
  float    expt,exptime[20],mlm,dla,dltaslt,dltlt,min_lambda,delta_lambda,
           exprat2,exprat,a,siglimit,dev,dlmt,maxdev,crpix1;
  int      bitpix,bpx,cnt,sln,nxs,anynl,slitnum,cntrline,naxis,size,nspc,nspec,
           cntrln[20],nslpx,nstk,clean,errvl,onedim,status,i,j,k,ninfile,nulvl,
           diag,nsp,shfld,shfl,nfl,cndwn,cnup,nsize,offst,nline[20],bth,*undef,
           niter,l,mxd,INT,indx,ngrow,ig,kg,andx,sltype,dcflag,dispaxis,nod,
	         Pord,ncard,dim2;
  long     firstel, firstelem[3],naxes[3],noxes[3];
  int      *sze;
  float    *ipoint,*epoint,**instack,*outstack,*outerr,FLOAT,*FLOATP,*exptot,
           *devtn;
  fitsfile *fptr[20],*outptr;

  errvl=clean=diag=0;
  dcflag = sltype = 0;
  crpix1 = 1.;
  dispaxis = 1;

  firstelem[0]=firstelem[1]=firstelem[2]=1;
  firstel=1;
  undef=NULL;
  //directories
  COS_HOME=malloc(sizeof(CHAR)*80);
  HOME=malloc(sizeof(CHAR)*80);
  DATA_DIR=malloc(sizeof(CHAR)*80);
  DATA_DIR=getenv("COSMOS_IMAGE_DIR");
  if(DATA_DIR==NULL) die("COSMOS_IMAGE_DIR undefined!");
  strcat(DATA_DIR,"/");
  COS_HOME=getenv("COSMOS_HOME");
  if(COS_HOME==NULL) die("COSMOS_HOME undefined!");
  HOME=getenv("HOME");

  //get input files

  ninfile=0;
  if(argc<5){
    printf("Enter output file name: ");
    fgets(string,sizeof(string),stdin);
    sscanf(string,"%s",outfile);
    fflush(stdin);
    while(1){
      printf("Enter input files. CR=>end: ");
      fgets(string,sizeof(string),stdin);
      if(strlen(string)==1) break;
      sscanf(string,"%s",&infile[ninfile][0]);
      ninfile++;}
    if(ninfile==0) return 1;
    }
  else{
    i=1;
    while(i<argc){
      if(!strcmp(argv[i],"-o")){
	strcpy(outfile,argv[i+1]);
	i+=2;
	continue;}
      if(!strcmp(argv[i],"-d")){
	diag=1;;
	i+=2;
	continue;}
      strcpy(&infile[ninfile][0],argv[i]);
      i+=1;
      if(ninfile==20){
	printf("Maximum number of input files reached\n");
	break;}
      ninfile++;}
    if(!strcmp(outfile,"")) die("Usage sumspec -o outputfile inputfiles");}

  //get parameters
  strcpy(file,"sumspec");
  if(OpenCosParm(file)!=0) die("Cannot open sumspec parameter file");
  if(ReadParm_r("siglimit",&siglimit)==1) die("parameter file error");
  if(ReadParm_b("clean",&clean)==1) die("parameter file error");
  if(ReadParm_b("both",&bth)==1) die("parameter file error");
  if(ReadParm_i("iterations",&niter)==1) die("parameter file error");
  if(ReadParm_i("grow",&ngrow)==1) die("parameter file error");

  //open all files
  dim2=1;
  append=strstr(outfile,"_2spec");
  if(append==NULL){
     append=strstr(outfile,"_1spec");
     dim2=0;
     clean=0;}
  if(append==NULL) die("output file name must end with _1spec or _2spec");
  status=0;
  for(i=0;i<ninfile;i++){
    strcpy(flnm,DATA_DIR);
    strcat(flnm,infile[i]);
    strcat(flnm,append);
    strcat(flnm,".fits");
    fits_open_file(&fptr[i],flnm,READONLY,&status);
    if(status) fits_die(flnm,status);}
  strcpy(flnm,"!");
  strcat(flnm,DATA_DIR);
  strcat(flnm,outfile);
  strcat(flnm,".fits");
  fits_create_file(&outptr,flnm,&status);
  if(status) fits_die(flnm,status);
  instack=malloc(ninfile*sizeof(FLOATP));
  devtn=malloc(ninfile*sizeof(FLOAT));
  sze=malloc(ninfile*sizeof(INT));

  //read file, sum spectra

  //first hdu
  fits_read_key(fptr[0],TLOGICAL,"SHUFFLED",&shfld,string,&status);
  if(status){
    shfld=0;
    status=0;}
  fits_read_key(fptr[0],TINT,"N_SLITS",&nspec,string,&status);
  if(status)fits_die(infile[0],status);
  if(dim2){
    fits_read_key(fptr[0],TFLOAT,"D_SLIT",&dltaslt,string,&status);
    if(status)fits_die(infile[0],status);}
  fits_read_key(fptr[0],TFLOAT,"CRVAL1",&min_lambda,string,&status);
  if(status) fits_die(infile[0],status);
  fits_read_key(fptr[0],TFLOAT,"CDELT1",&delta_lambda,string,&status);
  if(status) fits_die(infile[0],status);
  fits_read_key(fptr[0],TFLOAT,"EXPTIME",&exptime[0],string,&status);
  if(status) fits_die(infile[0],status);
  expt=exptime[0];
  fits_read_key(fptr[0],TINT,"NOD",&nod,string,&status);
  status = 0;
  fits_get_img_param(fptr[0],3,&bitpix,&naxis,naxes,&status);
  if(naxis==3) errvl=1;
  if(clean && !errvl) die("Cannot clean CR's: no error matrix");
  if(status) fits_die(infile[0],status);
  //check other files for consistency
  for(j=1;j<ninfile;j++){
    fits_read_key(fptr[j],TLOGICAL,"SHUFFLED",&shfl,string,&status);
    if(status){
      shfl=0;
      status=0;}
    if(shfld != shfl) die("Input files have inconsistent shuffling");
    fits_read_key(fptr[j],TINT,"N_SLITS",&nspc,string,&status);
    if(status){
      strcpy(line,infile[j]);
      strcat(line," N_SLITS");
      fits_die(line,status);}
    if(nspc != nspec) die("Input files have inconsistent number of spectra");
    if(dim2){
      fits_read_key(fptr[j],TFLOAT,"D_SLIT",&dltlt,string,&status);
      if(status){
        strcpy(line,infile[j]);
        strcat(line," D_SLIT");
        fits_die(line,status);}
      if(dltlt != dltaslt)die("Input files have inconsistent slit sampling");}
    fits_read_key(fptr[j],TFLOAT,"CRVAL1",&mlm,string,&status);
    if(status){
      strcpy(line,infile[j]);
      strcat(line," CRVAL1");
      fits_die(line,status);}
    if(mlm != min_lambda) die("Input files have inconsistent lambda scales");
    fits_read_key(fptr[j],TFLOAT,"CDELT1",&dla,string,&status);
    if(status){
      strcpy(line,infile[j]);
      strcat(line," CDELT1");
      fits_die(line,status);}
    if(dla != delta_lambda)die("Input files have inconsistent lambda scales");
    fits_read_key(fptr[j],TFLOAT,"EXPTIME",&exptime[j],string,&status);
    if(status){
      strcpy(line,infile[j]);
      strcat(line," EXPTIME");
      fits_die(line,status);}
    expt+=exptime[j];
    }

  //each hdu
    for(nsp=0;nsp<nspec;nsp++){
    for(nfl=0;nfl<ninfile;nfl++){
      strcpy(lane,infile[nfl]);
      strcat(lane,":");
      fits_get_img_param(fptr[nfl],3,&bitpix,&naxis,naxes,&status);
      strcat(lane,"1");
      if(status) fits_die(lane,status);
      size=naxes[0]*naxes[1]*(1+errvl);
      nline[nfl]=naxes[1];
      *(sze+nfl)=size;
      instack[nfl]=malloc(size*sizeof(FLOAT));
      fits_read_key(fptr[nfl],TINT,"SLITNUM",&slitnum,string,&status);
      strcat(lane,"2");
      if(status) fits_die(lane,status);
      fits_read_key(fptr[nfl],TSTRING,"OBJECT",&name,string,&status);
      strcat(lane,"3");
      if(status) fits_die(lane,status);
      fits_read_key(fptr[nfl],TINT,"DISPAXIS",&dispaxis,string,&status);
      fits_read_key(fptr[nfl],TFLOAT,"CRPIX1",&crpix1,string,&status);
      fits_read_key(fptr[nfl],TINT,"DC-FLAG",&dcflag,string,&status);
      fits_read_key(fptr[nfl],TINT,"SLITTYPE",&sltype,string,&status);
      status = 0;

      //check file consistency
      if(!nfl){
        bpx=bitpix;
        nxs=naxis;
        noxes[0]=naxes[0];
        sln=slitnum;
        strcpy(line,name);}
      else{
        if(bpx!=bitpix || nxs!=naxis || noxes[0]!=naxes[0]){
          die("Input files have inconsistent parameters");}
        if(strcmp(line,name)) die("Input files have inconsistent names");}
      //read data
      if(dim2){
        fits_read_key(fptr[0],TINT,"SLITLEN",&nslpx,string,&status);
        if(status){
          nslpx=naxes[1];
          status=0;}
        fits_read_key(fptr[nfl],TINT,"CNTRLINE",&cntrln[nfl],string,&status);
        strcat(lane,"4");
        if(status) fits_die(lane,status);}
      fits_read_img(fptr[nfl],TFLOAT,firstel,size,&nulvl,instack[nfl],&anynl,
                    &status);
      strcat(lane,"5");
      if(status) fits_die(lane,status);
      if(errvl){                                  //calculate variance
        for(j=size/2+1;j<size;j++){
          ipoint=instack[nfl]+j;
          if(*ipoint>0) *(ipoint) *= *ipoint;}
      }
    }
    //align spectra

    if(dim2){
      cndwn=1000;
      cnup=1000;
      for(j=0;j<ninfile;j++){
        if(nline[j]-cntrln[j]<cnup) cnup=nline[j]-cntrln[j];
        if(cntrln[j]<cndwn) cndwn=cntrln[j];}
      noxes[1]=cnup+cndwn;}
    else{
      noxes[1]=1;}
    nsize=noxes[1]*noxes[0];
    outstack=calloc((1+errvl)*nsize,sizeof(FLOAT));
    if(clean || errvl){
      exptot=calloc((1)*nsize,sizeof(FLOAT));
      for(i=0;i<nsize;i++) *(exptot+i)=expt;}                  //total exp time

    //sum spectra,errors

    for(j=0;j<ninfile;j++){
      offst= dim2 ? naxes[0]*(cntrln[j]-cndwn) : 0;
      for(k=0;k<nsize;k++) *(outstack+k) += *(instack[j]+k+offst);
      if(errvl){
        for(k=0;k<nsize;k++){
          a=*(instack[j]+k+offst+*(sze+j)/2);
          if(a>0){                                       //is this a bad point?
            *(outstack+k+nsize) += a;}         //if not, add variance to total
          else{
            *(exptot+k)-=exptime[j];}    //if is, decrease total exp time
        }
      }
    }

    //clean CR's
    if(clean){
      for(l=0;l<niter;l++){
          for(k=0;k<nsize;k++){
            mxd=-1;
            maxdev=0;
            for(j=0;j<ninfile;j++){
              ipoint=instack[j]+naxes[0]*(cntrln[j]-cndwn)+k;
              epoint=ipoint+*(sze+j)/2;
              exprat=exptime[j]/(*(exptot+k));
              exprat2=exprat*exprat;
              if(*epoint<0) continue;             //bad point, already excluded
              dlmt=siglimit*sqrt(*epoint+(*(outstack+k+nsize)*exprat2));
              dev=(*ipoint-exprat*(*(outstack+k)))/dlmt;
              if(dev>1.){
                if(dev>maxdev){
                  maxdev=dev;
                  mxd=j;}
                continue;}
              if(bth && dev<-1.){
                if(-dev>maxdev){
                  maxdev=-dev;
                  mxd=j;}
                continue;}
            }
            if(mxd>-1){
              if(diag){
                *(outstack+k) = 1000*mxd;
                continue;}
              for(ig=-ngrow;ig<=ngrow;ig++){
                for(kg=-ngrow;kg<=ngrow;kg++){
                  andx=k+kg+naxes[0]*ig;
                  if(andx<0 || andx>nsize-1) continue;
                  indx=andx+naxes[0]*(cntrln[mxd]-cndwn);
                  if(indx<0 || indx>nsize-1) continue;
                  ipoint=instack[mxd]+indx;
                  if(*(ipoint+*(sze+mxd)/2)<0) continue;
                  *(outstack+andx) -= *ipoint;
                  *(outstack+andx+nsize) -= *(ipoint+*(sze+mxd)/2);
                  *(ipoint+*(sze+mxd)/2)=-1.;
                  *(exptot+andx) -= exptime[mxd];}
              }
            }
          }
        }
      }

    if(errvl){
      for(k=nsize;k<2*nsize;k++)*(outstack+k) = sqrt(*(outstack+k));
      noxes[2]=2;}
    if(clean){
      //renormalize for removed data
      for(k=0;k<nsize;k++){
        if(*(exptot+k)){
          exprat=expt/(*(exptot+k));
          *(outstack+k) *= exprat;
          *(outstack+k+nsize) *= exprat;}
        else{
          *(outstack+k)=0;
          *(outstack+k+nsize) = -9999;}
      }
    }
    //write summed spectrum

    fits_create_img(outptr,FLOAT_IMG,naxis,noxes,&status);
    if(status) fits_die("Error writing output file",status);
    //writing original header information to output file first
    status=0;
    ncard=0;
    while(1){
      ncard++;
      if(ncard>108) break;
      fits_read_record(fptr[0],ncard,cardline,&status);
      if(status){fits_die("Error reading file header",status);}
      if(strstr(cardline,"END")==cardline) break;
      if(strstr(cardline,"NAXIS")!=NULL) continue;
      if(strstr(cardline,"SIMPLE")!=NULL) continue;
      if(strstr(cardline,"BITPIX")!=NULL) continue;
      if(strstr(cardline,"BSCALE")!=NULL) continue;
      if(strstr(cardline,"BZERO")!=NULL) continue;
      if(strstr(cardline,"EXPTIME")!=NULL) continue;
      status=0;
      fits_write_record(outptr,cardline,&status);
      if(status) die("Error writing output file");}

    // New header
    if(!nsp){
      if(dim2){
        fits_write_key(outptr,TLOGICAL,"SHUFFLED",&shfld,
                       "does data contain shuffled region",&status);
        fits_write_key(outptr,TINT,"NOD",&nod,"nod distance in pixels",&status);
        fits_write_key(outptr,TFLOAT,"D_SLIT",&dltaslt,
                       "slit interval in arcsec",&status);}
      fits_write_key(outptr,TINT,"N_SLITS",&nspec,"Number of spectra",
                     &status);
      fits_write_key(outptr,TFLOAT,"EXPTIME",&expt,
                     "Exposure time",&status);
      if(status) fits_die("Error writing output file",status);
    }

    fits_write_key(outptr,TINT,"DISPAXIS",&dispaxis,"",&status);
    fits_write_key(outptr,TSTRING,"CTYPE1","LINEAR","",&status);
    if(dim2) fits_write_key(outptr,TSTRING,"CTYPE2","LINEAR","",&status);
    fits_write_key(outptr,TSTRING,"WAT0_001","system=world","",&status);
    fits_write_key(outptr,TSTRING,"WAT1_001",
                "wtype=linear label=Wavelength units=Angstroms","",&status);
    if(dim2) fits_write_key(outptr,TSTRING,"WAT2_001","wtype=linear","",&status);

    fits_write_key(outptr,TFLOAT,"CRVAL1",&min_lambda,"",&status);
    fits_write_key(outptr,TFLOAT,"CDELT1",&delta_lambda,"",&status);
    fits_write_key(outptr,TFLOAT,"CD1_1",&delta_lambda,"",&status);
    fits_write_key(outptr,TFLOAT,"CRPIX1",&crpix1,"",&status);
    fits_write_key(outptr,TINT,"DC-FLAG",&dcflag,"",&status);
    if(dim2){
      fits_write_key(outptr,TINT,"SLITLEN",&noxes[1],"slit length in pixels",
                     &status);
      fits_write_key(outptr,TINT,"CNTRLINE",&cndwn,"Spectrm central row",
                    &status);}
    fits_write_key(outptr,TINT,"SLITTYPE",&sltype,"Type of aperture",&status);
    fits_write_key(outptr,TINT,"SLITNUM",&slitnum,"Slit number",&status);
    fits_write_key(outptr,TSTRING,"OBJECT",&name,"Name of object",&status);



    // New header

    if(status) fits_die("Error writing output file",status);
    nsize*=(1+errvl);
    fits_write_pix(outptr,TFLOAT,firstelem,nsize,outstack,&status);
    if(status) fits_die(outfile,status);
    printf("Processing slit %d\r",nsp);
    fflush(stdout);
    ffflus(outptr,&status);
    for(k=0;k<ninfile;k++){
      free(instack[k]);
      fits_movrel_hdu(fptr[k],1,undef,&status);}
    if(clean) free(exptot);
    free(outstack);}
  ffflus(outptr,&status);
  fits_close_file(outptr,&status);
  return 0;}
