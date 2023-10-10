/*****************************************************************************\
 *                                                                            *
 *   1DSPEC  extracts 1-dimensional spectra from 2-dimensional spectra        *
 *           using either simple summation or optimal extraction              *                                       *
 *                                                                            *
 *   VERSION  20 Feb 2018                                                     *
 *                                                                            *
\*****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fitsio.h"
#include "cosmos.h"

int main(int argc,char *argv[]){

  int       optimal,clean,align,hwidth,edge,shfld,nspec,bitpix,naxis,errvl,INT,
            status,o_value,nsp,i,j,k,size,nline,slitnum,noff,isloff,
            nulvl,frow,lrow,slght,arr_len,loop,ncosmic,kmax,anynl,Pord,ncard;
  int       *undef,*INTP,*irnk;
  char      name[80],infile[80],outfile[80],flnm[1000],string[133],line[80],
            file[80],CHAR,lane[80],cardline[80],object[80];
  char      *DATA_DIR,*HOME,*COS_HOME,*append;
  long     firstel, firstelem[3],naxes[3],noxes[3],n_lampix;
  float     min_lambda,delta_lambda,FLOAT,pp,wtot,ftot,var,vmax,gain,sumy,omax,
            yval;
  float     *FLOATP,*stack,*estack,**spectrum,**espectrum,*Pstack,**Parray,
            *Vstack,**Varray,*lamind,**Pcoef,*Pcstack,*Tflux,*instack,
            *Wstack,**Warray,*CRstack,**CRarray,*Tsigma,*rank;
  double    *xlist,*ylist,**elist,DOUBLE,*DOUBLEP;
  fitsfile  *fptr,*outptr;

  firstelem[0]=firstelem[1]=firstelem[2]=1;
  firstel=1;
  undef=NULL;
  Pord=2;
  //directories
  DATA_DIR=malloc(80*sizeof(CHAR));
  HOME=malloc(80*sizeof(CHAR));
  COS_HOME=malloc(80*sizeof(CHAR));
  DATA_DIR=getenv("COSMOS_IMAGE_DIR");
  if(DATA_DIR==NULL) die("COSMOS_IMAGE_DIR undefined!");
  strcat(DATA_DIR,"/");
  COS_HOME=getenv("COSMOS_HOME");
  if(COS_HOME==NULL) die("COSMOS_HOME undefined!");
  HOME=getenv("HOME");

  //input file
  strcpy(infile,DATA_DIR);
  strcat(infile,argv[1]);
  strcpy(outfile,"!");
  strcat(outfile,infile);
  strcat(infile,"_2spec.fits");
  strcat(outfile,"_1spec.fits");

  //get parameters
  strcpy(file,"1dspec");
  if(OpenCosParm(file)!=0) die("Cannot open 1dspec parameter file");
  if(ReadParm_b("optimal",&optimal)==1) die("parameter file error");
  if(ReadParm_b("clean",&clean)==1) die("parameter file error");
  if(ReadParm_b("align",&align)==1) die("parameter file error");
  if(ReadParm_i("hwidth",&hwidth)==1) die("parameter file error");
  if(ReadParm_i("edge",&edge)==1) die("parameter file error");
  if(ReadParm_r("gain",&gain)==1 && optimal==1) die("parameter file error: gain");


  //open files
  status=0;
  fits_open_file(&fptr,infile,READONLY,&status);
  if(status) fits_die(&infile,status);
  fits_read_key(fptr,TLOGICAL,"SHUFFLED",&shfld,string,&status);
  if(status){
    shfld=0;
    status=0;}
  fits_read_key(fptr,TINT,"N_SLITS",&nspec,string,&status);
  if(status)fits_die(&infile[0],status);
  fits_read_key(fptr,TFLOAT,"CRVAL1",&min_lambda,string,&status);
  if(status) fits_die(&infile[0],status);
  fits_read_key(fptr,TFLOAT,"CDELT1",&delta_lambda,string,&status);
  if(status) fits_die(&infile[0],status);
  if(status) fits_die(&infile[0],status);
  fits_get_img_param(fptr,3,&bitpix,&naxis,naxes,&status);
  if(naxis==3) errvl=1;
  if(clean && !errvl) die("Cannot clean CR's: no error matrix");
  if(status) fits_die(&infile,status);
  //output file
  fits_create_file(&outptr,outfile,&status);
  if(status) fits_die(&outfile,status);
  n_lampix=naxes[0];
  lamind = malloc(n_lampix*sizeof(FLOAT));
  for(i=0;i<n_lampix;i++) lamind[i]=(float) i;
  Tflux = calloc(n_lampix,sizeof(FLOAT));
  Tsigma=  calloc(n_lampix,sizeof(FLOAT));
  if(align){
     rank= calloc(n_lampix,sizeof(FLOAT));
     irnk=calloc(n_lampix,sizeof(INT));}
  if(optimal){
    o_value = Pord <3 ? 5 : Pord+2;
    xlist=malloc(sizeof(DOUBLE)*2*(o_value));
    ylist=malloc(sizeof(DOUBLE)*(o_value));
    elist=malloc(sizeof(DOUBLEP)*(o_value));
    for(i=0;i<o_value;i++){
      *(elist+i)=malloc(sizeof(DOUBLE)*(o_value));}
    }

  //loop through spectra, extract

  for(nsp=0;nsp<nspec;nsp++){
    fits_get_img_param(fptr,3,&bitpix,&naxis,naxes,&status);
    strcat(lane,"1");
    if(status) fits_die(lane,status);
    size=naxes[0]*naxes[1]*(1+errvl);
    nline=naxes[1];
    noxes[0]=naxes[0];
    noxes[1]=1;
    noxes[2]=1+errvl;
    instack=malloc(size*sizeof(FLOAT));
    fits_read_key(fptr,TINT,"SLITNUM",&slitnum,string,&status);
    strcat(lane,"2");
    if(status) fits_die(lane,status);
    fits_read_key(fptr,TSTRING,"OBJECT",&name,string,&status);
    strcat(lane,"3");
    if(status) fits_die(lane,status);
    //read data
    fits_read_key(fptr,TINT,"SLITLEN",&noff,string,&status);
    if(status){
      noff=naxes[1];
      status=0;}
    fits_read_key(fptr,TINT,"CNTRLINE",&isloff,string,&status);
    strcat(lane,"4");
    if(status) fits_die(lane,status);
    fits_read_img(fptr,TFLOAT,firstel,size,&nulvl,instack,&anynl,
                  &status);
    strcat(lane,"5");
    if(status) fits_die(lane,status);
    //setup arrays
    spectrum=malloc(naxes[1]*sizeof(FLOATP));
    for(i=0;i<naxes[1];i++) spectrum[i]=instack+i*naxes[0];
    if(errvl){                                  //calculate variance
      espectrum=malloc(naxes[1]*sizeof(FLOATP));
      for(i=0;i<naxes[1];i++) espectrum[i]=instack+(i+naxes[1])*naxes[0];}

    //extract

    if(optimal){
      frow= hwidth>isloff-edge ? edge : isloff-hwidth;
      lrow= isloff+hwidth>noff-edge-1 ? noff-edge-1: isloff+hwidth;
      slght=lrow-frow+1;
      arr_len=n_lampix*slght;
      // S,W array
      Vstack=(float *) calloc(arr_len,sizeof(FLOAT));
      Varray=(float **) calloc(slght,sizeof(FLOATP));
      for(i=0;i<slght;i++) Varray[i]=Vstack+i*n_lampix;
      Wstack=(float *) calloc(arr_len,sizeof(FLOAT));
      Warray=(float **) calloc(slght,sizeof(FLOATP));
      for(i=0;i<slght;i++) Warray[i]=Wstack+i*n_lampix;
      // CR array
      CRstack=(int *) calloc(arr_len,sizeof(INT));
      CRarray=(int **) calloc(slght,sizeof(INTP));
      for(i=0;i<slght;i++) CRarray[i]=CRstack+i*n_lampix;
      // P array
      Pstack=(float *) calloc(arr_len,sizeof(FLOAT));
      Parray=(float **) calloc(slght,sizeof(FLOATP));
      for(i=0;i<slght;i++) Parray[i]=Pstack+i*n_lampix;
      // coef array
      Pcstack= (float *) calloc(slght*(Pord+1),sizeof(FLOAT));
      Pcoef=(float **) calloc(slght,sizeof(FLOATP));
      for(i=0;i<slght;i++) Pcoef[i]=Pcstack+i*(Pord+1);

      //initial values of total flux, variance, pixel variance, weights,mask
      for(j=0;j<n_lampix;j++){
        Tflux[j]=0;
        for(k=0;k<slght;k++){
          *(CRarray[k]+j) = (*(espectrum[k+frow]+j) <=0) ? 0 : 1;
          *(Varray[k]+j)= pow(*(espectrum[k+frow]+j),2);
          Tflux[j]+= *(spectrum[k+frow]+j)*(*(CRarray[k]+j));}
        }
      //loop entire fit

      for(loop=0;loop<4;loop++){
        // clean CRs fit line profiles
        for(j=0;j<n_lampix;j++){
            //reset CR array
          for(k=0;k<slght;k++){
            *(CRarray[k]+j) = (*(espectrum[k+frow]+j) <=0 || Tflux[k]==0.) ? 0 : 1;}
          }
        while(1){
          ncosmic=0;
          for(j=0;j<n_lampix;j++){
            for(k=0;k<slght;k++){
              if(*(CRarray[k]+j)==1 && *(Varray[k]+j)>0. &&Tflux[j]!=0.){
                *(Parray[k]+j)= *(spectrum[k+frow]+j)/Tflux[j];
                *(Warray[k]+j)=pow(Tflux[j],2)/(*(Varray[k]+j));}
              else{
                *(Parray[k]+j)=0.0;
                *(Warray[k]+j)=0.0;}
              }
            //renormalize
            pp=0;
            for(k=0;k<slght;k++){
               pp+= *(Parray[k]+j);
               if(isnan(*Parray[k]+j)) printf("nan %d %d %d %d\n",slitnum,loop,k,j);}

            if(pp>0) for(k=0;k<slght;k++) *(Parray[k]+j)/=pp;
            }
          for(k=0;k<slght;k++){
            plyfit_w(lamind,Parray[k],Warray[k],n_lampix,Pord,Pcoef[k],xlist,
                     ylist,elist);
            if(clean){
              for(j=0;j<n_lampix;j++){
                //if subtracted n&s data, reject high or low pixels
                if(shfld){
                  if(*(CRarray[k]+j) && pow(*(spectrum[k+frow]+j)-polyvalue(
                     (float)j,Pcoef[k],Pord)*Tflux[j],2)>16.*(*(Varray[k]+j))){
                       ncosmic++;
                       *(CRarray[k]+j)=0;}
                     }
                //if  normal data only reject high pixels
                else{
                  if(*(CRarray[k]+j) && -polyvalue((float)j,Pcoef[k],Pord)*
                     Tflux[j]+(*(spectrum[k+frow]+j))>4.*sqrt(*(Varray[k]+j))){
                       ncosmic++;
                       *(CRarray[k]+j)=0;}
                     }
                   }
                 }
               }
               break;
          if(ncosmic==0) break;}
        //construct cleaned profile, total spectrum, new variances
        for(j=0;j<n_lampix;j++){
            //reset CR array
            for(k=0;k<slght;k++){
              *(CRarray[k]+j) = (*(espectrum[k+frow]+j) <=0) ? 0 : 1;}
            while(1){
              wtot=ftot=0;
              for(k=0;k<slght;k++){
                pp=polyvalue((float)j,Pcoef[k],Pord);
                *(Varray[k]+j)= pow(*(espectrum[k+frow]+j),2)+(*(CRarray[k]+j)*
                                 (pp*Tflux[j]-(*(spectrum[k+frow]+j))))/gain;
                if(*(Varray[k]+j)>0.){
                  ftot+=pp*(*(spectrum[k+frow]+j))*(*(CRarray[k]+j))/(*(Varray[k]+j));
                  wtot+=pp*pp*(*(CRarray[k]+j))/(*(Varray[k]+j));}
              }
              if(wtot>0){
                Tflux[j]= ftot/wtot;
                Tsigma[j]=sqrt(1/wtot);}
              else{
                Tsigma[j]=0.;
                for(k=0;k<slght;k++){
                  *(CRarray[k]+j)=0;}
                }
              if(clean){
                vmax=0;
                for(k=0;k<slght;k++){
                  var=*(CRarray[k]+j)*(*(spectrum[k+frow]+j)-polyvalue(
                      (float)j,Pcoef[k],Pord)*Tflux[j])/sqrt(*(Varray[k]+j));
                  if((shfld && fabs(var)>vmax) || (!shfld && var>vmax)){
                    vmax=var;
                    kmax=k;}
                  }
                if(vmax>5.){
                  *(CRarray[kmax]+j)=0;
                  continue;}
                }
              break;}
            }
          }//end of iteration
      free(Wstack);
      free(Warray);
      free(CRstack);
      free(CRarray);
      free(Pcstack);
      free(Pcoef);
      free(Pstack);
      free(Parray);
      free(Vstack);
      free(Varray);}

    //simple extraction
    else{
      if(align){
        k=0;
        for(i=0;i<n_lampix;i++){
          yval=0;
          for(j=edge;j<noff-edge;j++){
            if(*(spectrum[j]+i) > yval){
              yval = *(spectrum[j]+i);
              omax= (float) j;}
            }
          if(yval>0.){
            order(&k,&omax,irnk,rank);
            k++;}
          }
        isloff=(int) (rank[k/2]+0.5);}
      frow= hwidth>isloff-edge ? edge : isloff-hwidth;
      lrow= isloff+hwidth>noff-edge-1 ? noff-edge-1: isloff+hwidth;
      slght=lrow-frow+1;
      for(i=0;i<n_lampix;i++){
        Tflux[i]=Tsigma[i]=0.;
        for(j=frow;j<lrow+1;j++){
          Tflux[i]+=*(spectrum[j]+i);
          Tsigma[i]+=pow(*(espectrum[j]+i),2);}
        Tsigma[i]=sqrt(Tsigma[i]);}
      }
    //write output file
    firstelem[2]=1;
    fits_create_img(outptr,FLOAT_IMG,naxis,noxes,&status);
    if(status) fits_die("Error writing output file",status);
    //writing original header information to output file first
    status=0;
    ncard=0;
    while(1){
      ncard++;
      if(ncard>108) break;
      fits_read_record(fptr,ncard,cardline,&status);
      if(status){fits_die("Error reading file header",status);}
      if(strstr(cardline,"END")==cardline) break;
      if(strstr(cardline,"NAXIS")!=NULL) continue;
      if(strstr(cardline,"SIMPLE")!=NULL) continue;
      if(strstr(cardline,"BITPIX")!=NULL) continue;
      if(strstr(cardline,"BSCALE")!=NULL) continue;
      if(strstr(cardline,"BZERO")!=NULL) continue;
      if(strstr(cardline,"SLITLEN")!=NULL) continue;
      if(strstr(cardline,"CNTRLINE")!=NULL) continue;
      if(strstr(cardline,"SHUFFLE")!=NULL) continue;
      status=0;
      fits_write_record(outptr,cardline,&status);
      if(status) die("Error writing output file");}
    firstelem[2]=1;
    fits_write_pix(outptr,TFLOAT,firstelem,n_lampix,Tflux,&status);
    if(status) die("Error writing output file");
    if(errvl){
      firstelem[2]=2;
      fits_write_pix(outptr,TFLOAT,firstelem,n_lampix,Tsigma,&status);
      if(status) die("Error writing output file");}
    printf("Writing spectrum %d\r",slitnum);
    if(status){
      printf("Error writing spectrum file %d\n", status);
      return 1;}
    free(spectrum);
    free(espectrum);
    fits_movrel_hdu(fptr,1,undef,&status);}
   status=0;
   printf("\n");
   ffflus(outptr,&status);
   fits_close_file(outptr,&status);
   return 0;}
