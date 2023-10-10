/*
  spec2big   converts extracted spectra from extension format to equispec format

  questions:
  1) How to deal with shuffled data?

*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "fitsio.h"
#include "cosmos.h"
#include "cpgplot.h"

// main
int main(int argc, char *argv[]){

   // variables
   fitsfile *fptr, *fptr_o[2];
   int   anynull, bitpix, cntr, dcflag, hdupos, hdutype, i, ii, j, k, n_slits,
         naxis, nkeys, nod, nnaxis, nullval, status, wcsdim;
   long  fpixel[3], inc[3], lpixel[3], naxes[3], nnaxes[2], nnaxes1, nnaxes2,
         Tnelem;
   char  apnum[133], card[FLEN_CARD], file[133], key[133], line[133],
         shuffled[133], tmp[133];
   float cd, cdelt1, crpix1, crval1, dslit, exptime, ltm;

   // pointers
   long  *ca, *cb, *cl, *nelem;
   char  **file_o;
   float **spec, **noise, **image_o;

   // initialize
   nnaxis = 2;
   nullval = 0;
   status = 0;
   wcsdim = 3;
   dcflag = 0;
   crpix1 = cd = ltm = 1.;

   // begin
   file_o = (char **) malloc((2) * sizeof(char*));
   for(i=0;i<2;i++) file_o[i] = (char *) malloc((133) * sizeof(char));

   if(argc <= 1 || argc > 2){
      printf("proper invocation: spec2big fitsfile\n");
      exit(1);
   }
      
   strcpy(file,argv[1]);
   strcpy(tmp,"!");
   strcat(tmp,argv[1]);
   strcat(file,"_2spec.fits");
   printf("Reading %s\n",file);
   strcpy(file_o[0],tmp);
   strcpy(file_o[1],tmp);
   strcat(file_o[0],"_big.fits");
   strcat(file_o[1],"_sig.fits");

   // open fits file
   fits_open_file(&fptr,file,READONLY,&status);
   if(status) fits_die(file,status);

   fits_read_key(fptr,TINT,"N_SLITS",&n_slits,card,&status);
   fits_read_key(fptr,TFLOAT,"EXPTIME",&exptime,card,&status);
   fits_read_key(fptr,TFLOAT,"CRVAL1",&crval1,card,&status);
   fits_read_key(fptr,TFLOAT,"CDELT1",&cdelt1,card,&status);
   fits_read_key(fptr,TSTRING,"SHUFFLED",&shuffled,card,&status);
   fits_read_key(fptr,TINT,"NOD",&nod,card,&status);
   fits_read_key(fptr,TFLOAT,"D_SLIT",&dslit,card,&status);
   if(status){
      fits_get_errstatus(status,line);
      printf("%s\n",line);
      status = 0;
   }

   ca = (long *) malloc((n_slits) * sizeof(long));
   cb = (long *) malloc((n_slits) * sizeof(long));
   cl = (long *) malloc((n_slits) * sizeof(long));
   nelem = (long *) malloc((n_slits) * sizeof(long));

   Tnelem = 0;
   nnaxes2 = 0;
   for(i=0;i<n_slits;i++){
      fits_get_img_param(fptr,3,&bitpix,&naxis,naxes,&status);
      fits_read_key(fptr,TINT,"CNTRLINE",&cntr,card,&status);
      if(status){
         fits_get_errstatus(status,line);
         printf("%s\n",line);
         exit(1);
      }
      nelem[i] = naxes[0]*naxes[1];
      nnaxes1 = naxes[0];

      /*
      fits_get_hdrspace(fptr, &nkeys, NULL, &status);
      for (ii = 1; ii <= nkeys; ii++) {
         if (fits_read_record(fptr, ii, card, &status))break;
         printf("%s\n", card);
      }
      */

      ca[i] = nnaxes2 + 1;
      cl[i] = nnaxes2 + cntr + 1;
      cb[i] = nnaxes2 + naxes[1];

      // sum y-axis
      nnaxes2 += naxes[1];

      // sum number of elements in each image
      Tnelem += nelem[i];
      fits_movrel_hdu(fptr, 1, NULL, &status);
   }

   nnaxes[0] = nnaxes1;
   nnaxes[1] = nnaxes2;
   for(i=0;i<2;i++) {
      status = 0;
      fits_create_file(&fptr_o[i],file_o[i],&status);
      fits_create_img(fptr_o[i],bitpix,nnaxis,nnaxes,&status);

      fits_write_key(fptr_o[i],TINT,"N_SLITS",&n_slits,"",&status);
      fits_write_key(fptr_o[i],TFLOAT,"D_SLIT",&dslit,"",&status);
      fits_write_key(fptr_o[i],TFLOAT,"CRVAL1",&crval1,"",&status);
      fits_write_key(fptr_o[i],TFLOAT,"CDELT1",&cdelt1,"",&status);
      fits_write_key(fptr_o[i],TFLOAT,"CRPIX1",&crpix1,"",&status);
      fits_write_key(fptr_o[i],TFLOAT,"EXPTIME",&exptime,"",&status);
      fits_write_key(fptr_o[i],TSTRING,"SHUFFLED",&shuffled,"",&status);
      fits_write_key(fptr_o[i],TINT,"NOD",&nod,"",&status);

      fits_write_key(fptr_o[i],TSTRING,"WAT0_001","system=equispec","",&status);
      fits_write_key(fptr_o[i],TSTRING,"WAT1_001",
                    "wtype=linear label=Wavelength units=Angstroms","",&status);
      fits_write_key(fptr_o[i],TSTRING,"WAT2_001","wtype=linear","",&status);
      fits_write_key(fptr_o[i],TSTRING,"WAT3_001","wtype=linear","",&status);
      fits_write_key(fptr_o[i],TINT,"WCSDIM",&wcsdim,"",&status);
      fits_write_key(fptr_o[i],TINT,"DC-FLAG",&dcflag,"",&status);
      fits_write_key(fptr_o[i],TSTRING,"CTYPE1","LINEAR","",&status);
      fits_write_key(fptr_o[i],TFLOAT,"LTM1_1",&ltm,"",&status);
      fits_write_key(fptr_o[i],TFLOAT,"CD1_1",&cdelt1,"",&status);
      fits_write_key(fptr_o[i],TSTRING,"CTYPE2","LINEAR","",&status);
      fits_write_key(fptr_o[i],TFLOAT,"LTM2_2",&ltm,"",&status);
      fits_write_key(fptr_o[i],TFLOAT,"CD2_2",&cd,"",&status);
      fits_write_key(fptr_o[i],TSTRING,"CTYPE3","LINEAR","",&status);
      fits_write_key(fptr_o[i],TFLOAT,"LTM3_3",&ltm,"",&status);
      fits_write_key(fptr_o[i],TFLOAT,"CD3_3",&cd,"",&status);
      
      if(status){
         fits_get_errstatus(status,line);
         printf("%s\n",line);
         status = 0;
      }
   }

   // initialize arrays
   spec = (float **) malloc((n_slits) * sizeof(float*));
   noise = (float **) malloc((n_slits) * sizeof(float*));

   fits_movabs_hdu(fptr, 1, &hdutype, &status);
   for(i=0;i<n_slits;i++){
      // read fits file

      //printf("i = %i\n",i);
      fits_get_hdu_num(fptr, &hdupos);
      fits_get_img_param(fptr,3,&bitpix,&naxis,naxes,&status);
     
      // read and write header keywords
      fits_get_hdrspace(fptr, &nkeys, NULL, &status);

      fits_read_key(fptr,TSTRING,"OBJECT",&line,card,&status);
      sprintf(key,"OBJ%03d",i+1);
      fits_write_key(fptr_o[0],TSTRING,key,&line,"",&status);
      fits_write_key(fptr_o[1],TSTRING,key,&line,"",&status);

      sprintf(key,"APNUM%d",i+1);
      sprintf(apnum,"%i %s %.2f %.2f",i+1,line,(float) ca[i],(float) cb[i]);
      fits_write_key(fptr_o[0],TSTRING,key,&apnum,"",&status);
      fits_write_key(fptr_o[1],TSTRING,key,&apnum,"",&status);

      sprintf(key,"CNTRL%03d",i+1);
      fits_write_key(fptr_o[0],TINT,key,&cl[i],"",&status);
      fits_write_key(fptr_o[1],TINT,key,&cl[i],"",&status);

      sprintf(key,"CSECT%dA",i+1);
      fits_write_key(fptr_o[0],TINT,key,&ca[i],"",&status);
      fits_write_key(fptr_o[1],TINT,key,&ca[i],"",&status);


      sprintf(key,"CSECT%dB",i+1);
      fits_write_key(fptr_o[0],TINT,key,&cb[i],"",&status);
      fits_write_key(fptr_o[1],TINT,key,&cb[i],"",&status);

      // read into individual arrays
      spec[i] = (float *) malloc(nelem[i]*sizeof(float));
      noise[i] = (float *) malloc(nelem[i]*sizeof(float));

      // setup arrays for reading data
      inc[0] = inc[1] = inc[2] = 1;
      fpixel[0] = fpixel[1] = fpixel[2] = 1;
      lpixel[0] = naxes[0];
      lpixel[1] = naxes[1]; 
      lpixel[2] = naxes[2]; 

      // read spectra
      fits_read_subset(fptr, TFLOAT, fpixel, lpixel, inc, &nullval, spec[i],
                       &anynull, &status);

      // read noise
      fpixel[2] = 2;
      fits_read_subset(fptr, TFLOAT, fpixel, lpixel, inc, &nullval, noise[i],
                       &anynull, &status);

      fits_movrel_hdu(fptr, 1, NULL, &status);

      /*if(status){
         fits_get_errstatus(status,line);
         printf("%s\n",line);
         exit(1);
      }*/

   }
   fits_close_file(fptr,&status);

   free(ca);
   free(cb);
   free(cl);

   // read into big array
   image_o = (float **) malloc((2) * sizeof(float*));
   for(i=0;i<2;i++) image_o[i] = (float *) malloc((Tnelem) * sizeof(float));

   k = 0;
   for(i=0;i<n_slits;i++){
      for(j=0;j<nelem[i];j++){
         image_o[0][k] = spec[i][j];
         image_o[1][k] = noise[i][j];
         k += 1; 
      }
   }
   free(nelem);

   // write fits files
   for(i=0;i<2;i++) {
      status = 0;
      fits_write_img(fptr_o[i],TFLOAT,1,Tnelem,image_o[i],&status);
      fits_close_file(fptr_o[i],&status);
      /*
      if(status){
         fits_get_errstatus(status,line);
         printf("%s\n",line);
         exit(1);
      }*/
      printf("Writing %s\n",file_o[i]);
   }

   free(file_o);
   free(spec);
   free(noise);
   free(image_o);
}
