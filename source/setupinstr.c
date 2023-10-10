#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ioutils.h"
#include "kdcutil.h"
#include "mgutils.h"
#include "optutils.h"
#include "cosdat.h"


/*
 *  SetupInstr   inputs observational setup, and defines instrument
 *
 *  VERSION      4 Dec 03
 *
 */

int SetupInstr(obsdef *obsdata, element **instrument){

  extern char OPTICS_FILE[];
  char CHAR,flnm[80],*COS_HOME,INSTRMNT[80],GRATNG[80];
  int  imacs,lng,spec;
  extern element *Cur_grism,*Cur_grating;
  double tmp;
  stgrism *sg; //IMACS
  element *optics;

  COS_HOME=malloc(sizeof(CHAR)*80);
  COS_HOME=getenv("COSMOS_HOME");
  if(COS_HOME==NULL) die("COSMOS_HOME undefined!");

  strcpy(flnm,COS_HOME);
  strcat(flnm,"/sdata/");
  strcat(flnm,OPTICS_FILE);
  optics=get_optic_data(flnm,0);
//  optics=get_data(flnm,*optstore,WHITESPACE,"!#","!","#;","","");
  if(optics == NULL) die("Failure reading optics data file");

  imacs=lng=spec=0;
  GRangle=NULL;
  GRorder=NULL;
  strcpy(INSTRMNT,obsdata->instrument);
  if(!strncmp(obsdata->mode,"SPEC",4)) spec=1;
  if(!strcmp(INSTRMNT,"IMACS")){
    imacs=1;
    if(!strcmp(obsdata->camera,"LONG")){
      lng=1;
      GRangle = xalloc(double);
      GRorder = xalloc(int);
      *GRangle=(22.5+obsdata->gr_angle)*Degree;
      *GRorder=obsdata->gr_order;
      Cur_grating=find_element(optics,obsdata->grating);
      if(Cur_grating==NULL) die("Cannot find grating definition!");
      printf("*\n");
      if(Cur_grating->head && Cur_grating->head->data) {
          memcpy (&tmp, &(Cur_grating->head->data[4]), 8);
          tmp = tmp + Degree*obsdata->alignrot;
          //printf("l align %f\n",tmp);
          memcpy (&(Cur_grating->head->data[4]), &tmp, 8);}
       }
    }
  if(!lng){
    Cur_grism=find_element(optics,obsdata->grating);
    if(Cur_grism==NULL) die("Cannot find grism definition!");
    if(Cur_grism->head && Cur_grism->head->data) {
          memcpy (&tmp, &(Cur_grism->head->data[4]), 8);
          tmp = tmp + Degree*obsdata->alignrot;
          //printf("s align %f\n",tmp);
          memcpy (&(Cur_grism->head->data[4]), &tmp, 8);
       }
    if(spec){
      sg=Cur_grism->data;
      GRorder = xalloc(int);
      *GRorder=obsdata->gr_order;
      if(sg->angle < 0) *GRorder = -*GRorder;
  }
    }

  //camera name

  if(imacs){
    if(lng){
      strcat(INSTRMNT,"_l");}
    else{
      strcat(INSTRMNT,"_s");}
    if(spec){
      strcat(INSTRMNT,"c");}
    else{
      strcat(INSTRMNT,"d");}
    }

  if(!strcmp(INSTRMNT,"LDSS2")){
    if(!spec){
      strcat(INSTRMNT,"_d");}
    }

   if(!strcmp(INSTRMNT,"LDSS3")){
    if(!spec){
      strcat(INSTRMNT,"_D");}
   }

   *instrument=find_element(optics,INSTRMNT);
   return 0;}

/* ReadSMFfile reads in SMF file after defining optics
/
*/

int ReadSMFfile(char *SMF_FILE, Obs **obset){

  element *optics;
  extern char OPTICS_FILE[];
  char *COS_HOME,flnm[80],CHAR;

  COS_HOME=malloc(sizeof(CHAR)*80);
  COS_HOME=getenv("COSMOS_HOME");
  if(COS_HOME==NULL) die("COSMOS_HOME undefined!");
  strcpy(flnm,COS_HOME);
  strcat(flnm,"/sdata/");
  strcat(flnm,OPTICS_FILE);
  optics=get_optic_data(flnm,0);
  if(optics == NULL){
    printf("Failure reading optics data file");
    return 1;}
    *obset=read_smdf(SMF_FILE,optics,0);
  if(*obset == NULL){
    printf ("Error reading SMF file %s\n",SMF_FILE);
    return 1;}

  return 0;}
