 /*
 *  OpenCosParm- a variant of OpenParm that searches in $COSMOS_PAR_DIR and
  *              ~/Cospar directories for parameter file
 *
 *  Version 21 Sept 2016
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cosdat.h"
#include "cosmos.h"
#include <unistd.h>

int OpenCosParm(char prog[]){

  unsigned long    off;
  char   file[1000],*pardir;
  extern char **environ;

  //search for environment variable or truncated version
  pardir=*environ;
  while(1){
    if(strstr(pardir,"OSMOS_PAR_DIR") != NULL){
      pardir=index(pardir,'=')+1;
      break;}
    off=strcspn(pardir,"\0");
    if(off<2){
      printf("!\n! COSMOS_PAR_DIR is not defined\n!\n");
      return 1;}
    pardir=pardir+off+1;}



  //program parameters

  strcpy(file,pardir);
  if(*(file+strlen(file)-1) != '/') strcat(file,"/");
  strcat(file,prog);
  strcat(file,".json");
  if(OpenParm(file)==0) return 0;
  else return 1;

  return 0;}
