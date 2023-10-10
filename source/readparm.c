/*****************************************************************************\
*                                                                             *
*   READPARM   reads JSON par files. readparm_s reads string values           *
*                                    readparm_i reads interger values         *
*				                             readparm_r reads real values             *
*                                    readparm_b reads logical values          *
*                                                                             *
*                                    openparm must be called first            *
*                                                                             *
*                                    Parameters: parm: name of parameter      *
*                                                svalue,rvalue,ivalue: value  *
*                                                                             *
*                                    Return: 0 if successful, 1 if error      *
*                                                                             *
*   VERSION: 11 Mar 2018                                                      *
*                                                                             *
\*****************************************************************************/

#include  <stdlib.h>
#include  <string.h>
#include  <stdio.h>

int           i;
char          choice[6];
static int    npar;
static char   key[100][30],value[100][30];

int OpenParm(char filename[]){

  int   len;
  char  line[133],s1[30],s2[30],choice[6];
  FILE  *parfile;

  //read parameter file
  parfile=fopen(filename,"r");
  if(parfile==NULL){
    printf("Cannot open parameter file\n",filename);
    return 1;}

  //parse File
  if(!fgets(line,133,parfile)) die("Parameter file error");
  if(!strcmp(line,"{")) die("Parameter file error");
  npar=0;
  while(1){
    if(!fgets(line,133,parfile)) break;
    if((!strcmp(line,"}")) || (!strcmp(line,"}\n"))) break;
    sscanf(line,"%s %s",s1,s2);
    memset(key[npar],'\0',30);
    memset(value[npar],'\0',30);
    strncpy(key[npar],s1+1,strlen(s1)-3);
    len=strlen(s2);
    if(strchr(s2,',') !=NULL){
      memset(s2+len-1,'\0',1);
      len--;}
    if(strchr(s2,'"') !=NULL){
       strncpy(value[npar],s2+1,len-2);}
    else{
       strcpy(value[npar],s2);}
    npar++;}
  fclose(parfile);
  return 0;}


  int ReadParm_s(char parm[], char ss[]){
    for(i=0;i<npar;i++){
      if(!strcasecmp(parm,key[i])){
        if(sscanf(value[i],"%s",ss)){
          return 0;}
        }
        }
    return 1;
  }

  int ReadParm_b(char parm[], int *ii){
   for(i=0;i<npar;i++){
      if(!strcasecmp(parm,key[i])){
        if(!sscanf(value[i],"%s",choice)) return 1;
        if(!strcmp(choice,"true")){
          *ii=1;
          return 0;}
        *ii=0;
        if(!strcmp(choice,"false")){
          return 0;}
        return 1;}
      }
      return 1;}

  int ReadParm_i(char parm[], int *ii[]){
    for(i=0;i<npar;i++){
      if(!strcasecmp(parm,key[i])){
        if(sscanf(value[i],"%d",ii)) return 0;}
      }
    return 1;
  }

  int ReadParm_r(char parm[], float *rr){
    for(i=0;i<npar;i++){
      if(!strcasecmp(parm,key[i])){
        if(sscanf(value[i],"%f",rr)) return 0;}
      }
    return 1;
  }
