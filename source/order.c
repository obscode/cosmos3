/*
*   ORDER sorts a group of objects by the value of am. nob is the index
*   number of the object; irnk = ranked list of index numbers; rank =
*   ranked list of values of am. rank and irnk must be dimensioned in
*   calling program >= number of objects.
*   call ORDER once for each object.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void order(int *nob, float *am, int irnk[], float rank[]){

  int 	ifust,ilst,inow,j;

  switch(*nob){
  case 0:
    *irnk=*nob;
    *rank=*am;
    return;
  case 1:
    if(*am>=*rank){
      *(irnk+1)=*nob;
      *(rank+1)=*am;}
    else{
      *(irnk+1)=*irnk;
      *(rank+1)=*rank;
      *irnk=*nob;
      *rank=*am;}
    return;
  default:
    ifust=0;
    ilst=*nob-1;
    for(;;){
      inow=(ifust+ilst)/2;
      if((inow==ilst)||
	 (*(rank+inow)==*am)||
	 (*(rank+inow)<*am)&&((inow==*nob-1)||(*(rank+inow+1)>=*am))||
	 (*(rank+inow)>*am)&&((inow==0)||(*(rank+inow-1)<=*am))){
	if(*(rank+inow)<=*am) inow++;
	for(j=*nob-1;j>=inow;j--){
	  *(irnk+j+1)=*(irnk+j);
	  *(rank+j+1)=*(rank+j);}
	*(irnk+inow)=*nob;
	*(rank+inow)=*am;
	return;}
      else if(*(rank+inow)<*am) ifust=inow+1;
      else ilst=inow-1;}
    }
  }

void dorder(int *nob, double *am, int irnk[], double rank[]){

  int 	ifust,ilst,inow,j;

  switch(*nob){
  case 0:
    *irnk=*nob;
    *rank=*am;
    return;
  case 1:
    if(*am>=*rank){
      *(irnk+1)=*nob;
      *(rank+1)=*am;}
    else{
      *(irnk+1)=*irnk;
      *(rank+1)=*rank;
      *irnk=*nob;
      *rank=*am;}
    return;
  default:
    ifust=0;
    ilst=*nob-1;
    for(;;){
      inow=(ifust+ilst)/2;
      if((inow==ilst)||
	 (*(rank+inow)==*am)||
	 (*(rank+inow)<*am)&&((inow==*nob-1)||(*(rank+inow+1)>=*am))||
	 (*(rank+inow)>*am)&&((inow==0)||(*(rank+inow-1)<=*am))){
	if(*(rank+inow)<=*am) inow++;
	for(j=*nob-1;j>=inow;j--){
	  *(irnk+j+1)=*(irnk+j);
	  *(rank+j+1)=*(rank+j);}
	*(irnk+inow)=*nob;
	*(rank+inow)=*am;
	return;}
      else if(*(rank+inow)<*am) ifust=inow+1;
      else ilst=inow-1;}
    }
  }

