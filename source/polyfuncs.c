/******************************************************************************\
*                                                                              *
*  POLYFUNC  a collection of polynomial solution and evaluation routines       *
*                                                                              *
*  VERSION  28 Feb 2005                                                        *
*                                                                              *
\******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*
 *  POLYVALUE evaluates a polynomial
 *
*/


float polyvalue(float x,float coef[],int ord){

  float value;
  int i;
  double t,val;

  t=(double) x;
  val=(double) coef[0];
  for(i=1;i<=ord;i++){
    val+=t*coef[i];
    t*=x;}
  value=val;
  return value;}

void polyval(float x,float coef[],int ord,float *value){

  int i;
  double t,val;

  t=x;
  val=coef[0];
  for(i=1;i<=ord;i++){
    val+=t*coef[i];
    t*=x;}
  *value=(float)val;
  return;}


/*****************************************************************************\
*                                                                             *
*  POLYFIT fits a polynomial to data                                          *
*          PLYFIT is the version with externally-defined working arrays       *
*                                                                             *
*  USAGE: plyfit(xvalues,yvalues,nvalues,nterms,xlist,ylist,elist,coef)       *
*         buffers xlist,ylist, elist must be defined                          *
*              xlist=malloc(sizeof(double)*2*(order+1))                       *
*              ylist=malloc(sizeof(double*(order+1))                          *
*              elist=malloc(sizeob(double_pointer*(order+1))                  *
*             for(i=0;i<order+1;i++)elist[i]=malloc(sizeof(double)*(order+2)  *
*                                                                             *
*  VERSION 5 Juneq 2003                                                       *
*                                                                             *
\*****************************************************************************/

void plyfit(float *xvals,float *yvals,int nvals ,int order, float *coef,
             double *xlist, double *ylist,double **elist){

  double x,fact,sum;
  int   i,j,k,nterms;

  nterms=order+1;

  for(i=0;i<2*nterms;i++) *(xlist+i)=0;
  for(i=0;i<nterms;i++) *(ylist+i)=0;
  for(j=0;j<nvals;j++){
    x=1;
    for(i=0;i<2*nterms;i++){
      *(xlist+i)+=x;
      x*= *(xvals+j);}
    x=1;
    for(i=0;i<nterms;i++){
      *(ylist+i)+=*(yvals+j)*x;
      x*= *(xvals+j);}
  }

  //equation terms
  for(i=0;i<nterms;i++){
    for(j=0;j<nterms;j++){
      *(*(elist+i)+j)= *(xlist+i+j);}
    *(*(elist+i)+nterms)= *(ylist+i);}


  //solve simultaneous equations

  //reduce equations
  for(i=1;i<nterms;i++){
    for(j=i;j<nterms;j++){
      fact= *(*(elist+j)+i-1)/ *(*(elist+i-1)+i-1);
      for(k=i;k<nterms+1;k++){
	*(*(elist+j)+k)-= fact*(*(*(elist+i-1)+k));}
      *(*(elist+j)+i-1)=0;}
    }

  //solve for parameters
  for(i=nterms-1;i>=0;i--){
    sum=*(*(elist+i)+nterms);
    for(j=i+1;j<nterms;j++){
      sum-=*(coef+j)*(*(*(elist+i)+j));}
    *(coef+i)=sum/(*(*(elist+i)+i));}

  return;}



/*****************************************************************************\
*                                                                             *
*  POLYFIT fits a polynomial to data                                          *
*                                                                             *
*  USAGE: polyfit(xvalues,yvalues,nvalues,nterms,coef)                        *
*                                                                             *
*  VERSION 5 Juneq 2003                                                        *
*                                                                             *
\*****************************************************************************/

void polyfit(float *xvals,float *yvals,int nvals ,int order, double *coef){

  double **elist,*xlist,*ylist;
  double x,fact,sum;
  int   i,j,k,nterms;

  nterms=order+1;

  //xy factors
  xlist=malloc(sizeof(x)*2*nterms);
  ylist=malloc(sizeof(x)*nterms);
  for(i=0;i<2*nterms;i++) *(xlist+i)=0;
  for(i=0;i<nterms;i++) *(ylist+i)=0;

  for(j=0;j<nvals;j++){
    x=1;
    for(i=0;i<2*nterms;i++){
      *(xlist+i)+=x;
      x*= *(xvals+j);}
    x=1;
    for(i=0;i<nterms;i++){
      *(ylist+i)+=*(yvals+j)*x;
      x*= *(xvals+j);}
    }

  //equation terms
  elist=malloc(sizeof(xlist)*nterms);
  for(i=0;i<nterms;i++){
    *(elist+i)=malloc(sizeof(x)*(nterms+1));
   for(j=0;j<nterms;j++){
     *(*(elist+i)+j)= *(xlist+i+j);}
   *(*(elist+i)+nterms)= *(ylist+i);}


  //solve simultaneous equations

  //reduce equations
  for(i=1;i<nterms;i++){
    for(j=i;j<nterms;j++){
      fact= *(*(elist+j)+i-1)/ *(*(elist+i-1)+i-1);
      for(k=i;k<nterms+1;k++){
	*(*(elist+j)+k)-= fact*(*(*(elist+i-1)+k));}
      *(*(elist+j)+i-1)=0;}
    }

  //solve for parameters
  for(i=nterms-1;i>=0;i--){
    sum=*(*(elist+i)+nterms);
    for(j=i+1;j<nterms;j++){
      sum-=*(coef+j)*(*(*(elist+i)+j));}
    *(coef+i)=sum/(*(*(elist+i)+i));}

  //free memory
  for(i=0;i<nterms;i++) free(*(elist+i));
  free(elist);
  free(xlist);
  free(ylist);

  return;}

/*****************************************************************************\
*                                                                             *
*  POLYFIT_W fits a polynomial to data with weights                           *
*                                                                             *
*  USAGE: polyfit(xvalues,yvalues,weights,nvalues,nterms,coef)                *
*                                                                             *
*  VERSION 28 Aug 2006                                                        *
*                                                                             *
\*****************************************************************************/

void polyfit_w(float *xvals,float *yvals,float *wts, int nvals ,int order,
             double *coef){

  float  weight;
  double **elist,*xlist,*ylist;
  double x,fact,sum;
  int   i,j,k,nterms;

  nterms=order+1;

  //xy factors
  xlist=malloc(sizeof(x)*2*nterms);
  ylist=malloc(sizeof(x)*nterms);
  for(i=0;i<2*nterms;i++) *(xlist+i)=0;
  for(i=0;i<nterms;i++) *(ylist+i)=0;

  for(j=0;j<nvals;j++){
    weight=*(wts+i);
    x=1;
    for(i=0;i<2*nterms;i++){
      *(xlist+i)+=x*weight;
      x*= *(xvals+j);}
    x=1;
    for(i=0;i<nterms;i++){
      *(ylist+i)+=*(yvals+j)*x*weight;
      x*= *(xvals+j);}
  }

  //equation terms
  elist=malloc(sizeof(xlist)*nterms);
  for(i=0;i<nterms;i++){
    *(elist+i)=malloc(sizeof(x)*(nterms+1));
    for(j=0;j<nterms;j++){
      *(*(elist+i)+j)= *(xlist+i+j);}
    *(*(elist+i)+nterms)= *(ylist+i);}


  //solve simultaneous equations

  //reduce equations
  for(i=1;i<nterms;i++){
    for(j=i;j<nterms;j++){
      fact= *(*(elist+j)+i-1)/ *(*(elist+i-1)+i-1);
      for(k=i;k<nterms+1;k++){
        *(*(elist+j)+k)-= fact*(*(*(elist+i-1)+k));}
      *(*(elist+j)+i-1)=0;}
  }

  //solve for parameters
  for(i=nterms-1;i>=0;i--){
    sum=*(*(elist+i)+nterms);
    for(j=i+1;j<nterms;j++){
      sum-=*(coef+j)*(*(*(elist+i)+j));}
    *(coef+i)=sum/(*(*(elist+i)+i));}

  //free memory
  for(i=0;i<nterms;i++) free(*(elist+i));
  free(elist);
  free(xlist);
  free(ylist);

  return;}



/*****************************************************************************\
*                                                                             *
*  POLYFIT_w fits a polynomial to data with weighting                         *
*  PLYFIT_W is the version with externally-defined working arrays             *
*                                                                             *
*  USAGE: plyfit_w(xvalues,yvalues,weigths,nvalues,nterms,coef,xlist,ylist,   *
*                  elist)                                                     *
*         buffers xlist,ylist, elist must be defined                          *
*              xlist=malloc(sizeof(double)*2*(order+1))                       *
*              ylist=malloc(sizeof(double*(order+1))                          *
*              elist=malloc(sizeob(double_pointer*(order+1))                  *
*             for(i=0;i<order+1;i++)elist[i]=malloc(sizeof(double)*(order+2)  *
*                                                                             *
*  VERSION 27 Aug 2006                                                        *
*                                                                             *
\*****************************************************************************/

void plyfit_w(float *xvals,float *yvals,float *wts,int nvals ,int order,
              float *coef,double *xlist, double *ylist,double **elist){

  float weight;
  double x,fact,sum;
  int   i,j,k,nterms;

  nterms=order+1;

  for(i=0;i<2*nterms;i++) *(xlist+i)=0;
  for(i=0;i<nterms;i++) *(ylist+i)=0;
  for(j=0;j<nvals;j++){
    weight=*(wts+j);
    x=1;
    for(i=0;i<2*nterms;i++){
      *(xlist+i)+=x*weight;
      x*= *(xvals+j);}
    x=1;
    for(i=0;i<nterms;i++){
      *(ylist+i)+=*(yvals+j)*x*weight;
      x*= *(xvals+j);}
  }

  //equation terms
  for(i=0;i<nterms;i++){
    for(j=0;j<nterms;j++){
      *(*(elist+i)+j)= *(xlist+i+j);}
    *(*(elist+i)+nterms)= *(ylist+i);}


  //solve simultaneous equations

  //reduce equations
  for(i=1;i<nterms;i++){
    for(j=i;j<nterms;j++){
      fact= *(*(elist+j)+i-1)/ *(*(elist+i-1)+i-1);
      for(k=i;k<nterms+1;k++){
        *(*(elist+j)+k)-= fact*(*(*(elist+i-1)+k));}
      *(*(elist+j)+i-1)=0;}
  }

  //solve for parameters
  for(i=nterms-1;i>=0;i--){
    sum=*(*(elist+i)+nterms);
    for(j=i+1;j<nterms;j++){
      sum-=*(coef+j)*(*(*(elist+i)+j));}
    *(coef+i)=sum/(*(*(elist+i)+i));}

  return;}
