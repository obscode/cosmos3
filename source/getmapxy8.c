/*****************************************************************************\
*
*   GetMapXY8   calculates the postions of mask holes and slits on each ccd
*               from spectral map for dispersed images, and combines on
*               8-chip coord system used by mosaic
*
*   VERSION  7 Dec 2017
*
*   USAGE:   GetMapXY8(map_file,nlines,wavelength list, plotpos,slit_type,X8,Y8)
                     return value: number of slit images
\*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "cosmos.h"


int  GetMapXY8(char* mpfile, int nline, float* sline,int ibin, int plotpos,
              int which, int* X8, int* Y8){


    int      i,j,xdisper,ord_disp,ord_sag,ord_tilt,ord_sagit,ord_slen,slitnum,
             nslit,nlv,chip,ends,ixcd, iycd,n,ocen,prtl,npxl,xcd,ycd,numchip,sign;
    float    lambda,lambda1,slit_len,slmin,slmax,lmin,lmax,coef_disp[10],coef_xdisp[10],
             coef_sag[10],coef_tilt[10],coef_sagit[10],coef_slen[10],slenth,slnt,
             slmx,slmn,xcen,
                ycen,xcen2,ycen2,xccd,x8cd,y8cd,tilt,sagit,xcen1,yccd,lambda0;
    char      file[133],line[133],dewar[10],name[40],instr[10];
    char      *lin;
    FILE    *mapfile;

    ocen=ends=0;
    if(plotpos>0) ends=1;
    else if(plotpos<0) ocen=1;

    strcpy(file,mpfile);
    if(strstr(file,".map")==NULL)strcat(file,".map");
    mapfile=fopen(file,"r");
    if(mapfile==NULL){
        printf("mapfile %s does not exist!\n",file);
        return 0;}

    //Read mapping data

  fgets(line,133,mapfile);
  if(!sscanf(line,"Xdispersion = %d",&xdisper)){
    printf("Mapfile error\n");
    return 0;}
  fgets(line,133,mapfile);
  if(!sscanf(line,"Fit orders = %d %d %d %d %d",&ord_disp,&ord_sag,&ord_tilt,
	     &ord_sagit,&ord_slen)){
    printf("Mapfile error\n");
    return 0;}
  if(ord_disp>9 || ord_sag>9 || ord_tilt>9 || ord_sagit>9 || ord_slen>9){
    printf("Warning: maximum fit order = 9\n");
    return 0;}
  fgets(line,133,mapfile);
  fgets(line,133,mapfile);
  if(!sscanf(line,"Lambda = %f %f",&lambda0,&lambda1)){
    printf("Mapfile error\n");
    return 0;}
  fgets(line,133,mapfile);
  if(!sscanf(line,"Dewar = %s %d",dewar,&numchip)) die("Mapfile error");
  if(Readcdf(dewar)) die("Error reading dewar definition file");
  fgets(line,133,mapfile);

  /*------------------------Find line ends-------------------------------*/

  nslit=-1;
  nlv=0;
  while(1){
    if(fgets(line,133,mapfile)==NULL)  die("Unexpected end to mapfile");
    if(sscanf(line,"SLIT %d %s",&slitnum, name)) break;}

  //New slit
  while(1){
    if(!strncmp(line,"END",3)) break;                           //end of data
    if(!sscanf(line,"SLIT %d %s",&slitnum, name)) die("Mapfile error");
    if((fgets(line,133,mapfile))==NULL) die("Unexpected end to mapfile");
    if(!strncmp(line,"END",3)) break;                           //end of data
    if(!sscanf(line,"LENGTH = %f",&slit_len)) continue;     //no data in slit

    //found a slit with data
    nslit++;
    //New Chip
    while(2){
      if((fgets(line,133,mapfile))==NULL) die("Unexpected end to mapfile");
      if(!strncmp(line,"END",3)) break;
      if(!sscanf(line,"CHIP %d %f %f %f %f %d %d",&chip,&slmin,&slmax,&lmin,&lmax,
          &npxl,&prtl)) break;
      //coef_disp data
      if((fgets(line,133,mapfile))==NULL)die("Unexpected end to mapfile");
      lin=&line[0];
      for(i=0;i<=ord_disp;i++){
	if((lin=strpbrk(lin,"1234567890.-"))==NULL)
                                               die("Unexpected end to mapfile");
	sscanf(lin,"%f",&coef_disp[i]);
	if((lin=strpbrk(lin," "))==NULL) die("Unexpected end to mapfile");}
      //coef_xdisp dat
      if((fgets(line,133,mapfile))==NULL)die("Unexpected end to mapfile");
      lin=&line[0];
      for(i=0;i<=ord_disp;i++){
	if((lin=strpbrk(lin,"1234567890.-"))==NULL)
                                               die("Unexpected end to mapfile");
	sscanf(lin,"%f",&coef_xdisp[i]);
	if((lin=strpbrk(lin," "))==NULL)die("Unexpected end to mapfile");}
      //coef_sag dat
      if((fgets(line,133,mapfile))==NULL)die("Unexpected end to mapfile");
      lin=&line[0];
      for(i=0;i<=ord_sag;i++){
	if((lin=strpbrk(lin,"1234567890.-"))==NULL)
                                               die("Unexpected end to mapfile");
	sscanf(lin,"%f",&coef_sag[i]);
	if((lin=strpbrk(lin," "))==NULL)die("Unexpected end to mapfile");}
      //coef_tilt data
      if((fgets(line,133,mapfile))==NULL)die("Unexpected end to mapfile");
      lin=&line[0];
      for(i=0;i<=ord_tilt;i++){
	if((lin=strpbrk(lin,"1234567890.-"))==NULL)
                                              die("Unexpected end to mapfile");
	sscanf(lin,"%f",&coef_tilt[i]);
	if((lin=strpbrk(lin," "))==NULL)die("Unexpected end to mapfile");}
      //coef_sagit data
      if((fgets(line,133,mapfile))==NULL)die("Unexpected end to mapfile");
      lin=&line[0];
      for(i=0;i<=ord_sagit;i++){
	if((lin=strpbrk(lin,"1234567890.-"))==NULL)
                                               die("Unexpected end to mapfile");
	sscanf(lin,"%f",&coef_sagit[i]);
	if((lin=strpbrk(lin," "))==NULL)die("Unexpected end to mapfile");}
      //coef_slen data
      if((fgets(line,133,mapfile))==NULL)die("Unexpected end to mapfile");
      lin=&line[0];
      for(i=0;i<=ord_slen;i++){
	if((lin=strpbrk(lin,"1234567890.-"))==NULL)
                                               die("Unexpected end to mapfile");
	sscanf(lin,"%f",&coef_slen[i]);
	if((lin=strpbrk(lin," "))==NULL)die("Unexpected end to mapfile");}
  if((prtl && which>0)  || (!prtl && which<0)) continue;

      //find ends of each line

      slenth=slmax-slmin;
//      sign = slmax>slmin ? 1 : -1;
      sign=1;
      for(n=0;n<nline;n++){
        lambda=sline[n];
        if(lambda<lmin || lambda>lmax) continue;
        xcen=polyvalue(lambda,coef_disp,ord_disp);;
        ycen=polyvalue(lambda,coef_sag,ord_sag);
        tilt=polyvalue(lambda,coef_tilt,ord_tilt);
        if(slit_len>.1 && !ocen){
            slnt=polyvalue(lambda,coef_slen,ord_slen);
            slmx=slmax*slnt/slenth;
            slmn=slmin*slnt/slenth;
            sagit=polyvalue(lambda,coef_sagit,ord_sagit);
            xcen1=xcen+tilt*slmn/(slmx-slmn);
            xcen2=xcen1+tilt;
            if(!ends){
                ycen=ycen+(slmx+slmn)/2.;
                xcen=0.5*(xcen1+xcen2);}
            else{
                xcen=xcen1-sagit;
                xcen2=xcen+sign*tilt;
                ycen=ycen+slmn;
                ycen2=ycen+(slmx-slmn);}
            }
        for(i=0;i<=ends;i++){
            xccd=xcen;
            if(!xdisper){
            xccd=ycen;
            ycen=xcen;}
            mspos(chip,xccd,ycen,&x8cd,&y8cd);
            xcd=xccd/ibin+1.;
            ycd=(ycen)/ibin+1.;
            ixcd=floor((x8cd)/ibin+1.0);
            iycd=floor((y8cd)/ibin+1.0);
            X8[nlv]= numchip==8 ? ixcd:xcd;
            Y8[nlv]= numchip==8 ? iycd:ycd;
            nlv++;
            ycen=ycen2;
            xcen=xcen2;}
      }
    }
}
return nlv;}
