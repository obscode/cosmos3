/*****************************************************************************\
*                                                                             *
*  SPECTRAL-MAP uses a map file to predict the location of spectral lines     *
*               output in format to be used by xxxx                           *
*                                                                             *
*  VERSION  18 Mar 2005                                                       *
*                                                                             *
\*****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fitsio.h"
#include "cosmos.h"
#include "cpgplot.h"

int main(int argc,char *argv[]){

  char  CHAR,file[80],flnm[80],ifile[80],dfile[80],*DATA_DIR,line[133],*lin,
        name[40],mfile[80],oldline[133],lnfile[80],*COS_HOME,*HOME,
        dewar[10],objname[40],alphabet[53];  
  int   INT,status,bitpix,naxis,nelem,i,j,here,xdisper,ord_disp,ord_sag,ocen,
        ix8,ord_tilt,ord_sagit,ord_slen,nslit,chip,slitnum,ixpos,iypos,isloff,
        noff,lvmin,lvmax,search_height,search_width,iswmn,iswmx,nsw,nline,n,
        iy8,islmx,islmn,nsl,il,ii,jj,kk,ixtop,iytop,slit_width,nlv,iyp,ndline,
        ixcd,iycd,ibin,magfac,ysign,narg,maxord,ord_ddisp,ord_dsag,minlines,
        ends,ybin;
  float FLOAT,telscale,lambda0,lambda1,max_slit,slen,slit_len,slmin,slmax,
        xa[2],ya[2],x8,y8,lmin,lmax,lambda,xcen,ycen,slnt,ypos,xpos,dix,diy,
        x8cd,y8cd,xccd,yccd,zz,sumax,sum,coef_disp[10],ynew,coef_sag[10],
        coef_tilt[10],coef_sagit[10],coef_slen[10],ymid,slnth,slenth,tilt,
        slmx,slmn,aa,bb,bmx,bmn,iyd,d_coef[10],t,s_coef[10],coef_xdisp[10],
        x_coef[10],xnew,ycen2,xcen2;
  long  naxes[3],firstel,stack_len;
  int   *INTP,*image[9],*im,*anynl,**array[9],*isx,*isy;
  float *FLOATP,*stack,**spectrum,*wt,**weight,*sarray,**sbox,*ssx,*ssy;
  int   nulvl=0;
  double *xlist,*ylist,**elist,DOUBLE,*DOUBLEP;
  fitsfile *fptr[9];
  FILE *mapfile,*linefile,*outfile;
  fitsdef fitsinfo;
  int   xsize=2048;
  int   msize=8192;


  // size of linelist
  int nlist = 5000;
  float xshift[nlist],yshift[nlist],sline[nlist],lval[nlist],xval[nlist],
        lamshift[nlist];

  nulvl=0;
  firstel=1;
  nline=ends=ocen=0;
  strcpy(alphabet,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");

   //get input files

  if(argc<7){
    while(1){
      printf("Enter map file name: ");
      scanf("%s",mfile);
      strcpy(file,mfile);
      if(strstr(file,".map")==NULL)strcat(file,".map");
      mapfile=fopen(file,"r");
      if(mapfile!=NULL) break;
      printf("mapfile %s does not exist!\n",file);
      }
    while(1){
      printf("Enter line file name:   ");
      scanf("%s",dfile);
      linefile=fopen(dfile,"r");
      if(linefile!=NULL) break;
      printf("File % cannot be found/n");}
    printf("Enter binning: ");
    scanf("%d",&ibin);}
  else{
    if(!(strcmp(argv[1],"-m"))){
      narg=2;}
    else{
      if(!(strcmp(argv[3],"-m"))){
	narg=4;}
      else{
	if(!(strcmp(argv[5],"-m"))){
	  narg=6;}
	else{
	  printf("proper invocation: spectral-map -l linefile -m mapfile\n");
	  return 1;}
        }
      }
    strcpy(mfile,argv[narg]);
    strcpy(file,mfile);
    if(strstr(file,".map")==NULL)strcat(file,".map");
    mapfile=fopen(file,"r");
    if(mapfile==NULL){
      printf("map file %s does not exist!\n",file);
      return 1;}
    if(!(strcmp(argv[1],"-l"))){
      narg=2;}
    else{
      if(!(strcmp(argv[3],"-l"))){
	narg=4;}
      else{
	if(!(strcmp(argv[5],"-l"))){
	  narg=6;}
	else{
	  printf("proper invocation: spectral-map -l linefile -m mapfile -b binning\n");
	  return 1;}
        }
      }

    //is linelist a file or wavelength?
    if(strpbrk(argv[narg],alphabet)==NULL){
      nline=1;
      sscanf(argv[narg],"%f",&sline[0]);}
    else{
      sscanf(argv[narg],"%s",dfile);
      strcpy(dfile,argv[narg]);
      linefile=fopen(dfile,"r");
      if(linefile==NULL) die("cannot find line file");}
     if(!(strcmp(argv[1],"-b"))){
      narg=2;}
    else{
      if(!(strcmp(argv[3],"-b"))){
	narg=4;}
      else{
	if(!(strcmp(argv[5],"-b"))){
	  narg=6;}
	else{
	  printf("proper invocation: spectral-map  -m mapfile -b binning\n");
	  return 1;}
        }
      }
     sscanf(argv[narg],"%d",&ibin);
    if(argc==8){
      if(!strcmp(argv[7],"-c")) ocen=1;
      if(!strcmp(argv[7],"-e")) ends=1;
    }
  }
  
  strcpy(file,mfile);
  strcat(file,".xy");
  outfile=fopen(file,"w");
  if(outfile==NULL) die("cannot open output file");

  //binning
  if(ibin<10){
    ybin=ibin;}
  else{
    i=(int) ibin/10;
    ybin=ibin-10*i;
    ibin=i;}

  //Read mapping data

  fgets(line,133,mapfile);
  if(!sscanf(line,"Xdispersion = %d",&xdisper)){
    printf("Mapfile error\n");
    return 1;}
  fgets(line,133,mapfile);
  if(!sscanf(line,"Fit orders = %d %d %d %d %d",&ord_disp,&ord_sag,&ord_tilt,
	     &ord_sagit,&ord_slen)){
    printf("Mapfile error\n");
    return 1;}
  if(ord_disp>9 || ord_sag>9 || ord_tilt>9 || ord_sagit>9 || ord_slen>9){
    printf("Warning: maximum fit order = 9\n");
    return 1;}
  fgets(line,133,mapfile);
  fgets(line,133,mapfile);
  if(!sscanf(line,"Lambda = %f %f",&lambda0,&lambda1)){
    printf("Mapfile error\n");
    return 1;}
  fgets(line,133,mapfile);
  if(!sscanf(line,"Dewar = %s",dewar)) die("Mapfile error");
  if(Readcdf(dewar)) die("Error reading dewar definition file");
  fgets(line,133,mapfile);

  //Read line list

  if(!nline){
  while((fgets(line,133,linefile))!=NULL){
    sscanf(line,"%f",&sline[nline]);
    nline++;}
  }

   /*------------------------Find line ends-------------------------------*/

  nslit=-1;
  while(1){
    if(fgets(line,133,mapfile)==NULL){
      printf("Unexpected end of mapfile at slit %d\n",nslit);
      return 1;}
    if(sscanf(line,"SLIT %d %s",&slitnum, name)) break;}

  //New slit
  while(1){
    if(!strncmp(line,"END",3)) break;                           //end of data
    if(!sscanf(line,"SLIT %d %s",&slitnum, name)) die("Mapfile error");
    if(fgets(line,133,mapfile)==NULL){
      printf("Unexpected end of mapfile at slit %d\n",nslit);
      return 1;}
    if(!strncmp(line,"END",3)) break;                           //end of data
    if(!sscanf(line,"LENGTH = %f",&slit_len)) continue;     //no data in slit

    //found a slit with data
    nslit++;
    nlv=0;
    //New Chip
    while(2){
    if(fgets(line,133,mapfile)==NULL){
      printf("Unexpected end of mapfile at slit %d\n",nslit);
      return 1;}
      if(!strncmp(line,"END",3)) break;                                
      if(!sscanf(line,"CHIP %d %f %f %f %f",&chip,&slmin,&slmax,&lmin,&lmax))
                                                                        break;
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
    
   
      //find ends of each line

      slenth=slmax-slmin;
      for(n=0;n<nline;n++){

	//extraction box centered on line
	
	lambda=sline[n];
	if(lambda<lmin || lambda>lmax) continue;
	xcen=polyvalue(lambda,coef_disp,ord_disp);;	 
	ycen=polyvalue(lambda,coef_sag,ord_sag);
  tilt=polyvalue(lambda,coef_tilt,ord_tilt);
	if(slit_len>.1 && !ocen){
	  slnt=polyvalue(lambda,coef_slen,ord_slen);
	  slmx=slmax*slnt/slenth;
	  slmn=slmin*slnt/slenth;
	  if(!ends){
      ycen=ycen+(slmx+slmn)/2.;}
    else{
      xcen=xcen+tilt*slmn/(slmx-slmn);
      xcen2=xcen+tilt;
      ycen=ycen+slmn;
      ycen2=ycen+(slmx-slmn);}
  }
  for(i=0;i<=ends;i++){
    xccd=xcen;
    if(!xdisper){
      xccd=ycen;
      ycen=xcen;}
    mspos(chip,xccd,ycen,&x8cd,&y8cd);
    xccd=xccd/ibin+1.;
    yccd=(ycen)/ybin+1.;
    ixcd=floor((x8cd)/ibin+1.5);
    iycd=floor((y8cd)/ybin+1.5);
    fprintf(outfile,"%4d %4d %5d %s %d %6.1f %6.1f %6.1f %s %f\n",
                ixcd,iycd,slitnum,name,chip,xccd,yccd,lambda,objname,tilt);
    ycen=ycen2;
    xcen=xcen2;}
      }
    }
    }
  fclose(outfile);
  return 0;}


	
	
	
