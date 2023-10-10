/*****************************************************************************\
*                                                                             *
*  SUBSKY2      subrtracts sky from 2-dim spectra using Kelson algorith and   *
*               information in map files an calculates errors                 *
*               2d bspline version                                            *
*                                                                             *
*  VERSION  08 Sept 2005                                                      *
*                                                                             *
\*****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fitsio.h"
#include "cosmos.h"
//#include "bspline.h"
#include "cpgplot.h"
//#include "debugs.h"

int comp_nums(const long *, const long *);
double *lamtemp;

int main(int argc,char *argv[]){

  char     file[80],ofile[80],flnm[1000],ifile[80],dfile[80],*datadir,line[133],*lin,c[7],
           cardline[80],objname[40],bdfile[80],answer[80],dewar[10],CHAR,ccur,
           empty,*homedir,*objsub,*parfile;
  int      INT,status,nelem,i,j,here,xdisper,ord_disp,ord_sag,ord_tilt,nchip,
           ord_sagit,ord_slen,nslit,n_slit,chip,diag,breakpoint,iii,ni,slitnum,
           dsign,sagsn,exclude,ssign,is,ix,ixmx,iymx,ixx,iyy,nint,iytop,iybot,
           pnt,i1,i2,any,ii,size,ibin,knot0,knotn,n_knots,lwork,nest,mbx2,lpix,
           fpix,medthrsh,medbox,mbx1,nelem2,ncard,two,three,nulvl,narg,prob,
           iopt,curnt,nmed,badpix[9],x0,x1,y0,y1,l,k,nbpx,splorder,ier=0,uu,vv,
           bx,by,morder,b1,b2,kwork,nimx,twodim,coefsz,s_splorder,nargz,nargf,
           n_knots_s,slorder,wlorder,nestm,nest_s,lwork1,lwork2,edge,bpx,nchp,
           nshuffle,inshuffle,ibiny,xbin,ybin,bbin,bxbin,bybin,sub_sky,ixend,
           year,month,day,edge2,read_old,onelem,tmp;
  int      *INTP,anynl,*iwork,*ixval,*iyval,*irnk,*badx0[9],*badx1[9],
           *stack,*bady0[9],*bady1[9];
  long     *index,*asort;
  float    FLOAT,min_lambda,max_lambda,lambda0,lambda1,max_slit,slit_len,slmin,
           slmax,lmin,lmax,lambda,xcen,ycen,slnt,sagt,rag,wid,fmin,fmax,tilt,
           medval,stdev,thrsh,siglimit,noise,gain,xxx,scen,objshift,
           diag1,tthrsh,bthrsh,coef_xdisp[10],coef_disp[10],coef_sag[10],curv,
           diag0,coef_tilt[10],coef_sagit[10],xcur,ycur,ypixmid,sagmid,
           delknot,ki,di,coef_slen[10],s_delknot,dog,dig,dag,yyy,fdate,tmp2;
  float    *FLOATP,*image[9],*im,**array[9],*error[9],*err,**earray[9],*flptr,
           *oimage[9],*oim,**oarray[9];
  double   flux,*knot,*goodpix,*work,*lamval,*rnk,*slval,*coefs,aa,bb,cc,smfac,
           *sltemp,*knot_s,*work1,*work2,zeroval,*goodtemp,*flxval,*spl_vals;
  float    *dlamval,*dsplpix;
  long     naxes[3],firstel,firstelem[3];
  fitsfile *fptr[9], *fptr_o[9], *oldfptr[9];
  fitsdef  fitsinfo,ofitsinfo;
  FILE     *badfile,*mapfile;
  int      tp;
  double   v1,lval1,lvaln,slval1,slvaln,loff,resid;
  int      resort,nsort,soff,osort;
  int      one=1;

  work1=work2=coefs=NULL;
  iwork=NULL;
  nbpx=100;                      //maximum number of bad pixel groups per chip
  two=2;
  three=3;
  breakpoint=prob=0;
  nargz=0;
  zeroval=0.000001;                           //value for weights set to zero
  //spline fitting parameters
  smfac=0.;
  nulvl=0;
  bpx=diag=0;
  read_old=0;
  firstelem[0]=firstelem[1]=firstelem[2]=1;
  firstel=1;
  objsub=malloc(sizeof(CHAR)*80);
  strcpy(objsub,"");
  parfile=malloc(sizeof(CHAR)*80);
  strcpy(parfile,"");
  datadir=malloc(sizeof(CHAR)*80);
  datadir=getenv("COSMOS_IMAGE_DIR");
  if(datadir==NULL){
    printf("COSMOS_IMAGE_DIR undefined!\n");
    return 1;}
  strcat(datadir,"/");
  homedir=malloc(sizeof(CHAR)*80);
  homedir=getenv("COSMOS_HOME");
  if(homedir==NULL){
    printf("COSMOS_HOME undefined!\n");
    return 1;}
  strcat(homedir,"/");
  bbin=bxbin=bybin=1;

  //get input files

  if(argc<5){
    while(1){
      printf("Enter map file name: ");
      scanf("%s",dfile);
      strcpy(file,dfile);
      strcat(file,".map");
      mapfile=fopen(file,"r");
      if(mapfile!=NULL) break;
      printf("File %s does not exist!\n",file);}
    printf("Enter image file name:   ");
    scanf("%s",dfile);
    printf("Enter bad pixel file name:   ");
    fflush(stdout);
    fgets(line,133,stdin);
    if(strlen(line)>=2){
      nargz=1;
      sscanf(line,"%s",bdfile);}
  }
  else{
    nargf=narg=nargz=0;
    for(i=1;i<argc;i++){
      if(!(strcmp(argv[i],"-m"))){
        narg=i+1;
        continue;}
      if(!(strcmp(argv[i],"-z"))){
        nargz=i+1;
        continue;}
      if(!(strcmp(argv[i],"-d"))){
        diag=1;
        continue;}
      if(!(strcmp(argv[i],"-f"))){
        nargf=i+1;
        continue;}
      if(!(strcmp(argv[i],"-o"))){
        strcpy(objsub,argv[i+1]);
        continue;}
      if(!(strcmp(argv[i],"-r"))){
        read_old=1;
        continue;}
      if(!(strcmp(argv[i],"-p"))){
        strcpy(parfile,argv[i+1]);
        continue;}
    }
    if(narg==0 || nargf==0){
      printf("proper invocation: subsky -f framename -m mapfile\n");
      return 1;}
    strcpy(dfile,argv[narg]);
    strcpy(file,dfile);
    strcat(file,".map");
    mapfile=fopen(file,"r");
    if(mapfile==NULL){
      printf("File %s does not exist!\n",file);
      return 1;}
    strcpy(dfile,argv[nargf]);
    if(nargz)  strcpy(bdfile,argv[nargz]);
  }

  //get parameters
  strcpy(file,"subsky");
  if(strcmp(parfile,"")) strcpy(file,parfile);
  if(OpenCosParm(file)!=0) die("Cannot open subsky parameter file");
  if(ReadParm_r("MINLAMBDA",&min_lambda)==1){
    printf("parameterfile error 0!\n");
    return 1;}
  if(ReadParm_r("MAXLAMBDA",&max_lambda)==1){
    printf("parameterfile error 0!\n");
    return 1;}
  if(ReadParm_r("SIGLIMIT",&siglimit)==1){
    printf("parameterfile error 0!\n");
    return 1;}
  if(ReadParm_i("EXCLUDE",&exclude)==1) die("parameter file error");
  if(ReadParm_r("NOISE",&noise)==1){
    printf("parameterfile error 0!\n");
    return 1;}
  if(ReadParm_r("GAIN",&gain)==1){
    printf("parameterfile error 0!\n");
    return 1;}
  if(ReadParm_i("MEDBOX",&medbox)==1){
    printf("parameterfile error 0!\n");
    return 1;}
  if(ReadParm_r("DIAG_0",&diag0)==1){
    printf("parameterfile error 0!\n");
    return 1;}
  if(ReadParm_r("DIAG_1",&diag1)==1){
    printf("parameterfile error 0!\n");
    return 1;}
  if(ReadParm_i("EDGE",&edge)==1){
    printf("parameterfile error 0!\n");
    return 1;}
  if(ReadParm_i("EDGE2",&edge2)==1) edge2=edge;
  if(ReadParm_r("OBJSHIFT",&objshift)==1) objshift=0.0;
  if(ReadParm_r("DELTAKNOT",&delknot)==1) die("parameter file error");
  if(ReadParm_i("SPL_ORDER",&splorder)==1) die("parameter file error");
  if(ReadParm_i("SPACEORDER",&slorder)==1) slorder=1;
  if(ReadParm_b("2D_SPLINE",&twodim)==1) die("parameter file error");
  if(ReadParm_b("SUB_SKY",&sub_sky)==1) sub_sky=1;
//  if(ReadParm_r("s_deltaknot",&s_delknot)==1) die("parameter file error");
//  if(ReadParm_i("s_splineorder",&s_splorder)==1) die("parameter file error");

  medbox/=2;
  medbox=2*medbox+1;
  rnk= (double *) malloc(sizeof(double)*medbox);
  irnk= (int *) malloc(sizeof(int)*medbox);
  stack= (int *) malloc(sizeof(int)*medbox);
  medthrsh=(int)(medbox*0.3);
  mbx1=medbox-1;
  mbx2=medbox/2;
  //slit knots
  if(twodim){
    n_knots_s=4;
    nest_s=n_knots_s+4;
    knot_s= (double *) malloc(sizeof(double)*nest_s);}

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
  if(min_lambda==0 && max_lambda==0){
    min_lambda=lambda0;
    max_lambda=lambda1;}
  else{
    if(lambda0>min_lambda){
      printf("Warning: minimum lambda set to %f\n",lambda0);
      min_lambda=lambda0;}
    if(lambda1<max_lambda){
      printf("Warning: maximum lambda set to %f\n",lambda1);
      max_lambda=lambda1;}
  }
  fgets(line,133,mapfile);
  if(!sscanf(line,"Dewar = %s %d",dewar,&nchip)) die("Mapfile error");
  max_slit=0.;
  n_slit=0;
  rewind(mapfile);

  //temporary fix for non-intuitive cases
  if((strstr(dewar,"E2V") != NULL)||(strstr(dewar,"LDSS3") != NULL)){
    tmp = edge;
    edge = edge2;
    edge2 = tmp;
    tmp2 = -objshift;
    objshift = tmp2;
    }

  //get badpixel data

  //count badpixel groups
  for(i=1;i<=nchip;i++) badpix[i]=0;
  if(!nargz){
    strcpy(bdfile,homedir);
    strcat(bdfile,"sdata/badpix/");
    strcat(bdfile,dewar);}
  strcat(bdfile,".badpix");
  while(1){
    if((badfile=fopen(bdfile,"r"))==NULL){
      printf("Cannot find bad pixel file!\n");
      break;}
    else{
      printf("Found bad pixel file %s\n",bdfile);}
    bpx=1;
    if(!fgets(line,133,badfile)){
      printf("Warning: bad pixel file is empty\n");}
    if(sscanf(line,"%d %d %d %d %d %d",&nchp,&x0,&x1,&y0,&y1,&bbin)<6){
      bbin=bxbin=bybin=1;}
    else{
      if(bbin<10){
        bxbin=bybin=bbin;}
      else{
        i=(int) bbin/10;
        bybin=bbin-10*i;
        bxbin=i;}
    }
    badpix[nchp]++;
    while(1){
      if(!fgets(line,133,badfile)) break;
      sscanf(line,"%d %d %d %d %d",&nchp,&x0,&x1,&y0,&y1);
      badpix[nchp]++;}
    rewind(badfile);
    break;}

  strcpy(flnm,datadir);
  strcat(flnm,dfile);

  //if read_old, check for previous subsky output
  //and open file
  if(read_old){
    for(i=1;i<=nchip;i++){
        strcpy(ofile,flnm);
        addbar(ofile);
        subbars(ofile);
        strcat(ofile,"_s");
        if(nchip>1){
          strcat(ofile,"_c");
          sprintf(c,"%d",i);
          strcat(ofile,c);}
        strcat(ofile,".fits");
        status=OpenFitsFile(ofile,&oldfptr[i],&ofitsinfo);
        if(status) die("cannot find output file");
	}
    }

  //open input files, make output files

  for(i=1;i<=nchip;i++){
    strcpy(file,flnm);
    if(nchip>1){
      addbar(file);
      strcat(file,"c");
      sprintf(c,"%d",i);
      strcat(file,c);}
    strcat(file,".fits");
    status=0;
    status=OpenFitsFile(file,&fptr[i],&fitsinfo);
    if(status) die("cannot find input file");
    if(!diag){
      strcpy(file,"!");
      strcat(file,flnm);
      addbar(file);
      subbars(file);
      strcat(file,"_s");
      if(nchip>1){
        strcat(file,"_c");
        sprintf(c,"%d",i);
        strcat(file,c);}
      strcat(file,".fits");
      fits_create_file(&fptr_o[i],file,&status);}
  }

  //Read ccd data

  naxes[0]=fitsinfo.naxes[0];
  naxes[1]=fitsinfo.naxes[1];
  naxes[2]=2;
  ibin=ibiny=fitsinfo.binning;
  if(fitsinfo.ybinning) ibiny=fitsinfo.ybinning;
  xbin=ibin;
  ybin=ibiny;
  if(!xdisper){
    xbin=ibiny;
    ybin=ibin;}
  nelem=naxes[0]*naxes[1];
  status=0;
  nshuffle=fitsinfo.nshuffle;
  //fixup for new N&S mask design screwup
  sscanf(fitsinfo.date,"%d-%d-%d",&year,&month,&day);
  fdate=year+ (float) month/12+(float) day/365;
  if(fdate > 2007.455) nshuffle=-nshuffle;

  if(bxbin>ibin || bybin>ibiny){
    printf("Warning: bad pixel file has higher binning than data\n");}

  //create badpixel arrays
  if(bpx){
    for(i=1;i<=nchip;i++){
      badx0[i]=malloc(sizeof(INT)*badpix[i]);
      badx1[i]=malloc(sizeof(INT)*badpix[i]);
      bady0[i]=malloc(sizeof(INT)*badpix[i]);
      bady1[i]=malloc(sizeof(INT)*badpix[i]);
      badpix[i]=0;}
    //read badpixel maps
    while(1){
      if(!fgets(line,133,badfile)) break;
      sscanf(line,"%d %d %d %d %d",&nchp,&x0,&x1,&y0,&y1);
      if(x0<1 || y0<1 || x1>naxes[0]*ibin || y1>naxes[1]*ibiny){
        printf("Bad pixel values in badpixel file: %d %d %d %d %d\n",nchp,x0,x1,y0,y1);
        return(1);}
      if(nchp>nchip) break;
      *(badx0[nchp]+badpix[nchp])=(x0-1)/(ibin/bxbin); // iraf to c counting conversion
      *(badx1[nchp]+badpix[nchp])=(x1-1)/(ibin/bxbin);
      *(bady0[nchp]+badpix[nchp])=(y0-1)/(ibiny/bybin);
      *(bady1[nchp]+badpix[nchp])=(y1-1)/(ibiny/bybin);
      badpix[nchp]++;}
  }

  for(i=1;i<=nchip;i++){
    //image arrays
    image[i]=malloc(sizeof(float)*naxes[0]*naxes[1]*2);
        im=image[i];
    array[i]=malloc(sizeof(float *)*naxes[1]);
    for(j=0;j<naxes[1];j++) *(array[i]+j)=im+naxes[0]*j;
    //error arrays follow image arrays
    error[i]=image[i]+naxes[0]*naxes[1];
    err=error[i];
    earray[i]=malloc(sizeof(float *)*naxes[1]);
    for(j=0;j<naxes[1];j++) *(earray[i]+j)=err+naxes[0]*j;
    fits_read_img(fptr[i],TFLOAT,firstel,nelem,&nulvl,im,&anynl,&status);
    if(status)fits_die("Data file error",status);
    //flag bad pixels
    for(l=0;l<badpix[i];l++){
      for(j=*(badx0[i]+l);j<=*(badx1[i]+l);j++){
	for(k=*(bady0[i]+l);k<=*(bady1[i]+l);k++){
//	  *(*(array[i]+k)+j)=nan(&empty);}    //bad pixels set to nan
	  *(*(array[i]+k)+j)=0.0;}    //bad pixels set to zero
        }
      }
    //    fits_close_file(fptr[i],&status);
    printf("Read ccd chip %d\r",i);
    fflush(stdout);}
  printf("\n");

  //read in old output if read_old
  if(read_old){
    naxes[0]=ofitsinfo.naxes[0];
    naxes[1]=ofitsinfo.naxes[1];
    naxes[2]=ofitsinfo.naxes[2];
    onelem=naxes[0]*naxes[1]*naxes[2];
    for(i=1;i<=nchip;i++){
       oimage[i]=malloc(sizeof(float)*naxes[0]*naxes[1]*naxes[2]);
           oim=oimage[i];
       oarray[i]=malloc(sizeof(float *)*naxes[1]);
       for(j=0;j<naxes[1];j++) *(oarray[i]+j)=oim+naxes[0]*j;
       fits_read_img(oldfptr[i],TFLOAT,firstel,onelem,&nulvl,oim,&anynl,&status);
       if(status)fits_die("Data file error",status);
       //fits_close_file(oldfptr[i],&status);
       printf("Read old ccd chip %d\r",i);
       fflush(stdout);}
     printf("\n");
     }

  //if diag output

  if(diag){
    cpgopen("/xwindow");
    wid=6.;
    cpgpap(wid,1.);}


   /*---------------------------Extract Spectra-------------------------------*/

  nslit=-1;
  while((fgets(line,133,mapfile))!=NULL){
    if(sscanf(line,"SLIT %d %s",&slitnum,objname)) continue;
    if(!read_old){
      if(strcmp(objsub,"")){
        if(!(strcmp(objname,objsub))){
          breakpoint=1;}
        else{
          continue;}
        }
      }
    if(sscanf(line,"LENGTH = %f",&slit_len)){
      //found a slit with data
      nslit++;
      if(!(strcmp(objname,objsub))){
        printf("Reducing slit %d \n",slitnum);}
      else{
        printf("Reducing slit %d \r",slitnum);}
      fflush(stdout);
      //for each chip of slit
      while(1){
	if((fgets(line,133,mapfile))==NULL){
	  breakpoint=1;
	  break;}
        if(!sscanf(line,"CHIP %d %f %f %f %f",&chip,&slmin,&slmax,&lmin,&lmax))
	  {//end of slit data
	    sscanf(line,"SLIT %d %s",&slitnum,objname);
	    break;}
	//coef_disp data
	if((fgets(line,133,mapfile))==NULL){
	  printf("Unexpected end to mapfile slit %d b\n",slitnum);
	  return 1;}
	lin=&line[0];
	for(i=0;i<=ord_disp;i++){
	  if((lin=strpbrk(lin,"1234567890.-"))==NULL){
	    printf("Unexpected end to mapfile slit %d c\n",slitnum);
	    return 1;}
	  sscanf(lin,"%f",&coef_disp[i]);
	  if((lin=strpbrk(lin," "))==NULL){
	    printf("Unexpected end to mapfile slit %d d\n",slitnum);
	    return 1;}}
	//coef_xdisp data
	if((fgets(line,133,mapfile))==NULL){
	  printf("Unexpected end to mapfile slit %d b\n",slitnum);
	  return 1;}
	lin=&line[0];
	for(i=0;i<=ord_disp;i++){
	  if((lin=strpbrk(lin,"1234567890.-"))==NULL){
	    printf("Unexpected end to mapfile slit %d c\n",slitnum);
	    return 1;}
	  sscanf(lin,"%f",&coef_xdisp[i]);
	  if((lin=strpbrk(lin," "))==NULL){
	    printf("Unexpected end to mapfile slit %d d\n",slitnum);
	    return 1;}}
	//coef_sag data
	if((fgets(line,133,mapfile))==NULL){
	  printf("Unexpected end to mapfile slit %d e\n",slitnum);
	  return 1;}
	lin=&line[0];
	for(i=0;i<=ord_sag;i++){
	  if((lin=strpbrk(lin,"1234567890.-"))==NULL){
	    printf("Unexpected end to mapfile slit %d f\n",slitnum);
	    return 1;}
	  sscanf(lin,"%f",&coef_sag[i]);
	  if((lin=strpbrk(lin," "))==NULL){
	    printf("Unexpected end to mapfile slit %d g\n",slitnum);
	    return 1;}
	}
	//coef_tilt data
	if((fgets(line,133,mapfile))==NULL){
	  printf("Unexpected end to mapfile slit %d h\n",slitnum);
	  return 1;}
	lin=&line[0];
	for(i=0;i<=ord_tilt;i++){
	  if((lin=strpbrk(lin,"1234567890.-"))==NULL){
	    printf("Unexpected end to mapfile slit %d i\n",slitnum);
	    return 1;}
	  sscanf(lin,"%f",&coef_tilt[i]);
	  if((lin=strpbrk(lin," "))==NULL){
	    printf("Unexpected end to mapfile slit %d j\n",slitnum);
	    return 1;}
	}
	//coef_sagit data
	if((fgets(line,133,mapfile))==NULL){
	  printf("Unexpected end to mapfile slit %d k\n",slitnum);
	  return 1;}
	lin=&line[0];
	for(i=0;i<=ord_sagit;i++){
	  if((lin=strpbrk(lin,"1234567890.-"))==NULL){
	    printf("Unexpected end to mapfile slit %d l\n",slitnum);
	    return 1;}
	  sscanf(lin,"%f",&coef_sagit[i]);
	  if((lin=strpbrk(lin," "))==NULL){
	    printf("Unexpected end to mapfile slit %d m\n",slitnum);
	    return 1;}
	}
	//coef_slen data
	if((fgets(line,133,mapfile))==NULL){
	  printf("Unexpected end to mapfile slit %d n\n",slitnum);
	  return 1;}
	lin=&line[0];
	for(i=0;i<=ord_slen;i++){
	  if((lin=strpbrk(lin,"1234567890.-"))==NULL){
	    printf("Unexpected end to mapfile slit %d o\n",slitnum);
	    return 1;}
	  sscanf(lin,"%f",&coef_slen[i]);
	  if((lin=strpbrk(lin," "))==NULL){
	    printf("Unexpected end to mapfile slit %d p\n",slitnum);
	    return 1;}
	}

	//do trace, get Flux(lambda)
	if(nshuffle) inshuffle = (slmax>slmin) ? nshuffle : -nshuffle;
	if(lmin<min_lambda) lmin=min_lambda;
	if(lmax>max_lambda) lmax=max_lambda;
	if(lmax<=lmin+10.) continue;
	lambda=lmin;
	xcen=polyvalue(lambda,coef_disp,ord_disp)/xbin;
	ix=(int) xcen;
	ixmx=naxes[0];
	iymx=naxes[1];
	if(!xdisper){
	  ixmx=iymx;
	  iymx=naxes[0];}
	if(ix<0)ix=0;
	if(ix>ixmx-1) ix=ixmx-1;
	//direction of dispersion
	dsign = (polyvalue(lmax,coef_disp,ord_disp)/xbin>xcen) ? 1 : -1;
	ssign = (slmax>slmin) ? 1: -1;
        //ixend = (int) polyvalue(lmax,coef_disp,ord_disp)/xbin;

	//allocate arrays
        ixend = (int) fabs(polyvalue(lmax,coef_disp,ord_disp)/xbin-ix);
        ixend = (ixend>ixmx) ? ixmx: ixend;
	//size=((int)(2*ssign*fabs(slmax-slmin+1)/ybin))*((int)dsign*(ixend));
	size=((int)(2*fabs(slmax-slmin+1)/ybin))*(ixend);
  if(nshuffle) size*=2;
	lamtemp=malloc(sizeof(double)*size);
	lamval=malloc(sizeof(double)*size);
	dlamval=malloc(sizeof(float)*size);
	slval=malloc(sizeof(double)*size);
	sltemp=malloc(sizeof(double)*size);
	if(diag) dsplpix=malloc(sizeof(double)*size);
	flxval=malloc(sizeof(double)*size);
	ixval=malloc(sizeof(INT)*size);
	iyval=malloc(sizeof(INT)*size);
	index=malloc(sizeof(long)*size);
	goodpix=malloc(sizeof(double)*size);
	goodtemp=malloc(sizeof(double)*size);
	for(i=0;i<size;i++){
	  *(index+i)=i;
	  *(lamval+i)=0.0;
	  *(dlamval+i)=0.0;
	  *(goodtemp+i)=1.;}

	nint=0;
        nimx=0;

	//loop over all x values for which lambda<max_lambda
	while(ix<ixmx && ix>=0){
	  ycen=polyvalue(lambda,coef_sag,ord_sag)/ybin;
	  slnt=polyvalue(lambda,coef_slen,ord_slen)/ybin;
	  tilt=polyvalue(lambda,coef_tilt,ord_tilt)/xbin;
	  curv=polyvalue(lambda,coef_sagit,ord_sagit)/xbin;

 	  //shift object by 'objshift' pixels for exclude
	  if (chip <= 4) scen = ycen - objshift;
	  else scen = ycen + objshift;
	  iytop=(int)(ycen+slmax*(slnt/(slmax-slmin)));
	  iybot=(int)(ycen+slmin*(slnt/(slmax-slmin)))+1;

	  //move in direction of increasing lambda
	  if(dsign*tilt<0){
	    i1=iybot;
	    i2=iytop;
	    di=ssign;}
	  else{
	    i1=iytop;
	    i2=iybot;
	    di=-ssign;}
	  any=0;
	  ni=abs(i2-i1)+1;
    if(ni>nimx) nimx=ni;
	  i=i1-di;
	  for(iii=0;iii<ni;iii++){
	    i+=di;
	    if(i<0 || i>=iymx) continue;
      yyy = (i-ycen)/slnt;
      ypixmid=(iytop+iybot)/2;
      sagmid=curv*4.*(powf((ypixmid-i)/slnt,2)-powf((ypixmid-ycen)/slnt,2));
      //note that tilt and sagitta signs are reversed compared within
      //etract-spec, because xxx value we want is the xxx that gives the
      //wavelength that is the same as that at object slit position
//      sagmid=-sagmid;
      xxx = ix - (tilt*yyy*ssign*dsign - sagmid);
      lambda=polyvalue(xxx*xbin,coef_xdisp,ord_disp);
      if(isnan(lambda) || lambda>lmax) continue;
	    any=1;

       if (lambda>= lmin && lambda <= lmax) {
             //printf(" %.3f  %.3f %.3f\n",lambda,ix,i);
	   *(lamtemp+nint)=lambda;
           *(sltemp+nint)=iii;
	    //ignore edge and excluded strip
	    if(iii<edge || iii>=ni-edge2) *(goodtemp+nint)=zeroval;
	    if(exclude){
	      if(abs(i-scen)<exclude)  *(goodtemp+nint)=zeroval;}
	    if(xdisper){
	      ixx=ix;
	      iyy=i;}
	    else{
	      ixx=i;
	      iyy=ix;}
	    *(ixval+nint)=ixx;
	    *(iyval+nint)=iyy;
	    nint++;

      if(nshuffle && i+inshuffle>=0 && i+inshuffle<iymx){
        *(ixval+nint)=ixx;
        *(iyval+nint)=iyy+inshuffle;
        nint++;}
      if(nint>=size){
        printf("Bad mapping on slit %d\n",slitnum);
        any=0;
        break;}
    }
    }
    if(!any) break;
    ix+=dsign;}

    //for all other slits, copy old results into array
    if(read_old){
      if(strcmp(objsub,"") && strcmp(objname,objsub)){
	for(i=1;i<nint;i++){
           flptr=*(array[chip]+*(iyval+*(index+i)))+*(ixval+*(index+i));
           if(*flptr != 0.0 && !prob){
             if (*flptr > 0) {
                *(*(earray[chip]+*(iyval+*(index+i)))+*(ixval+*(index+i)))=
                 sqrt(gain*(*flptr)+noise*noise)/gain;
              } else {
                *(*(earray[chip]+*(iyval+*(index+i)))+*(ixval+*(index+i)))= noise/gain;
              }
           } else {
             *flptr=0;
             if(prob){
               *(*(earray[chip]+*(iyval+*(index+i)))+*(ixval+*(index+i)))=noise/gain;}
             else{
               *(*(earray[chip]+*(iyval+*(index+i)))+*(ixval+*(index+i)))=-999999.;}
           }
           *flptr=(float)*(*(oarray[chip]+*(iyval+*(index+i)))+*(ixval+*(index+i)));
	   }
        printf("Copying slit %d \r",slitnum);
        fflush(stdout);
        free(ixval);
        free(iyval);
        free(lamval);
        free(slval);
        free(index);
        free(goodpix);
        free(dlamval);
        free(sltemp);
        free(flxval);
        free(goodtemp);
        if(diag) free(dsplpix);
        continue;
        }
      }

	//sort by lambda


  if(!nshuffle){
    resort = 1;
    nsort = 0;
    osort = 32767;

    while (resort) {
      nsort += 1;
     asort = (long *) malloc(nint*sizeof(long));
     for(k=0;k<nint;k++) asort[k]=index[k];
     qsort(asort, nint, sizeof(long), (void *)comp_nums);
     for(i=0;i<nint;i++) index[i]=asort[i];
     free(asort);

      resort = 0;

      for(i=0;i<nint;i++){
        *(lamval+i)=*(lamtemp+*(index+i));
        *(dlamval+i)=(float) *(lamtemp+*(index+i));
        *(slval+i)=*(sltemp+*(index+i));
        *(goodpix+i)=*(goodtemp+*(index+i));
        //if 2 values identical tweak so curfit doesn't complain
        loff = ((double) rand()/RAND_MAX);
        if(i && *(lamval+i) == *(lamval+i-1)) {
            aa = *(lamval+i);
            bb = (double) (loff * 1e-4);
            cc = aa + bb;
            *(lamtemp+*(index+i)) = cc;
            resort += 1;
            }
        }
      osort = resort;
    }


    }
    for(i=0;i<nint;i++){
      *(flxval+i)=(float)*(*(array[chip]+*(iyval+*(index+i)))+*(ixval+
                                                                *(index+i)));
      //set weights to zero for bad pixels
      //	  if(isnan(*(flxval+i))) *(goodpix+i)=zeroval;}
      if(*(flxval+i)==0.0) *(goodpix+i)=zeroval;}
    free(lamtemp);

    //if nod&shuffle, skip sky fit
    if(nshuffle || !sub_sky) goto A;
    //if(nshuffle) goto A;

    	//DIAG OUTPUT

    if(diag && lmin<diag0 && lmax>diag1){
	  fmin=100000.;
	  fmax=-100000.;
	  for(i=1;i<nint;i++){
//	    if(*(lamval+i)<diag0 || *(lamval+i)>diag1 || isnan(*(flxval+i))) continue;
	    if(*(lamval+i)<diag0 || *(lamval+i)>diag1 || *(flxval+i)==0.0) continue;
	    *(dlamval+i)=(float) (*(lamval+i));
	    *(dsplpix+i)=(*(flxval+i))*(*(goodpix+i));
	    if(*(flxval+i)<fmin)fmin=*(flxval+i);
	    if(*(flxval+i)>fmax)fmax=*(flxval+i);}
	  fmax=fmax+0.1*(fmax-fmin);
	  fmin=fmin-0.1*(fmax-fmin);
	  cpgpage;
	  cpgsch(1);
	  cpgscr(0,1,1,1);
	  cpgscr(1,0,0,0);
	  cpgask(0);
	  cpgsci(1);
	  cpgenv(diag0,diag1,fmin,fmax,0,0);
	  cpglab("Lambda","Flux",objname);
	  cpgpt(nint,dlamval,dsplpix,-1);}

	//take running median and reject CR's

	if(!nint){
	  printf("Unable to fit sky of slit %d, chip %d %d (1)\n",slitnum,chip);
	  continue;}
	//first group
	nmed=0;
	lpix=-1;
	while(nmed<medbox){
	  lpix++;
	  flux=*(flxval+lpix);
	  if(flux != 0.0){
	    dorder(&nmed,&flux,irnk,rnk);
	    nmed++;}
	  if(nmed==mbx2) curnt=lpix;}
	//sigma calculation
	//kelson's way 30th percentile; too sensitive on line wings
	//better way 50th percentile
	medval=*(rnk+mbx2);
	stdev=sqrt(gain*medval+noise*noise)/gain;
	thrsh=medval+siglimit*stdev;
	for(i=0;i<medbox;i++){
	  if(*(irnk+i)==mbx2){
	    if(*(rnk+i)>thrsh) *(goodpix+curnt)=zeroval;
	    break;}
	  }
	//now move down list one (good) pixel at a time
	while(curnt<nint-mbx2){
	  //next center pixel
	  curnt++;
//	  if(isnan(*(flxval+curnt)))continue;
	  if(*(flxval+curnt) == 0.0)continue;
	  //pitch 1st pixel
	  for(i=0;i<medbox;i++){
	    if(!(*(irnk+i))){
	      j=i;
	      break;}
	    }
	  for(i=j+1;i<medbox;i++){
	    *(irnk+i-1)= *(irnk+i);
	    *(rnk+i-1)= *(rnk+i);}
	  for(i=0;i<medbox-1;i++) *(irnk+i)=*(irnk+i)-1;
	  //add new pixel
	  while(lpix<nint-1){
	    lpix++;
	    flux=*(flxval+lpix);
	    if(flux != 0.0){
	      dorder(&mbx1,&flux,irnk,rnk);
	      break;}
	    }
	  medval=*(rnk+mbx2);
	  stdev=sqrt(gain*medval+noise*noise)/gain;
	  tthrsh=medval+siglimit*stdev;
	  //lower cutoff also
	  bthrsh=medval-siglimit*stdev;
	  flux=*(flxval+curnt);
	  if(flux>tthrsh || flux<bthrsh) *(goodpix+curnt)=zeroval;}

	//setup knots

  //lambda knots
if  (xdisper){
	  knotn=*(ixval+nint-1);
	  knot0=*(ixval);}
	else{
	  knotn=*(iyval+nint-1);
	  knot0=*(iyval);}
	n_knots=(int)((knotn-knot0+1)/delknot);
  if(n_knots<0){
	  n_knots=-n_knots;
	  i=knotn;
	  knotn=knot0;
	  knot0=i;}
  nest=nestm=n_knots+2*splorder+2;
	knot=malloc(sizeof(double)*nest);
	ki=knot0-delknot;
	di=delknot;
	if(polyvalue(xbin*ki,coef_xdisp,ord_disp)>polyvalue(xbin*(ki+di),coef_xdisp,
						       ord_disp)){
	  ki=knotn+delknot;
	  di=-delknot;}
        i = 0;
	for(i1=splorder+1;i1<n_knots+splorder+1;i1++){
	  ki+=di;
          loff=polyvalue(xbin*ki,coef_xdisp,ord_disp);
          if ((loff >= *lamval) && (loff <= *(lamval+nint-1))) {
            ii = i+splorder+1;
	    *(knot+ii) = loff;
            //printf("%d %f %f %f\n", i,*(knot+ii),*lamval,*(lamval+nint-1));
            i += 1;
           }
         }
        n_knots = i;
	lval1=*lamval;
	lvaln=*(lamval+nint-1);
  //check for correct positioning of knots

        if (*(knot+splorder+1)<= lval1) {
           *(knot+splorder+1)+= 0.5*((*knot+splorder+2)-(*knot+splorder+1));
           }
        if (*(knot+n_knots+splorder)>= lvaln) {
           *(knot+n_knots+splorder) -= 0.5*((*knot+n_knots+splorder)- (*knot+n_knots+splorder-1));
           }
        if (lval1>= *(knot+splorder+1)) {
           lval1= *(knot+splorder+1)-1;
           }
        if (lvaln<= *(knot+n_knots+splorder)) {
           lvaln= *(knot+n_knots+splorder)+1.;
           }

//  slorder=1;
  wlorder=splorder;
  uu=nest-wlorder-1;
  vv=nest_s-slorder-1;
  morder = (wlorder>slorder) ? wlorder+1: slorder+1;
  bx=wlorder*vv+slorder+1;
  by=slorder*uu+wlorder+1;
  if(bx<=by){
    b1=bx;
    b2=b1+vv-slorder;}
  else{
    b1=by;
    b2=b1+uu-wlorder;}
  if(twodim){
    lwork1=uu*vv*(2+b1+b2)+2*(uu+vv+morder*(nint+nestm)+nestm-
                              wlorder-slorder)+b2+1;}
  else{
    lwork1=nint*(splorder+1)+nest*(7+3*splorder);}
  if(!(work1=realloc(work1,sizeof(double)*lwork1)))
                                    die("\nMemory failure work1");
  if(twodim){
    slval1=-.3;
    slvaln=nimx+.3;
    *(knot_s+2)=-0.2;
    *(knot_s+3)=-0.1;
    *(knot_s+4)=nimx+0.1;
    *(knot_s+5)=nimx+0.2;
    lwork2=uu*vv*(6*vv+2)+6*vv+1;
    if(!(work2=realloc(work2,sizeof(double)*lwork2)))
      die("\nMemory failure work2");
    kwork=nint+(nest-2*wlorder-1)*(nest_s-2*slorder-1);
    coefsz=nest*nest_s;}
  else{
    kwork=nest;
    coefsz=nest;}
  if(!(iwork=realloc(iwork,sizeof(INT)*kwork)))
    die("\nMemory failure iwork");
  if(!(coefs=realloc(coefs,sizeof(double)*coefsz)))
    die("\nMemory failure coef");
  iopt=-1;
  ier=prob=0;
  v1=1.0e-10;
  //fit spline to spectrum
//  for(i=0;i<nint;i++) if(isnan(*(flxval+i))) *(flxval+i)=0;
  for(i=0;i<nint;i++) if(*(flxval+i) == 0.0) *(flxval+i)=0;
  //printf("Fitting slit %d\r",slitnum);
  //fflush(stdout);
  if(twodim){
    surfit_(&iopt,&nint,lamval,slval,flxval,goodpix,&lval1,&lvaln,&slval1,
            &slvaln,&wlorder,&slorder,&smfac,&nest,&nest_s,&nestm,&v1,
            &n_knots,knot,&n_knots_s,knot_s,coefs,&resid,work1,&lwork1,
            work2,&lwork2,iwork,&kwork,&ier);
    if(ier>0 || isnan(*(coefs+n_knots_s*n_knots/2))){
      printf("Problem fiting sky of slit %d, chip %d %d (1)\n",slitnum,chip,
             ier);
      prob=1;}
    else{
      bispev2_(knot,&n_knots,knot_s,&n_knots_s,coefs,&wlorder,&slorder,
         lamval,slval,&nint,flxval,&ier);
    }
  }
  else{
    curfit_(&iopt,&nint,lamval,flxval,goodpix,&lval1,&lvaln,&wlorder,&smfac,
            &nest,&n_knots,knot,coefs,&resid,work1,&lwork1,iwork,&ier);
    if(ier || isnan(*(coefs+n_knots/2))){
      printf("Problem fiting sky of slit %d, chip %d %d (1)\n",slitnum,chip,
             ier);
      prob=1;}
    else{
      splev_(knot,&n_knots,coefs,&splorder,lamval,flxval,&nint,&ier);}
      if(diag && lmin<diag0 && lmax>diag1){
         for(i=0;i<nint;i++){
          *(dsplpix+i)=(float) *(flxval+i);
        }
      }
    if(ier){
      printf("Problem fiting sky of slit %d, chip %d %d (1)\n",slitnum,chip,
             ier);
      prob=1;}
  }


     //DIAG OUTPUT

  if(diag && lmin<diag0 && lmax>diag1){
    cpgsci(2);
    cpgline(nint,dlamval,dsplpix);
    printf("\nHit CR to continue, Q to quit: ");
    fflush(stdout);
    fgets(answer,80,stdin);
    if(!strncasecmp(answer,"q",1)) return 0;
//    cpgclos();
    }
  free(knot);

  //Subtract sky, Calculate errors
A: for(i=1;i<nint;i++){
     flptr=*(array[chip]+*(iyval+*(index+i)))+*(ixval+*(index+i));
     if(*flptr != 0.0 && !prob){
       if (*flptr > 0) {
          *(*(earray[chip]+*(iyval+*(index+i)))+*(ixval+*(index+i)))=
           sqrt(gain*(*flptr)+noise*noise)/gain;
        } else {
          *(*(earray[chip]+*(iyval+*(index+i)))+*(ixval+*(index+i)))= noise/gain;
        }
       if(!nshuffle && sub_sky) *flptr-= (float) *(flxval+i);
     } else {
       *flptr=0;
       if(prob){
         *(*(earray[chip]+*(iyval+*(index+i)))+*(ixval+*(index+i)))=noise/gain;}
       else{
         *(*(earray[chip]+*(iyval+*(index+i)))+*(ixval+*(index+i)))=-999999.;}
     }
   }

  free(ixval);
  free(iyval);
  free(lamval);
  free(slval);
  free(index);
  free(goodpix);
  free(dlamval);
  free(sltemp);
  free(flxval);
  free(goodtemp);
  if(diag) free(dsplpix);
      }
    }
    if(breakpoint) break;}

  //Write out image files

  if(diag) return 0;
  printf("\nWriting output files\n");
  for(i=1;i<=nchip;i++){
    //printf("Writing ccd chip %d\r",i);
    im=image[i];
    err=error[i];
    status=0;
    //write image array
    firstelem[2]=1;
    nelem2=nelem*2;
    fits_create_img(fptr_o[i],FLOAT_IMG,three,naxes,&status);
    if(read_old){
      fits_write_pix(fptr_o[i],TFLOAT,firstelem,nelem2,oim,&status);
      }
    fits_write_pix(fptr_o[i],TFLOAT,firstelem,nelem2,im,&status);
    //write error array
    status=0;
    //firstelem[2]=2;
    //fits_write_pix(fptr_o[i],TFLOAT,firstelem,nelem,err,&status);
    //fix up header
    fits_delete_record(fptr_o[i],8,&status);
    fits_delete_record(fptr_o[i],8,&status);
    //fits_delete_record(fptr_o[i],6,&status);
    status=0;
    ncard=0;
    while(1){
      ncard++;
      if(ncard>108) break;
      fits_read_record(fptr[i],ncard,cardline,&status);
      if(status){fits_die("Error reading file header",status);}
      if(strstr(cardline,"END")==cardline) break;
      if(strstr(cardline,"NAXIS")!=NULL) continue;
      if(strstr(cardline,"SIMPLE")!=NULL) continue;
      if(strstr(cardline,"BITPIX")!=NULL) continue;
      if(strstr(cardline,"BSCALE")!=NULL) continue;
      if(strstr(cardline,"BZERO")!=NULL) continue;
      status=0;
      fits_write_record(fptr_o[i],cardline,&status);
      if(status) die("Error writing output file");}
    //if(status){
    //   printf("Error writing data for chip %d %d\n",i,status);
    //return 1;}
    printf("Wrote ccd chip %d\r",i);
    fflush(stdout);
    fits_close_file(fptr_o[i],&status);
    status=0;
    free(image[i]);
    free(array[i]);
    if(read_old){
      free(oimage[i]);
      free(oarray[i]);
      }
    }
  printf("\n");

  return 0;}

int comp_nums(const long *num1, const long *num2)
{
   if (*(lamtemp+*num1) <  *(lamtemp+*num2)) return -1;
   if (*(lamtemp+*num1) == *(lamtemp+*num2)) return  0;
   if (*(lamtemp+*num1) >  *(lamtemp+*num2)) return  1;
}
