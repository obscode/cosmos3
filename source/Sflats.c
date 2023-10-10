/*****************************************************************************\
*                                                                             *
*  Sflats sums and normalizes spectroscopic flats, with CR rejection          *
*                                                                             *
*  VERSION  08 Aug  2005                                                      *
*                                                                             *
\*****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fitsio.h"
#include "cosmos.h"


int main(int argc,char *argv[]){

  char  file[20][80],flnm[1000],ifile[80],dfile[80],*datadir,line[133],*lin,c[7],
        ftype[6],jfile[1000],ofile[80],answer[80],CHAR,*ch,*name,dewar[10],
        bfile[80],bdfile[80],*homedir,date[12];
  int   pfit,up,INT,status,bitpix,naxis,i,j,k,here,xdisper,ord_disp,ord_sag,
        ntot,medbox,dwn,nbias,nolm,ord_tilt,ord_sagit,ord_slen,nslit,n_slit,
        ord_flat,chip,slitnum,max_slit,n80,breakpoint,iii,ni,ixpos,iypos,lvmin,
        lvmax,dsign,ij,anynl,dimen,clean,proff,search,strsz,strip,b1,b2,uu,vv,
        splf,ssign,mdbx1,is,ix,ixmx,iymx,ixx,iyy,nint,nexty,iytop,iybot,i1,i2,
        any,jj,ii,size,lpix,minpixnum,iwl0,iwl1,nsv,isv,nwl,ilv,xsz,ysz,biasx0,
        biasx1,biasy0,biasy1,mdbx2,xsize,ysize,nchp,nlv,nshuffle,skysub,nfile,
        srcsz,iycen,n200,n1,ibin,nxtpx,overscan,biaslins,filen,iwl,iwlmx,nbpx,
        nchip,islmx,islcen,nax0,nax1,inshuffle,optimal,nk1,nms,narg,badpix[9],
        x0,x1,y0,y1,bpx,l,nargz,ngood,shfld,xshuffle,yshuffle,ibiny,lbin,sbin,
        bbin,bxbin,bybin,year,month,day;
  float FLOAT,delta_lambda,delta_slit,telscale,min_lambda,max_lambda,lambda0,
        flax,awl,SCALE,v1,lambda1,slen,slit_len,slmin,slmax,lmin,lmax,sloff,
        dltaslt,lv,lambda,xcen,ycen,slnt,scatl,slitend,slitpos,ypos,yp,xpos,
        wlcen,yval1,yvaln,wid,fmin,fmax,dix,diy,value,d11,d21,d12,d22,tilt,
        coef_disp[10],coef_sag[10],coef_tilt[10],coef_sagit[10],coef_slen[10],
        coef_flat[10],coef_xdisp[10],resid,lval1,lvaln,flux,medval,stdev,thrsh,
        siglimit,noise,gain,xxx,exptime[20],exptot,flx,slcen,llow,lhigh,tthrsh,
        bthrsh,smfac,delknotl,delknoty,ki,trprev,di,fdate;
  long  newaxes[2],naxes[3],noxes[2],firstel,firstelem[3],stack_len,fpixel[2],
        nelem,noutlm,flen,lpixel[2],nuxes[2];
  int   *INTP,**good,*gd,*gimage,*sltnum,*isloff,*noff,
        *upper,*lower,*ixval,*iyval,*badx0[9],*badx1[9],*bady0[9],*bady1[9],
        *irnk;
  float *FLOATP,**stack,**spectrum,*sm,**sumspc,*image[20],*im,*oimage,*oim,
        **oarray,**array[20],*s_image[20],*s_im,**s_array[20],*bimage,**barray,
        *trace,*medtrc,*ytemp,*prfl,*index,*noefs,*nootl,*knotl,*lamval,*rnk,
        *lamtemp,*ptl,*ptx,*nooty,*work2,*knoty,*yval,*goodtemp,*gdpxg;
  int   nulvl=0;
  double   DOUBLE,*xlist,*ylist,**elist,a1;
  fitsfile *fptr_obj[20][9],*outptr,*bptr[9];
  fitsdef  fitsinfo;
  FILE     *mapfile,*datafile,*badfile;
  float    llimit;

  //parameters
  breakpoint=bpx=0;
  exptot=0.;
  nbpx=100;                      //maximum number of bad pixel groups per chip
  smfac=0.;
  nulvl=0;
  firstelem[0]=firstelem[1]=firstelem[2]=1;
  firstel=1;
  ch=NULL;
  inshuffle=xshuffle=yshuffle=shfld=0;
  SCALE=0.015; //pixel size, in mm
  minpixnum=50;//minimum # pixels on a chip
  ofile[0]='.';
  bbin=bxbin=bybin=1;

   //directories
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

  //get parameters
  optimal=clean=0;
  if(OpenCosParm("Sflats")!=0) die("Cannot open parameter file");
  if(ReadParm_r("MINLAMBDA",&min_lambda)==1)die("parameter file error");
  if(ReadParm_r("MAXLAMBDA",&max_lambda)==1)die("parameter file error");
  if(ReadParm_i("ORDER",&ord_flat)==1) die("parameter file error");
  nbias=0;
  if(ReadParm_i("MED_BOX",&medbox)==1) die("parameter file error");
  mdbx1=medbox-1;
  mdbx2=medbox/2;
  if(ReadParm_r("FLOOR",&llimit)==1) llimit = 0.1;
  if(ReadParm_s("BIAS_FILE",bfile)==0 && strcasecmp(bfile,"none")) nbias=1;
  if(ReadParm_s("FIT_MODE",ftype)==1) die("parameter file error");
  pfit=1;
  if(!strcasecmp(ftype,"median")) pfit=0;
  if(ReadParm_i("SHUFFLED",&shfld)==1)die("parameter file error");

 //get input files
  if(argc<5){
    while(1){
      printf("Enter map file name: ");
      scanf("%s",dfile);
      strcpy(ifile,dfile);
      strcat(ifile,".map");
      mapfile=fopen(ifile,"r");
      if(mapfile!=NULL) break;
      printf("File %s does not exist!\n",file);}
    nfile=0;
    while(1){
      printf("Enter flat file names CR=>end:   ");
      fflush(stdout);
      fgets(line,133,stdin);
      if(strlen(line)<2){
	if(nfile>0) break;
	else continue;}
      addbar(line);
      strcpy(file[nfile],line);
      nfile++;}
    if(nfile>2){
      printf("Enter output file name: ");
      scanf("%s",ofile);}
    printf("Enter bias file name:   ");
    fflush(stdout);
    fgets(line,133,stdin);
    if(strlen(line)>=2){
      nbias=1;
      sscanf(line,"%s",bfile);}
    printf("Enter bad pixel file name:   ");
    fflush(stdout);
    fgets(line,133,stdin);
    if(strlen(line)>=2){
      nargz=1;
      sscanf(line,"%s",bdfile);}
  }
  else{
    narg=0;
    for(i=1;i<argc;i++){
      if(!(strcmp(argv[i],"-b"))){
	narg=i+1;
	break;}
      }
    if(narg){
      if(strcasecmp(argv[narg],"none")){
	nbias=1;
	strcpy(bfile,argv[narg]);}
      }
    narg=nargz=0;
    for(i=1;i<argc;i++){
      if(!(strcmp(argv[i],"-m"))){
        narg=i+1;
        continue;}
      if(!(strcmp(argv[i],"-z"))){
        nargz=i+1;
        continue;}
            }
    if(!narg)die("proper invocation:Sflats -f frame[s] -m mapfile -b biasfile [-o output file]\n");
    strcpy(dfile,argv[narg]);
    strcpy(flnm,dfile);
    strcat(flnm,".map");
    mapfile=fopen(flnm,"r");
    if(mapfile==NULL){
      printf("File %s does not exist!\n",flnm);
      return 1;}
    if(nargz) strcpy(bdfile,argv[nargz]);
    if(!(strcmp(argv[1],"-f"))){
      narg=2;}
    else{
      if(argc>3 && !(strcmp(argv[3],"-f"))){
	narg=4;}
      else{
        if(argc>5 && !(strcmp(argv[5],"-f"))){
          narg=6;}
        else{
          if(argc>7 && !(strcmp(argv[7],"-f"))){
            narg=8;}
          else{
            printf("proper invocation: Sflats -f framename[s] -m mapfile -b biasfile [-o output file]\n");
            return 1;}
        }
      }
    }
    nfile=0;
    while(1){
      strcpy(file[nfile],argv[narg]);
      nfile++;
      if(narg==argc-1  || !strcmp(argv[narg+1],"-b") || !strcmp(argv[narg+1],"-m") || !strcmp(argv[narg+1],"-o") || !strcmp(argv[narg+1],"-z")) break;
      narg++;}
    narg=0;
    for(i=1;i<argc;i++){
      if(!(strcmp(argv[i],"-o"))){
        narg=i+1;
        break;}
        }
      if(!narg && nfile>2)die("proper invocation: Sflats -f frame[s] -m mapfile -b biasfile -o output file\n");
    if(narg) strcpy(ofile,argv[narg]);
    }

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
  if(!sscanf(line,"Scale ~ %f",&telscale)){
    printf("Mapfile error\n");
    return 1;}
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
  fgets(line,133,mapfile);
  xxx=0.;
  n_slit=0;

  //find #of mapped slits, width of widest
  while((fgets(line,133,mapfile))!=NULL){
    if(!sscanf(line,"LENGTH = %f",&slen)) continue;
    n_slit++;
    if(slen>xxx) xxx=slen;}
  max_slit=1+(int)(xxx/SCALE);
  isloff=malloc(n_slit*sizeof(INT));
  sltnum=malloc(n_slit*sizeof(INT));
  noff=malloc(n_slit*sizeof(INT));
  n80=80*n_slit;
  name=malloc(n80*sizeof(CHAR));
  //buffers
  nsv=max_slit;
  if(nsv<medbox+5) nsv=medbox+5;
  rnk=malloc(sizeof(FLOAT)*nsv);
  irnk=malloc(sizeof(INT)*nsv);

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

  //open all data files

  exptot=0;
  for(i=0;i<nfile;i++){
    strcpy(flnm,datadir);
    strcat(flnm,file[i]);
    addbar(flnm);
    for(j=1;j<=nchip;j++){
      strcpy(jfile,flnm);
      if(nchip>1){
        strcat(jfile,"c");
        sprintf(c,"%d",j);
        strcat(jfile,c);}
      strcat(jfile,".fits");
      status=0;
      status=OpenFitsFile(jfile,&fptr_obj[i][j],&fitsinfo);
      if(status){
        printf("file %s does not exist!",jfile);
        return 1;}
      if(!strcmp(dewar,"LDSS3") && nchip==1){
        fits_read_key(fptr_obj[i][j],TINT,"NMOSAIC",&nms,line,&status);
        if(status) die("flat file was not properly stitched");}
    }
    exptot+=fitsinfo.exptime;
    exptime[i]=fitsinfo.exptime;}
  naxes[0]=fitsinfo.naxes[0];
  naxes[1]=fitsinfo.naxes[1];
  overscan=fitsinfo.overscan;
  biaslins=fitsinfo.biaslins;
  bitpix=fitsinfo.bitpix;
  nshuffle=fitsinfo.nshuffle;
  if(shfld){
    if(nshuffle){
      die("Cannot shuffle already shuffled images!");}
    else{
      nshuffle=shfld;
      //fixup for new N&S mask design screwup
      sscanf(fitsinfo.date,"%d-%d-%d",&year,&month,&day);
      fdate=year+ (float) month/12+(float) day/365;
      if(fdate > 2007.455) nshuffle=-nshuffle;
    }
  }

  bitpix=abs(bitpix);
  nelem=naxes[0]*naxes[1];
  xsize=naxes[0]-overscan;
  ysize=naxes[1]-biaslins;
  noutlm=xsize*ysize;
  nuxes[0]=xsize;
  nuxes[1]=ysize;
  ibin=ibiny=fitsinfo.binning;
  if(fitsinfo.ybinning) ibiny=fitsinfo.ybinning;
  lbin=ibin;
  sbin=ibiny;
  if(!xdisper){
    lbin=ibiny;
    sbin=ibin;}
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
      x0 = (x0-1)/(ibin/bxbin); // iraf to c counting conversion
      x1 = (x1-1)/(ibin/bxbin);
      y0 = (y0-1)/(ibiny/bybin);
      y1 = (y1-1)/(ibiny/bybin);
      if(x0<0 || y0<0 || x1>nuxes[0]-1 || y1>nuxes[1]-1){
        printf("Bad pixel values in badpixel file: %d %d %d %d %d %d %d\n",nchp,x0,x1,nuxes[0],y0,y1,nuxes[1]);
        return(1);}
      if(nchp>nchip) break;
      *(badx0[nchp]+badpix[nchp])=x0;
      *(badx1[nchp]+badpix[nchp])=x1;
      *(bady0[nchp]+badpix[nchp])=y0;
      *(bady1[nchp]+badpix[nchp])=y1;
      badpix[nchp]++;}
  }

  //constuct output file name

  if(ofile[0]=='.'){
    if(nfile==1){
      strcpy(ofile,file[0]);}
    else{
      if(nfile==2){
        i=strlen(file[0]);
        if(strlen(file[1])<i) i=strlen(file[1]);
        for(j=0;j<i;j++){
          if(file[0][j] != file[1][j]) break;}
        strcpy(ofile,file[0]);
        strcat(ofile,"-");
        strcat(ofile,file[1]+j);}
    }
  }
  strcat(ofile,"_flat");

  //bias frame
  if(nbias){
    bimage=malloc(sizeof(FLOAT)*nelem);
    barray=malloc(sizeof(FLOATP)*naxes[1]);
    for(j=0;j<naxes[1];j++) *(barray+j)=bimage+naxes[0]*j;
    for(j=1;j<=nchip;j++){
      strcpy(flnm,datadir);
      strcat(flnm,bfile);
      if(nchip>1){
        strcat(flnm,"c");
        sprintf(line,"%d",j);
        strcat(flnm,line);}
      strcat(flnm,".fits");
      status=0;
      fits_open_file(&bptr[j],flnm,READONLY,&status);
      if(status) die("error opening bias file");}
    }

  //bounds, trace buffers
  ilv = (xdisper) ? xsize: ysize;
  upper=calloc(ilv,sizeof(INT));
  lower=calloc(ilv,sizeof(INT));
  trace=calloc(ilv,sizeof(FLOAT));
  medtrc=calloc(ilv,sizeof(FLOAT));
  index=calloc(ilv,sizeof(FLOAT));
  for(i=0;i<ilv;i++) index[i]=i;

  //input buffers for each data frame
  for(i=0;i<nfile;i++){
    image[i]=calloc(nelem,sizeof(FLOAT));
    array[i]=malloc(sizeof(FLOATP)*naxes[1]);
    for(j=0;j<naxes[1];j++) *(array[i]+j)=image[i]+j*naxes[0];}
  //bad pixel array
  gimage=calloc(nelem,sizeof(INT));
  for(j=0;j<nelem;j++) *(gimage+j)=1;
  good=malloc(sizeof(INTP)*naxes[1]);
  for(j=0;j<naxes[1];j++) *(good+j)=gimage+j*naxes[0];
  //output buffer
  oimage=calloc(noutlm,sizeof(FLOAT));
  oarray=malloc(sizeof(FLOATP)*ysize);
  for(j=0;j<ysize;j++) *(oarray+j)=oimage+j*xsize;
  //fit buffers
  xlist=calloc(2*(ord_flat+1),sizeof(DOUBLE));
  ylist=calloc(ord_flat+1,sizeof(DOUBLE));
  elist=calloc(ord_flat+1,sizeof(xlist));
  for(i=0;i<ord_flat+1;i++){
    *(elist+i)=malloc(sizeof(DOUBLE)*(ord_flat+2));}


  /*------------process spectral flats, one chip at a time-------------------*/

  for(nchp=1;nchp<=nchip;nchp++){

    //flag bad pixels
    for(l=0;l<badpix[nchp];l++){
      for(j=*(badx0[nchp]+l);j<=*(badx1[nchp]+l);j++){
        for(k=*(bady0[nchp]+l);k<=*(bady1[nchp]+l);k++){
          *(*(good+k)+j)=0;}
      }
    }
    //output image file
    strcpy(ifile,"!");
    strcat(ifile,datadir);
    strcat(ifile,ofile);
    if(nchip>1){
      strcat(ifile,"_");
      strcat(ifile,"c");
      sprintf(c,"%d",nchp);
      strcat(ifile,c);}
    strcat(ifile,".fits");
    status=0;
    fits_create_file(&outptr,ifile,&status);
    bitpix=-32;
    naxis=2;
    fits_create_img(outptr,bitpix,naxis,nuxes,&status);
    //fits_update_key(outptr,TINT,"BIASLINS",0,ch,&status);
    //fits_update_key(outptr,TINT,"OVERSCAN",0,ch,&status);
    if(!strcmp(dewar,"LDSS3") && nchip==1){
      fits_write_key(outptr,TINT,"NMOSAIC",&nchip,"Number of chips mosaiced",&status);}
    if(status) fits_die("Unable to create output spectrum file",status);

    //clear output array
    for(i=0;i<xsize*ysize;i++) *(oimage+i)=1.0;
    //bias file
    if(nbias){
      fits_read_img(bptr[nchp],TFLOAT,firstel,nelem,&nulvl,bimage,&anynl,
		    &status);
      if(status) fits_die("Error readidng bias file",status);
      fits_close_file(bptr[nchp],&status);}
    //read in chip data
    status=0;
    for(i=0;i<nfile;i++){
      fits_read_img(fptr_obj[i][nchp],TFLOAT,firstel,
                    nelem,&nulvl,image[i],&anynl,&status);
      if(status) fits_die("Error reading data file",status);
      fits_close_file(fptr_obj[i][nchp],&status);
      //subtract bias
      if(nbias){
        for(j=0;j<naxes[1];j++){
          for(k=0;k<naxes[0];k++) *(*(array[i]+j)+k)-= *(*(barray+j)+k);}
      }
      if(overscan && biaslins){
        xsz=xsize;
        ysz=ysize;
        biasx0=xsz+3;
        biasx1=(xsize+overscan)-3;
        biasy0=ysz+3;
        biasy1=(ysize+biaslins)-3;
        subbias(array[i],naxes,biasx0,biasx1,biasy0,biasy1,xsz,ysz);}
    }

      rewind(mapfile);
      nslit=-1;
      //coef_disp data
      while(1){
        if((fgets(line,133,mapfile))==NULL)  die("Unexpected end to mapfile");
        if(sscanf(line,"SLIT %d %s",&slitnum,name+80*(nslit+1))) break;}

      //New slit
      while(2){
        if(!strncmp(line,"END",3)) break;                         //end of data
        if(!sscanf(line,"SLIT %d %s",&slitnum,name)) die("Mapfile error");
        sltnum[nslit+1]=slitnum;
        if((fgets(line,133,mapfile))==NULL) die("Unexpected end to mapfile");
        if(!strncmp(line,"END",3)) break;                        //end of data
        if(!sscanf(line,"LENGTH = %f",&slit_len)) continue;   //no data in slit
      slitnum+=1;
      //found a slit with data
      nslit++;

      //New Chip
      while(3){
	if((fgets(line,133,mapfile))==NULL) die("Unexpected end to mapfile");
	if(!strncmp(line,"END",3)) break;
	if(!sscanf(line,"CHIP %d %f %f %f %f",&chip,&slmin,&slmax,&lmin,&lmax))
	  break;
	if(nshuffle){
    inshuffle = (slmax>slmin) ? nshuffle : -nshuffle;
    if(!shfld){
      if(xdisper){
        yshuffle=inshuffle;}
      else{
        xshuffle=inshuffle;}
    }
  }
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
	if(chip != nchp || min_lambda>lmax || max_lambda<lmin) continue;
	/*nb notation l=>lambda
                      wl=>unbinned pixels in lambda direction
		      sl=>unbinned pixels in slit direction
		      x=> binned x pixels
		      y=> binned y pixels
	*/

	//do trace, find limits of spectrum, flux vs x

	//direction of dispersion
	fpixel[0]=fpixel[1]=10000.;
	lpixel[0]=lpixel[1]=0;
	lambda=lmin;
	if(lambda<min_lambda) lambda=min_lambda;
	wlcen=polyvalue(lambda,coef_disp,ord_disp);
	if(wlcen < 0) wlcen=0;
	if(xdisper){
	  if(wlcen>=ibin*xsize) wlcen=ibin*xsize-1;}
	else{
	  if(wlcen>=ibiny*ysize) wlcen=ibiny*ysize-1;}
	awl=wlcen;
	iwl0=(int)awl;
	dsign = (polyvalue(lmax,coef_disp,ord_disp)>wlcen) ? 1 : -1;
	ssign = (slmax>slmin) ? 1: -1;
	awl-=dsign*lbin*0.5;
	iwl=(int)awl;
	iwlmx=xsize*ibin-1;
	islmx=ysize*ibiny-1;
	if(!xdisper){
	  iwlmx=ysize*ibiny-1;
	  islmx=xsize*ibin-1;}
	if(wlcen>iwlmx) wlcen=iwlmx;

	k=0;
	while(1){
	  awl+=dsign*lbin*0.5;
	  iwl=(int)awl;
	  lambda=polyvalue(iwl,coef_xdisp,ord_disp);
	  if(lambda>max_lambda || lambda>lmax || iwl<0 || iwl>iwlmx){
	    iwl1=iwl-dsign;
	    break;}
	  slnt=polyvalue(lambda,coef_slen,ord_slen);
	  slcen=polyvalue(lambda,coef_sag,ord_sag);
	  ilv=iwl/lbin;
	  isv=slcen/sbin;
	  if(xdisper) {
	    ixx=(iwl-ssign*slnt/2.)/ibin;
	    iyy= (ssign>0) ? slcen+slnt*slmin/(slmax-slmin):
	                     slcen+slnt*slmax/(slmax-slmin);
	    iyy/=ibiny;
	    lower[ilv]=iyy;
	    if(lower[ilv]<0)lower[ilv]=0;
	    upper[ilv]=iyy+ssign*slnt/ibiny+1;
	    if(upper[ilv]>ysize-1)upper[ilv]=ysize-1;
	    dwn=isv-(isv-iyy)/1.5;
	    if(dwn<0) dwn=0;
	    up=isv+(upper[ilv]-isv)/1.5;
	    if(up>=naxes[1])up=naxes[1]-1;
	    nsv=(up-dwn+1);
            ngood=0;
	    for(k=0;k<nsv;k++){
               if(!(*(*(good+k+dwn)+ilv))) continue;
	         order(&ngood,*(array[0]+k+dwn)+ilv,irnk,rnk);
               ngood++;}
	    ngood/=2;
            if(ngood){
              trprev=trace[ilv]=*(rnk+ngood)/exptime[0];
            } else {
              trace[ilv]=trprev;}
          }
    else {
	    iyy=(iwl-ssign*slnt/2.)/ibiny;
            ixx= (ssign>0) ? slcen+slnt*slmin/(slmax-slmin):
            slcen+slnt*slmax/(slmax-slmin);
            ixx/=ibin;
            lower[ilv]=ixx;
            upper[ilv]=ixx+ssign*slnt/ibin+1;
            if(lower[ilv]<0)lower[ilv]=0;
            if(upper[ilv]>xsize-1)upper[ilv]=xsize-1;
            dwn=isv-(isv-ixx)/1.5;
            if(dwn<0) dwn=0;
            up=isv+(upper[ilv]-isv)/1.5;
            if(up>=naxes[0]) up=naxes[0]-1;
            nsv=(up-dwn+1);
            ngood=0;
            for(k=0;k<nsv;k++){
              if(!(*(*(good+ilv)+k+dwn))) continue;
      	      order(&ngood,*(array[0]+ilv)+k+dwn,irnk,rnk);
              ngood++;}
            ngood/=2;
            if(ngood){
              trprev=trace[ilv]=*(rnk+ngood)/exptime[0];
            } else {
              trace[ilv]=trprev;}
           }
	  if(fpixel[0]>ixx+xshuffle) fpixel[0]=ixx+xshuffle;
    if(fpixel[0]>ixx) fpixel[0]=ixx;
	  if(fpixel[1]>iyy+yshuffle) fpixel[1]=iyy+yshuffle;
    if(fpixel[1]>iyy) fpixel[1]=iyy;
	  ixx+=ssign*slnt/ibin;
	  iyy+=ssign*slnt/ibiny;
	  if(lpixel[0]<ixx+xshuffle) lpixel[0]=ixx+xshuffle;
    if(lpixel[0]<ixx) lpixel[0]=ixx;
	  if(lpixel[1]<iyy+yshuffle) lpixel[1]=iyy+yshuffle;
    if(lpixel[1]<iyy) lpixel[1]=iyy;}
	if(!k){
	  printf("failed to find slitbox %d %d\n",chip,slitnum);
	  continue;}
	fpixel[0]-=3;
	fpixel[1]-=3;
	lpixel[0]+=4;
	lpixel[1]+=4;
	if(fpixel[0]<0) fpixel[0]=0;
	if(lpixel[0]>=xsize) lpixel[0]=xsize-1;
	if(fpixel[1]<0) fpixel[1]=0;
	if(lpixel[1]>=ysize) lpixel[1]=ysize-1;
	noxes[0]=nax0=lpixel[0]-fpixel[0]+1;
	noxes[1]=nax1=lpixel[1]-fpixel[1]+1;
/*	printf("noxes = %d %d\n",noxes[0],noxes[1]);
	printf("fpixel = %d %d\n",fpixel[0],fpixel[1]);
	printf("lpixel = %d %d\n\n",lpixel[0],lpixel[1]); */
	nolm=nax0*nax1;

	iwl=(int) wlcen;
	proff=0;

	/*----------------------------average--------------------------------*/

	sm=calloc(nolm,sizeof(FLOAT));
	sumspc=malloc(sizeof(FLOATP)*nax1);
	for(j=0;j<nax1;j++) *(sumspc+j)=sm+nax0*j;
  for(filen=0;filen<nfile;filen++){
    for(j=0;j<noxes[1];j++){
      for(i=0;i<noxes[0];i++){*(*(sumspc+j)+i) +=
        *(*(array[filen]+j+fpixel[1])+i+fpixel[0]);}
    }
  }
	  for(j=0;j<nolm;j++) *(sm+j) /= exptot;

	/*-----------------fit spectrum, normalize---------------------------*/

	nwl=abs(iwl1-iwl0)/lbin+1;
	if(nwl<minpixnum) continue;
	iwl=iwl0;
	if(iwl>iwl1) iwl=iwl1;
	iwl/=lbin;

	//polynomial fit
	if(pfit){
	  plyfit(&index[iwl],&trace[iwl],nwl,ord_flat,coef_flat,xlist,ylist,
                 elist);}
	//median fit
	else{
	  nsv = (nwl>medbox) ? medbox: nwl;
	  for(i=0;i<nsv;i++) order(&i,trace+i+iwl,irnk,rnk);
	  flx=*(rnk+nsv/2);
	  if(nwl>medbox)nsv/=2;
	  for(i=0;i<nsv;i++) *(medtrc+i)=flx;
	  if(nwl>medbox){
	    //now move down list one pixel at a time
	    nxtpx=nsv-1;
	    lpix=mdbx1;
	    while(lpix<nwl-1){
	      //pitch 1st pixel
	      for(i=0;i<=medbox;i++){
          if(*(irnk+i) != 0) continue;
          j=i;
          break;}
	      for(i=j+1;i<medbox;i++){
          *(irnk+i-1)= *(irnk+i);
          *(rnk+i-1)= *(rnk+i);}
	      for(i=0;i<mdbx1;i++) *(irnk+i)=*(irnk+i)-1;
	      lpix++;
	      nxtpx++;
	      //add new pixel
	      order(&mdbx1,trace+lpix+iwl,irnk,rnk);
	      *(medtrc+nxtpx)=*(rnk+mdbx2);}
	    awl=*(rnk+mdbx2);
	    for(i=nxtpx+1;i<nwl;i++) *(medtrc+i)=awl;}
	  }

  //flatten spectral region

	ij=iwl-1;
	for(j=0;j<nwl;j++){
	  ij++;;
	  if(pfit){
	    flx=polyvalue((float)ij,coef_flat,ord_flat);}
	  else{
	    flx=*(medtrc+j);}
	  for(i=lower[ij];i<upper[ij];i++){
	    if(xdisper){
	      *(*(oarray+i)+ij) = 0.0;
              if (flx > 0) *(*(oarray+i)+ij) = *(*(sumspc+i-fpixel[1])+ij-fpixel[0])/flx;
              if(*(*(oarray+i)+ij)<llimit) *(*(oarray+i)+ij)=0.;
        if(inshuffle && inshuffle+i<ysize && inshuffle+i>=0){
          if(!shfld){
            *(*(oarray+i+inshuffle)+ij) = 0.0;
            if (flx>0) *(*(oarray+i+inshuffle)+ij) = *(*(sumspc+i+inshuffle-fpixel[1])+ij-
                                           fpixel[0])/flx;}
          else{
            *(*(oarray+i+inshuffle)+ij) = 0.0;
            if (flx>0) *(*(oarray+i+inshuffle)+ij) = *(*(sumspc+i-fpixel[1])+ij-
                                           fpixel[0])/flx;}
        }
      }
	    else{
	      *(*(oarray+ij)+i)  = 0.0;
              if (flx>0) *(*(oarray+ij)+i)  = *(*(sumspc+ij-fpixel[1])+i-fpixel[0])/flx;
             if(*(*(oarray+ij)+i)<llimit) *(*(oarray+ij)+i)=0.;
        if(nshuffle && inshuffle+i<xsize && inshuffle+i>=0){
          if(!shfld){
            *(*(oarray+ij)+i+inshuffle)  = 0.0;
            if (flx>0) *(*(oarray+ij)+i+inshuffle)  = *(*(sumspc+ij-fpixel[1])+i+inshuffle-
                                            fpixel[0])/flx;}
          else{
            *(*(oarray+ij)+i+inshuffle)  = 0.0;
            if (flx>0) *(*(oarray+ij)+i+inshuffle)  = *(*(sumspc+ij-fpixel[1])+i-
                                            fpixel[0])/flx;}
        }
	    }
	  }
  }
	  free(sm);
	  free(sumspc);
    printf("chip %d object %d\r",nchp,slitnum);
	  if(breakpoint) break;}
    }


    /*********************write out chip**************************************/

    fits_write_pix(outptr,TFLOAT,firstelem,noutlm,oimage,&status);
    if(status){
       printf("Error writing data for chip %d %d\n",nchp,status);
    }
    printf("Wrote ccd chip %d           \r",nchp);
    fflush(stdout);
    fits_close_file(outptr,&status);}

  return 0;}
