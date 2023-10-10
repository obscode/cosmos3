/*****************************************************************************\
*                                                                             *
*  ADJUST-MAP uses comparison arc frame[s] to make adjustments to a map file  *
*                                                                             *
*  VERSION  28 Mar 2018                                                       *
*                                                                             *
\*****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fitsio.h"
#include "cosmos.h"
#include "cpgplot.h"

void wcalc(int,float*,float*,int,float*);

int main(int argc,char *argv[]){

  char     CHAR,file[80],flnm[80],ifile[80],dfile[80],*DATA_DIR,line[140],
           *lin,name[40],dataline[6384],mfile[80],lnfile[80],lane[140],
           dewar[10],dfile2[80],tfile[80],tfile2[80],answer[80],ofile[80],
           sfile[80],cfile[80],c2file[80],objstr[80],objstr1[80];
  int      INT,status,nelem,i,j,xdisper,ord_disp,ord_sag,ord_tilt,ord_sagit,
           ord_slen,nslit,chip,slitnum,search_height,search_width,iswmn,f_tilt,
           iswmx,dmax0,dmax1,dmax2,di,nargo,ne,ord_dsagu,ord_ddispu,stype,
           nsw,nline,n,islmx,islmn,nsl,il,ii,jj,oldixtop,oldiytop,slit_width,
           nlv,ndline,breakpoint,ibin,ysign,narg, maxord,ord_ddisp,ybin,xbin,
           ord_dsag,nulvl,minlines,ixxtop2,oldixtop2,oldiytop2,interp,nbad,
           ixxtop,iyytop,odsp,osag,nargd,diag,nlpts,niter,ibiny,nlv_orig,nclump,
           ord_dtilt,nsl1,ixxtop0,ixxtop1,otlt,nchip,anynl,debug,chipx,chipy,
           flx,nsl3,nsl7,nsl9,frint,prtl,tedge,k,da,dad,year,month,day,ord_dtiltu,
           iii,jjj,pgdev1,pgdev2,pgdev3,wtpt,pw,h_ord_ddisp,h_ord_dsag;
  float    FLOAT,lambda0,lambda1,slit_len,slmin,slmax,x8,y8,lmin,lmax,lambda,
           xcen,ycen,slnt,ypos,xpos,x88,y88,sumax,sum,botval,botval2,sagmid,
           fac,factr[20],coef_disp[10],ynew,coef_sag[10],coef_tilt[10],
           coef_sagit[10],yp,ymid,edge,yedge,coef_slen[10],slnth,slenth,
           tilt,slmx,slmn,aa,bb,bmx,bmn,d_coef[10],t,s_coef[10],curv,span,
           coef_xdisp[10],xnew,sum2,sumax2,sumax0,residual,ystdev,xcut,ycut,
           sumx0,sum0,sum1,sumx1,sumx2,invfrnt,sumax1,sumx,sumy0,cc,dd,
           yresidual,ydemin,ydemax,ixtop,iytop,shifti,shiftj,ixtop2,iytop2,
           sumy,sumy1,ypos2,t_coef[10],demin,demax,siglimit,wid,sumy2,stdev,
           d,d2,xdf,fmax0,fmax1,fmax2,value,ww,tdemin,tdemax;
  size_t   arsize;
  long     naxes[3],firstel;
  int      *INTP,*image[9],*im,**array[9],*isx,*isy,*image2[9],*im2,
           **array2[9],*iresiduals,*iyresiduals;
  float    *FLOATP,*ssx,*ssy,*laval,*deval,*residuals,*yresiduals,*ydeval,
          *tdeval;
  double   *xlist,*ylist,**elist,DOUBLE,*DOUBLEP,pweight;
  fitsfile *fptr[9],*sptr[9];
  FILE     *mapfile,*linefile,*newfile,*posfile,*rmsfile;
  fitsdef  fitsinfo;
  dewdat    dewinfo;
  //linelist related
  float     *xshift_orig,*yshift_orig,*tshift_orig,*lval_orig,*cchp,*xshift,
            *yshift,*sline,*lval,*xval,*lamshift,*tshift,*lval_b,*xshift_b,
            *yshift_b,*weights,*lweights,*wval,*wv,*wval_orig,*tshift_b,
            *yweights,*tweights,*xweights;
  int       *throw_out;
  float     histo_y[300], histo_min, histo_max, histo_delt, hy;
  int       histo_n[300], nhy, histo_nmax, histo;

  debug=wtpt=0;
  edge=1.0;
  nulvl=interp=0;
  firstel=1;
  if(debug)posfile=fopen("adjust-map.dat","w");//debug
  minlines=2;                         //fit-order<=n_lines+minlines
  frint=1;

  //directories
  DATA_DIR=malloc(sizeof(CHAR)*80);
  DATA_DIR=getenv("COSMOS_IMAGE_DIR");
  if(DATA_DIR==NULL) die("COSMOS_IMAGE_DIR undefined!");
  strcat(DATA_DIR,"/");

  //get parameters

  strcpy(file,"adjust-map");
  if(OpenCosParm(file)!=0) die("cannot open parameter file");
  if(ReadParm_i("SEARCH_HGT",&search_height)==1)die("Parameter file error");
  if(ReadParm_i("SEARCH_WID",&search_width)==1)die("Parameter file error");
  if(ReadParm_i("SLIT_WIDTH",&slit_width)==1)die("Parameter file error");
  iswmn=-edge-(int)slit_width/2;
  iswmx=edge+(int)slit_width/2;
  for(i=0;i<iswmx+2;i++) factr[i]=exp(-i*i/(slit_width*slit_width*0.5));
  //for(i=0;i<iswmx+2;i++) factr[i]=0;
  factr[0]=1;
  nsw=iswmx-iswmn+1;
  if(ReadParm_i("ORD_DISP",&ord_ddisp)==1) die("Parameter file error");
  if(ReadParm_i("ORD_SAG",&ord_dsag)==1) die("Parameter file error");
  if(ReadParm_i("ORD_TILT",&ord_dtilt)==1) ord_dtilt=-1;
  ord_tilt=-1;
  if(ReadParm_b("FIT_TILT",&f_tilt)==1) f_tilt=0;
  if(!f_tilt) ord_dtilt=-1;
  if(ReadParm_s("LINELIST",lnfile)==1) die("Parameter file error");
  if(ReadParm_i("EDGE",&tedge)==1) die("Parameter file error");
  if(ReadParm_i("ITERATIONS",&niter)==1) die("Parameter file error");
  if(ReadParm_r("SIGLIMIT",&siglimit)==1) die("Parameter file error");
  if(ReadParm_b("OUTLIER-WT",&pw)==1) die("Parameter file error");
  if(ReadParm_i("NCLUMP",&nclump)==1) die("Parameter file error");
  if(ReadParm_b("HISTOGRAM",&histo)==1) die("Parameter file error");
  if(ReadParm_i("H_ORD_DISP",&h_ord_ddisp)==1) h_ord_ddisp = ord_ddisp;
  if(ReadParm_i("H_ORD_SAG",&h_ord_dsag)==1) h_ord_dsag = ord_dsag;

  linefile=fopen(lnfile,"r");
  if(linefile==NULL) die("cannot open line list file");


  //Read line list, set up line-related arrays

  nline=0;
  //count lines
  while((fgets(line,140,linefile))!=NULL){
    nline++;}
  rewind(linefile);
  //setup arrays
  arsize=sizeof(FLOAT)*nline;
  xshift_orig=malloc(arsize);
  yshift_orig=malloc(arsize);
  tshift_orig=malloc(arsize);
  weights=malloc(arsize);
  xweights=malloc(arsize);
  yweights=malloc(arsize);
  tweights=malloc(arsize);
  lweights=malloc(arsize);
  lval_orig=malloc(arsize);
  wv=malloc(arsize);
  wval=malloc(arsize);
  wval_orig=malloc(arsize);
  cchp=malloc(arsize);
  xshift=malloc(arsize);
  yshift=malloc(arsize);
  sline=malloc(arsize);
  lval=malloc(arsize);
  xval=malloc(arsize);
  lamshift=malloc(arsize);
  tshift=malloc(arsize);
  lval_b=malloc(arsize);
  xshift_b=malloc(arsize);
  yshift_b=malloc(arsize);
  tshift_b=malloc(arsize);
  throw_out=malloc(sizeof(INT)*nline);
  for(i=0;i<arsize;i++) *(xweights+i)=*(yweights+i)=*(tweights+i)=*(weights+i)=
                        *(lweights+i)=1;
  //read line list
  for(i=0;i<nline;i++){
    fgets(line,140,linefile);
    ne=sscanf(line,"%f %f",&sline[i],&ww);
    if(ne==2) *(lweights+i)=ww;}

  //setup arrays
  maxord=ord_ddisp;
  maxord = (ord_dsag>maxord) ? ord_dsag: maxord;
  maxord = (ord_dtilt>maxord) ? ord_dtilt: maxord;
  maxord+=1;
  xlist=malloc(sizeof(DOUBLE)*2*maxord);
  ylist=malloc(sizeof(DOUBLE)*maxord);
  elist=malloc(sizeof(DOUBLEP)*maxord);
  for(i=0;i<maxord;i++){
    *(elist+i)=malloc(sizeof(DOUBLE)*(maxord+1));}

  //get input files
  if(argc<5){
    while(1){
      printf("Enter map file name: ");
      scanf("%s",mfile);
      strcpy(file,mfile);
      if(strstr(file,".map")==NULL) strcat(file,".map");
      mapfile=fopen(file,"r");
      if(mapfile!=NULL) break;
      printf("mapfile %s does not exist!\n",file);
      }
    printf("Enter comparison arc frame name[s] CR=>end: ");
    scanf("%s",dfile);
    fgets(line,140,stdin);
    if(strncmp(line,"\n",1)){
      sscanf(line,"%s",&dfile2);
      interp=1;}
    }
  else{
    narg=nargd=diag=nargo=0;
    for(i=1;i<argc-1;i++){
      if(!(strcmp(argv[i],"-m"))){
        narg=i+1;
        continue;}
      if(!(strcmp(argv[i],"-d"))){
        nargd=i+1;
        continue;}
      if(!(strcmp(argv[i],"-o"))){
        nargo=i+1;
        continue;}
    }
    if(!narg) die("proper invocation: adjust-map -f framename -m mapfile\n");
    strcpy(mfile,argv[narg]);
    strcpy(file,mfile);
    if(strstr(file,".map")==NULL) strcat(file,".map");
    fprintf(stderr, "%s\n", file);
    mapfile=fopen(file,"r");
    if(mapfile==NULL) die("Cannot open map file");
    narg=0;
    for(i=1;i<argc-1;i++){
      if(!(strcmp(argv[i],"-f"))){
	narg=i+1;
	break;}
      }
    strcpy(dfile,argv[narg]);
    if(argc>narg+1 && strcmp(argv[narg+1],"-m")&& strcmp(argv[narg+1],"-d") && strcmp(argv[narg+1],"-o")){
      interp=1;
      strcpy(dfile2,argv[narg+1]);
      if(!strcmp(dfile,dfile2)) die("Two data files must be different");}
    }
  if(nargd) sscanf(argv[nargd],"%d",&diag);
  if(nargo) strcpy(ofile,argv[nargo]);
  strcpy(tfile,dfile);
  if(interp) strcpy(tfile2,dfile2);

  //output file
//  if(!diag){
  if (1){
    if(interp){
      i=strlen(dfile);
      j=strlen(dfile2);
      if(!strcmp(&dfile[i-2],&dfile2[j-2]) && (!strcmp(&dfile[i-2],"_b") ||
                                             !strcmp(&dfile[i-2],"f"))){
        dfile[i-2]=dfile[i];
        i-=2;}
      if(strlen(dfile2)<i)i=strlen(dfile2);
      for(j=0;j<i;j++){
        if(dfile[j] != dfile2[j]) break;}
      strcat(dfile,"-");
      strcat(dfile,dfile2+j);
      }
    if(!strcmp(mfile,dfile)){
      printf("spectral map %s.map moved to %s.map.old\n",mfile,mfile);
      strcat(mfile,".map");
      fclose(mapfile);
      strcpy(ifile,"mv ");
      strcat(ifile,mfile);
      strcat(ifile," ");
      strcat(mfile,".old");
      strcat(ifile,mfile);
      system(ifile);
      mapfile=fopen(mfile,"r");}
    lin=strrchr(dfile,'/');
    if(lin==NULL){
      lin=dfile;}
    else{
      lin+=1;}
    strcpy(file,lin);
    if (!diag) {
      if(!nargo) strcpy(dfile,file);
      else strcpy(dfile,ofile);
      strcat(dfile,".map");
      if((newfile=fopen(dfile,"w"))==NULL){
        printf("error opening output file %s\n",dfile);
        return 1;}
    }
    //  strcpy(file,argv[3]);
    strcpy(sfile,file);
    strcat(sfile,".rms");
    strcpy(cfile,file);
    strcat(cfile,".ps/cps");
    strcpy(c2file,file);
    strcat(c2file,"y.ps/cps");
    if((rmsfile=fopen(sfile,"w"))==NULL){
      printf("error opening output file %s\n",sfile);
      return 1;}
    fprintf(rmsfile,"#%4s %15s %10s\n", "N", "ID", "STDEV");
    }

  //Read mapping data

  fgets(line,140,mapfile);
  if(!sscanf(line,"Xdispersion = %d",&xdisper)) die("Mapfile error");
  if(!diag) fprintf(newfile,"%s",line);
  fgets(line,140,mapfile);
  if(!sscanf(line,"Fit orders = %d %d %d %d %d",&ord_disp,&ord_sag,&ord_tilt,
             &ord_sagit,&ord_slen)) die("Mapfile error");
  if(ord_disp>9 || ord_sag>9 || ord_tilt>9 || ord_sagit>9 || ord_slen>9){
    printf("Warning: maximum fit order = 9\n");
    return 1;}
  if(!diag) fprintf(newfile,"%s",line);
  fgets(line,140,mapfile);
  if(!diag) fprintf(newfile,"%s",line);
  fgets(line,140,mapfile);
  if(!sscanf(line,"Lambda = %f %f",&lambda0,&lambda1)) die("Mapfile error");
  if(!diag) fprintf(newfile,"%s",line);
  fgets(line,140,mapfile);
  if(!sscanf(line,"Dewar = %s %d",dewar,&nchip)){
    strcpy(dewar,"SITE");                         //default for old map files
    if(!diag) fprintf(newfile,"Dewar = SITE 8\n");
    nchip=8;}
  else{
    if(!diag) fprintf(newfile,"%s",line);}
  if(Readcdf(dewar)==1){
    printf("Error reading dewar definition file %s\n",dewar);
    return 1;}
  Getchipdat(&dewinfo);
  chipx=dewinfo.xchip;
  chipy=dewinfo.ychip;
  //open files
  strcpy(flnm,DATA_DIR);
  strcat(flnm,tfile);
  strcpy(ifile,flnm);
  for(i=1;i<=nchip;i++){
    strcpy(file,ifile);
    if(nchip==1){
      strcat(file,".fits");}
    else{
      addbar(file);
      strcat(file,"c");
      sprintf(line,"%d.fits",i);
      strcat(file,line);}
    status=OpenFitsFile(file,&fptr[i],&fitsinfo);
    if(status) fits_die(file,status);}
  //if 2 image files...
  if(interp){
    strcpy(flnm,DATA_DIR);
    strcat(flnm,tfile2);
    strcpy(ifile,flnm);
    for(i=1;i<=nchip;i++){
      strcpy(file,ifile);
      if(nchip==1){
        strcat(file,".fits");}
      else{
        addbar(file);
        strcat(file,"c");
        sprintf(line,"%d.fits",i);
        strcat(file,line);}
      status=(OpenFitsFile(file,&sptr[i],&fitsinfo));
      if(status) fits_die(file,status);}
    }


  //Check for correct SITE dewar
  sscanf(fitsinfo.date,"%d-%d-%d",&year,&month,&day);
  if(year>=2005 && month>=8){
    if(strstr(dewar,"SITE") != NULL && strstr(dewar,"SITE2")==NULL){
printf("\n***************************************************************\nWARNING: for IMACS data after Aug 1, 2005, you must specify the\n instrument as \"IMACS2\" in defineobs, or the dewar as \"SITE2\"\n in the dewar definition file. See documentation for details\n****************************************************************\n\n");
    return 1;}
  }
  //Read ccd data

  naxes[0]=fitsinfo.naxes[0];
  naxes[1]=fitsinfo.naxes[1];
  ibin=ibiny=fitsinfo.binning;
  if(fitsinfo.ybinning) ibiny=fitsinfo.ybinning;
  xbin=ibin;
  ybin=ibiny;
  if(!xdisper){
    xbin=ibiny;
    ybin=ibin;}
  nelem=naxes[0]*naxes[1];
  status=0;
  for(i=1;i<=nchip;i++){
    image[i]=malloc(sizeof(INT)*naxes[0]*naxes[1]);
    im=image[i];
    array[i]=malloc(sizeof(INTP)*naxes[1]);
    for(j=0;j<naxes[1];j++) *(array[i]+j)=im+naxes[0]*j;
    status=ReadFitsFile(&fptr[i],&fitsinfo,TINT,im);
    if(status) fits_die("Error readinng image files",status);
    if(interp){
      image2[i]=malloc(sizeof(INT)*naxes[0]*naxes[1]);
      im2=image2[i];
      array2[i]=malloc(sizeof(INTP)*naxes[1]);
      for(j=0;j<naxes[1];j++) *(array2[i]+j)=im2+naxes[0]*j;
      status=ReadFitsFile(&sptr[i],&fitsinfo,TINT,im2);
      if(status) fits_die("Error readinng image files",status);}
    printf("Read ccd chip %d\r",i);
    fflush(stdout);}
  printf("\n");

    pgdev2=cpgopen(cfile);
    cpgslct(pgdev2);
    if(f_tilt){cpgsubp(3,3);}
    else{cpgsubp(4,3);}
    cpgpap(10.,0.85);

    nlpts=(lambda1-lambda0+1);
    laval=malloc(sizeof(FLOAT)*nlpts);
    deval=malloc(sizeof(FLOAT)*nlpts);
    ydeval=malloc(sizeof(FLOAT)*nlpts);
    tdeval=malloc(sizeof(FLOAT)*nlpts);
    for(i=0;i<nlpts;i++) *(laval+i)=lambda0+i;

/*    pgdev3=cpgopen(c2file);
    cpgslct(pgdev3);
    cpgsubp(3,3);
    cpgpap(10.,0.85);*/

  //if diag output
    if(diag){
    pgdev1=cpgopen("/xwindow");
    cpgslct(pgdev1);
    wid=8.;
    cpgsubp(1,2);
    cpgpap(wid,1.);
    }

    //pgplot postscript

   while(1){
    if(fgets(line,140,mapfile)==NULL)  die("Unexpected end to mapfile a");
    if(sscanf(line,"SLIT %d %s %d",&slitnum, name,&stype)) break;
    if(!diag) fprintf(newfile,"%s",line);}

  //New slit
  while(1){
    if(diag==0) fprintf(newfile,"%s",line);
    if(!strncmp(line,"END",3)) break;                             //end of data
    if(!sscanf(line,"SLIT %d %s %d",&slitnum, name,&stype)) die(line);
    if((fgets(line,140,mapfile))==NULL) die("Unexpected end to mapfile b");
    if(diag==0) fprintf(newfile,"%s",line);
    if(!strncmp(line,"END",3)) break;                             //end of data
    if(!sscanf(line,"LENGTH = %f",&slit_len)) continue;       //no data in slit

    printf("Analyzing slit  %d\r",slitnum);
    fflush(stdout);

    if(histo){
      histo_min = 999;
      histo_max = -999;}

    if (stype == 1) {
      ord_ddispu = h_ord_ddisp;
      ord_dtiltu = -1;
      ord_dsagu  = h_ord_dsag;
    } else {
      ord_ddispu = ord_ddisp;
      ord_dtiltu = ord_dtilt;
      ord_dsagu  = ord_dsag;
    }

    //found a slit with data
    nslit++;
    nlv=0;
    ndline=0;
    //process all chips for this slit
    while(2){
      if((fgets(line,140,mapfile))==NULL) die("Unexpected end to mapfile c");
      if(!sscanf(line,"CHIP %d %f %f %f %f %d %d",&chip,&slmin,&slmax,&lmin,
                 &lmax,&i,&prtl)) break;
      strcpy(&dataline[ndline*140],line);
      ndline++;
      //coef_disp data
      if((fgets(line,140,mapfile))==NULL)die("Unexpected end to mapfile d");
      strcpy(&dataline[ndline*140],line);
      ndline++;
      lin=&line[0];
      for(i=0;i<=ord_disp;i++){
        if((lin=strpbrk(lin,"1234567890.-"))==NULL)
                                            die("Unexpected end to mapfile e");
        sscanf(lin,"%f",&coef_disp[i]);
        if((lin=strpbrk(lin," "))==NULL) die("Unexpected end to mapfile f");}

      //coef_xdisp dat
      if((fgets(line,140,mapfile))==NULL)die("Unexpected end to mapfile g");
      strcpy(&dataline[ndline*140],line);
      ndline++;
      lin=&line[0];
      for(i=0;i<=ord_disp;i++){
        if((lin=strpbrk(lin,"1234567890.-"))==NULL)
                                            die("Unexpected end to mapfile h");
        sscanf(lin,"%f",&coef_xdisp[i]);
        if((lin=strpbrk(lin," "))==NULL)die("Unexpected end to mapfile i");}
      //coef_sag dat
      if((fgets(line,140,mapfile))==NULL)die("Unexpected end to mapfile j");
      strcpy(&dataline[ndline*140],line);
      ndline++;
      lin=&line[0];
      for(i=0;i<=ord_sag;i++){
        if((lin=strpbrk(lin,"1234567890.-"))==NULL)
                                          die("Unexpected end to mapfile k");
        sscanf(lin,"%f",&coef_sag[i]);
        if((lin=strpbrk(lin," "))==NULL)die("Unexpected end to mapfile l");}
      //coef_tilt data
      if((fgets(line,140,mapfile))==NULL)die("Unexpected end to mapfile m");
      strcpy(&dataline[ndline*140],line);
      ndline++;
      lin=&line[0];
      for(i=0;i<=ord_tilt;i++){
        if((lin=strpbrk(lin,"1234567890.-"))==NULL)
                                            die("Unexpected end to mapfile n");
        sscanf(lin,"%f",&coef_tilt[i]);
        if((lin=strpbrk(lin," "))==NULL)die("Unexpected end to mapfile o ");}
      //coef_sagit data
      if((fgets(line,140,mapfile))==NULL)die("Unexpected end to mapfile p");
      strcpy(&dataline[ndline*140],line);
      ndline++;
      lin=&line[0];
      for(i=0;i<=ord_sagit;i++){
        if((lin=strpbrk(lin,"1234567890.-"))==NULL)
          die("Unexpected end to mapfile");
        sscanf(lin,"%f",&coef_sagit[i]);
        if((lin=strpbrk(lin," "))==NULL)die("Unexpected end to mapfile q");}
      //coef_slen data
      if((fgets(line,140,mapfile))==NULL)die("Unexpected end to mapfile r");
      strcpy(&dataline[ndline*140],line);
      ndline++;
      lin=&line[0];
      for(i=0;i<=ord_slen;i++){
        if((lin=strpbrk(lin,"1234567890.-"))==NULL)
          die("Unexpected end to mapfile");
        sscanf(lin,"%f",&coef_slen[i]);
        if((lin=strpbrk(lin," "))==NULL)die("Unexpected end to mapfile s");}


  /*---------------------------Tweak Map-------------------------------*/

      breakpoint=0;
      nslit=-1;

      //center each line

      slenth=slmax-slmin;
      for(n=0;n<nline;n++){

	//extraction box centered on line

	lambda=sline[n];
	if(lambda<lmin || lambda>lmax) continue;
	xcen=polyvalue(lambda,coef_disp,ord_disp)/xbin;;
	ycen=polyvalue(lambda,coef_sag,ord_sag)/ybin;
	slnt=polyvalue(lambda,coef_slen,ord_slen)/ybin;
	tilt=polyvalue(lambda,coef_tilt,ord_tilt)/xbin;
  curv=polyvalue(lambda,coef_sagit,ord_sagit)/xbin;
	slmx=slmax*slnt/slenth;
	slmn=slmin*slnt/slenth;
	ymid=ycen+(slmx+slmn)/2.;
	ysign=1;
	if(slmx<slmn){
    tilt=-tilt;
	  t=slmx;
	  slmx=slmn;
	  slmn=t;
	  ysign=-1;}
	islmx=(int)slmx+edge;
	islmn=(int)slmn-edge;
	slnth=slmx-slmn;
	nsl=islmx-islmn+1;
  nsl1=0.5*nsl;
  nsl3=0.3*nsl;
  nsl7=0.7*nsl;
  nsl9=0.9*nsl;
	//is search box completely on chip?
	aa=tilt*slmx/(slmx-slmn);
	bb=tilt*slmn/(slmx-slmn);
	bmx = (tilt>0) ? aa+search_width+xcen: bb+search_width+xcen;
	bmn= (tilt>0) ? bb-search_width+xcen: aa-search_width+xcen;
	if(xdisper){
	  if(bmx>chipx || bmn<0 || ycen+slmx+search_height> chipy  ||
	     ycen+slmn-search_height< 0 ) continue;}
	else{
	  if(bmx>chipy || bmn<0 || ycen+slmx+search_height> chipx  ||
	     ycen+slmn-search_height< 0 ) continue;}
	lval[nlv]=lambda;
  wval[nlv]=lweights[n];
	xval[nlv]=xcen;
	//interpolation arrays
	isx=malloc(sizeof(INT)*nsl);
	isy=malloc(sizeof(INT)*nsl);
	ssx=malloc(sizeof(FLOAT)*nsl);
	ssy=malloc(sizeof(FLOAT)*nsl);
	for(i=0;i<nsl;i++){
	  il=i+islmn;
    sagmid=curv*4.*(powf(il/slnt,2)-powf((ymid-ycen)/slnt,2));
	  aa=xcen+tilt*il/(slmx-slmn)-sagmid;
//   aa=xcen;
	  *(isx+i)=(int)aa;
	  *(ssx+i)=aa-*(isx+i);
	  aa=ycen+il;
	  *(isy+i)=(int)aa;
	  *(ssy+i)=aa-*(isy+i);}

	//peak up on line

	sumax=sumax2=0;
	breakpoint=0;
  botval=botval2=999999.;
	for(shifti=-search_width;shifti<=search_width;shifti+=1.0){
	  for(shiftj=-search_height;shiftj<=search_height;shiftj+=0.25){
	    sum=sum2=0;
	    //sum over extraction box
	    for(jj=0;jj<nsl;jj++){
	      for(ii=iswmn;ii<=iswmx;ii++){
		xpos=*(isx+jj)+shifti+ii+ *(ssx+jj);
		ypos=*(isy+jj)+shiftj+ *(ssy+jj);
		if(!xdisper){
		  yp=xpos;
		  xpos=ypos;
		  ypos=yp;}
		if(xpos<0 || xpos>=naxes[0] || ypos<0 || ypos>=naxes[1]){
		  breakpoint=1;
		  break;}
		flx=i_interpol(array[chip],naxes,xpos,ypos);
    sum+=flx;
    if(flx<botval) botval=flx;
		if(interp){
      flx=i_interpol(array2[chip],naxes,xpos,ypos);
      sum2+=flx;
      if(flx<botval2) botval2=flx;}
        }
      }
	    if(sum>sumax){
	      sumax=sum;
	      ixtop=shifti;
	      iytop=shiftj;}
	    if(interp && sum2>sumax2){
	      sumax2=sum2;
	      ixtop2=shifti;
	      iytop2=shiftj;}
	    }
	  if(breakpoint) break;}
	if(breakpoint) continue;


  // find center of gravity of line

  sum=sum0=sumx0=sum1=sumy=sumy0=sumy1=0;
  //sum over extraction box
  for(jj=tedge;jj<nsl-tedge;jj++){
    for(ii=iswmn-1;ii<=iswmx+1;ii++){
      xpos=ixtop+ii+*(isx+jj)+ *(ssx+jj);
      ypos=iytop+*(isy+jj)+ *(ssy+jj);
      if(!xdisper){
        yp=xpos;
        xpos=ypos;
        ypos=yp;}
      flx=i_interpol(array[chip],naxes,xpos,ypos)-botval;
      if(jj<=nsl1){
        sum0+=flx;
        sumy0+=flx*jj;}
      else{
        sum1+=flx;
        sumy1+=flx*jj;}
      if(interp){
        xpos=ixtop2+ii+*(isx+jj)+ *(ssx+jj);
        ypos=iytop2+*(isy+jj)+ *(ssy+jj);
        if(!xdisper){
          yp=xpos;
          xpos=ypos;
          ypos=yp;}
        flx=i_interpol(array[chip],naxes,xpos,ypos)-botval2;
        sum2+=flx;
        sumy2+=flx*jj;}
    }
  }
  sum=sum0+sum1;
  sumy=sumy0+sumy1;
  xdf=sumx/sum;
  cchp[nlv]=chip;
  yshift[nlv]=ybin*(iytop);
  if(!tedge) yshift[nlv]+=sumy/sum-nsl/2.;
  if(interp) yshift[nlv]=0.5*(yshift[nlv]+ybin*iytop2);

  //convolve in wavelength direction to find peak

  fmax0=fmax1=fmax2=0;
  for(di=-20;di<=20;di++){
    sumx0=sumx1=sum2=0;
    for(jj=tedge;jj<nsl-tedge;jj++){
      for(ii=0;ii<=iswmx+2;ii++){
        xpos=*(isx+jj)+*(ssx+jj)+ixtop+di/10.+ii;
        ypos=iytop+*(isy+jj)+*(ssy+jj);
        if(!xdisper){
          yp=xpos;
          xpos=ypos;
          ypos=yp;}
        flx=i_interpol(array[chip],naxes,xpos,ypos);
        if(jj<=nsl1){
          sumx0+=flx*factr[ii];}
        else{
          sumx1+=flx*factr[ii];}
        xpos=*(isx+jj)+*(ssx+jj)+ixtop+di/10.-ii;
        ypos=iytop+*(isy+jj)+*(ssy+jj);
        if(!xdisper){
          yp=xpos;
          xpos=ypos;
          ypos=yp;}
        flx=i_interpol(array[chip],naxes,xpos,ypos);
        if(jj<=nsl1){
          sumx0+=flx*factr[ii];}
        else{
          sumx1+=flx*factr[ii];}
        if(!interp) continue;
        //interpolate between 2 frames
        xpos=ixtop2+ii+di/10.+*(isx+jj)+ *(ssx+jj);
        ypos=iytop2+*(isy+jj)+ *(ssy+jj);
        if(!xdisper){
          yp=xpos;
          xpos=ypos;
          ypos=yp;}
        flx=i_interpol(array[chip],naxes,xpos,ypos);
        sum2+=flx*factr[ii];
        xpos=ixtop2-ii+di/10.+*(isx+jj)+ *(ssx+jj);
        ypos=iytop2+*(isy+jj)+ *(ssy+jj);
        if(!xdisper){
          yp=xpos;
          xpos=ypos;
          ypos=yp;}
        flx=i_interpol(array[chip],naxes,xpos,ypos);
        sum2+=flx*factr[ii];
          }
      }
    if(sumx0>fmax0){
      fmax0=sumx0;
      dmax0=di;}
    if(sumx1>fmax1){
      fmax1=sumx1;
      dmax1=di;}
    if(sum2>fmax2){
      fmax2=sum2;
      dmax2=di;}
    }

  tshift[nlv]=-0.1*(dmax1-dmax0);
  tshift[nlv]*=nsl/((sumy1/sum1)-(sumy0/sum0));
  xshift[nlv]=xbin*(ixtop+0.05*(dmax0+dmax1));
  if(interp) xshift[nlv]=0.5*(xshift[nlv]+xbin*(ixtop2+0.1*dmax2));
  if(chip>4){                                             //SITE specific!!
    if(!xdisper){
	    xshift[nlv]=-xshift[nlv];
	    yshift[nlv]=-yshift[nlv];}
    }
  if(xdisper){
    xnew=xbin*xcen;
    ynew=ybin*ymid;}
  else{
    xnew=xbin*ycen;
    ynew=ybin*xcen;}
	mspos(chip,xnew,ynew,&x8,&y8);
	if(debug)fprintf(posfile,"%d %s %8.2f %7.2f %7.2f %d %7.2f %7.2f %6.2f%6.2f %d %d %8.4f\n",slitnum,name,lambda,x8,y8,chip,xnew,ynew,xshift[nlv],
                   yshift[nlv],dmax0,dmax1,tilt);
  free(isx);
	free(isy);
	free(ssx);
	free(ssy);
  if(histo){
    if (yshift[nlv] < histo_min) histo_min = yshift[nlv];
    if (yshift[nlv] > histo_max) histo_max = yshift[nlv];}
	nlv++;}//end of line loop
        }//end of chip loop

    if(diag!=0){
      da=slitnum/diag;
      dad=da*diag;}
      if(histo){
        histo_delt = 0.25*ybin;
        for(i=0;i<(int) ((histo_max-histo_min)/histo_delt+1);i++) {
           histo_y[i] = histo_min + i*histo_delt;
           histo_n[i] = 0;
        }
        for(i=0;i<nlv;i++) {
           nhy = (int) ((yshift[i]-histo_min)/histo_delt);
           histo_n[nhy] += 1;}
        }

    //fit position corrections
    for(i=0;i<=ord_ddispu;i++) d_coef[i]=0.;
    for(i=0;i<=ord_dsagu;i++) s_coef[i]=0;
    for(i=0;i<=ord_dtiltu;i++) t_coef[i]=0;
    odsp=otlt=osag=-1;
    stdev=0;
    ystdev=0;
    if(nlv>=minlines){
      nlv_orig = nlv;
      for(i=0;i<nlv;i++) {
        xshift_orig[i] = xshift[i];
        yshift_orig[i] = yshift[i];
        tshift_orig[i] = tshift[i];
        lval_orig[i] = lval[i];
        wval_orig[i]=wval[i];}
      for(j=0;j<niter;j++){
        nlv = 0;
        nbad = 0;
        if(histo){
          for(i=0;i<nlv_orig;i++) {
            nhy = (int) ((yshift_orig[i]-histo_min)/histo_delt);
            if (histo_n[nhy] > 2 &&
               (yshift_orig[i] < -ybin*(search_height-0.25) ||
               yshift_orig[i] > ybin*(search_height-0.25))) {
                    throw_out[i] = 1;
                    nbad++;
            }
          }
        }

        for(i=0;i<nlv_orig;i++) {
          if (!(throw_out[i])) {
            xshift[nlv] = xshift_orig[i];
            yshift[nlv] = yshift_orig[i];
            tshift[nlv] = tshift_orig[i];
            lval[nlv] = lval_orig[i];
            wval[nlv]=wval_orig[i];
            nlv++;}
          else{
            nbad++;}

        }
        if(ord_ddispu>=0){
          odsp = (nlv>=minlines+ord_ddispu) ? ord_ddispu : nlv-minlines;
          odsp = (j==0 && odsp>1 && niter>1) ? odsp-1 : odsp;
          if(pw) wcalc(nclump,lval,xshift,nlv,xweights);
          for(k=0;k<nlv;k++) wv[k]=wval[k]*xweights[k];
          plyfit_w(lval,xshift,wv,nlv,odsp,d_coef,xlist,ylist,elist);}
        if(ord_dsagu>=0){
          osag = (nlv>=minlines+ord_dsagu) ? ord_dsagu : nlv-minlines;
          osag = (j==0 && osag>1) ? osag-1 : osag;
          if(pw) wcalc(nclump,lval,yshift,nlv,yweights);
          for(k=0;k<nlv;k++) wv[k]*=wval[k]*yweights[k];
          plyfit_w(lval,yshift,wv,nlv,osag,s_coef,xlist,ylist,elist);}
        if(ord_dtiltu>=0){
          otlt = (nlv>=minlines+ord_dtiltu) ? ord_dtiltu : nlv-minlines;
          if(pw) wcalc(nclump,lval,tshift,nlv,tweights);
          for(k=0;k<nlv;k++) wv[k]*=wval[k]*tweights[k];
          plyfit_w(lval,tshift,wv,nlv,otlt,t_coef,xlist,ylist,elist);}
        if(niter>1 && j<niter-1){
          if (nlv > 5) {
            residuals = malloc((nlv)*sizeof(FLOAT));
            iresiduals = malloc((nlv)*sizeof(INT));
            yresiduals = malloc((nlv)*sizeof(FLOAT));
            iyresiduals = malloc((nlv)*sizeof(INT));
            }
          for(i=0;i<nlv;i++) {
            residual = fabs( xshift[i]-polyvalue(lval[i],d_coef,odsp));
            yresidual = fabs( yshift[i]-polyvalue(lval[i],s_coef,osag));
            if(isnan(residual) || isnan(yresidual)){
              printf("numerical error %d\n",slitnum);
              continue;}
            stdev += pow(residual,2.0);
            ystdev += pow(yresidual,2.0);
            if (nlv>5) {
              order(&i,&residual,iresiduals,residuals);
              order(&i,&yresidual,iyresiduals,yresiduals);
            }
          }
          stdev=sqrt(stdev/nlv);
          ystdev=sqrt(ystdev/nlv);
          if (nlv>5) {
             if(dad==slitnum) printf("%5d orders= %d %d  stdevs = %9.3f %9.3f ",slitnum,odsp,osag,stdev,ystdev);
             stdev=1.49*residuals[(int) (nlv/2)+1];
             ystdev=1.49*yresiduals[(int) (nlv/2)+1];
             if(dad==slitnum) printf("%9.3f %9.3f\n",stdev,ystdev);
             free(residuals);
             free(iresiduals);
             free(yresiduals);
             free(iyresiduals);
             }
             else {
             if(dad==slitnum) printf("stdevs = %9.3f %9.3f\n",stdev,ystdev);
          }
          i=0;
          nbad=0;
          xcut = (stdev > 0.2) ? siglimit*stdev : siglimit*0.2;
          ycut = (ystdev > 0.2) ? siglimit*ystdev : siglimit*0.2;
          while(i<nlv_orig){
            d=fabs(xshift_orig[i]-polyvalue(lval_orig[i],d_coef,odsp));
            d2=fabs(yshift_orig[i]-polyvalue(lval_orig[i],s_coef,osag));
//            printf("%f %f %f %f\n",lval[i],d,stdev,xshift[i]);
            throw_out[i] = 0;
            if(d>xcut || d2 > ycut){
              xshift_b[nbad]=xshift_orig[i];
              yshift_b[nbad]=yshift_orig[i];
              tshift_b[nbad]=tshift_orig[i];
              lval_b[nbad]=lval_orig[i];
              throw_out[i] = 1;
              nbad++;
//              for(k=i+1;k<nlv;k++){
//                xshift[i]=xshift[i+1];
//                lval[i]=lval[i+1];}
//              nlv--;
            }
            i++;}
        }
      }
    }
    else{
        printf("Insuffient data to correct slit %d\n",slitnum);
      fflush(stdout);}

    //plot to postscript file 1
    cpgslct(pgdev2);
    cpgpage;
    cpgsch(2);
//    cpgscr(0,1,1,1);
//    cpgscr(1,0,0,0);
    cpgask(0);
    cpgsci(1);
    demin=10000;
    demax=-10000;
    for(iii=0;iii<nlv;iii++){
      if(xshift[iii]<demin) demin=xshift[iii];
      if(xshift[iii]>demax) demax=xshift[iii];}
    for(iii=0;iii<nbad;iii++){
      if(xshift_b[iii]<demin) demin=xshift_b[iii];
      if(xshift_b[iii]>demax) demax=xshift_b[iii];}

    for(iii=0;iii<nlpts;iii++){
      *(deval+iii)=polyvalue(lambda0+iii,d_coef,odsp);
//      if(*(deval+iii)<demin) demin=*(deval+iii);
//      if(*(deval+iii)>demax) demax=*(deval+iii);
    }
    demin-=1;
    demax+=1;
    cpgenv(lambda0,lambda1,demin,demax,0,0);
    //cpglab("Lambda","D Lambda",name);
    cpglab("\\gl (\\A)","\\gD\\gl (pixels)",name);
    if (stype == 1) sprintf(objstr,"%i - HOLE",slitnum);
    if (stype == 2) sprintf(objstr,"%i - SLIT",slitnum);
    cpgsch(1.75);
    sprintf(objstr1,"%-s","DISP");
    cpgmtxt("T",0.5,0,0.0,objstr);
    cpgmtxt("T",0.5,0.85,0.0,objstr1);
    if (pw) {
      cpgsci(4);
      for(k=0;k<nlv;k++) {
       cpgsch(4*sqrt(0.25*log10f(wval[k]+2)*xweights[k]));
       cpgpt1(lval[k],xshift[k],17);
      }
    } else {
      cpgsch(10);
      cpgsci(4);
      cpgpt(nlv,lval,xshift,1);}
    cpgsci(2);
    cpgpt(nbad,lval_b,xshift_b,1);
    cpgsci(1);
    cpgline(nlpts,laval,deval);
    cpgsci(1);
    cpgupdt();

    //plot to postscript file 2
    cpgpage;
    cpgsch(2);
//    cpgscr(0,1,1,1);
//    cpgscr(1,0,0,0);
    cpgask(0);
    cpgsci(1);
    ydemin=10000;
    ydemax=-10000;
    for(iii=0;iii<nlv;iii++){
//      printf("%d %7.1f %3.1f\n",cchp[iii],lval[iii],yshift[iii]);
      if(yshift[iii]<ydemin) ydemin=yshift[iii];
      if(yshift[iii]>ydemax) ydemax=yshift[iii];}
    for(iii=0;iii<nbad;iii++){
      if(yshift_b[iii]<ydemin) ydemin=yshift_b[iii];
      if(yshift_b[iii]>ydemax) ydemax=yshift_b[iii];}
    for(iii=0;iii<nlpts;iii++){
      *(ydeval+iii)=polyvalue(lambda0+iii,s_coef,osag);
 //     if(*(ydeval+iii)<ydemin) ydemin=*(ydeval+iii);
 //     if(*(ydeval+iii)>ydemax) ydemax=*(ydeval+iii);
    }
    ydemin-=1;
    ydemax+=1;
    cpgenv(lambda0,lambda1,ydemin,ydemax,0,0);
    //cpglab("Lambda","D Lambda",name);
    cpglab("\\gl (\\A)","\\gDY (pixels)",name);
    cpgsch(1.75);
    sprintf(objstr1,"%-s","SAG");
    cpgmtxt("T",0.5,0,0.0,objstr);
    cpgmtxt("T",0.5,0.85,0.0,objstr1);
    if (pw) {
      cpgsci(4);
      for(k=0;k<nlv;k++) {
       cpgsch(4*sqrt(0.25*log10f(wval[k]+2)*yweights[k]));
       cpgpt1(lval[k],yshift[k],17);
      }
    } else {
      cpgsch(10);
      cpgsci(4);
      cpgpt(nlv,lval,xshift,1);}
    cpgsci(2);
    cpgpt(nbad,lval_b,yshift_b,1);
    cpgsci(1);
    cpgline(nlpts,laval,ydeval);
    cpgsci(1);
    cpgupdt();

        //plot postscript 3
        cpgpage;
        cpgsch(2);
        cpgask(0);
        cpgsci(1);
        tdemin=10000;
        tdemax=-10000;
        for(iii=0;iii<nlv;iii++){
          if(tshift[iii]<tdemin) tdemin=tshift[iii];
          if(tshift[iii]>tdemax) tdemax=tshift[iii];}
        for(iii=0;iii<nlpts;iii++){
          *(tdeval+iii)=polyvalue(lambda0+iii,t_coef,otlt);}
        span = tdemax-tdemin+1;
        tdemin-=0.1*span;
        tdemax+=0.1*span;
        cpgenv(lambda0,lambda1,tdemin,tdemax,0,0);
        cpglab("\\gl (\\A)","\\gD\\gl (pixels)",name);
        cpgsch(1.75);
        sprintf(objstr1,"%-s","TILT");
        cpgmtxt("T",0.5,0,0.0,objstr);
        cpgmtxt("T",0.5,0.85,0.0,objstr1);
        cpgsch(2);
        if (pw) {
          cpgsci(4);
          for(k=0;k<nlv;k++) {
           cpgsch(4*sqrt(0.25*log10f(wval[k]+2)*tweights[k]));
           cpgpt1(lval[k],tshift[k],17);
          }
        } else {
          cpgsch(10);
          cpgsci(4);
          cpgpt(nlv,lval,tshift,1);
        }
        cpgsch(10);
        cpgsci(2);
        cpgpt(nbad,lval_b,tshift_b,1);
        for(k=0;k<nlv_orig;k++) {
         if (throw_out[k]) {
         cpgsch(2.0);
         cpgpt1(lval_orig[k],tshift_orig[k],17);
         }
        }
        cpgsci(1);
        cpgline(nlpts,laval,tdeval);
        cpgsci(1);
        cpgupdt();



    //plot diagnostics
    if(diag > 0 && dad==slitnum){
      if (niter>1) printf("stdev = %f\n",stdev/siglimit);
      cpgslct(pgdev1);

      cpgpage;
      cpgsch(2);
      cpgscr(0,1,1,1);
      cpgscr(1,0,0,0);
      cpgask(0);
      cpgsci(1);
      cpgenv(lambda0,lambda1,demin,demax,0,0);
      //cpglab("Lambda","D Lambda",name);
      cpglab("\\gl (\\A)","\\gD\\gl (pixels)",name);
      cpgsch(15);
      cpgsci(4);
      cpgpt(nlv,lval,xshift,1);
      cpgsci(2);
      cpgpt(nbad,lval_b,xshift_b,1);
      cpgsci(1);
      cpgline(nlpts,laval,deval);
      cpgsci(1);

      cpgpage;
      cpgsch(2);
//      cpgscr(0,1,1,1);
//      cpgscr(1,0,0,0);
      cpgask(0);
      cpgsci(1);
      cpgenv(lambda0,lambda1,ydemin,ydemax,0,0);
      //cpglab("Lambda","D Lambda",name);
      cpglab("\\gl (\\A)","\\gDY (pixels)",name);
      cpgsch(15);
      cpgsci(4);
      cpgpt(nlv,lval,yshift,1);
      cpgsci(2);
      cpgpt(nbad,lval_b,yshift_b,1);
      cpgsci(1);
      cpgline(nlpts,laval,ydeval);
      cpgsci(1);


      printf("Hit CR to continue, Q to quit, G to go: ");
      fflush(stdout);
      fgets(answer,80,stdin);
      if(!strncasecmp(answer,"q",1)) return 0;
      if(!strncasecmp(answer,"g",1)) diag=-1;

      }
    //or apply corrections


    if(diag==0){
      i=0;
      printf("Correcting slit %d\r",slitnum);
      fflush(stdout);
      while(i<ndline){
      //chip data
      strncpy(lane,&dataline[i*140],140);
      i++;
      sscanf(lane,"CHIP %d",&chip);
      fprintf(newfile,"%s",lane);
      //x=f(lambda)
      strncpy(lane,&dataline[i*140],140);
      i++;
      lin=&lane[0];
      for(ii=0;ii<=ord_disp;ii++){
        if((lin=strpbrk(lin,"1234567890.-"))==NULL){
          printf("Unexpected end to mapfile slit %d\n line = %s",slitnum,lin);
          return 1;}
        sscanf(lin,"%f",&coef_disp[ii]);
        if((lin=strpbrk(lin," "))==NULL){
          printf("Unexpected end to mapfile slit %d d\n",slitnum);
          return 1;}}
      for(j=0;j<=ord_ddispu;j++){
        if(isnan(d_coef[0])) break;
        if(xdisper || chip<=4){                           //SITE specific!!
          coef_disp[j]+=d_coef[j];}
        else{
          coef_disp[j]-=d_coef[j];}
	    }
      for(ii=0;ii<=ord_disp;ii++) fprintf(newfile,"%13.6e ",coef_disp[ii]);
          fprintf(newfile,"\n");
      //lambda=f(x)
      strncpy(lane,&dataline[i*140],140);
      i++;
      lin=&lane[0];
      for(ii=0;ii<=ord_disp;ii++){
        if((lin=strpbrk(lin,"1234567890.-"))==NULL){
          printf("Unexpected end to mapfile slit %d c\n",slitnum);
          return 1;}
        sscanf(lin,"%f",&coef_xdisp[ii]);
        if((lin=strpbrk(lin," "))==NULL){
          printf("Unexpected end to mapfile slit %d d\n",slitnum);
          return 1;}}
      //	for(j=0;j<=odsp;j++){
      //	  coef_xdisp[j]+=x_coef[j];}
      for(ii=0;ii<=ord_disp;ii++) fprintf(newfile,"%13.6e ",coef_xdisp[ii]);
      fprintf(newfile,"\n");
      //y=f(lambda)
      strncpy(lane,&dataline[i*140],140);
      i++;
      lin=&lane[0];
      for(ii=0;ii<=ord_sag;ii++){
        if((lin=strpbrk(lin,"1234567890.-"))==NULL){
          printf("Unexpected end to mapfile slit %d f\n",slitnum);
          return 1;}
        sscanf(lin,"%f",&coef_sag[ii]);
        if((lin=strpbrk(lin," "))==NULL){
          printf("Unexpected end to mapfile slit %d g\n",slitnum);
          return 1;}
	    }
      for(j=0;j<=osag;j++){
        if(isnan(s_coef[0]) || prtl) break;
        if(xdisper || chip<=4){                           //SITE specific!!
          coef_sag[j]+=s_coef[j];}
        else{
          coef_sag[j]-=s_coef[j];}
	    }
      for(ii=0;ii<=ord_sag;ii++) fprintf(newfile,"%13.6e ",coef_sag[ii]);
          fprintf(newfile,"\n");
      //tilt
      strncpy(lane,&dataline[i*140],140);
      i++;
      lin=&lane[0];
      for(ii=0;ii<=ord_tilt;ii++){
        if((lin=strpbrk(lin,"1234567890.-"))==NULL){
          printf("Unexpected end to mapfile slit %d f\n",slitnum);
          return 1;}
        sscanf(lin,"%f",&coef_tilt[ii]);
        if((lin=strpbrk(lin," "))==NULL){
          printf("Unexpected end to mapfile slit %d g\n",slitnum);
          return 1;}
      }
      for(j=0;j<=otlt;j++){
        if(isnan(t_coef[0])) break;
        coef_tilt[j]+=t_coef[j];}
      for(ii=0;ii<=ord_tilt;ii++) fprintf(newfile,"%13.6e ",coef_tilt[ii]);
      fprintf(newfile,"\n");
      //other functions
      for(j=0;j<2;j++){
        strncpy(lane,&dataline[i*140],140);
        fprintf(newfile,"%s",lane);
        i++;}
    }
  }
  fprintf(rmsfile,"%5d %15s %10.5f %10.5f\n", slitnum,name,stdev,ystdev);
  fflush(rmsfile);
  }
  cpgslct(pgdev2);
  cpgclos();
/*  cpgslct(pgdev3);
  cpgclos();*/
  fclose(rmsfile);
  if (diag==0) fclose(newfile);
  return 0;}

/*****************************************************************************\
*                                                                             *
*   WCALC  calculates weights for points in lsf calculations, based on        *
*          clumpiness                                                         *
*                                                                             *
\ ****************************************************************************/

void wcalc2(float *xvals, float *yvals, int nvals, float *wts){

  int   INT,i,j,k,nclump;
  float FLOAT,xmin, xmax,window,medist,wndw,d,dist;
  float *rank;
  int   *irnk;

  irnk=malloc(sizeof(INT)*nvals*nvals/2);
  rank=malloc(sizeof(FLOAT)*nvals*nvals/2);


  //calculate x value window
  nclump=5;
  nclump = (nclump>sqrt(nvals)) ? sqrt(nvals) : nclump;//#points in window
  xmin=9e+99;
  xmax=-9e+99;
  for(i=0;i<nvals;i++){
    xmin = (*(xvals+i)<xmin) ? *(xvals+i) : xmin;
    xmax = (*(xvals+i)>xmax) ? *(xvals+i) : xmax;}
  window=0.5*nclump*(xmax-xmin)/nvals;

  //caluclate median distance between all points
    k=0;
  for(i=0;i<nvals;i++){
    for(j=i+1;j<nvals;j++){
      d=fabsf(*(yvals+i)-*(yvals+j));
      order(&k,&d,irnk,rank);
      k++;}
  }
  medist=(*(rank+k/2));

  //for each point, find median distance to other points within window,
  //and from that calculate weight

  for(i=0;i<nvals;i++){
    wndw=window;
    if(i<nclump/2 || i> nvals-1-nclump/2) wndw=2*window;
    k=0;
    for(j=0;j<nvals;j++){
      if(i==j) continue;
      if(fabsf(*(xvals+i)-*(xvals+j))>wndw) continue;
      d=fabsf(*(yvals+i)-*(yvals+j));
      order(&k,&d,irnk,rank);
      k++;}
    if(k>=2){
      j=k/2;
      dist=(*(rank+j));
      if(j*2==k)dist=0.5*(dist+*(rank+j-1));
      //weight equals ratio of median distance between all points to median
      //distance to points in window
      *(wts+i)=medist/dist;}
    else{
      if(k==1){
        d = (d>0.1) ? d : 0.1;
        *(wts+i)=1/d;}
      else{
        *(wts+i)=1.0;}
    }
  }
  free(irnk);
  free(rank);
  return;}

  /*****************************************************************************\
  *                                                                             *
  *   WCALC  calculates weights for points in lsf calculations, based on        *
  *          clumpiness                                                         *
  *                                                                             *
  \ ****************************************************************************/

  void wcalc(int nclump, float *xvals, float *yvals, int nvals, float *wts){

    int   INT,i,j,k;
    float FLOAT,xmin, xmax,ymin,ymax,medist,wndw,d,dist,xfac,yfac;
    float  yvali,xvali;
    float  yvalj,xvalj;
    float  mn,m1;
    float *rank;
    int   *irnk;

    irnk=malloc(sizeof(INT)*nvals*nvals/2);
    rank=malloc(sizeof(FLOAT)*nvals*nvals/2);


    //calculate x,y value window
    if (!nclump) nclump=4;
    xmin=ymin=9e+99;
    xmax=ymax=-9e+99;

    for(i=0;i<nvals;i++){
      xmin = (*(xvals+i)<xmin) ? *(xvals+i) : xmin;
      xmax = (*(xvals+i)>xmax) ? *(xvals+i) : xmax;
      ymin = (*(yvals+i)<ymin) ? *(yvals+i) : ymin;
      ymax = (*(yvals+i)>ymax) ? *(yvals+i) : ymax;}

      if (xmax > xmin) {
        xfac=pow(1/(xmax-xmin),2);
      } else {
        xfac=4.0;
      }
      if (ymax > ymin) {
        yfac=pow(1/(ymax-ymin),2);
      } else {
        yfac=4.0;
      }

      //caluclate median distance between all points
      k=0;
      for(i=0;i<nvals;i++){
        for(j=i+1;j<nvals;j++){
          xvali = *(xvals+i);
          yvali = *(yvals+i);
          xvalj = *(xvals+j);
          yvalj = *(yvals+j);
          if ( isnan(xvali) || isnan(xvalj) || isnan(yvali) || isnan(yvalj)) {
            d = 0.0;
            } else {
            d=fabs(sqrt(pow(yvali-yvalj,2)*yfac+pow(xvali-xvalj,2)*xfac));
            order(&k,&d,irnk,rank);
            k++;}
        }
      }
      medist=(*(rank+k/2));

      //for each point, find distance to nth nearest neighbor, weight is inverse

      mn = 0;
      for(i=0;i<nvals;i++){
        k=0;
        for(j=0;j<nvals;j++){
          if(i==j) continue;
          xvali = *(xvals+i);
          yvali = *(yvals+i);
          xvalj = *(xvals+j);
          yvalj = *(yvals+j);
          if ( isnan(xvali) || isnan(xvalj) || isnan(yvali) || isnan(yvalj)) {
            d = 0.0;
           } else {
            d=fabs(sqrt(pow(yvali-yvalj,2)*yfac+pow(xvali-xvalj,2)*xfac));
            order(&k,&d,irnk,rank);
            k++;}
          }
        m1 = *(rank+nclump-1);
        if (!m1 || !medist) {
          *(wts+i)=1.0;
        } else {
           if (m1 == 0.0) {
           *(wts+i)=0.0;
          } else {
           *(wts+i)=pow(medist/m1,1);
          }
        }
        mn += *(wts+i);
        }

      mn = mn/nvals;
      if (mn != 0.0) {
       for(i=0;i<nvals;i++){
         *(wts+i)= *(wts+i)/mn;
       //  printf("weight for %d is %f and mean is %f\n", i+1, *(wts+i),mn);
       }
      }

      free(irnk);
      free(rank);
      return;}

  int comp_nums(const float *num1, const float *num2)
  {
     if (*num1 <  *num2) return -1;
     if (*num1 == *num2) return  0;
     if (*num1 >  *num2) return  1;
  }
