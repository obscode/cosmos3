/*****************************************************************************\
*                                                                             *
*  EXTRACT-SPEC  extracts 1 or 2-dim spectra using information in map files   *
*                  each spectrum is a separate 2-d fits extension             *
*                                                                             *
*  VERSION  01 Feb 2018                                                       *
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

  char     CHAR,file[80],flnm[1000],ifile[1000],dfile[80],line[133],string[133],
           c[7],name[40],dewar[10],ameans[12],answer[80],ch,trc[7],cardline[80],
           maskobj[80];
  char     *DATA_DIR,*HOME,*COS_HOME,*lin;
  int      nulvl,inshuffle,INT,status,nelem,i,j,here,xdisper,ord_disp,ord_sag,
           ord_tilt,ord_sagit,ord_slen,nslit,n_slit,n_slitpix,n_lampix,chip,
           nn,nshuffle,anynl,ibin,slitnum,ixpos,iypos,isloff,noff,lvmin,lvmax,
           nchip,nprv,n1,n2,n3,search,srcsz,nitr,one,edata,narg,shfld,k,dimen,
           ipv,nssub,n_hole,ipmax,sltype,useholes,nout,slt,logsam,nnn,nod,open,
           strace,star_ave,obj_ave,obj_each,ord_curv,ntrcpt,spcwdth,ord_trace,
           totpt,nstar,both_ave,nobj,pxintvl,npxnt,kpmx,kpmn,fobj,pchp,pislff,
           ibiny,otrace,edge,od1,year,month,day,ncard,hwidth,frow,lrow,slght,
           arr_len,o_value,loop,ncosmic,nprev,kmax,ncr,Pord,slnn,notilt;
  float    FLOAT,delta_lambda,delta_slit,telscale,min_lambda,max_lambda,ff,
           lambda0,lambda1,max_slit,slen,slit_len,slmin,slmax,lmin,lmax,sloff,
           curv,nlv,lambda,xcen,ycen,slnt,scatl,slitpos,ypos,yp,xpos,dltaslt,o,
           tilt,coef_disp[10],coef_sag[10],coef_tilt[10],coef_sagit[10],wshft,
           apv,prmax,pv,prav,offset,f1,f2,f3,coef_slen[10],slct,curmin,curmax,
           coef_curv[10],coef_trace[10],offmax,offmin,wid,sumpos,sumflx,lamean,
           xmin,xmax,ymin,ymax,xps,yps,theta,srot,wrot,camscale,prmin,ef,Pav,
           initoff,asloff,xmm,ymm,dscale,wt,pcirc,pline,slmmm,dsss,dssr,dsst,sign,
           coef_ddsp[9],flxng,fdate,pav,ww,pp,gain,vmax,var,wtot,ftot,a,b,Ncosm[1200],
           pixscale,ypixmid,sagmid;
  long     naxes[3],newaxes[3],firstel,firstelem[3],stack_len;
  int      *INTP,*iprtot,*prfl,*irnk,*irnk2,*type,*sltnm,*CRstack,**CRarray;
  float   *FLOATP,*stack,*estack,**spectrum,**espectrum,*image[9],*error[9],
          *rnk,*prvl,**array[9],**earray[9],*im,*err,*pnt,*offst,*trace,*wavl,
          *prv,*rnk2,*ftrace,*mtrace,*npoff,*mwavl,*xposs,*yposs,*Pstack,
          **Parray,*Vstack,**Varray,*lamind,**Pcoef,*Pcstack,*Tflux,*Tsigma,
          *Wstack,**Warray;
  double   *xlist,*ylist,**elist,DOUBLE,*DOUBLEP;
  fitsfile *fptr[9], *outptr;
  fitsdef  fitsdata;
  FILE     *mapfile,*datfile,*dat0file,*dat1file,*dat2file;

  one=1;
  pcirc=30.;
  pline=100.;
  dscale=0;
  nulvl=nprv=nout=totpt=nstar=strace=otrace=edge=0;
  firstelem[0]=firstelem[1]=firstelem[2]=1;
  firstel=1;
  edata=1;
  estack=stack=trace=mtrace=ftrace=wavl=NULL;
  espectrum=spectrum=NULL;
  prav=offset=theta =0;
  curmin=curmax=0;
  coef_trace[0]=0;
  ord_trace=0;
  xmin=ymin=10000;
  xmax=ymax=-10000;
  spcwdth=5; //half width used to get spectrum trace center of gravity
  pxintvl=10; //wavelength pixel binning of trace
  Pord=2; //order of spectrum trace for 1-d exraction

  //directories
  DATA_DIR=malloc(sizeof(CHAR)*80);
  DATA_DIR=getenv("COSMOS_IMAGE_DIR");
  if(DATA_DIR==NULL)die ("COSMOS_IMAGE_DIR undefined!\n");
  strcat(DATA_DIR,"/");

  //get parameters
  if(OpenCosParm("extract-spec")) die("Cannot find parameter file ");
  if(ReadParm_i("DIMENSION",&dimen)==1) die("parameter file error: dimension");
  if(ReadParm_s("SAMPLING",line)==1) die("parameter file error: sampling");
  logsam=1;
  if(strstr(line,"LOG")==NULL && strstr(line,"LOG")==NULL) logsam=0;
  if(ReadParm_r("DELTALAM",&delta_lambda)==1) die("parameterfile error:deltalam");
  star_ave=obj_ave=obj_each=both_ave=0;
  if(ReadParm_i("SEARCH",&search)==1) die("parameter file error: search");
  if(search){
    if(ReadParm_s("ALIGN",ameans)==1) die("parameterfile error: align");
    if(!strcasecmp(ameans,"none")) search=0;
    if(!strcasecmp(ameans,"star_ave")) star_ave=1 ;
    if(!strcasecmp(ameans,"obj_ave"))obj_ave =1 ;
    if(!strcasecmp(ameans,"both_ave"))both_ave =1 ;
    if(!strcasecmp(ameans,"obj_each"))obj_each=1;}
  if(search){
    if(ReadParm_s("TRACE",trc)==1) die("parameterfile error: trace");
    if(!strcasecmp(trc,"object")){
      otrace=1;}
    else{
      if(!strcasecmp(trc,"star")){
        strace=1;}
      }
    if(strace||otrace)if(ReadParm_i("TRACE_ORDR",&ord_curv)==1)
      die("parameterfile error: trace_order");
    if(obj_ave || both_ave){
      if(ReadParm_r("OBJ_FRAC",&slct)==1) die("parameterfile error: obj_frac");}
  }
  if(ReadParm_r("INIT_OFF",&initoff)==1) initoff=0;;
  if(ReadParm_r("DELTASLIT",&dltaslt)==1) die("parameterfile error: deltaslit");
  if(ReadParm_r("MINLAMBDA",&min_lambda)==1) die("parameterfile error: minlambda");
  if(ReadParm_r("MAXLAMBDA",&max_lambda)==1) die("parameterfile error: maxlambda");
  if(ReadParm_b("SUB_NS",&nssub)==1 || dimen==1) nssub=1;
  if(ReadParm_i("EDGE",&edge)==1) die("parameter file error edge");
  if(ReadParm_b("USE_HOLES",&useholes)==1) useholes=1;
  if(dimen==1) if(ReadParm_i("HWIDTH",&hwidth)==1) die("parameter file error: hwidth");
  if(ReadParm_r("GAIN",&gain)==1 && dimen==1) die("parameter file error: gain");
  if(search>0){
    srcsz=2*search+1;
    prvl=malloc(srcsz*sizeof(FLOAT));
    prv=malloc(srcsz*sizeof(FLOAT));
    prfl=malloc(srcsz*sizeof(INT));}

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
    scanf("%s",file);}
  else{
    if(!(strcmp(argv[1],"-m"))){
      narg=2;}
    else{
      if(!(strcmp(argv[3],"-m"))){
        narg=4;}
      else{
        printf("proper invocation: extract-spec -f framename -m mapfile\n");
        return 1;}
    }
    strcpy(dfile,argv[narg]);
    strcpy(file,dfile);
    strcat(file,".map");
    mapfile=fopen(file,"r");
    if(mapfile==NULL){
      printf("File %s does not exist!\n",file);
      return 1;}
    if(!(strcmp(argv[1],"-f"))){
      narg=2;}
    else{
      if(!(strcmp(argv[3],"-f"))){
        narg=4;}
      else{
        printf("proper invocation: extract-spec -f framename -m mapfile\n");
        return 1;}
    }
//    if(argc==6) slnn= atoi(argv[5]);
  notilt=0;
  if(argc==6 && !(strcmp(argv[5],"-t")))notilt=1;
    strcpy(file,argv[narg]);}

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
  od1=ord_disp-1;
  if(ord_disp>9 || ord_sag>9 || ord_tilt>9 || ord_sagit>9 || ord_slen>9){
    printf("Warning: maximum fit order = 9\n");
    return 1;}
  fgets(line,133,mapfile);
  if(sscanf(line,"Scale ~ %f Camscale = %f",&telscale,&camscale)<2){
    printf("Mapfile error: outdated mapfile?\n");
    return 1;}
  delta_slit=dltaslt/telscale;
  fgets(line,133,mapfile);
  if(!sscanf(line,"Lambda = %f %f",&lambda0,&lambda1)){
    printf("Mapfile error\n");
    return 1;}
  if(min_lambda==0){
    min_lambda=lambda0;}
  else{
    if(lambda0>min_lambda){
      printf("Warning: minimum lambda set to %f\n",lambda0);
      min_lambda=lambda0;}
    }
    if(max_lambda==0){
      max_lambda=lambda1;}
    else{
      if(lambda1<max_lambda){
        printf("Warning: maximum lambda set to %f\n",lambda1);
        max_lambda=lambda1;}
      }
  lamean=(max_lambda+min_lambda)/2.;
  if(logsam){
    n_lampix=(max_lambda-min_lambda)/delta_lambda+2;
    delta_lambda=log10(max_lambda/min_lambda)/n_lampix;
    n_lampix=log10(max_lambda/min_lambda)/delta_lambda+2;}
  else{
    n_lampix=(max_lambda-min_lambda)/delta_lambda+2;}
  if(dimen==1){
     lamind=malloc(sizeof(FLOAT)*n_lampix);
     for(i=0;i<n_lampix;i++) lamind[i]=(float) i;}

  fgets(line,133,mapfile);
  if(!sscanf(line,"Dewar = %s %d",dewar,&nchip)) die("Mapfile error");
  if(!strncmp(dewar,"SITE",4) || !strncmp(dewar,"E2V",3)){
    pcirc=50.;
    pline=550.;}
  fgets(line,133,mapfile);
  max_slit=0.;
  n_slit=n_hole=0;

  //find #of mapped slits, width of widest, average pixel scale
  pixscale=0.;
  while((fgets(line,133,mapfile))!=NULL){
    if(sscanf(line,"SLIT %d %s %d",&i,line,&sltype)) continue;
    if(!sscanf(line,"LENGTH = %f",&slen)) continue;
    if(sltype!=2){
      n_hole++;
      continue;}
    else{
      n_slit++;
      if(slen>max_slit) max_slit=slen;
      if(fgets(line,133,mapfile)==NULL) die("Mapfile error");
      if(!sscanf(line,"CHIP %d %f %f",&i,&slmin,&slmax)) die("Mapfile error");
      pixscale+=slen/(slmax-slmin);
      //compiler bug fix
      i=n_slit;}
    }
  if(!n_slit) die("No slits found");
  pixscale = pixscale/n_slit;
  if((strace || star_ave || both_ave) && !n_hole){
     die("Your map is missing the needed alignment holes");}
  rnk=malloc(sizeof(FLOAT)*4*(n_slit+n_hole));
  irnk=malloc(sizeof(INT)*4*(n_slit+n_hole));
  rnk2=malloc(sizeof(FLOAT)*n_lampix/pxintvl);
  irnk2=malloc(sizeof(INT)*n_lampix/pxintvl);
  iprtot=malloc(sizeof(INT)*(n_slit+n_hole));
  npoff=malloc(sizeof(FLOAT)*(n_slit+n_hole));
  xposs=malloc(sizeof(FLOAT)*(n_slit+n_hole));
  yposs=malloc(sizeof(FLOAT)*(n_slit+n_hole));
  type=malloc(sizeof(INT)*(n_slit+n_hole));
  sltnm=malloc(sizeof(INT)*n_slit);
  offst=malloc((n_slit+n_hole)*sizeof(FLOAT));
  for(i=0;i<n_slit+n_hole;i++) *(offst+i)=initoff;

  if(strace || otrace){
    mtrace=malloc(sizeof(FLOAT)*n_lampix*(n_slit+n_hole));
    ftrace=malloc(sizeof(FLOAT)*n_lampix);
    wavl=malloc(sizeof(FLOAT)*n_lampix);
    mwavl=malloc(sizeof(FLOAT)*n_lampix*(n_slit+n_hole));}
  if(search) trace=malloc(sizeof(FLOAT)*n_lampix);
  if(dimen==1){
    Tflux = calloc(n_lampix,sizeof(FLOAT));
    Tsigma=  calloc(n_lampix,sizeof(FLOAT));}

 //read data files

  strcpy(flnm,DATA_DIR);
  strcat(flnm,file);
  strcpy(ifile,flnm);
  for(i=1;i<=nchip;i++){
    strcpy(file,flnm);
    if(nchip>1){
      addbar(file);
      strcat(file,"c");
      sprintf(c,"%d",i);
      strcat(file,c);}
    strcat(file,".fits");
    status=0;
    if(OpenFitsFile(file,&fptr[i],&fitsdata)){
      printf("file %s does not exist!",file);
      return 1;}
      }

  //Read ccd data

  ibin=ibiny=fitsdata.binning;
  if(fitsdata.ybinning) ibiny=fitsdata.ybinning;
  naxes[0]=fitsdata.naxes[0];
  naxes[1]=fitsdata.naxes[1];
  if(fitsdata.naxis==2) edata=0;
  nelem=naxes[0]*naxes[1]*(1+edata);
  nshuffle=fitsdata.nshuffle;
  //fixup for new N&S mask design screwup
  sscanf(fitsdata.date,"%d-%d-%d",&year,&month,&day);
  fdate=year+ (float) month/12+(float) day/365;
  if(fdate > 2007.455) nshuffle=-nshuffle;
  nod=(int) (0.5+sqrt(fitsdata.ranod*fitsdata.ranod+fitsdata.decnod*
       fitsdata.decnod)/dltaslt);
  shfld = (nshuffle && !nssub) ? 1 : 0;

  for(i=1;i<=nchip;i++){
    image[i]=malloc(sizeof(FLOAT)*naxes[0]*naxes[1]*(1+edata));
    im=image[i];
    array[i]=malloc(sizeof(FLOATP)*naxes[1]);
    for(j=0;j<naxes[1];j++) *(array[i]+j)=im+naxes[0]*j;
    fits_read_img(fptr[i],TFLOAT,firstel,nelem,&nulvl,im,&anynl,&status);
    if(status)fits_die("Data file error",status);
    if(edata){
      error[i]=image[i]+naxes[0]*naxes[1];
      err=error[i];
      earray[i]=malloc(sizeof(FLOATP)*naxes[1]);
      for(j=0;j<naxes[1];j++) *(earray[i]+j)=err+naxes[0]*j;}
    printf("Read ccd chip %d\r",i);
    fflush(stdout);}
  printf("\n");

//output image file
  strcpy(file,"!");
  strcat(file,ifile);
  addbar(file);
  subbars(file);
  if(dimen==1){
    strcat(file,"_1spec.fits");}
  else{
    strcat(file,"_2spec.fits");}
  status=0;
  fits_create_file(&outptr,file,&status);
  if(status){
    printf("Unable to create output spectrum file (%d)\n",file);
    return 1;}

  //setup for trace determination and 1-d spatial image determination
  if(strace || otrace || dimen==1){
    o_value = ord_curv<3 ? 5 : ord_curv+2;
    xlist=malloc(sizeof(DOUBLE)*2*(o_value));
    ylist=malloc(sizeof(DOUBLE)*(o_value));
    elist=malloc(sizeof(DOUBLEP)*(o_value));
    for(i=0;i<o_value;i++){
      *(elist+i)=malloc(sizeof(DOUBLE)*(o_value));}
  }

  open=0;
  wid=8.;


  /*--------Find slit offsets and Extract Spectra----------------------------*/

  nn = (search) ? 0:1;
  if(!nn) nn = (strace || otrace) ? -1: 0;
  for(nitr=nn;nitr<=1;nitr++){
    //if nitr==0 and spectra traced, calculate mean spectral trace
    if(nitr==0 && (strace||otrace)){
      ord_trace=ord_curv;
      if(!open){cpgopen("/xwindow");open=1;}
      cpgpage;
      cpgpap(wid,1.);
      cpgsch(1);
      cpgscr(0,1,1,1);
      cpgscr(1,0,0,0);
      cpgask(0);
      cpgsci(1);
      cpgenv(min_lambda,max_lambda,-search,search,0,0);
      cpglab("Lambda","Delta","Apply mean correction? [Y/N]");
      cpgpt(totpt,mwavl,mtrace,16);
      //fit trace
      plyfit(mwavl,mtrace,totpt,ord_trace,coef_trace,xlist,ylist,elist);
      for(j=0;j<n_lampix;j++){
        *(wavl+j)=min_lambda+j*delta_lambda;
        *(ftrace+j)=polyvalue(*(wavl+j),coef_trace,ord_trace);}
       cpgsci(2);
      cpgslw(2);
      cpgline(n_lampix,wavl,ftrace);
      cpgslw(1);
      printf("Correct for spectrum curvature? [Y/N]: ");
      fflush(stdout);
//      fgets(answer,80,stdin);
//      if(!strncasecmp(answer,"q",1)) return 0;}
      cpgcurs(&xps,&yps,&ch);
      if(ch!='Y' && ch!='y'){
      ord_trace=0;
      coef_trace[0]=0;}
  }

    //if nitr=1, ready to extract, calculate various offsets
    if(nitr>0){
      curmax=curmin=0;
      //find min,max offsets if trace
      if(strace||otrace){
        curmin=1000.;
        curmax=-1000;
        for(lambda=min_lambda;lambda<=max_lambda;lambda+=delta_lambda){
          o=polyvalue(lambda,coef_trace,ord_trace);
          if(o>curmax) curmax=o;
          if(o<curmin) curmin=o;}
        }
      //calculate average offset if used
      if(star_ave || obj_ave || both_ave){
        prav=0;
        xmin-=0.1*(xmax-xmin);
        xmax+=0.09*(xmax-xmin);
        ymin-=0.1*(ymax-ymin);
        ymax+=0.09*(ymax-ymin);
        xmm=xmin+.03*(xmax-xmin);
        ymm=ymax+.03*(ymax-ymin);
        for(j=0;j<2;j++){
          nnn=0;
          if(!open){cpgopen("/xwindow");open=1;}
          cpgpage;
          cpgpap(wid,1.);
          cpgsch(1);
          cpgscr(0,1,1,1);
          cpgscr(1,0,0,0);
          cpgask(0);
          cpgsci(1);
          cpgenv(xmin,xmax,ymin,ymax,0,0);
          cpglab("X","Y","");
          cpgslw(5);
          cpgsfs(1);
          if(!star_ave){
            cpgsci(1);
            for(i=(1-slct)*nobj;i<=nobj;i++){
              nnn++;
              apv= *(npoff+*(sltnm+irnk[i]));
              xps=*(xposs+*(sltnm+irnk[i]));
              yps=*(yposs+*(sltnm+irnk[i]));
              cpgcirc(xps,yps,pcirc);
              cpgmove(xps,yps);
              xps+=pline*apv;
              cpgdraw(xps,yps);
              prav+=apv;}
            }
          if(!obj_ave){
            for(i=0;i<nslit;i++){
              if(*(type+i)==2) continue;
              nnn++;
              apv= *(npoff+i);
              xps=*(xposs+i);
              yps=*(yposs+i);
              cpgsci(2);
              cpgcirc(xps,yps,pcirc);
              cpgsci(1);
              cpgmove(xps,yps);
              xps+=pline*apv;
              cpgdraw(xps,yps);
              prav+=apv;}
          }
          cpgslw(1);
          if(j==0){
            offset=initoff+prav/nnn;//average slit offset, in pixels
            sprintf(answer,"Average slit offset = %3.1f pixels. Apply offset? (Y/N) "
                   ,offset);
            offset-=initoff;
            cpgtext(xmm,ymm,answer);
            cpgcurs(&xps,&yps,&ch);
            if(ch != 'y' && ch != 'Y'){
              cpgsci(0);
              cpgtext(xmm,ymm,answer);
              cpgsci(1);
              break;}
            for(i=0;i<nslit;i++) *(offst+i)+=offset;
            for(i=0;i<nslit;i++) *(npoff+i)-=offset;
            continue;}
        }
        //calculate rotation
        srot=wrot=0;
        if(!star_ave){
          for(i=(1-slct)*nobj;i<=nobj;i++){
            apv= *(npoff+*(sltnm+irnk[i]));
            xps=*(xposs+*(sltnm+irnk[i]));
            yps=*(yposs+*(sltnm+irnk[i]));
            wt=(yps>0)?yps:-yps;
            srot+=apv*wt/yps;
            wrot+=wt;}
        }
        if(!obj_ave){
          for(i=0;i<nslit;i++){
            if(*(type+i)==2) continue;
            apv= *(npoff+i);
            xps=*(xposs+i);
            yps=*(yposs+i);
            wt=(yps>0)?yps:-yps;
            srot+=apv*wt/yps;
            wrot+=wt;}
        }
          theta=srot/wrot;
          theta*=57.2958;
          sprintf(answer,"Average rotation = %4.2f degrees. Apply rotation? (Y/N) ", theta);
          theta/=57.2958;
          cpgtext(xmm,ymm,answer);
          cpgcurs(&xps,&yps,&ch);
          if(ch != 'y' && ch != 'Y') theta = 0.;

          //calculate scale
          if(!open){cpgopen("/xwindow");open=1;}
          cpgpage;
          cpgpap(wid,1.);
          cpgsch(1);
          cpgscr(0,1,1,1);
          cpgscr(1,0,0,0);
          cpgask(0);
          cpgsci(1);
          cpgenv(xmin,xmax,ymin,ymax,0,0);
          cpglab("X","Y","");
          cpgslw(5);
          cpgsfs(1);
          srot=0;
          wrot=0;
          if(!star_ave){
            cpgsci(1);
            for(i=(1-slct)*nobj;i<=nobj;i++){
              apv= *(npoff+*(sltnm+irnk[i]));
              xps=*(xposs+*(sltnm+irnk[i]));
              yps=*(yposs+*(sltnm+irnk[i]));
              wt=(xps>0)?xps:-xps;
              apv-=yps*theta;
              srot+=apv*wt/xps;
              wrot+=wt;
              cpgcirc(xps,yps,pcirc);
              cpgmove(xps,yps);
              xps+=pline*apv;
              cpgdraw(xps,yps);}
            }
          if(!obj_ave){
            for(i=0;i<nslit;i++){
              if(*(type+i)==2) continue;
              nn++;
              apv= *(npoff+i);
              xps=*(xposs+i);
              yps=*(yposs+i);
              wt=(xps>0)?xps:-xps;
              apv-=yps*theta;
              srot+=apv*wt/xps;
              wrot+=wt;
              cpgsci(2);
              cpgcirc(xps,yps,pcirc);
              cpgsci(1);
              cpgmove(xps,yps);
              xps+=pline*apv;
              cpgdraw(xps,yps);}
          }
          cpgsci(0);
          cpgtext(xmm,ymm,answer);
          cpgsci(1);
          dscale=srot/wrot;
          sprintf(answer,"Average scale change = %f. Apply? (Y/N) ", dscale);
          cpgtext(xmm,ymm,answer);
          cpgcurs(&xps,&yps,&ch);
          if(ch != 'y' && ch != 'Y') dscale=0;
          if(!open){cpgopen("/xwindow");open=1;}
          cpgpage;
          cpgpap(wid,1.);
          cpgsch(1);
          cpgscr(0,1,1,1);
          cpgscr(1,0,0,0);
          cpgask(0);
          cpgsci(1);
          cpgenv(xmin,xmax,ymin,ymax,0,0);
          cpglab("X","Y","");
          cpgslw(5);
          cpgsfs(1);
          wrot=srot=0;
          if(!star_ave){
            cpgsci(1);
            for(i=(1-slct)*nobj;i<=nobj;i++){
              apv= *(npoff+*(sltnm+irnk[i]));
              xps=*(xposs+*(sltnm+irnk[i]));
              yps=*(yposs+*(sltnm+irnk[i]));
              apv=apv-yps*theta-xps*dscale;
              srot+=apv*apv;
              wrot+=1.;
              cpgcirc(xps,yps,pcirc);
              cpgmove(xps,yps);
              xps+=pline*apv;
              cpgdraw(xps,yps);}
          }
          if(!obj_ave){
            for(i=0;i<nslit;i++){
              if(*(type+i)==2) continue;
              nn++;
              apv= *(npoff+i);
              xps=*(xposs+i);
              yps=*(yposs+i);
              apv=apv-yps*theta-xps*dscale;
              srot+=apv*apv;
              wrot+=1.;
              cpgsci(2);
              cpgcirc(xps,yps,pcirc);
              cpgsci(1);
              cpgmove(xps,yps);
              xps+=pline*apv;
              cpgdraw(xps,yps);}
          }
          srot=sqrt(srot/wrot);
          cpgsci(0);
          cpgtext(xmm,ymm,answer);
          cpgsci(1);
          sprintf(answer,"offset RMS = %4.2f pixels. Hit CR to continue\n",srot);
          cpgtext(xmm,ymm,answer);
          cpgcurs(&xps,&yps,&ch);
          cpgclos();
       }
          else{dscale=0;}

      //or use individual offsets

      if(obj_each){
        if(!open){cpgopen("/xwindow");open=1;}
        cpgpage;
        cpgpap(wid,1.);
        cpgsch(1);
        cpgscr(0,1,1,1);
        cpgscr(1,0,0,0);
        cpgask(0);
        cpgsci(1);
        xmax=nslit;
        ymin=-search;
        ymax=search;
        cpgenv(0,xmax,ymin,ymax,0,0);
        cpglab("n","delta","");
        cpgslw(5);
//        cpgsfs(1);
//        cpgsci(1);
        for(i=0;i<=nobj;i++){
          xps=*(xposs+irnk[i]);
          yps=*(yposs+irnk[i]);
          apv= *(npoff+*(sltnm+irnk[i]))-yps*theta-xps*dscale;
          xps=i;
          cpgpt1(xps,apv,4);}
        printf("Enter index of first object to use individual offset: ");
        scanf("%d",&fobj);
        prav=0;
        for(i=fobj;i<nobj;i++) prav+=*(npoff+*(sltnm+irnk[i]));
        prav/=nslit-fobj;
        for(i=0;i<fobj;i++) *(offst+*(sltnm+irnk[i]))=+prav;
        for(i=fobj;i<nobj;i++) *(offst+*(sltnm+irnk[i]))+=*(npoff+
                                                           *(sltnm+irnk[i]));}
      }
    //process slits
    if(nitr==0) printf("Finding slit offsets\n");
    if(nitr==1) printf("Extracting spectra\n");
    if(nitr<0) printf("Finding spectrum curvature\n");
    rewind(mapfile);
    nslit=-1;
    nobj=0;
    while(1){
      if((fgets(line,133,mapfile))==NULL)  die("Unexpected end to mapfile");
      if(sscanf(line,"SLIT %d %s",&slitnum, name)) break;}

    //New slit
    while(2){
      if(!strncmp(line,"END",3)) break;                           //end of data
      if(!sscanf(line,"SLIT %d %s %d",&slitnum, name,&sltype))
        die("Mapfile error");
      if((fgets(line,133,mapfile))==NULL) die("Unexpected end to mapfile");
      if(!strncmp(line,"END",3)) break;                           //end of data
      if(!sscanf(line,"LENGTH = %f POS = %f %f",&slit_len,&xps,&yps)) continue;     //no data in slit
      //found a slit with data
      nslit++;
      pchp=0;
      *(type+nslit)=sltype;
      if(!nitr && (sltype==1 || obj_ave || both_ave)){
        xps*=camscale;
        yps*=fabs((double)camscale);
        xmin = (xps<xmin) ? xps: xmin;
        ymin = (yps<ymin) ? yps: ymin;
        xmax = (xps>xmax) ? xps: xmax;
        ymax = (yps>ymax) ? yps: ymax;
        *(xposs+nslit)=xps;
        *(yposs+nslit)=yps;}
      fflush(stdout);
      if(nitr<1){
        for(i=0;i<srcsz;i++) *(prfl+i)=0;
        for(i=0;i<srcsz;i++) *(prvl+i)=0;}
      //extracted spectrum
      //bugfix
      noff=n_slitpix=(int)(slit_len/delta_slit)-2*edge*pixscale/delta_slit;
      if(shfld) n_slitpix*=2;
      ntrcpt=0;
      stack_len=n_slitpix*n_lampix;
      if(nitr){
        stack=(float *) realloc(stack,sizeof(float)*stack_len);
        for(i=0;i<stack_len;i++) *(stack+i)=0.;
        spectrum=(float **) realloc(spectrum,sizeof(float *)*n_slitpix);
        for(i=0;i<n_slitpix;i++) *(spectrum+i)=stack+i*n_lampix;
        if(edata){
          estack=(float *) realloc(estack,stack_len*sizeof(float));
          espectrum=(float **) realloc(espectrum,sizeof(float *)*n_slitpix);
          for(i=0;i<n_slitpix;i++) *(espectrum+i)=estack+i*n_lampix;}
      }
      //New Chip
      while(3){
        if((fgets(line,133,mapfile))==NULL) die("Unexpected end to mapfile");
        if(!strncmp(line,"END",3)) break;
        if(!sscanf(line,"CHIP %d %f %f %f %f",&chip,&slmin,&slmax,&lmin,&lmax))
          break;
        if(nshuffle) inshuffle = (slmax>slmin) ? nshuffle : -nshuffle;
        //dnod=nod*(slmax-slmin)/slit_len;
        sign = (slmax>slmin) ? 1: -1;
        //coef_disp data
        if((fgets(line,133,mapfile))==NULL) die("Unexpected end to mapfile");
        lin=&line[0];
        for(i=0;i<=ord_disp;i++){
          if((lin=strpbrk(lin,"1234567890.-"))==NULL)
            die("Unexpected end to mapfile");
          sscanf(lin,"%f",&coef_disp[i]);
          if((lin=strpbrk(lin," "))==NULL) die("Unexpected end to mapfile");}
        //dispersion derivative
        for(i=1;i<=ord_disp;i++){
          coef_ddsp[i-1]=i*coef_disp[i];}
        //skip inverse disp data
        if((fgets(line,133,mapfile))==NULL) die("Unexpected end to mapfile");
        //coef_sag dat
        if((fgets(line,133,mapfile))==NULL) die("Unexpected end to mapfile");
        lin=&line[0];
        for(i=0;i<=ord_sag;i++){
          if((lin=strpbrk(lin,"1234567890.-"))==NULL)
            die("Unexpected end to mapfile");
          sscanf(lin,"%f",&coef_sag[i]);
          if((lin=strpbrk(lin," "))==NULL) die("Unexpected end to mapfile");}
        //coef_tilt data
        if((fgets(line,133,mapfile))==NULL) die("Unexpected end to mapfile");
        lin=&line[0];
        for(i=0;i<=ord_tilt;i++){
          if((lin=strpbrk(lin,"1234567890.-"))==NULL)
            die("Unexpected end to mapfile");
          sscanf(lin,"%f",&coef_tilt[i]);
          if((lin=strpbrk(lin," "))==NULL) die("Unexpected end to mapfile");}
        //coef_sagit data
        if((fgets(line,133,mapfile))==NULL) die("Unexpected end to mapfile");
        lin=&line[0];
        for(i=0;i<=ord_sagit;i++){
          if((lin=strpbrk(lin,"1234567890.-"))==NULL)
            die("Unexpected end to mapfile");
          sscanf(lin,"%f",&coef_sagit[i]);
          if((lin=strpbrk(lin," "))==NULL) die("Unexpected end to mapfile");}
        //coef_slen data
        if((fgets(line,133,mapfile))==NULL) die("Unexpected end to mapfile");
        lin=&line[0];
        for(i=0;i<=ord_slen;i++){
          if((lin=strpbrk(lin,"1234567890.-"))==NULL)
            die("Unexpected end to mapfile");
          sscanf(lin,"%f",&coef_slen[i]);
          if((lin=strpbrk(lin," "))==NULL) die("Unexpected end to mapfile");}

        //do extraction

        //if not a star and nitr=-1 or (nitr=0 and star_ave method used,skip
        if(sltype==2 && ((nitr<0 && !otrace) || (nitr==0 && star_ave))) continue;
        slmmm=sign*(offst[nslit]+yps*theta+xps*dscale);
        slmin=slmin-slmmm;
        slmax=slmax-slmmm;
        //bugfix
        asloff=((slmin*slit_len/(slmax-slmin))/delta_slit)-edge*pixscale/delta_slit;
        isloff=(int)asloff;
        isloff+=sign;
        if(pchp &&isloff != pislff){
          isloff=pislff;}
        pchp=chip;
        pislff=isloff;
        sloff=isloff*delta_slit;
        if(logsam){
          lvmin=(int)(log10(lmin/min_lambda)/delta_lambda);
          lvmax=(int)(log10(lmax/min_lambda)/delta_lambda);}
        else{
          lvmin=((int)((lmin-min_lambda)/delta_lambda));
          lvmax=((int)((lmax-min_lambda)/delta_lambda));}
/*        printf("delta_lambda= %f\n",delta_lambda);
        printf("lmin = %f\n",lmin);
        printf("lmax = %f\n",lmax);
        printf("lvmin = %i\n",lvmin);
        printf("lvmax = %i\n\n",lvmax);*/
        if(lvmin<1) lvmin=1;
        if(lvmax>=n_lampix-1) lvmax=n_lampix-2;
        nlv=(int)lvmax-lvmin + 1;
        npxnt=0;
        sumflx=sumpos=0;
        for(j=lvmin;j<=lvmax;j++){
          if(logsam){
            // ln(10) = 2.30259
            lambda = min_lambda*exp(2.30259*(j*delta_lambda));}
          else{
            lambda = min_lambda+j*delta_lambda;}
          //printf("lambda = %f\n",lambda);
          npxnt++;
          xcen=polyvalue(lambda,coef_disp,ord_disp);
          //bug workaround
          if(isnan(xcen)){
            lambda+=.01;
            xcen=polyvalue(lambda,coef_disp,ord_disp);}
          ycen=polyvalue(lambda,coef_sag,ord_sag)+sign*(offst[nslit]
                +yps*theta+xps*dscale+polyvalue(lambda,coef_trace,ord_trace));
          if(isnan(ycen)){
            lambda+=.01;
            ycen=polyvalue(lambda,coef_sag,ord_sag)+sign*(offst[nslit]
                +yps*theta+xps*dscale+polyvalue(lambda,coef_trace,ord_trace));}
          if(isnan(xcen) || isnan(ycen) || isnan(sloff)){
            printf("numerical error slitnum %d %d\n",slitnum,j);
            return 1;}
          slnt=polyvalue(lambda,coef_slen,ord_slen);
          tilt=polyvalue(lambda,coef_tilt,ord_tilt);
          if(notilt) tilt=0;
          curv=polyvalue(lambda,coef_sagit,ord_sagit);
          scatl=slnt/slit_len;
          prmax=-1000000.;
          prmin=100000.;
          ipmax=0;
    //      tilt=-tilt;
          if(nitr<1) for(k=0;k<srcsz;k++) *(prv+k)=0;
         // if(nitr==1) printf("%d %f %f %f\n",slitnum,xcen,ycen,sloff);
          for(i=0;i<noff;i++){
            slitpos=sloff+i*delta_slit;
            slitpos*=scatl;
            ypos=(ycen+slitpos)/ibiny;
            //bugfix:
            ypixmid=ycen+slmin+slnt/2;
            sagmid=curv*4.*(powf((ypixmid-(ycen+slitpos))/slnt,2)-powf((ypixmid
                  -ycen)/slnt,2));
            xpos=(xcen+tilt*slitpos/(slmax-slmin)-sagmid)/ibin;
            if(!xdisper){
              xpos=(ycen+slitpos)/ibin;
              ypos=(xcen+tilt*slitpos/slnt-sagmid)/ibiny;}

            if(xpos<0 || ypos<0 || xpos>naxes[0]-2 || ypos>naxes[1]-2)continue;

            if(i>n_slitpix || j>n_lampix){
              printf("!! %d %d %d %d\n",ixpos,iypos,i,j);
              fflush(stdout);}

            //extract spectrum along slit at this wavelength

            //photometric calibrations to flux/angstrom/arcsec-binning
            flxng=(polyvalue(lambda,coef_ddsp,od1)/ibin)*delta_slit*scatl/ibiny;
            if(isnan(flxng)){
              lambda+=.01;
              flxng=(polyvalue(lambda,coef_ddsp,od1)/ibin)*delta_slit*
                scatl/ibiny;}
            if(flxng<0) flxng=-flxng;
            if(nitr==1){
              if(edata){
                *(*(espectrum+i)+j)=flxng*e_interpol(earray[chip],earray[chip],
                                               naxes,xpos,ypos);
                *(*(spectrum+i)+j)=flxng*e_interpol(array[chip],earray[chip],
                                                  naxes,xpos,ypos);}
              else{
                *(*(spectrum+i)+j)=flxng*f_interpol(array[chip],naxes,xpos,ypos);}

              //deal with shuffled region, if n&s
              if(nshuffle){
                //if desired, subtract shuffled spectrum
                if(nssub){
                  if(edata){
                    *(*(spectrum+i)+j)-=flxng*e_interpol(array[chip],earray[chip],
                                        naxes,xpos,ypos+inshuffle);
                    ff=flxng*e_interpol(earray[chip],earray[chip],naxes,xpos,
                                  ypos+inshuffle);
                    pnt= *(espectrum+i)+j;
                    if(ff<=0 || *pnt<=0){
                      *pnt=-1000.;
                      *(*(spectrum+i)+j)=0;}
                    else{
                      *pnt=sqrt((*pnt)*(*pnt)+ff*ff);}
                  }
                  else{
                    *(*(spectrum+i)+j) -= flxng*f_interpol(array[chip],naxes,
                                                           xpos,ypos+inshuffle);}
                  }
                //else, extract seperately


                else{
                  if(edata){
                   *(*(espectrum+i+noff)+j)=flxng*e_interpol(earray[chip],
                                          earray[chip],naxes,xpos,ypos+inshuffle);
                    *(*(spectrum+i+noff)+j) = flxng*e_interpol(array[chip],
                                          earray[chip],naxes,xpos,ypos+inshuffle);
                  }
                  else{
                    *(*(spectrum+i+noff)+j) = flxng*f_interpol(array[chip],
                                                naxes,xpos,ypos+inshuffle);}

                  if(isnan(*(*(spectrum+i+noff)+j))
                     ||(edata &&isnan(*(*(espectrum+i+noff)+j)))){
                    *(*(spectrum+i+noff)+j)=0;
                    if(edata) *(*(espectrum+i+noff)+j)=-999999;}
                }
              }
              if(isnan(*(*(spectrum+i)+j))||
                 (edata&&isnan(*(*(espectrum+i)+j)))){
                *(*(spectrum+i)+j)=0;
                if(edata) *(*(espectrum+i)+j)=-999999;}
            }

            //or check if this is peak of spectrum trace at this wavelength
            else{
              pv = (xdisper) ? sign*(ypos-ycen/ibiny): sign*(xpos-ycen/ibin);
              pv+=0.5;
              ipv=(int)pv;
              if(ipv<=search && ipv>=-search){
                if(edata){
                  ef=e_interpol(earray[chip],earray[chip],naxes,xpos,ypos);
                  ff=e_interpol(array[chip],earray[chip],naxes,xpos,ypos);}
                else{
                  ff=f_interpol(array[chip],naxes,xpos,ypos);}
                if(!isnan(ff) && (!edata || ef>0)){
                  *(prvl+ipv+search)+=ff;
                  *(prv+ipv+search)=ff;
                  prmin = (ff<prmin) ? ff: prmin;
                  if(prmax<ff){
                    prmax=ff;
                    ipmax=ipv+search;}
                  }
                else{
                  *(prv+ipv+search)=-99999.;}
              }
            }
          }
/**********************end of slit position loop******************************/

          //determine center of gravity of spectrum trace at this wavelength
          if((nitr<0 && ((strace && sltype==1) || (otrace && sltype==2)))
             ||nitr==0){
            kpmn=ipmax-spcwdth;
            kpmx=ipmax+spcwdth;
            kpmn = (kpmn<0) ? 0:kpmn;
            kpmx = (kpmx>srcsz-1) ? srcsz-1:kpmx;
            if(kpmx-ipmax < ipmax-kpmn) kpmn=2*ipmax-kpmx;
            if(ipmax-kpmn < kpmx-ipmax) kpmx = 2*ipmax-kpmn;
            for(k=kpmn;k<=kpmx;k++){
              if(*(prv+k)<-90000.) continue;
              sumflx+=(*(prv+k)-prmin)*k;
              sumpos+=*(prv+k)-prmin;}
            if(npxnt==pxintvl){
              if((strace||otrace)&&nitr<0)*(wavl+ntrcpt)=lambda-
                delta_lambda*pxintvl/2.;
              if(sumflx>0){
                *(trace+ntrcpt)=(sumflx/sumpos)-search;
//               *(trace+ntrcpt)=ipmax-search;
                ntrcpt++;}
              sumflx=sumpos=0;
              npxnt=0;}
          }
          if(nitr==0) *(prfl+ipmax)= *(prfl+ipmax)+1;


/***************************end of lambda loop********************************/
      }
/********************end of while(3) all chip analysis************************/
      }

      //if nitr==0 and right type object calculate offset for object
      if(nitr==0 && (both_ave || (sltype==1 && star_ave) || (sltype==2 &&
         (obj_ave||obj_each)))) {
        //find median profile offset
        for(i=0;i<ntrcpt;i++) order(&i,trace+i,irnk2,rnk2);
        *(npoff+nslit)=*(trace+irnk2[ntrcpt/2]);
//        printf("offset %f object %s\n",*(npoff+nslit),name);

        ipmax=0;
        prmin=1000000.;
        for(i=0;i<srcsz;i++){
          if(*(prvl+i)>ipmax){
            ipmax= *(prvl+i);
            ipv=i;}//*(prvl+ipv) is height of profile
          prmin=(prmin>*(prvl+i)) ? *(prvl+i):prmin;
        }
        /*alternative, find center of gravity of summed profile
        sumflx=sumpos=0;
        for(k=ipv-spcwdth;k<ipv+spcwdth;k++){
          sumflx+=(*(prv+k)-prmin)*k;
          sumpos+=*(prv+k);}
        *(npoff+nslit)=(sumflx/sumpos)-search;*/

        if(sltype==2 && (both_ave || obj_ave || obj_each)){//rank object strengths
          order(&nobj,prvl+ipv,irnk,rnk);
          *(sltnm+nobj)=nslit;//slit # of nth object
//            printf("%d %d %f %f\n",nobj,nslit,*(npoff+nslit),*(prvl+ipv));
          nobj++;}
        }

      //if nitr<0 fit spectrum traces
      if(nitr<0 && ((strace && sltype==1) || (otrace && sltype==2))){
        //plot trace
        if(!open){cpgopen("/xwindow");open=1;}
        cpgpage;
        cpgpap(wid,1.);
        cpgsch(1);
        cpgscr(0,1,1,1);
        cpgscr(1,0,0,0);
        cpgask(0);
        cpgsci(1);
        cpgenv(min_lambda,max_lambda,-search,search,0,0);
        strcpy(answer,name);
        strcat(answer," Hit A to accept, R to reject, Q to quit\n");
        cpglab("Lambda","Delta",answer);
        cpgpt(ntrcpt,wavl,trace,16);
        //fit trace
        plyfit(wavl,trace,ntrcpt,ord_curv,coef_curv,xlist,ylist,elist);
        for(j=0;j<ntrcpt;j++){
          *(ftrace+j)=polyvalue(*(wavl+j),coef_curv,ord_curv);}
        cpgsci(2);
        cpgline(ntrcpt,wavl,ftrace);
        printf("Hit A to accept, R to reject, Q to quit\n");
        fflush(stdout);
        cpgcurs(&xps,&yps,&ch);
        if(ch=='q' || ch=='Q') return 0;
//        fgets(answer,80,stdin);
//        if(!strncasecmp(answer,"q",1)) return 0;
        nstar++;
        //add shifted trace to composite
        if(ch=='a' || ch=='A'){
//        if(!strncasecmp(answer,"a",1)){
          wshft=polyvalue(lamean,coef_curv,ord_curv);
          for(j=0;j<ntrcpt;j++){
            *(mwavl+totpt)=*(wavl+j);
            *(mtrace+totpt)=*(trace+j)-wshft;
            totpt++;}
          }
        cpgclos;}

/****************************end of offset calculations*********************/

      /*--------------------write spectrum hdu-----------------------------*/

      //if nitr==1 and right kind of object, and 2-d extraction, write data
      if(nitr<1 || (!useholes && sltype!=2)) continue;

      // 2-dimensional Data
      if(dimen==2){
        newaxes[1]=noff;
        if(shfld) newaxes[1]*=2;}
      else{
        newaxes[1]=1;}
      newaxes[0]=n_lampix;
      newaxes[2]=2;
      //output hdu
      if(edata){
        fits_create_img(outptr,FLOAT_IMG,3,newaxes,&status);}
      else{
        fits_create_img(outptr,FLOAT_IMG,2,newaxes,&status);}
      if(status){
        printf("Error creating spectrum file 1023\n");
        return 1;}
//      fits_delete_record(outptr,7,&status);
//      fits_delete_record(outptr,7,&status);
      if(status){
        printf("Error creating spectrum file\n");
        return 1;}
      //first hdu? write info
      if(!nout){
        nout=1;
      //writing original header information to output file first
        status=0;
        ncard=0;
        while(1){
          ncard++;
          if(ncard>108) break;
          fits_read_record(fptr[1],ncard,cardline,&status);
          if(status){fits_die("Error reading file header",status);}
          if(strstr(cardline,"END")==cardline) break;
          if(strstr(cardline,"NAXIS")!=NULL) continue;
          if(strstr(cardline,"SIMPLE")!=NULL) continue;
          if(strstr(cardline,"BITPIX")!=NULL) continue;
          if(strstr(cardline,"BSCALE")!=NULL) continue;
          if(strstr(cardline,"BZERO")!=NULL) continue;
          if(strstr(cardline,"OBJECT")!=NULL) continue;
          if(strstr(cardline,"COMMENT")!=NULL && strlen(cardline)==7) continue;
          status=0;
          fits_write_record(outptr,cardline,&status);
          if(status) die("Error writing output file");}
          if(dimen==2){
            fits_write_key(outptr,TLOGICAL,"SHUFFLED",&shfld,
                           "does data contain shuffled region",&status);
            if(status){
                printf("Error 1 writing output spectrum file (%d)\n",status);
                return 1;} 
            fits_write_key(outptr,TINT,"NOD",&nod,"nod distance in pixels",&status);
            if(status){
                printf("Error 2 writing output spectrum file (%d)\n",status);
                return 1;}
            fits_write_key(outptr,TFLOAT,"D_SLIT",&dltaslt,
                           "slit interval in arcsec",&status);
            if(status){
                printf("Error 3 writing output spectrum file (%d)\n",status);
                return 1;}}
          if(useholes) n_slit+=n_hole;
          fits_write_key(outptr,TINT,"N_SLITS",&n_slit,"Number of spectra",
                         &status);
          if(status){
            printf("Error 4 writing output spectrum file (%d)\n",status);
            return 1;}
        }

        if(status){
          printf("Unable to create output spectrum file (%d)\n",status);
          return 1;}
        isloff=-isloff+1;
        if(isloff<0) die("Bad spectrum location");

        float crpix1,log_min_lambda,log_delta_lambda;
        int dcflag,dispaxis;

        crpix1 = 1.;
        dispaxis = 1;

        // Standard headers
        fits_write_key(outptr,TINT,"DISPAXIS",&dispaxis,"",&status);
        fits_write_key(outptr,TSTRING,"CTYPE1","LINEAR","",&status);
        fits_write_key(outptr,TSTRING,"CTYPE2","LINEAR","",&status);
        fits_write_key(outptr,TSTRING,"WAT0_001","system=world","",&status);
        fits_write_key(outptr,TSTRING,"WAT1_001",
                    "wtype=linear label=Wavelength units=Angstroms","",&status);
        fits_write_key(outptr,TSTRING,"WAT2_001","wtype=linear","",&status);

        if(logsam){
          // Log
          log_min_lambda = log10(min_lambda);
          log_delta_lambda = log10(max_lambda/min_lambda)/newaxes[0];
          fits_write_key(outptr,TFLOAT,"CRVAL1",&log_min_lambda,"",&status);
          fits_write_key(outptr,TFLOAT,"CDELT1",&log_delta_lambda,"",&status);
          fits_write_key(outptr,TFLOAT,"CD1_1",&log_delta_lambda,"",&status);
          fits_write_key(outptr,TFLOAT,"CRPIX1",&crpix1,"",&status);
          dcflag = 1;
          fits_write_key(outptr,TINT,"DC-FLAG",&dcflag,"",&status);
        } else {
          // Linear
          fits_write_key(outptr,TFLOAT,"CRVAL1",&min_lambda,"",&status);
          fits_write_key(outptr,TFLOAT,"CDELT1",&delta_lambda,"",&status);
          fits_write_key(outptr,TFLOAT,"CD1_1",&delta_lambda,"",&status);
          fits_write_key(outptr,TFLOAT,"CRPIX1",&crpix1,"",&status);
          dcflag = 0;
          fits_write_key(outptr,TINT,"DC-FLAG",&dcflag,"",&status);
        }

        // COSMOS Specific headers
          fits_write_key(outptr,TINT,"SLITNUM",&slitnum,"Slit number",&status);
          fits_write_key(outptr,TSTRING,"OBJECT",name,"Name of object",&status);
          fits_write_key(outptr,TINT,"SLITTYPE",&sltype,"Type of aperture",&status);
          if(dimen==2){
          fits_write_key(outptr,TINT,"SLITLEN",&noff,"slit length in pixels",
                         &status);
          fits_write_key(outptr,TINT,"CNTRLINE",&isloff,"Spectrm central row",
                         &status);}
          if(status){
            printf("Error 5 writing output spectrum file (%d)\n",status);
            return 1;}
        fits_read_key(fptr[1],TSTRING,"OBJECT",&maskobj,string,&status);
        if(status){
            printf("Error reading input file header (%d)\n",status);
            return 1;}
        fits_write_key(outptr,TSTRING,"MASKOBJ",&maskobj,
                       "original value of OBJECT",&status);

        if(status){
          printf("Unable to create output spectrum file (%d)\n",status);
          return 1;}
        if(dimen==2){
          firstelem[2]=1;
          fits_write_pix(outptr,TFLOAT,firstelem,stack_len,stack,&status);
          if(edata){
            firstelem[2]=2;
            fits_write_pix(outptr,TFLOAT,firstelem,stack_len,estack,&status);}
          printf("Writing spectrum %d\r",slitnum);
          if(status){
            printf("Error writing spectrum file %d\n", status);
            return 1;}
//          if(edata) free(estack);
          }

      // 1-dimension data. do optimal extraction

      else{
        frow= hwidth>isloff ? 0 : isloff-hwidth;
        lrow= isloff+hwidth>noff-1 ? noff-1 : isloff+hwidth;
        slght=lrow-frow+1;
        arr_len=n_lampix*slght;
        // S,W array
        Vstack=(float *) calloc(arr_len,sizeof(FLOAT));
        Varray=(float **) calloc(slght,sizeof(FLOATP));
        for(i=0;i<slght;i++) Varray[i]=Vstack+i*n_lampix;
        Wstack=(float *) calloc(arr_len,sizeof(FLOAT));
        Warray=(float **) calloc(slght,sizeof(FLOATP));
        for(i=0;i<slght;i++) Warray[i]=Wstack+i*n_lampix;
        // CR array
        CRstack=(int *) calloc(arr_len,sizeof(INT));
        CRarray=(int **) calloc(slght,sizeof(INTP));
        for(i=0;i<slght;i++) CRarray[i]=CRstack+i*n_lampix;
        // P array
        Pstack=(float *) calloc(arr_len,sizeof(FLOAT));
        Parray=(float **) calloc(slght,sizeof(FLOATP));
        for(i=0;i<slght;i++) Parray[i]=Pstack+i*n_lampix;
        // coef array
        Pcstack= (float *) calloc(slght*(Pord+1),sizeof(FLOAT));
        Pcoef=(float **) calloc(slght,sizeof(FLOATP));
        for(i=0;i<slght;i++) Pcoef[i]=Pcstack+i*(Pord+1);

        //initial values of total flux, variance, pixel variance, weights,mask
        for(j=0;j<n_lampix;j++){
          Tflux[j]=0;
          for(k=0;k<slght;k++){
            *(CRarray[k]+j) = (*(espectrum[k+frow]+j) <=0) ? 0 : 1;
            *(Varray[k]+j)= pow(*(espectrum[k+frow]+j),2);
            Tflux[j]+= *(spectrum[k+frow]+j)*(*(CRarray[k]+j));}
          }
        //loop entire fit

        for(loop=0;loop<4;loop++){
          // clean CRs fit line profiles
          for(j=0;j<n_lampix;j++){
              //reset CR array
            for(k=0;k<slght;k++){
              *(CRarray[k]+j) = (*(espectrum[k+frow]+j) <=0 || Tflux[k]==0.) ? 0 : 1;}
            }
          while(1){
            ncosmic=0;
            for(j=0;j<n_lampix;j++){
              for(k=0;k<slght;k++){
                if(*(CRarray[k]+j)==1 && *(Varray[k]+j)>0. &&Tflux[j]!=0.){
                  *(Parray[k]+j)= *(spectrum[k+frow]+j)/Tflux[j];
                  *(Warray[k]+j)=pow(Tflux[j],2)/(*(Varray[k]+j));}
                else{
                  *(Parray[k]+j)=0.0;
                  *(Warray[k]+j)=0.0;}
                }
              //renormalize
              pp=0;
              for(k=0;k<slght;k++){
                 pp+= *(Parray[k]+j);
                 if(isnan(*Parray[k]+j)) printf("nan %d %d %d %d\n",slitnum,loop,k,j);}

              if(pp>0) for(k=0;k<slght;k++) *(Parray[k]+j)/=pp;
              }
            for(k=0;k<slght;k++){
              plyfit_w(lamind,Parray[k],Warray[k],n_lampix,Pord,Pcoef[k],xlist,
                       ylist,elist);
              for(j=0;j<n_lampix;j++){
                //if subtracted n&s data, reject high or low pixels
                if(nshuffle){
                  if(*(CRarray[k]+j) && pow(*(spectrum[k+frow]+j)-polyvalue(
                     (float)j,Pcoef[k],Pord)*Tflux[j],2)>16.*(*(Varray[k]+j))){
                       ncosmic++;
                       *(CRarray[k]+j)=0;}
                     }
                //if  normal data only reject high pixels
                else{
                  if(*(CRarray[k]+j) && -polyvalue((float)j,Pcoef[k],Pord)*
                     Tflux[j]+(*(spectrum[k+frow]+j))>4.*sqrt(*(Varray[k]+j))){
                       ncosmic++;
                       *(CRarray[k]+j)=0;}
                     }
                }
                 }
                 break;
            if(ncosmic==0) break;}
          //construct cleaned profile, total spectrum, new variances
          for(j=0;j<n_lampix;j++){
            Ncosm[j]=0;
              ncr=0;

              //reset CR array
              for(k=0;k<slght;k++){
                *(CRarray[k]+j) = (*(espectrum[k+frow]+j) <=0) ? 0 : 1;}
              while(1){
                wtot=ftot=0;
                for(k=0;k<slght;k++){
                  pp=polyvalue((float)j,Pcoef[k],Pord);
                  *(Varray[k]+j)= pow(*(espectrum[k+frow]+j),2)+(*(CRarray[k]+j)*
                                   (pp*Tflux[j]-(*(spectrum[k+frow]+j))))/gain;
                  if(*(Varray[k]+j)>0.){
                    ftot+=pp*(*(spectrum[k+frow]+j))*(*(CRarray[k]+j))/(*(Varray[k]+j));
                    wtot+=pp*pp*(*(CRarray[k]+j))/(*(Varray[k]+j));}
                }
                if(wtot>0){
                  Tflux[j]= ftot/wtot;
                  Tsigma[j]=sqrt(1/wtot);}
                else{
                  Tsigma[j]=0.;
                  for(k=0;k<slght;k++){
                    *(CRarray[k]+j)=0;}
                  }
                vmax=0;
                for(k=0;k<slght;k++){
                  var=*(CRarray[k]+j)*(*(spectrum[k+frow]+j)-polyvalue(
                      (float)j,Pcoef[k],Pord)*Tflux[j])/sqrt(*(Varray[k]+j));
                  if((nshuffle && fabs(var)>vmax) || (!nshuffle && var>vmax)){
                    vmax=var;
                    kmax=k;}
                  }
                if(vmax>5.){
                  *(CRarray[kmax]+j)=0;
                  Ncosm[j]++;
                  ncr++;
                  continue;}
                break;}
              }
            }//end of iteration

        free(Wstack);
        free(Warray);
        free(CRstack);
        free(CRarray);
        free(Pcstack);
        free(Pcoef);
        free(Pstack);
        free(Parray);
        free(Vstack);
        free(Varray);
        //write output file
        firstelem[2]=1;
        fits_write_pix(outptr,TFLOAT,firstelem,n_lampix,Tflux,&status);
        if(edata){
          firstelem[2]=2;
          fits_write_pix(outptr,TFLOAT,firstelem,n_lampix,Tsigma,&status);}
        printf("Writing spectrum %d\r",slitnum);
        if(status){
          printf("Error writing spectrum file %d\n", status);
          return 1;}

      }
    } //end of slit
  }
  // Free some memory will you!!!
  if(search>0){
    free(prvl);
    free(prv);
    free(prfl);
  }

  for(i=1;i<=nchip;i++){
    free(image[i]);
    free(array[i]);
    if(edata){
      free(earray[i]);
    }
  }

  if(strace || otrace || dimen==1){
    free(xlist);
    free(ylist);
    for(i=0;i<ord_curv+2;i++){
      free(*(elist+i));
    }
    free(elist);
  }

  free(rnk);
  free(irnk);
  free(rnk2);
  free(irnk2);
  free(iprtot);
  free(npoff);
  free(xposs);
  free(yposs);
  free(type);
  free(sltnm);
  free(offst);

  if(strace || otrace){
    free(mtrace);
    free(ftrace);
    free(wavl);
    free(mwavl);
  }
  if(search) free(trace);

  printf("\n\n%d spectra extracted\n",n_slit);
  fits_close_file(outptr,&status);
  return 0;}
