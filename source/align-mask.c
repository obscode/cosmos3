/*****************************************************************************\
*                                                                                                                                          *
*  ALIGN-MASK    compares predicted  and observed positions of a direct mask                       *
*                image to compute corrections to a dewar offset file                                                 *
*                                                                                                                                          *
*  VERSION       8 Jan 2018                                                                                             *
*                                                                                                                                          *
*  hidden features; if extra input parameter -c n, do only chip n                                    *
*                   if extra input parameter -p, write file with  x,y,dx,dy
*                   if extra parameter -l input one spectral line to use                           *
*                   q in answers quits program                                                                  *
\*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>

#include "cpgplot.h"
#include "ioutils.h"
#include "kdcutil.h"
#include "mgutils.h"
#include "optutils.h"
#include "cosmos.h"
#include "clardy.h"

int comp_nums(const float *, const float *);

int  main(int argc,char *argv[]){

  int     indx[10000][9],bitpix,naxis,nob[9],narg,status,here,i,chip,n,ind,boxwdth,boxht,misline[2],
			  nulvl,halfht,halfwd,nobj,srchbox,nelem,chip1,j,INT,update,nmiss,misalign,not,
			  magfac,np1,np2,chipnum,anynl,nneg,npos,s1,s2,s3,s4,npr,ntot,n_pairs,nimag,sign,
			  nchip,contig,k,xcorn,ycorn,bflag,ixp,iyp,ii,jj,ixtop,iytop,ypar,ndb,changes,
			  fxpxl,fypxl,i1,i2,i3,i4,nx0,nx1,isspec,out,nlines,aver,ibin,ybin,autol,nit,
        holes,niter,rsign;
  int      **array,*image,*INTP,**row,*sbox,asize;
  float   FLOAT,lmb[100],xccd,yccd,xcfp,ycfp,xlfp,ylfp,xlcd,ylcd,lmin,lmax,
			  xcd,ycd,x8,y8,dj,xcen,xmid,ymid,dx,dy,sumax,sum,xnew,ynew,scal,xxx,dav,
			  yyy,deltax,deltay,avneg,avpos,dtheta,xv,yv,pscale,edge,theta,thetas,th,dalign[100],
			  thetav,xavv,yavv,divid2,sc,scale,xm8,ym8,sumx,sumy,lines[100],bx,by,misal,delal,
			  ddx,mnx,mxx,mny,mxy,aspct,wid,xa[2],ya[2],a2,a3,slope,xdd,ydd,sigma,
			  t,flmx,flmn,flx,flxrat,lamb,sigma1,xrfp,yrfp,lcutor,rat,totthet,
			  totdx,totdy,totmis,totsc;
  double relat[2];
  int       *irnk,*irnk_y,*chp,*aptype,*lind,*lnd,NS;
  float    *tlist,*daxis,*mslope,*x,*y,*xs,*ys,*xd,*yd,*rank,*rank_y,*apwidth,
           *oaxis,*xfp,*yfp,*lambda;
  double   DOUBLE,cur_wavl, cur_temp,sc2,x1,x2,y1,y2,xs1,xs2,ys1,ys2,r_1,rs_1;
  double *xl,*yl,*xr,*yr;
  char     flnm[80],file[1000],maskfile[80],imacsdir[133],line[133],SMF_FILE[80],
			  camdef[80],camoff[80],linfile[80],binning[10],answer[80],dfile[1000],lane[80],
			  ifile[1000],CHAR,oldfile[200],odfile[200],command[500],parm[80],value[80];
  char     *rind,*DATA_DIR,*COS_HOME;
  long     naxes[3],firstel,firste[2],nn;
  vect2    mpleft,mpright,vect;
  vect3    mskpos,ccdpos,campos;
  objq     *oq, *oq1;
  slit      slitdat;
  obsdef   obsdata;
  Obs     * obset = NULL;
  element  *instrument;
  FILE     *outfile,*infile,*linefile,*datfile;
  fitsfile *fptr[9];
  dewdat   dewinfo;
  fitsdef  fitsinfo;
  float     pi=3.1415926;

  float    dxmed,dymed,xmad,ymad,lcut,ceiling;
  int      l,loop,nread,maxpairs;
  float    *csort;

  edge=1.;
  update=0;
  nulvl=0;
  out=0;
  firstel=1;
  chipnum=0;
  nlines=0;
  lmax=0;
  lmin=10000.;
  NS=0;
  maxpairs=1000;
  totdx=totdy=totmis=totthet=0;
  totsc=1.0;
  niter=4;//number of iterations in auto mode
  //  flxrat=1.20; //required minimum ratio of max to min pixel for real aperture

  //data directory

  DATA_DIR=malloc(sizeof(CHAR)*80);
  DATA_DIR=getenv("COSMOS_IMAGE_DIR");
  if(DATA_DIR==NULL){
	 printf("COSMOS_IMAGE_DIR undefined!\n");
	 return 1;}
  strcat(DATA_DIR,"/");

  //get parameters

  strcpy(file,"align-mask");
  if(OpenCosParm(file)!=0) die("Cannot open align-mask parameter file!");
  if(ReadParm_i("SEARCHBOX",&srchbox)==1) die("parameterfile error a");
  if(ReadParm_i("MAGFACTOR",&magfac)==1) die("parameterfile error b");
  if(ReadParm_s("LAMFILE",linfile)==1) die("parameterfile error c");
  if(ReadParm_i("NAVER",&aver)==1) die("parameter file error d");
  if(ReadParm_r("THRESHOLD",&flxrat)==1) die("parameter file error e");
  if(ReadParm_r("SIGLIMIT",&lcutor)==1) lcutor=10;
  if(ReadParm_r("MAXFLUX",&ceiling)==1) ceiling=50000;
  if(ReadParm_r("WINDOWSIZE",&wid)==1) die("parameterfile error f");
  if(ReadParm_b("MISALIGN",&misalign)==1) misalign=0;
  if(ReadParm_b("AUTO",&autol)==1)autol=0;
  if(ReadParm_b("USE_HOLES",&holes)==1) die("parameterfile error g");
  if(!autol) niter=20;
  /*------------------------get input data----------------------------------*/

  //prompt for data

  if(argc<5){
	  printf("Usage: align-mask -o obserset -f imagefile\n");
	  return 0;}

  //data on command line

  else{
	 if(!(strcmp(argv[1],"-o"))){
		narg=2;}
	 else{
		if(!(strcmp(argv[3],"-o"))){
	narg=4;}
		else{
	printf("proper invocation: align-mask -o obserset -f imagefile\n");
	return 1;}
	 }
	 strcpy(maskfile,argv[narg]);
	 if(ReadObsDef(maskfile,&obsdata)!=0){
		printf("Error reading observation definition file %s!\n",maskfile);
		return 1;}
   //rotation handedness
   rsign=1;
   if(!strcmp("IMACS",obsdata.instrument) && !strcmp("SHORT",obsdata.camera)){
     rsign=-1;}
	 if(!(strcmp(argv[1],"-f"))){
		narg=2;}
	 else{
		if(!(strcmp(argv[3],"-f"))){
	narg=4;}
		else{
	printf("proper invocation: align-mask -o obserset -f imagefile\n");
	return 1;}
		}
	 strcpy(dfile,argv[narg]);
	 if(argc>=6){
		if(!(strcmp(argv[5],"-c"))){
		  sscanf(argv[6],"%d",&chipnum);}
		if(!(strcmp(argv[5],"-p"))){
		  datfile=fopen("align-mask.dat","w");
		  out=1;}
	 if(!(strcmp(argv[5],"-l"))){
		  sscanf(argv[6],"%f",lines);
		  nlines=1;}
		}
	 }

  //get observing setup

  //mask data
  strcpy(flnm,obsdata.mask);
  strcat(flnm,".SMF");
  if(ReadSMFfile(flnm,&obset)){
		i=1;
	return 1;}

	 if(SetupCamera(&obsdata)==1) return 1;
  //chip data
  Getchipdat(&dewinfo);

  /********************* get data file info **********************************/

  strcpy(ifile,DATA_DIR);
  strcat(ifile,dfile);
  here=1;
  nchip=dewinfo.nchip;
  nread=0;
  for(i=1;i<=nchip;i++){
	 strcpy(file,ifile);
	 //printf("%s\n",file);
	 if(nchip==1){
		strcat(file,".fits");}
	 else{
     if(strchr(dfile,'_')) addbar(file);
		strcat(file,"c");
		sprintf(line,"%d.fits",i);
		strcat(file,line);}
//	 printf("\rReading %s",file);
	 status=OpenFitsFile(file,&fptr[i],&fitsinfo);
	 if(status) fits_die("Error opening data file" ,status);
	 //if(status==104){printf("Could not read file %s\n" ,file);}
	 nread++;}
  if(nread==0) fits_die("Error opening data file" ,status);
  bitpix=fitsinfo.bitpix;
  naxes[0]=fitsinfo.naxes[0];
  naxes[1]=fitsinfo.naxes[1];
  ibin=fitsinfo.binning;
  ybin=fitsinfo.ybinning;
  if(!ybin) ybin=ibin;
  nelem=naxes[0]*naxes[1];
  image=malloc(sizeof(INT)*nelem);
  array=malloc(sizeof(INTP)*naxes[1]);
  for(j=0;j<naxes[1];j++) *(array+j)=image+naxes[0]*j;

  //Additional info

  cur_temp=obset->temp;
  if(!strcmp(obsdata.grating,"NS")||!strcmp(obsdata.instrument,"LDSS3")){
		NS=1;}
  if(strcmp(obsdata.mode,"SPEC")){
	 lines[0]=obset->cw;
	 nlines=1;
	 isspec=0;}
  else{
	 isspec=1;
	 //read spectral line file
	 if(!nlines){
		  linefile=fopen(linfile,"r");
		  if(linefile==NULL) die("Cannot open line file");
		  nlines=0;
		  nmiss=0;
		  while(nlines<100){
				if(fgets(line,133,linefile)==NULL) break;
				sscanf(line,"%f",&lines[nlines]);
				if(lines[nlines]<lmin) lmin=lines[nlines];
				if(lines[nlines]>lmax) lmax=lines[nlines];
				nlines++;}
		  mslope=malloc(sizeof(FLOAT)*nlines*nlines);}
	 }

	 status=cpgopen("/xwindow");
	 if(status <= 0){
		  printf("plot error %d",status);
		  return 1;}
	 //count objects
	 nobj=0;
	 oq=obset->ob;
	 oq1=oq;
	 while(1){
		  oq=oq->next;
		  nobj++;
		  if(oq==oq1) break;}
	 nimag=nobj*nlines;
	 //create arrays
	 lind=malloc(sizeof(INT)*nimag);
	 apwidth=malloc(sizeof(FLOAT)*nimag);
	 xl=malloc(sizeof(DOUBLE)*nimag);
	 yl=malloc(sizeof(DOUBLE)*nimag);
	 xfp=malloc(sizeof(FLOAT)*nimag);
	 yfp=malloc(sizeof(FLOAT)*nimag);
	 xr=malloc(sizeof(DOUBLE)*nimag);
	 yr=malloc(sizeof(DOUBLE)*nimag);
	 aptype=malloc(sizeof(INT)*nimag);
	 x=malloc(sizeof(FLOAT)*nimag);
	 y=malloc(sizeof(FLOAT)*nimag);
	 xs=malloc(sizeof(FLOAT)*nimag);
	 ys=malloc(sizeof(FLOAT)*nimag);
	 xd=malloc(sizeof(FLOAT)*nimag);
	 yd=malloc(sizeof(FLOAT)*nimag);
	 chp=malloc(sizeof(INT)*nimag);
	 lnd=malloc(sizeof(INT)*nimag);
   lambda=malloc(sizeof(FLOAT)*nimag);
	 asize = maxpairs>2*nimag+1 ? maxpairs : 2*nimag+1;
	 irnk=malloc(sizeof(INT)*asize);
	 irnk_y=malloc(sizeof(INT)*2*nimag+1);
	 rank=malloc(sizeof(FLOAT)*asize);
	 rank_y=malloc(sizeof(FLOAT)*2*nimag+1);

	 if(SetupCamera(&obsdata)==1) return 1;

    //iterate fit

	 for(nit=0;nit<=niter;nit++){
      changes=0;
		//Instrument setup
		SetupInstr(&obsdata,&instrument);
		lcut=lcutor;
		//find predicted positions of all slits
		for(i=1;i<=nchip;i++) nob[i]=0;
		nobj=0;
		sigma=0.;
		oq=obset->ob;
		oq1=oq;

		//read in slit positions
		while(1){
			if(oq==NULL) die("SMF file error");
			slitdat=oq->slit;
			for(i=0;i<nlines;i++){
        if(!holes && slitdat.shape !=2) break;
				cur_wavl=lamb=lines[i];
				while(1){
					lind[nobj]=i;
					//center
					mskpos = get3de2v(oq->smpos,0.0);
					ccdpos = Op_transform(mskpos,instrument,cur_wavl,cur_temp);
					xcfp=ccdpos.x;
					ycfp=ccdpos.y;
					chip=fp2ccd(xcfp,ycfp,&xccd,&yccd,cur_wavl);
					//is image on a chip?
					if(chip<1) break;
					//just use one chip?
					if(chipnum && chip != chipnum) break;
					//width
					vect.x=slitdat.width;
					vect.y=0.;
					mpleft=sum2vect(oq->smpos,vect);
					mskpos = get3de2v(mpleft,0.0);
					campos = Op_transform(mskpos,instrument,cur_wavl,cur_temp);
					xlfp=campos.x;
					ylfp=campos.y;
					chip1=fp2ccd(xlfp,ylfp,&xlcd,&ylcd,cur_wavl);
					if(chip!=chip1)break;
					apwidth[nobj]=(sqrt((xlcd-xccd)*(xlcd-xccd)/(ibin*ibin)+(ylcd-yccd)*(ylcd-
											 yccd)/(ybin*ybin)));
					if(apwidth[nobj]<1.) apwidth[nobj]=1.;
					//"left"=minvalue slit end
					mpleft=sum2vect(oq->smpos,lslit(oq->slit));
					mskpos = get3de2v(mpleft,0.0);
					campos = Op_transform(mskpos,instrument,cur_wavl,cur_temp);
					xlfp=campos.x;
					ylfp=campos.y;
					chip1=fp2ccd(xlfp,ylfp,&xlcd,&ylcd,cur_wavl);
					//is left end on same chip?
					if(chip!=chip1)break;
          lambda[nobj]=cur_wavl;
					xl[nobj]=xlcd/ibin;
					yl[nobj]=ylcd/ybin;
					//"right"=maxvalue slit end
					mpright=sum2vect(oq->smpos,rslit(oq->slit));
					mskpos = get3de2v(mpright,0.0);
					campos = Op_transform(mskpos,instrument,cur_wavl,cur_temp);
					xrfp=campos.x;
					yrfp=campos.y;
					xfp[nobj]=0.5*(xlfp+xrfp);
					yfp[nobj]=0.5*(ylfp+yrfp);
					chip1=fp2ccd(xrfp,yrfp,&xcd,&ycd,cur_wavl);
					//is right end on same chip?
					if(chip !=chip1) break;
					xr[nobj]=xcd/ibin;
					yr[nobj]=ycd/ybin;
					aptype[nobj]=slitdat.shape;
					indx[nob[chip]][chip]=nobj;
					nob[chip]++;
					nobj++;
					break;}
					}
				oq=oq->next;
				if(oq==oq1) break;}
		ntot=0;
		mxx=mxy=-100000.;
		mnx=mny=100000.;

		//for each chip, locate each slit image

		//read in image
		for(chip=1;chip<=nchip;chip++){
			//printf("Analysing chip %d\r",chip);
			//fflush(stdout);
			if(nob[chip]==0) continue;
			status=0;
			status=ReadFitsFile(&fptr[chip],&fitsinfo,TINT,image);
			if(status) fits_die("Error reading data file",status);
			for(n=0;n<nob[chip];n++){
				ind=indx[n][chip];

				//extraction box centered on aperture

				//if hole
				if(aptype[ind] != 2){
					xmid=(xr[ind]+xl[ind])/2.;                    //xmid,ymid = predicted centroid
					ymid=(yr[ind]+yl[ind])/2.;
					halfht=halfwd=(int)(apwidth[ind]+0.5);
					boxwdth=boxht=2*halfht+1;
					sbox=malloc(sizeof(INT)*boxwdth*boxwdth);
					row=malloc(sizeof(INTP)*boxwdth);
					for(i=0;i<boxht;i++){
						row[i]=sbox+i*boxwdth;
						for(j=0;j<boxwdth;j++){
							 *(row[i]+j)=1;}
						}
					}

				//if slit
				else{
					bx=dx=(xr[ind]-xl[ind])/2.;
					xmid=(xr[ind]+xl[ind])/2.;
					if(bx<0) bx=-bx;
					boxwdth=2*((int)(bx)+edge);
					by=dy=(yr[ind]-yl[ind])/2.;
					ymid=(yr[ind]+yl[ind])/2.;
					if(by<0) by=-by;
					boxht=2*((int)(by)+edge);
					slope=dx/dy;
					halfwd=boxwdth/2;
					halfht=boxht/2;
					sbox=malloc(sizeof(INT)*boxwdth*boxht);
					row=malloc(sizeof(INTP)*boxht);
					for(i=0;i<boxht;i++){
						row[i]=sbox+i*boxwdth;
						xcen=halfwd+(i-halfht)*slope;
						for(j=0;j<boxwdth;j++){
							dj=j-xcen;
							if(dj<0)dj=-dj;
							if(dj<halfwd){
								 *(row[i]+j)=1;}
							else{
								 *(row[i]+j)=0;}
							}
						}
					 }

					 //peak up on aperture
				sumax=0.;
				flmn=999999999.;
				flmx=0.;
				bflag=0;
				for(i=-srchbox;i<=srchbox;i++){
					xcorn=floor(xmid-halfwd)+i;
					for(j=-srchbox;j<=srchbox;j++){
						ycorn=floor(ymid-halfht)+j;
						sum=0.;
						//sum over offset extraction box
						for(jj=0;jj<boxht;jj++){
							iyp=ycorn+jj;
							if(iyp<0 || iyp>=naxes[1]){
							  bflag=1;
							  break;}
							for(ii=0;ii<boxwdth;ii++){
							  ixp=xcorn+ii;
							  if(ixp<0 || ixp>=naxes[0]){
								  bflag=1;
								  break;}
								flx= *(array[iyp]+ixp);
								if(flx<flmn) flmn=flx;
								if(flx>flmx)flmx=flx;
								sum+=flx*(*(row[jj]+ii));
								if(flx>flmx)flmx=flx;
								if(flx<flmn)flmn=flx;}
							if(bflag) break;}
						if(bflag) continue;
						if(sum>sumax){
							sumax=sum;
							ixtop=i;
							iytop=j;}
					}
					 }
				if(flmn<0) flmn=1;
				if((flmx/flmn<flxrat || flmx > ceiling) && nit==0){
//					printf("missed image at %d %f %f\n",chip,xmid,ymid);
					continue;}
				//Take center of gravity within (enlarged) aperture
				sumx=0;
				sumy=0;
				sum=0;
				bflag=0;
				xcorn=floor(xmid-halfwd+ixtop);
				ycorn=floor(ymid-halfht+iytop);
				for(jj=-3;jj<boxht+3;jj++){
					iyp=ycorn+jj;
					if(iyp<0 || iyp>=naxes[1]){
						bflag=1;
						break;}
					for(ii=-3;ii<boxwdth+3;ii++){
					ixp=xcorn+ii;
					if(ixp<0 || ixp>=naxes[0]){
						bflag=1;
						break;}
					flx= *(array[iyp]+ixp);
					sumx+=flx*ixp;
					sumy+=flx*iyp;
					sum+=flx;}
					if(bflag) break;}
				if(bflag) continue;
				xnew=sumx/sum;                    //xnew, ynew = observed centroid
				ynew=sumy/sum;
				xmid*=ibin;
				ymid*=ybin;
				xnew*=ibin;
				ynew*=ybin;
				ccd8(chip,xnew,ynew,&x8,&y8);   //x8,y8 = 8 chip observed centroid
				ccd8(chip,xmid,ymid,&xm8,&ym8); //xm8, ym8 = 8chip predicted centroid
				x[ntot]=xm8;//predicted
				y[ntot]=ym8;
				xs[ntot]=x8;//observed
				ys[ntot]=y8;
				xd[ntot]=x8-xm8;//observed-predicted
				yd[ntot]=y8-ym8;
				chp[ntot]=chip;
				lnd[ntot]=lind[ind];
				// fix for lack of clipping issue (4-30-08)
				sigma+=((x8-xm8)*(x8-xm8)+(y8-ym8)*(y8-ym8));
				if(mnx>x8)mnx=x8;
				if(mny>y8)mny=y8;
				if(mxx<x8)mxx=x8;
				if(mxy<y8)mxy=y8;
				if(out) fprintf(datfile,"%d %f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
										    chip,lambda[ntot],xfp[ind],yfp[ind],xmid,ymid,xnew,ynew,
                        x[ntot],y[ntot],xs[ntot],ys[ntot]);
				ntot++;
				free(row);
				free(sbox);}
				}
		  out=0;

		  for(loop=0;loop<3;loop++) {
				csort = (float *) malloc(ntot*sizeof(float));
				for(k=0;k<ntot;k++) csort[k]=xd[k];
				qsort(csort, ntot, sizeof(float), (void *)comp_nums);
				dxmed = csort[ntot/2];
				for(k=0;k<ntot;k++) csort[k]=yd[k];
				qsort(csort, ntot, sizeof(float), (void *)comp_nums);
				dymed = csort[ntot/2];
				for(k=0;k<ntot;k++) csort[k]=fabs(xd[k]-dxmed);
				qsort(csort, ntot, sizeof(float), (void *)comp_nums);
				xmad = 1.49*csort[ntot/2];
				for(k=0;k<ntot;k++) csort[k]=fabs(yd[k]-dymed);
				qsort(csort, ntot, sizeof(float), (void *)comp_nums);
				ymad = 1.49*csort[ntot/2];
 //           printf("Delta=(%.2f %.2f) Sigma=(%.2f %.2f)\n", dxmed,dymed,xmad,ymad);
//				xmad = (xmad > 1) ? xmad : 1.0;
//				ymad = (ymad > 1) ? ymad : 1.0;
				free(csort);
				for(k=0;k<ntot;k++) {
					 if ((fabs(yd[k]-dymed) > lcut*ymad) || (fabs(xd[k]-dxmed) > lcut*xmad)) {
 //                   printf("Clipping out (%.2f,%.2f)\n",x[k],y[k]);
						  for(l=k;l<ntot-1;l++) {
								x[l] = x[l+1];
								y[l] = y[l+1];
								xs[l] = xs[l+1];
								ys[l] = ys[l+1];
								xd[l] = xd[l+1];
								yd[l] = yd[l+1];
								lnd[l]=lnd[l+1];
								chp[l] = chp[l+1];}
						  ntot = ntot-1;
						  k -= 1;}
					 }
				lcut /= 1.25;}
		  sigma = 0.0;
				for(k=0;k<ntot;k++) {
					 sigma1=xd[k]*xd[k]+yd[k]*yd[k];
					 sigma+=sigma1;}


		  /*   PLOT RESIDUALS                                                      */

		  if(ntot<2) die("\ninsufficient matches to do a fit\n");
		  ddx=mxx-mnx;
		  if(mxy-mny > ddx) ddx=mxy-mny;
		  mxx+=0.1*ddx;
		  mnx-=0.1*ddx;
		  mxy+=0.1*ddx;
		  mny-=0.1*ddx;
        if(!autol || nit==0 || nit==niter){
   		  aspct=(mxy-mny)/(mxx-mnx);
   		  //wid=9.0;
   		  if(!wid){wid=5.0;}
   		  if(aspct*wid>10.){wid=10./aspct;}
   		  cpgpap(wid,aspct);
   		  cpgsch(1);
   		  cpgscr(0,1,1,1);
   		  cpgscr(1,0,0,0);
   		  cpgask(0);
   		  cpgsci(1);
   		  cpgenv(mnx,mxx,mny,mxy,0,0);
   		  cpglab("X","Y",maskfile);
   		  cpgsch(1);
   		  cpgslw(3);
   		  cpgpt(ntot,x,y,-1);
   		  for(i=0;i<ntot;i++){
   				cpgslw(10);
   				cpgpt(1,x+i,y+i,-1);
   				cpgslw(1);
   				xa[0]=x[i];
   				ya[0]=y[i];
   				if(aver>1){
   				for(j=0;j<ntot;j++){
   				sc=pow(x[i]-x[j],2)+pow(y[i]-y[j],2);
   				order(&j,&sc,irnk,rank);}
   				sc2=sc=0;
   				for(j=0;j<aver;j++){
   				sc2+=xd[irnk[j]];
   				sc+=yd[irnk[j]];}
   				xdd=sc2/aver;
   				ydd=sc/aver;}
   				else{
   				xdd=xd[i];
   				ydd=yd[i];}
   				xa[1]=x[i]-magfac*xdd;
   				ya[1]=y[i]-magfac*ydd;
   				cpgline(2,xa,ya);}
	/*		printf("\nIterate solution? ");
			fflush(stdout);
			fflush(stdin);
			fgets(answer,80,stdin);
			if(!strncasecmp(answer,"n",1)) break;*/

		  /*   CALCULATE RELATIVE SCALES                                           */

		  sigma=sqrt(sigma/ntot);
        printf("\n\n%d matches found. sigma= %6.2f pixels  ",ntot,sigma);}
//		  if(!autol || nit==niter){
        if(nit==niter) break;
        if(!autol || nit==0){
   		  printf("continue? ");
   	  	  fgets(answer,80,stdin);
   		  if(!strncasecmp(answer,"q",1)) return 0;
   	  	  if(!strncasecmp(answer,"n",1)){
   			  if(nit>0) break;
   			  return 0;}
        }
		 divid2=ceil(ntot/sqrt((float)maxpairs));
		  if(divid2<1) divid2=1;
		  n_pairs=0;
		  for(i=0;i<ntot;i+=divid2){
				for(j=i+1;j<ntot;j+=divid2){
					 sc2=(pow(x[i]-x[j],2.)+pow(y[i]-y[j],2.))/
					 (pow(xs[i]-xs[j],2.)+pow(ys[i]-ys[j],2.));
					 sc=sqrt(sc2);
					 order(&n_pairs,&sc,irnk,rank);
					 n_pairs++;}
				}
		  //take average of central 60% of values
		  np1=0.2*n_pairs;
		  np2=0.8*n_pairs;
		  scale=0.;
		  for(i=np1;i<np2;i++) scale+=rank[i];
		  scale/=(np2-np1);
        not=1;
        if(!autol){
   		  printf("\nAverage scale =          %7.5f    Reset scale? ",scale);
   		  fflush(stdout);
   		  fgets(answer,80,stdin);
   		  if(!strncasecmp(answer,"q",1)) return 0;
   		  if(!strncasecmp(answer,"n",1)){
   				scale=1;
               not=0;}
            }
		  if(autol || not){
				changes=1;
            totsc*=scale;
		  /*    RESCALE, FIND CENTROIDS                                              */

				  xavv=yavv=0.;
				  for(i=0;i<ntot;i++){
						x[i]/=scale;
	//					xs[i]/=scale;
						xd[i]=xs[i]-x[i];
						xavv+=x[i];
						y[i]/=scale;
	//					ys[i]/=scale;
						yd[i]=ys[i]-y[i];
						yavv+=y[i];}
				  xavv/=ntot;
				  yavv/=ntot;
				  if(!autol){
					  cpgeras();
					  cpgenv(mnx,mxx,mny,mxy,0,0);
					  cpglab("X","Y",maskfile);
					  cpgsch(1);
					  cpgslw(3);
					  cpgpt(ntot,x,y,-1);
					  for(i=0;i<ntot;i++){
							cpgslw(10);
							cpgpt(1,x+i,y+i,-1);
							cpgslw(1);
							xa[0]=x[i];
							ya[0]=y[i];
							if(aver>1){
							for(j=0;j<ntot;j++){
							sc=pow(x[i]-x[j],2)+pow(y[i]-y[j],2);
							order(&j,&sc,irnk,rank);}
							sc2=sc=0;
							for(j=0;j<aver;j++){
							sc2+=xd[irnk[j]];
							sc+=yd[irnk[j]];}
							xdd=sc2/aver;
							ydd=sc/aver;}
							else{
							xdd=xd[i];
							ydd=yd[i];}
							xa[1]=x[i]-magfac*xdd;
							ya[1]=y[i]-magfac*ydd;
							cpgline(2,xa,ya);}
					}
				}

		  /*      CALCULATE ROTATION ANGLE                                           */

		  /*use robust method if chance of disperser misalignment                */

		  if(misalign && isspec) {
				if(!strcmp(fitsinfo.dewarori,"NS") || NS){
					 daxis=xd;
				 	 oaxis=y;
           sign=1;
				 	 ddx=mxy-mny;}
				else{
					 daxis=yd;
				 	 oaxis=x;
           sign=-1;
				 	 ddx=mxx-mnx;}
		  //get rotation of each slit image find median

		  		 ndb=0;
				 for(j=0;j<ntot;j++){
					 if(fabs(*(oaxis+j))<0.1*ddx) continue;
					 ndb++;
					 rat=*(daxis+j)/(*(oaxis+j));
					 order(&ndb,&rat,irnk,rank);}
				 thetav=-rank[ndb/2];}

// otherwise, use general method

		else{
		  thetav=0;
		  tlist=malloc(sizeof(thetas)*ntot*ntot/2+1);
		  npr=0;
		  for(i=0;i<ntot;i+=divid2){
				x2=x[i];
				y2=y[i];
				xs2=xs[i];
				ys2=ys[i];
				for(j=i+1;j<ntot;j+=divid2){
				x1=x2-x[j];
				y1=y2-y[j];
				xs1=xs2-xs[j];
				ys1=ys2-ys[j];
				r_1=sqrt(x1*x1+y1*y1);
				theta=acos((double)x1/r_1);
				if(y1<0) theta=-theta;
				rs_1=sqrt(xs1*xs1+ys1*ys1);
				thetas=acos((double)xs1/rs_1);
				if(ys1<0) thetas=-thetas;
				thetas-=theta;
				if(thetas>pi) thetas-=2.*pi;
				if(thetas<-pi) thetas+=2.*pi;
				*(tlist+npr)=thetas;
				rank[npr]=r_1;
				npr++;}
		  }

		  /*check for angles straddling 180 degrees*/

		  nneg=npos=avneg=avpos=0;
		  for(i=0;i<npr;i++){
				if(*(tlist+i)<0){
					 nneg++;
					 avneg+=*(tlist+i);}
				else{
					 npos++;
					 avpos+=*(tlist+i);}
				}
		  if(npos>0 && nneg>0){
				avneg/=nneg;
				avpos/=npos;
				if(avneg<-pi/2. && avpos>pi/2.){
					 for(i=0;i<npr;i++){
						  if(*(tlist+i)<0) *(tlist+i)+=2.*pi;}
					 }
				}
		  dtheta=0;
		  thetas=0;
		  for(i=0;i<npr;i++){
			  thetas+=*(tlist+i)*rank[i];
				dtheta+=rank[i];}
		  thetav=-thetas/dtheta;
		  free(tlist);}

		  dtheta=thetav*57.2958;
        not=0;
		  if(!autol){
			  printf("\nAverage rotation angle = %7.3f    Reset angle? ",dtheta);
			  fflush(stdout);
			  fflush(stdin);
			  fgets(answer,80,stdin);
			  if(!strncasecmp(answer,"q",1)) return 0;
			  if(!strncasecmp(answer,"n",1)){
					thetav=0.;
               not=1;}
			  }
		  /*   ROTATE                                                                */
		  if(!not){
            totthet+=dtheta;
				changes=1;
				xv=xavv;
				yv=yavv;
				for(i=0;i<ntot;i++){
				x1=x[i]-xv;
				y1=y[i]-yv;
				r_1=sqrt(x1*x1+y1*y1);
				theta=acos((double)x1/r_1);
				if(y1<0) theta=-theta;
				theta-=thetav;
				x[i]=r_1*cos(theta)+xv;
				xd[i]=xs[i]-x[i];
				y[i]=r_1*sin(theta)+yv;
				yd[i]=ys[i]-y[i];}
				if(!autol){
					cpgeras();
					cpgenv(mnx,mxx,mny,mxy,0,0);
					cpglab("X","Y",maskfile);
					cpgsch(1);
					cpgslw(3);
					cpgpt(ntot,x,y,-1);
					for(i=0;i<ntot;i++){
						 cpgslw(10);
						 cpgpt(1,x+i,y+i,-1);
						 cpgslw(1);
						 xa[0]=x[i];
						 ya[0]=y[i];
						 if(aver>1){
						 for(j=0;j<ntot;j++){
						 sc=pow(x[i]-x[j],2)+pow(y[i]-y[j],2);
						 order(&j,&sc,irnk,rank);}
						 sc2=sc=0;
						 for(j=0;j<aver;j++){
						 sc2+=xd[irnk[j]];
						 sc+=yd[irnk[j]];}
						 xdd=sc2/aver;
						 ydd=sc/aver;}
						 else{
						 xdd=xd[i];
						 ydd=yd[i];}
						 xa[1]=x[i]-magfac*xdd;
						 ya[1]=y[i]-magfac*ydd;
						 cpgline(2,xa,ya);}
				 }
}

		  /*   CALCULATE X,Y OFFSETS                                                 */

		  deltax=0;
		  deltay=0;
		  n_pairs=0;
		  for(i=0;i<ntot;i++){
//				xxx=x[i]-xs[i];
//				yyy=y[i]-ys[i];
				order(&n_pairs,&xd[i],irnk,rank);
				order(&n_pairs,&yd[i],irnk_y,rank_y);
				n_pairs++;}
		  np1=0.2*n_pairs;
		  np2=0.8*n_pairs;
		  deltax=deltay=0.;
		  for(i=np1;i<np2;i++){
				deltax+=rank[i];
				deltay+=rank_y[i];}
		  deltax=deltax/(np2-np1);
		  deltay=deltay/(np2-np1);
        not=0;
		  if(!autol){
			  printf("\nAverage X,Y Offsets =  %6.1f %6.1f pixels. Reset positions? "
			         ,deltax,deltay);
			  fflush(stdout);
			  fflush(stdin);
			  fgets(answer,80,stdin);
			  if(!strncasecmp(answer,"q",1)) return 0;
			  if(!strncasecmp(answer,"n",1)){
					deltax=0.;
					deltay=0.;
               not=1;}
				else{changes=1;}
		     }
		  if(!not){
           totdx+=deltax;
           totdy+=deltay;
				for(i=0;i<ntot;i++){
					x[i]+=deltax;
					yd[i]-=deltay;
					y[i]+=deltay;
					xd[i]-=deltax;}
					if(!autol){
						cpgeras();
						cpgenv(mnx,mxx,mny,mxy,0,0);
						cpglab("X","Y",maskfile);
						cpgsch(1);
						cpgslw(3);
						cpgpt(ntot,x,y,-1);
						for(i=0;i<ntot;i++){
							 cpgslw(10);
							 cpgpt(1,x+i,y+i,-1);
							 cpgslw(1);
							 xa[0]=x[i];
							 ya[0]=y[i];
							 if(aver>1){
							 for(j=0;j<ntot;j++){
							 sc=pow(x[i]-x[j],2)+pow(y[i]-y[j],2);
							 order(&j,&sc,irnk,rank);}
							 sc2=sc=0;
							 for(j=0;j<aver;j++){
							 sc2+=xd[irnk[j]];
							 sc+=yd[irnk[j]];}
							 xdd=sc2/aver;
							 ydd=sc/aver;}
							 else{
							 xdd=xd[i];
							 ydd=yd[i];}
							 xa[1]=x[i]-magfac*xdd;
							 ya[1]=y[i]-magfac*ydd;
							 cpgline(2,xa,ya);}
					 }
				deltax/=dewinfo.scale;
				deltay/=dewinfo.scale;}
		  th=theta*57.2958;
		  //printf("%10.5f %10.5f %10.5f %10.5f\n",th, scal,xxx,yyy);

		  /*DETERMINE DISPERSER MISALIGNMENT                               */

		  if(misalign && isspec) {
				if(!strcmp(fitsinfo.dewarori,"NS") || NS){
					 daxis=yd;}
				else{
					 daxis=xd;}
		  //get mean shift for each line
				for(i=0;i<nlines;i++){
					 dav=0;
					 ndb=0;
					 for(j=0;j<ntot;j++){
						  if(lnd[j]==i){
								dav+=*(daxis+j);
								ndb++;}
						  }
					 dalign[i]= (ndb > 0) ? dav/ndb : 0;}

/*				cpgpap(wid,0.6);
				cpgsch(1);
				cpgscr(0,1,1,1);
				cpgscr(1,0,0,0);
				cpgask(0);
				cpgsci(1);
				cpgenv(lmin-100.,lmax+100.,-10.,10.,0,0);
				cpglab("Lambda","Delta",maskfile);
				cpgsch(1);
				cpgslw(3);
				cpgpt(nlines,lines,dalign,25);
				//fit relation*/
				polyfit(lines,dalign,nlines,1,relat);
				//printf("%f %f\n",relat[0],relat[1]);
/*				xa[1]=lines[nlines-1];
				ya[0]=relat[0]+relat[1]*lines[0];
				ya[1]=relat[0]+relat[1]*lines[nlines-1];
				cpgline(2,xa,ya);*/
//				delal = strcmp(obsdata.instrument,"LDSS3") ? -141.*relat[1] : -70*relat[1];
        delal=-141*relat[1];
  			if(!autol){
					printf("Disperser rotation = %f. Apply rotation? ",delal);
					fflush(stdout);
					fflush(stdin);
					fgets(answer,80,stdin);
					if(!strncasecmp(answer,"q",1)) return 0;
					if(!strncasecmp(answer,"n",1)){
						 delal=0;}
					else{
						 changes=1;}
                }
            totmis+=delal;}
		  if(changes){
   		  GetCof(&theta,&scal,&xxx,&yyy);
   		  theta=(theta-thetav*rsign);
   		  scal*=scale;
   		  xxx-=deltax;
   		  yyy-=deltay;
   		  obsdata.alignrot+=delal;
   		  SetCof(theta,scal,xxx,yyy);
   		  cpgpage();}
        else{
           break;}
     }

	  /* CORRECT COF FILE AND REWRITE                                            */

	 printf("Total corrections:\n\nscale = %f\nrotation = %f\ndelta(x) = %f\ndelta(y) = %f\ndisperser rotation = %f\n",
           totsc,totthet,totdx,totdy,totmis);
    printf("Save changes? ");
	 fflush(stdout);
	 fflush(stdin);
	 fgets(answer,80,stdin);
	 if(!strncasecmp(answer,"n",1)) return 0;

	 strcpy(file,obsdata.dewoff);
	 strcat(file,".dewoff");
	 infile=fopen(file,"r");
    if((infile=fopen(file,"r"))==NULL){
        //look in sdata/dewoff
        COS_HOME=malloc(sizeof(CHAR)*80);
        COS_HOME=getenv("COSMOS_HOME");
        if(COS_HOME==NULL){
            printf("COSMOS_HOME undefined!\n");
            return 1;}
        strcpy(lane,COS_HOME);
        strcat(lane,"/sdata/dewoff/");
        strcat(lane,file);
        if((infile=fopen(lane,"r"))==NULL){
          printf("cant find %s\n",lane);
          return 1;}
        }


	 //if image file and dewoff file have same name, create backup copy

	 if(!(strcmp(obsdata.dewoff,dfile))){
		  strcat(file,"%");
		  outfile=fopen(file,"w");
		  while(1){
		  if(!fgets(line,133,infile)) break;
		  fprintf(outfile,"%s",line);}
		  fclose(outfile);
		  fclose(infile);
		  infile=fopen(file,"r");}

	 strcpy(file,dfile);
	 strcat(file,".dewoff");
	 outfile=fopen(file,"w");
	 if(outfile==NULL){
		  printf("Unable to create new dewar offset file %s\n",file);
		  return 1;}
	 printf("Creating new dewar offset file %s\n",file);
	 theta*=57.2958;
	 fprintf(outfile,"%8.4f %8.6f %9.4f %9.4f %d\n",theta,scal,xxx,yyy,ypar);
	 fclose(infile);
	 fclose(outfile);

	 //update original obsdef file with name of new dewar offset file
	 strcpy(file,maskfile);
	 strcat(file,".obsdef");
	 printf("Updating %s\n",file);
	 strcpy(oldfile,file);
	 strcat(oldfile,"%");
	 strcpy(command,"mv -f ");
	 strcat(command,file);
	 strcat(command," ");
	 strcat(command,oldfile);
	 system(command);

	 infile=fopen(oldfile,"r");
	 outfile=fopen(file,"w");
	 while(fgets(line,133,infile)){
		  sscanf(line,"%s",parm);
		  if(!strcmp(parm,"DEWOFF")){
				fprintf(outfile,"DEWOFF          %s\n",dfile);}
		  else if(!strcmp(parm,"D_ALIGNROT")){
				fprintf(outfile,"D_ALIGNROT      %f\n",obsdata.alignrot);}
		  else{
				fprintf(outfile,"%s",line);}
		  }
  return  0;}

int comp_nums(const float *num1, const float *num2)
{
	if (*num1 <  *num2) return -1;
	if (*num1 == *num2) return  0;
	if (*num1 >  *num2) return  1;
  return 0;
}
