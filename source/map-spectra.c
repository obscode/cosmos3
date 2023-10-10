/*****************************************************************************\
*
*  MAP_SPECTRA  produces parameter files to be used to map spectra produced
*               by a mask/observation combination onto the ccd's. Parameter
*               files are input to sky subtraction, wavelength fitting, and
*               spectrum extraction programs
*
*   VERSION  5 Jan 2018
*
*   USAGE:   map-spectra obsname
\*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ioutils.h"
#include "kdcutil.h"
#include "mgutils.h"
#include "optutils.h"
#include "cosmos.h"
#include "clardy.h"

int  main(int argc,char *argv[]){

  Obs       *obset = NULL;
  element   *optics, *instrument, *p;
  vect3     mskpos,campos,cenpos,objpos;
  vect2     mpl,mpleft,mpright,mpdif,mpcen,ldif,rdif,mpdelt;
  double    cur_wavl, cur_temp;
  int       i,j,gr_order,chip,curslit,minpixl,nchip,xchip,ychip,xdisper,npxl,
            npixels,npx,nbd,lastchip,chip1,order,ind,dorder,corder,*ip,
            ord_disper,nbad,mindist,ord_sag,nchp,ord_tilt,ord_sagit,ord_slen,
            maxord,mesage,fsln,nval,badd,badnum[1000],useaps,useprtl,prtl,
            chipnum,leftprtl,rightprtl,bdnum[1000],chip0;
  float     gr_angle,xccd,yccd,xfp,yfp,xlfp,ylfp,lambda,lambda0,lambda1,
            dlambda1,telscale,temp,slitobj,slittop,slitbot,sltp,curvx,meanl,
            xccd1,localscale,minrat,yccd1,disper,camscale,slitleft,slitright,
            lamin,lamax,mplft,mprt,FLOAT,lamstart,coefs[12],slitmin,slitmax,
            slitmean,leftpos,rightpos,midpos,centr,rat,edge,lamstop,xrcd,yrcd,
            xlcd,ylcd,slitln,xl,xr,yl,yr,stdev,dmean,xrfp,yrfp,xmin,xmax,ymin,
            ymax;
  float     *xval,*yval,*lamval,*tilt,*sagit,*slitlen,*frctn;
  double    *xlist,*ylist,**elist,maskang,DOUBLE;
  char       **badname,**bdname,*CHARP;
  char       cmra,flnm[80],file[80],maskfile[80],INSTRMENT[80],GRATNG[80],
             CHAR,SMF_FILE[80],camdef[80],camoff[80],imacsdir[80];
  slit       slitdat;
  objq       *oq, *oq1;
  obsdef     obsdata;
  objdat     objdata;
  dewdat     dewinfo;
  FILE       *infile,*outfile,*curvfile;

  float   pixint=3.;//fitting interval in pixels
  telscale=2.899;
  nbd=nbad=0;
  mindist=3;//minimum distance of object from chip edge
  //order of fits
  ord_disper=7;
  ord_sag=5;
  ord_tilt=5;
  ord_sagit=5;
  ord_slen=5;

  /*------------------------get observ data----------------------------------*/

  //prompt for parameters
  if(argc<2){
    while(1){
      printf("Enter observation name: ");
      scanf("%s",maskfile);
      if(ReadObsDef(maskfile,&obsdata)==0) break;
      printf("Cannot open observation definition file!\n");}
    }
  //obsdef file on commancd line
  else{
    sscanf(argv[1],"%s",maskfile);
    if(ReadObsDef(maskfile,&obsdata)!=0)
      {printf("%Cannot open observation definition file!\n");}
    }

  //parameters
  strcpy(file,"map-spectra");
  if(OpenCosParm(file)!=0) die("Cannot open map-spectra parameter file");
  if(ReadParm_r("MINLAMBDA",&lambda0)==1) die("parameter file error");
  if(ReadParm_r("MAXLAMBDA",&lambda1)==1)die("parameter file error");
  if(ReadParm_i("MINPIXL",&minpixl)==1) die("parameter file error");
  minpixl/=pixint;
  if(ReadParm_b("USE_HOLES",&useaps)==1) die("parameter file error");
  if(ReadParm_b("USE_PRTLS",&useprtl)==1) useprtl=0;

  /*------------------------- get mask data--------------------------------- */
  strcpy(flnm,obsdata.mask);
  strcpy(SMF_FILE,flnm);
  strcat(SMF_FILE,".SMF");
  if(ReadSMFfile(SMF_FILE,&obset)) return 1;
  cur_temp=obset->temp;

  //read obsservation data
  if(!strcmp(obsdata.mode,"DIRECT")){
    printf("map-spectra requires a spectroscopic image\n");
    return 1;}
  SetupInstr(&obsdata,&instrument);
  if(SetupCamera(&obsdata)==1) return 1;


  /*----------- get CCD parmaeters and mapping -----------------------------*/
   Getchipdat(&dewinfo);
   chipnum=dewinfo.nchip;

  /*--------------------------setup  arrays_______________________________*/
  maxord=ord_disper;
  maxord = (ord_sag>maxord) ? ord_sag: maxord;
  maxord = (ord_sagit>maxord) ? ord_sagit: maxord;
  maxord = (ord_tilt>maxord) ? ord_tilt: maxord;
  maxord+=1;
  xlist=malloc(sizeof(DOUBLE)*2*maxord);
  ylist=malloc(sizeof(DOUBLE)*maxord);
  elist=malloc(sizeof(xlist)*maxord);
  for(i=0;i<maxord;i++){
    *(elist+i)=malloc(sizeof(DOUBLE)*(maxord+1));}
  //  nval=(lambda1-lambda0)/disper+10;
  nval=10000;
  xval=malloc(sizeof(FLOAT)*nval);
  yval=malloc(sizeof(FLOAT)*nval);
  tilt=malloc(sizeof(FLOAT)*nval);
  sagit=malloc(sizeof(FLOAT)*nval);
  lamval=malloc(sizeof(FLOAT)*nval);
  slitlen=malloc(sizeof(FLOAT)*nval);
  frctn=malloc(sizeof(FLOAT)*nval);
  badname=malloc(sizeof(CHARP)*1000);
  for(i=0;i<1000;i++) *(badname+i)=malloc(sizeof(CHAR)*20);
  bdname=malloc(sizeof(CHARP)*1000);
  for(i=0;i<1000;i++) *(bdname+i)=malloc(sizeof(CHAR)*20);
  /*----------------------determine slit scale, dispersion___________________*/
  oq=obset->ob;
  oq1=oq;
  meanl=(lambda0+lambda1)/2.;
  dmean=(lambda1/lambda0)/10.;
  while(1){
    mskpos = get3de2v(oq->smpos,0.0);
    cur_wavl=meanl;
    campos = Op_transform(mskpos,instrument,cur_wavl,cur_temp);
    xfp=campos.x;
    yfp=campos.y;
    chip1=fp2ccd(xfp,yfp,&xccd1,&yccd1,cur_wavl);
    if(!chip1) goto C;
    cur_wavl+=1.;
    campos = Op_transform(mskpos,instrument,cur_wavl,cur_temp);
    xfp=campos.x;
    yfp=campos.y;
    chip=fp2ccd(xfp,yfp,&xccd,&yccd,cur_wavl);
    if(chip != chip1) goto C;
    //along which axis does dispersion run?
    xdisper=0;
    if(fabs((double)((xccd-xccd1)/(yccd-yccd1))) > 1.) xdisper=1;
    if(xdisper) disper=fabs(pixint/(double)(xccd-xccd1));
    else disper=fabs(pixint/(double)(yccd-yccd1));
    //what is camera demagnification,handedness?
    mskpos.x+=1;
    mskpos.y+=1.;
    campos = Op_transform(mskpos,instrument,cur_wavl,cur_temp);
    xfp=campos.x;
    yfp=campos.y;
    chip1=fp2ccd(xfp,yfp,&xccd1,&yccd1,cur_wavl);
    if(chip != chip1) goto C;
    if(xdisper) camscale=fabs((double)(yccd-yccd1));
    else camscale=fabs((double)(xccd-xccd1));
    if(!strcasecmp(obsdata.camera,"LONG")) camscale=-camscale;
    break;
C:  oq=oq->next;
    if(oq==oq1){                            //haven't found a line within image
      meanl+=dmean;                         //change wavelength
      if(meanl>lambda0 && meanl<lambda1) continue;
      if(meanl>lambda1){          //hit upper wavelength, go in other direction
        dmean=-dmean;
        meanl=(lambda0+lambda1)/2.+dmean;
        continue;}
      die("Cannot map slits");}
  }
  strcpy(file,maskfile);
  strcat(file,".map");
  outfile=fopen(file,"w");
  fprintf(outfile,"Xdispersion =  %d\nFit orders = %d %d %d %d %d\n",xdisper,
          ord_disper,ord_sag,ord_tilt,ord_sagit,ord_slen);
  fprintf(outfile,"Scale ~ %10.4f Camscale = %8.3f\nLambda  = %7.1f %7.1f\n",telscale,camscale,lambda0,
          lambda1);
//  fprintf(outfile,"Dewar = %s %d\n\n",obsdata.dewar,dewinfo.nchip);
  fprintf(outfile,"Dewar = %s %d\n\n",obsdata.dewar,chipnum);

  /*------------------------ map slits --------------------------------------*/

  oq=oq1;
  curslit=1;
  mesage=0;

  //loop over all slits
  while(1){
    badd=0;
    printf("Processing slit %d\r",curslit);
    nchp=0;
    slitdat=oq->slit;
    objdata=oq->dat;
    if(slitdat.shape !=2 && !useaps){
      oq=oq->next;
      mesage=0;
      printf("skipping aperture before slit %d\n",curslit);
      if(oq==oq1) break;
      continue;}
    fflush(stdout);
    strcpy(*(bdname+nbd),objdata.name);
    bdnum[nbd]=curslit;
    strcpy(*(badname+nbad),objdata.name);
    badnum[nbad]=curslit;
    objpos=get3de2v(oq->smpos,0.0); //object position, mask coordinates
    mpleft=sum2vect(oq->smpos,lslit(oq->slit));
	mpright=sum2vect(oq->smpos,rslit(oq->slit));
    ldif=lslit(oq->slit);
    rdif=rslit(oq->slit);
    mpdif=sum2vect(ldif,rdif);
    mpdelt=mul2vect(mpdif,0.5);
	mpcen=sum2vect(oq->smpos,mpdelt);
	cenpos=get3de2v(mpcen,0.0); //center position, mask coordinates
    prtl=leftprtl=rightprtl=0;

    //find partial and bad slits, determine maximum included length of partial slits

    lambda=lambda0-disper;;
    fsln=0;
    minrat=10000.;
    npixels=0;
    while(lambda<lambda1-disper){
      lambda+=disper;
      campos = Op_transform(objpos,instrument,lambda,cur_temp);  //object pos, camera coord
      xfp=campos.x;                                                //object x, camera coord
      yfp=campos.y;                                                //object y, camera coord
      chip=fp2ccd(xfp,yfp,&xccd,&yccd,lambda); //xccd,yccd = object position
      if(chip<1) continue;//object not on a chip
      if(xccd<mindist || xccd>dewinfo.xchip-mindist-1 || yccd<mindist || yccd>dewinfo.ychip-mindist-1)
         continue;//object too close to edge
      //get ends of slit
      //"left"=minvalue slit end
      xl=mpleft.x;
      yl=mpleft.y;
      mskpos = get3de2v(mpleft,0.0);
      campos = Op_transform(mskpos,instrument,lambda,cur_temp);
      xlfp=campos.x;                                          //left end of slit, camera coord
      ylfp=campos.y;                                          // "    "   "   "      "     "
      chip0=fp2ccd(xlfp,ylfp,&xlcd,&ylcd,lambda);
      //"right"=maxvalue slit end
      xr=mpright.x;
      yr=mpright.y;
      mskpos = get3de2v(mpright,0.0);
      campos = Op_transform(mskpos,instrument,lambda,cur_temp);
      xrfp=campos.x;
      yrfp=campos.y;
      chip1=fp2ccd(xrfp,yrfp,&xrcd,&yrcd,lambda);
      if(chip0==0 || chip1==0) continue; //end too far off chip
      if((chip0<1 || chip0 != chip) && (chip1<1 || chip1 != chip)){
          continue;} //both ends of curved slit off chip
      if(chip0<0){
        if(xdisper){
           if(xlcd<0 || xlcd>dewinfo.xchip-1) continue; //slit end falls off END of chip
           if (yccd<0.5*dewinfo.ychip-1) edge=1;
           else edge=dewinfo.ychip-2;
           rat=fabs((yccd-edge)/(yccd-yrcd));}
        else{
          if(ylcd<0 || ylcd>dewinfo.ychip-1) continue;
          if(xccd<0.5*dewinfo.xchip-1) edge=1;
          else edge=dewinfo.xchip-2;
          rat=fabs((xccd-edge)/(xccd-xrcd));}
        prtl=leftprtl=1;
        if(!useprtl){
          break;}
        if(rat<minrat) minrat=rat;}
        if(chip1<0){
          if(xdisper){
             if(xrcd<0 || xrcd>dewinfo.xchip-1) continue;
             if (yccd<0.5*dewinfo.ychip-1) edge=1;
             else edge=dewinfo.ychip-2;
             rat=fabs((yccd-edge)/(yccd-ylcd));
          }
          else{
            if(yrcd<0 || yrcd>dewinfo.ychip-1) continue;
            if(xccd<0.5*dewinfo.xchip-1) edge=1;
            else edge=dewinfo.xchip-2;
            rat=fabs((xccd-edge)/(xccd-xlcd));}
          prtl=rightprtl=1;
          if(!useprtl){
            break;}
          if(rat<minrat) minrat=rat;}
          npixels++;}

      //adjust edge
      if(npixels<minpixl){
         badd=1;
         nbad++;}
      else if(prtl){
        nbd++;
        if(!useprtl){
          badd=1;
         }
        else{
          if(leftprtl){
            slitdat.alen=minrat*slitdat.blen;}
          else{
            slitdat.blen=minrat*slitdat.alen;}
          ldif=lslit(slitdat);
          rdif=rslit(slitdat);
          mpleft=sum2vect(oq->smpos,ldif);
          mpright=sum2vect(oq->smpos,rdif);


          mpdif=sum2vect(ldif,rdif);
          mpdelt=mul2vect(mpdif,0.5);
          mpcen=sum2vect(oq->smpos,mpdelt);
          cenpos=get3de2v(mpcen,0.0);}
          }

      lambda=lambda0-disper;;
      fsln=0;
      while(lambda<lambda1){
        if(badd) break;
        npxl=0;
        lastchip=0;
        xmax=ymax=0;
        xmin=ymin=5000;
        //loop while on one chip
        while(1){
          lambda+=disper;
          if(lambda>lambda1) break;
          campos = Op_transform(objpos,instrument,lambda,cur_temp);  //object pos, camera coord
          xfp=campos.x;                                                //object x, camera coord
          yfp=campos.y;                                                //object y, camera coord
          chip=fp2ccd(xfp,yfp,&xccd,&yccd,lambda); //xccd,yccd = object position
          if(chip<1 || (chip!=lastchip && npxl)) break;
          if(xccd<mindist || xccd>dewinfo.xchip-mindist-1 || yccd<mindist || yccd>dewinfo.ychip-mindist-1)
             break;//object too close to edge
          if(!npxl) lamin=lambda;
          slitobj= (xdisper) ? yccd: xccd;
          //get ends of slit
          //"left"=minvalue slit end
          xl=mpleft.x;
          yl=mpleft.y;
          mskpos = get3de2v(mpleft,0.0);
          campos = Op_transform(mskpos,instrument,lambda,cur_temp);
          xlfp=campos.x;                                          //left end of slit, camera coord
          ylfp=campos.y;                                          // "    "   "   "      "     "
          chip0=fp2ccd(xlfp,ylfp,&xlcd,&ylcd,lambda);
          if(chip0==0) continue;
          //does spectrum fall off left edge of chip?
          if(chip0<=0){
            if((xdisper && (xlcd<0 || xlcd>dewinfo.xchip-1)) ||
            (!xdisper && (ylcd<0 || ylcd>dewinfo.ychip-1)))
                 break;
                 printf("mapping failure\nleft %d %f %f\n", curslit,lambda);
                 printf("%d %f %f %f %f\n",chip,xfp,yfp,xccd,yccd);
            printf("%d %f %f %f %f\n",chip0,xlfp,ylfp,xlcd,ylcd);
          return 1;}
          if(xlcd<xmin) xmin=xlcd;
          if(xlcd>xmax) xmax=xlcd;
          if(ylcd<ymin) ymin=ylcd;
          if(ylcd>ymax) ymax=ylcd;
          slitmin = (xdisper) ? ylcd-slitobj: xlcd-slitobj;
          //"right"=maxvalue slit end
          xr=mpright.x;
          yr=mpright.y;
          mskpos = get3de2v(mpright,0.0);
          campos = Op_transform(mskpos,instrument,lambda,cur_temp);
          xrfp=campos.x;
          yrfp=campos.y;
          chip1=fp2ccd(xrfp,yrfp,&xrcd,&yrcd,lambda);
          if(chip1==0) continue; //end too far off chip
          if((chip0<1 || chip0 != chip) && (chip1<1 || chip1 != chip)){
              continue;} //both ends of curved slit off chip
          //does spectrum fall off right end of chip?
          if(chip1<0){
            if(xdisper && (xrcd<0 || xrcd>dewinfo.xchip-1)) break;
            if(!xdisper && (yrcd<0 || yrcd>dewinfo.ychip-1)) break;
            printf("mapping failure\nright%d %f\n", curslit,lambda);
            printf("%d %f %f %f %f\n",chip,xfp,yfp,xccd,yccd);
            printf("%d %f %f %f %f\n",chip1,xrfp,yrfp,xrcd,yrcd);
            return 1;}
        if(xrcd<xmin) xmin=xrcd;
        if(xrcd>xmax) xmax=xrcd;
        if(yrcd<ymin) ymin=yrcd;
        if(yrcd>ymax) ymax=yrcd;
        slitmax = (xdisper) ? yrcd-slitobj: xrcd-slitobj;
         slitln=sqrt((xl-xr)*(xl-xr)+(yl-yr)*(yl-yr));

	  //spectrum x,y location,y slitlength vs lambda
	if(xdisper){
	  yval[npxl]=yccd;
	  xval[npxl]=xccd;
    slitlen[npxl]=yrcd-ylcd;}
  else{
    yval[npxl]=xccd;
    xval[npxl]=yccd;
    slitlen[npxl]=xrcd-xlcd;}
  lamval[npxl]=lambda;

	//curvature calculation

  //tilt & curvature

  campos = Op_transform(cenpos,instrument,lambda,cur_temp);
  xfp=campos.x;
  yfp=campos.y;
  chip=fp2ccd(xfp,yfp,&xccd,&yccd,lambda);
  if(xdisper){
    centr=xccd;
    leftpos=xlcd;
    rightpos=xrcd;}
  else{
    centr=yccd;
    leftpos=ylcd;
    rightpos=yrcd;}
  tilt[npxl]=rightpos-leftpos;
  midpos=0.5*(rightpos+leftpos);
  sagit[npxl]=centr-midpos;
  lambda+=disper;
  lastchip=chip;
  npxl++;}
      if(!lastchip) continue;
      //hit end of chip
      if(npxl < minpixl)continue;

      /*------------fit and write data for spectrum segment------------------*/

      lamax=lambda-disper;
      if(!fsln){
        xfp=oq->smpos.x;
        yfp=oq->smpos.y;
        fprintf(outfile,"SLIT %d %s %d\n",curslit,objdata.name,slitdat.shape);
        fprintf(outfile,"LENGTH = %f POS = %8.3f %8.3f\n",slitln,xfp,yfp);
          fsln=1;}
      if(rightprtl) prtl=-1;
      fprintf(outfile,"CHIP %d %7.5f %7.5f %8.2f %8.2f %d %d %7.1f %7.1f %7.1f %7.1f\n",
                lastchip,slitmin,slitmax,lamin,lamax,npxl,prtl,xmin,xmax,ymin,ymax);
      nchp++;
      //dispersion solution
      plyfit(lamval,xval,npxl,ord_disper,coefs,xlist,ylist,elist);
      for(i=0;i<=ord_disper;i++) fprintf(outfile,"%13.6e ",coefs[i]);
      fprintf(outfile,"\n");
      //inverse dispersion solution
      plyfit(xval,lamval,npxl,ord_disper,coefs,xlist,ylist,elist);
      for(i=0;i<=ord_disper;i++) fprintf(outfile,"%13.6e ",coefs[i]);
      fprintf(outfile,"\n");
      //spectrum sag
      plyfit(lamval,yval,npxl,ord_sag,coefs,xlist,ylist,elist);
      for(i=0;i<=ord_sag;i++) fprintf(outfile,"%13.6e ",coefs[i]);
      fprintf(outfile,"\n");
      //spectrum tilt
      plyfit(lamval,tilt,npxl,ord_tilt,coefs,xlist,ylist,elist);
      for(i=0;i<=ord_tilt;i++) fprintf(outfile,"%13.6e ",coefs[i]);
      fprintf(outfile,"\n");
      //spectrum saggitta
      plyfit(lamval,sagit,npxl,ord_sagit,coefs,xlist,ylist,elist);
      for(i=0;i<=ord_sagit;i++) fprintf(outfile,"%13.6e ",coefs[i]);
      fprintf(outfile,"\n");
      //slitlength
      plyfit(lamval,slitlen,npxl,ord_slen,coefs,xlist,ylist,elist);
      for(i=0;i<=ord_slen;i++) fprintf(outfile,"%13.6e ",coefs[i]);
      fprintf(outfile,"\n");
      }
    oq=oq->next;
    mesage=0;
    if(oq==oq1) break;
    curslit++;}
  fprintf(outfile,"END\n\n");

  if(nbd){
    if(useprtl){
      printf("\nThe following partial slits were mapped:\n\n");}
    else{
      printf("\nThe following partial slits were not mapped:\n\n");}
    for(i=0;i<nbd;i++){
      printf("Slit %3d   %s\n",bdnum[i],*(bdname+i));}
    }
  if(nbad){
    printf("\nThe following bad slits were not mapped:\n\n");
    for(i=0;i<nbad;i++){
      printf("Slit %3d   %s\n",badnum[i],*(badname+i));}
    }


  //free memory
  for(i=0;i<maxord;i++) free(*(elist+i));
  free(elist);
  free(xlist);
  free(ylist);
  free(xval);
  free(yval);
  free(tilt);
  free(sagit);
  free(lamval);
  free(slitlen);
  free(frctn);
  for(i=0;i<1000;i++) {
    free(*(badname+i));
    free(*(bdname+i));}
  free(badname);
  free(bdname);
   return 0;}
