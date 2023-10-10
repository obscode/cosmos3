/*****************************************************************************
*
*   FP2CCD uses the data provided in a CAMERA DEFINITION FILE,  a
*                 DEWAR OFFSET FILE, and a CAMERA DISTORTION FILE to transform
*                 coordinates from the spectrograph focal plane to pixel locations
*                 on the ccd camera. Returns chip # and on-chip coordinates
*
*   READCDF reads in the camera definition file data and must be called
*                     before the first call to fp2ccd
*
*   READDISTOR reads in the camera distortion data and must be called
*                           before the first call to fp2ccd
*
*   READCOF reads in the dewar offset  file data and must be called
*                      before the first call to fp2ccd
*
*   [GET/SET]COF  retrieve/reset dewar offset parameters
*
*   SETUPCAMERA inputs obsdef info for possible using in mapping corrections
*                               and calls READCDF, READDISTOR, and READCOF
*
*   GETCHIPDAT returns chip data (number,size)
*
*   CCD8 puts chip coordinates on an 8-chip scale
*
*   MOSPOS puts chip coordinates on mosaic scale
*
*   ZERNIKE takes list of zernike polynomials and their values and calculates
*                   corrections to R, THETA for input R, THETA
*
*   ZERNIKE_XY does the same, returning corrections in X,Y for input X,Y
*
*   VERSION 26 Nov 2016
*
\*****************************************************************************/

#include  <stdlib.h>
#include  <string.h>
#include  <stdio.h>
#include   <math.h>
#include  "cosmos.h"

static int    nchip,xchip,ychip,xz[8],sxz[8],yz[8],syz[8],sl_n,sts1_n,stt1_n,
              i,oldsite,sts2_n,stt2_n,srs1_n,srt1_n,srs2_n,srt2_n,srz_n,xarsize,
              yarsize,hand,xxn[8],yyn[8],distor,oldway,nrzern,zrterm[40],ndone,
              ntzern,ztterm[40];
static float  scale,xmin[8],xmax[8],ymin[8],ymax[8],dtheta[8],xx0[8],yy0[8],
              dcth,dcx,dcy,dscale,sts1_coef[4],sl_coef[4],zrval[40],ztval[40],
              dt,stt1_coef[4],sts2_coef[4],stt2_coef[4],srs1_coef[4],misal,
              drr,srt1_coef[4],srs2_coef[4],srt2_coef[4],srz_coef[4];
static char   distfunc[80],dth;
static obsdef obsinfo;

//
// SetupCamera   calls Readcdf and Readcof and inputs instrument setup infor
// for possible use in correcting mapping
//

int SetupCamera(obsdef* obsdat){

    obsinfo=*obsdat;
    //printf("dewar %s\n",obsinfo.dewar);
    if(Readcdf(obsinfo.dewar)==1){
        printf("Error reading dewar definition file %s\n",obsinfo.dewar);
        return 1;}
    if(strstr(obsinfo.instrument,"IMACS") != NULL && !strcasecmp(obsinfo.camera,"SHORT")){
        if(ReadDistor(obsinfo.distor)==1){
            printf("Error reading camera distortion file %s for %s\n",
            obsinfo.distor,obsinfo.camera);
            return 1;}
        }
    if(Readcof(obsinfo.dewoff)==1){
        printf("Error reading dewar offset file %s\n",obsinfo.dewoff);
        return 1;}
    return 0;}

//
// Readcdf
//

int Readcdf(char filename[]){

    FILE     *parmfile;
    char     CHAR,line[133],imacsdir[130];
    int      i,nopar;
    char     *COS_HOME;

    COS_HOME=malloc(sizeof(CHAR)*80);
    COS_HOME=getenv("COSMOS_HOME");
    if(COS_HOME==NULL){
        printf("COSMOS_HOME undefined!\n");
        return 1;}
    strcpy(imacsdir,COS_HOME);
    strcat(imacsdir,"/sdata/dewdef/");
    strcat(imacsdir,filename);
    strcat(imacsdir,".dewdef");
    dscale=1.;
    if((parmfile=fopen(imacsdir,"r"))==NULL){
        printf("Cannot find dewar definition file\n");
        return 1;}
    for(;;){
        if(!fgets(line,133,parmfile)){
            fclose(parmfile);
            printf("Error 1 reading dewar definition file\n");
            return 1;}
    if(*line=='#') continue;
    i=sscanf(line,"%d %d %d %f %d %d %d",&nchip,&xchip,&ychip,&scale,&xarsize,
                 &yarsize,&hand);
    //handling legacy IMACS SITE.dewdef files
    if(i==4){
        printf("\nError: you are using the old SITE.dewdef file\n\n");
        return 1;}

    xchip-=1;
    ychip-=1;
    break;}
    for(i=0;i<nchip;i++){
        for(;;){
            if(!fgets(line,133,parmfile)){
	            fclose(parmfile);
	            return 1;}
            if(*line=='#') continue;
            sscanf(line,"%f %f %f %f %d %d",&xmin[i],&xmax[i],&ymin[i],&ymax[i],
                      &xxn[i],&yyn[i]);
            break;}
        for(;;){
            if(!fgets(line,133,parmfile)){
                fclose(parmfile);
 	            return 1;}
            if(*line=='#') continue;
            sscanf(line,"%d %d %d %d",&xz[i],&sxz[i],&yz[i],&syz[i]);
            break;}
        for(;;){
            if(!fgets(line,133,parmfile)){
                fclose(parmfile);
            	return 1;}
            if(*line=='#') continue;
            sscanf(line,"%f",&dtheta[i]);
            dtheta[i]/=57.2958;
            break;}
        for(;;){
            if(!fgets(line,133,parmfile)){
	            fclose(parmfile);
                 return 1;}
            if(*line=='#') continue;
            sscanf(line,"%f %f",&xx0[i],&yy0[i]);
            break;}
        }
    return 0;}

//
//ReadDistor
//

int ReadDistor(char filename[]){

    FILE     *parmfile;
    char     CHAR,line[133],imacsdir[130];
    int      i,nopar;
    char     *COS_HOME;

    COS_HOME=malloc(sizeof(CHAR)*80);
    COS_HOME=getenv("COSMOS_HOME");
    if(COS_HOME==NULL){
        printf("COSMOS_HOME undefined!\n");
        return 1;}
    strcpy(imacsdir,COS_HOME);
    strcat(imacsdir,"/sdata/distor/");
    strcat(imacsdir,filename);
    strcat(imacsdir,".distor");
    dscale=1.;
    if((parmfile=fopen(imacsdir,"r"))==NULL){
        printf("Cannot find dewar distortion file\n");
        return 1;}
    distor=oldway=0;
    sl_n=nrzern=ntzern=0;

    //read lambda corrections

    for(;;){
         if(!fgets(line,133,parmfile)){
            printf("No wavelength or distortion corrections applied\n:");
            return 0;}
         if(*line=='#') continue;
        sscanf(line,"%d %e %e %e %e",&sl_n,&sl_coef[0],&sl_coef[1],&sl_coef[2],
                  &sl_coef[3]);
        break;}
    if(sl_n==0) printf("No wavelength corrections applied\n:");

    //read radial distortion corrections

    for(;;){
        if(!fgets(line,133,parmfile)){
            printf("No distortion corrections applied\n");
            return 0;}
        if(*line=='#')continue;
        ndone=0;
        i=sscanf(line,"%d %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f",&nrzern,
                    &zrterm[ndone],&zrval[ndone],&zrterm[ndone+1],&zrval[ndone+1],
                    &zrterm[ndone+2],&zrval[ndone+2],&zrterm[ndone+3],&zrval[ndone+3],
                    &zrterm[ndone+4],&zrval[ndone+4],&zrterm[ndone+5],&zrval[ndone+5],
                    &zrterm[ndone+6],&zrval[ndone+6],&zrterm[ndone+7],&zrval[ndone+7],
                    &zrterm[ndone+8],&zrval[ndone+8]);
        break;}
    ndone+=(i-1)/2;
    if(ndone<nrzern){
        while(1){
            if(!fgets(line,133,parmfile)) break;
            i=sscanf(line,"%d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f",
                        &zrterm[ndone],&zrval[ndone],&zrterm[ndone+1],&zrval[ndone+1],
                        &zrterm[ndone+2],&zrval[ndone+2],&zrterm[ndone+3],&zrval[ndone+3],
                        &zrterm[ndone+4],&zrval[ndone+4],&zrterm[ndone+5],&zrval[ndone+5],
                        &zrterm[ndone+6],&zrval[ndone+6],&zrterm[ndone+7],&zrval[ndone+7],
                        &zrterm[ndone+8],&zrval[ndone+8],&zrterm[ndone+9],&zrval[ndone+9]);
            ndone+=i/2;
            if(ndone==nrzern) break;}
        }
    if(nrzern==0) printf("No radial distortion corrections applied\n");
    if(ndone<nrzern){
        printf("Error reading distortion file: missing terms\n");
        return 1;}

    //read angular distortion corrections

    for(;;){
        if(!fgets(line,133,parmfile)) {
            printf("No angular distortion corrections applied\n");
            return 0;}
        if(*line=='#') continue;
        ndone=0;
        i=sscanf(line,"%d %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f",&ntzern,
                    &ztterm[ndone],&ztval[ndone],&ztterm[ndone+1],&ztval[ndone+1],
                    &ztterm[ndone+2],&ztval[ndone+2],&ztterm[ndone+3],&ztval[ndone+3],
                    &ztterm[ndone+4],&ztval[ndone+4],&ztterm[ndone+5],&ztval[ndone+5],
                    &ztterm[ndone+6],&ztval[ndone+6],&ztterm[ndone+7],&ztval[ndone+7],
                    &ztterm[ndone+8],&ztval[ndone+8]);
        break;}
    if(i==0){
        printf("No angular distortion correction applied\n");
        return 0;}
    ndone+=(i-1)/2;
    if(ndone<ntzern){
        while(1){
            if(!fgets(line,133,parmfile))break;
            i=sscanf(line,"%d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f",
                        &ztterm[ndone],&ztval[ndone],&ztterm[ndone+1],&ztval[ndone+1],
                        &ztterm[ndone+2],&ztval[ndone+2],&ztterm[ndone+3],&ztval[ndone+3],
                        &ztterm[ndone+4],&ztval[ndone+4],&ztterm[ndone+5],&ztval[ndone+5],
                        &ztterm[ndone+6],&ztval[ndone+6],&ztterm[ndone+7],&ztval[ndone+7],
                        &ztterm[ndone+8],&ztval[ndone+8],&ztterm[ndone+9],&ztval[ndone+9]);
              ndone+=i/2;
              if(ndone==ntzern) break;}
        }
    if(ndone<ntzern){
        printf("Error reading distortion file: missing terms\n");
        return 1;}
    fclose(parmfile);
    return 0;}

//
//Readcof
//

int Readcof(char filename[]){

    FILE     *parmfile;
    char     CHAR,line[133],lane[133];
    int      i;
    char     *COS_HOME;

    strcpy(line,filename);
    if(strstr(line,".dewoff")==NULL) strcat(line,".dewoff");
    if((parmfile=fopen(line,"r"))==NULL){
        //look in sdata/dewoff
        COS_HOME=malloc(sizeof(CHAR)*80);
        COS_HOME=getenv("COSMOS_HOME");
        if(COS_HOME==NULL){
            printf("COSMOS_HOME undefined!\n");
            return 1;}
        strcpy(lane,COS_HOME);
        strcat(lane,"/sdata/dewoff/");
        strcat(lane,line);
        if((parmfile=fopen(lane,"r"))==NULL){
          printf("cant find %s in either current directory or sdata/dewoff\n",
                 filename);
          return 1;}
        }
    for(;;){
        if(!fgets(line,133,parmfile)){
            fclose(parmfile);
            return 1;}
        if(*line=='#') continue;
        if(sscanf(line,"%f %f %f %f ",&dcth,&dscale,&dcx,&dcy)<4){
            return 1;};
        dcth/=57.2958;
        break;}
    fclose(parmfile);
    return 0;}

//
//SetCof
//
void SetCof(float dc, float ds, float dx, float dy){

    dcth=dc;
    dscale=ds;
    dcx=dx;
    dcy=dy;
    return ;}
    //
    //GetCof
    //
void GetCof(float* dc, float* ds, float* dx, float* dy){
    *dc=dcth;
    *ds=dscale;
    *dx=dcx;
    *dy=dcy;
    return ;}

//
//Getchipdat
//

void Getchipdat(dewdat* dewinfo){

    dewinfo->nchip=nchip; //#chips
    dewinfo->xchip=xchip+1; //#x pixels
    dewinfo->ychip=ychip+1; //#y pixels
    dewinfo->xarsize=xarsize; //array size
    dewinfo->yarsize=yarsize;
    dewinfo->hand=hand; //chips contiguous?
    dewinfo->scale=scale/dscale;
    dewinfo->theta=dcth;
    for(i=1;i<=nchip;i++){
        dewinfo->sx[i]=sxz[i-1];
        dewinfo->sy[i]=syz[i-1];
        dewinfo->xz[i]=xz[i-1];
        dewinfo->yz[i]=yz[i-1];
        dewinfo->xloc[i]=xxn[i-1];
        dewinfo->yloc[i]=yyn[i-1];}
    return;}

int Getchipdat_f(char* var,float* val){
    if(!(strcmp(var,"scale"))){
      *val=scale/dscale;
      return(0);}
    if(!(strcmp(var,"theta"))){
      *val=dcth;
      return(0);}
    return(1);}

  int Getchipdat_i(char* var,int* val){
      if(!(strcmp(var,"xarsize"))){
        *val=xarsize;
        return(0);}
      if(!(strcmp(var,"yarsize"))){
        *val=yarsize;
        return(0);}
      return(1);}

//
//fp2ccd
//

int fp2ccd(float xfp, float yfp, float *xccd, float *yccd, float lambda){

    int           i,chip;
    double     r,theta,cth,rr;
    float        dscl,xx,zz,tt;
    //moe parameters
    float     moecoef[6]={1.83521176e-01,-7.49109947e-04,-1.21785697e-04,
                                    -8.07591646e-07,1.38731509e-08,1.67187628e-10};

     /*-------------------correct for dewar offsets-----------------------------*/

    r=sqrt((double)xfp*(double)xfp+(double)yfp*(double)yfp);
    if(xfp==0){
        if(yfp==0){
        theta=0;}
        else{
            theta=1.5707963;
            if(yfp<0){
                theta=-theta;}
            }
        }
    else if(yfp==0){
        theta=0;
        if(xfp<0){
            theta=3.141592654;}
        }
    else{
        theta= atan(yfp/xfp);
        if(xfp<0){
            if(yfp>0){
                theta+=3.141592654;}
            else{
                theta-=3.141592654;}
            }
        }
    dscl=dscale;

    //IMACS scale,distortion corrections

    if(strstr(obsinfo.instrument,"IMACS") != NULL){
        //MOE
        if(strstr(obsinfo.grating,"MOE") != NULL){
            xx=polyvalue(xfp,moecoef,5);
            xfp-=xx-.07*(obsinfo.gr_order-6);
            r=sqrt((double)xfp*(double)xfp+(double)yfp*(double)yfp);
            cth=xfp/r;
            if(cth>-1. && cth<1.){
                theta= acos(cth);}
            else{
                theta = (cth>0) ? 0.: 1.570796;}
            if(yfp<0) theta=-theta;
            }
        //short camera
        if(strcmp(obsinfo.camera,"SHORT") ==0){
            //scale vs lambda
            if(sl_n){
                dscl=dscale*(1.0 + 0.001*polyvalue(lambda,sl_coef,sl_n));}
            //distortion map
            rr=r/66.0;//maximum radius in focal plane coordinates
            if(nrzern>0){
                r+=zernike(rr,theta,nrzern,zrterm, zrval)/scale;}
            if(ntzern>0){
                zz=zernike(rr,theta,ntzern,ztterm,ztval)/(r*scale);
                theta+=zz;}
            }
        //long camera

        }

    theta-=dcth;
    xfp=r*cos(theta);
    yfp=r*sin(theta);
    xfp/=dscl;
    yfp/=dscl;
    xfp-=dcx;
    yfp-=dcy;

    // one contiguous chip
    if(nchip==1){
        yfp*=hand;                         //parity vs imacs dewar
        *xccd=xfp*scale+xarsize/2;
        *yccd=yfp*scale+yarsize/2;
        if(*xccd>=0 && *xccd<xarsize && *yccd>=0 && *yccd<yarsize){
            chip=1;}
        else{
            chip=-1;}
        return chip;}

    //multiple chips
    else{

        //identify chip
        chip=-1;
        for(i=0;i<nchip;i++){
            if(xfp<xmin[i] || xfp>xmax[i] || yfp<ymin[i] || yfp>ymax[i]) continue;
            chip=i;
            break;}
        if(chip<0){
            chip=0;
            return chip;}

    //transform coordinates

    xfp-=xx0[chip];
    yfp-=yy0[chip];
    r=sqrt(xfp*xfp+yfp*yfp);
    theta=acos((double)xfp/r);
    if(yfp<0) theta=-theta;
    theta-=dtheta[chip];
    r*=scale;
    xfp=r*cos(theta);
    yfp=r*sin(theta);
    *xccd=xz[chip]+sxz[chip]*xfp;
    *yccd=yz[chip]+syz[chip]*yfp;
    if(*xccd>=0 && *xccd<=xchip && *yccd>=0 & *yccd<=ychip){
        chip+=1;}
    else{
        chip=-1-chip;}
        return chip;}
    }

int ccd8(int chip, float xccd, float yccd, float *x8, float *y8){

  int chp;

    if(nchip==1){
        *x8=xccd-xarsize/2;
        *y8=hand*(yccd-yarsize/2);}
    else{
        chp=chip-1;
        *y8=yy0[chp]*scale +yz[chp]+syz[chp]*yccd;
        *x8=xx0[chp]*scale +xz[chp]+sxz[chp]*xccd;}
    return 0;}


/*
 *
 *  mspos  converts ccd coordinates to mosaic coordinates
 *
 */

int mspos(int chip,float xccd,float yccd,float *x8cd,float *y8cd){

    int chop;

    if(chip<0) chip=-chip;
    chop=chip-1;
    if(nchip==1){
        *x8cd=xccd;
        *y8cd=yccd;
        if(sxz[chop]<0) *x8cd--;
        if(syz[chop]<0) *y8cd--;
    return 0;}
    *x8cd=xxn[chop]*(xchip+1)+sxz[chop]*xccd;
    if(sxz[chop]<0) *x8cd=*x8cd-1;
    *y8cd=yyn[chop]*(ychip+1)+syz[chop]*yccd;
    if(syz[chop]<0) *y8cd=*y8cd-1;
    return 0;}

/*
 *
 *  zernike returns distortion correction at location r,theta from specified set of zernike terms
 *
 */

float zernike(float r, float theta, int nzr, int *ztrm, float *val){

    int i;
    float corr;
    static void *label[37]={&&A,&&B,&&C,&&D,&&E,&&F,&&G,&&H,&&I,&&J,&&K,&&L,
                                      &&M,&&N,&&O,&&P,&&Q,&&R,&&S,&&T,&&U,&&V,&&W,&&X,
                                      &&Y,&&Z,&&a,&&b,&&c,&&d,&&e,&&f,&&g,&&h,&&i,&&j,&&k};

    corr=0;
    for(i=0;i<nzr;i++){
    goto *(label[ztrm[i]-1]);
  A:corr+=val[i]*1;
    continue;
  B:corr+=val[i]*2*r*cos(theta);
     continue;
  C:corr+=val[i]*2*r*sin(theta);
    continue;
  D:corr+=val[i]*sqrt(3)*(2*pow(r,2)-1);
    continue;
  E:corr+=val[i]*sqrt(6)*pow(r,2)*sin(2*theta);
    continue;
  F:corr+=val[i]*sqrt(6)*pow(r,2)*cos(2*theta);
    continue;
  G:corr+=val[i]*sqrt(8)*(3*pow(r,2)-2)*r*sin(theta);
    continue;
  H:corr+=val[i]*sqrt(8)*(3*pow(r,2)-2)*r*cos(theta);
    continue;
  I:corr+=val[i]*sqrt(8)*pow(r,3)*sin(3*theta);
    continue;
  J:corr+=val[i]*sqrt(8)*pow(r,3)*cos(3*theta);
    continue;
  K:corr+=val[i]*sqrt(5)*(1-6*pow(r,2)+6*pow(r,4));
    continue;
  L:corr+=val[i]*sqrt(10)*(4*pow(r,2)-3)*pow(r,2)*cos(2*theta);
    continue;
  M:corr+=val[i]*sqrt(10)*(4*pow(r,2)-3)*pow(r,2)*sin(2*theta);
    continue;
  N:corr+=val[i]*sqrt(10)*pow(r,4)*cos(4*theta);
    continue;
  O:corr+=val[i]*sqrt(10)*pow(r,4)*sin(4*theta);
    continue;
  P:corr+=val[i]*sqrt(12)*(10*pow(r,4)-12*pow(r,2)+3)*r*cos(theta);
    continue;
  Q:corr+=val[i]*sqrt(12)*(10*pow(r,4)-12*pow(r,2)+3)*r*sin(theta);
    continue;
  R:corr+=val[i]*sqrt(12)*(5*pow(r,2)-4)*pow(r,3)*cos(3*theta);
    continue;
  S:corr+=val[i]*sqrt(12)*(5*pow(r,2)-4)*pow(r,3)*sin(3*theta);
    continue;
  T:corr+=val[i]*sqrt(12)*pow(r,5)*cos(5*theta);
    continue;
  U:corr+=val[i]*sqrt(12)*pow(r,5)*sin(5*theta);
    continue;
  V:corr+=val[i]*sqrt(7)*(20*pow(r,6)-30*pow(r,4)+12*pow(r,2)-1);
    continue;
  W:corr+=val[i]*sqrt(14)*(15*pow(r,4)-20*pow(r,2)+6)*pow(r,2)*sin(2*theta);
    continue;
  X:corr+=val[i]*sqrt(14)*(15*pow(r,4)-20*pow(r,2)+6)*pow(r,2)*cos(2*theta);
    continue;
  Y:corr+=val[i]*sqrt(14)*(6*pow(r,2)-5)*pow(r,4)*sin(4*theta);
    continue;
  Z:corr+=val[i]*sqrt(14)*(6*pow(r,2)-5)*pow(r,4)*cos(4*theta);
    continue;
  a:corr+=val[i]*sqrt(14)*pow(r,6)*sin(6*theta);
    continue;
  b:corr+=val[i]*sqrt(14)*pow(r,6)*cos(6*theta);
    continue;
  c:corr+=val[i]*4*(35*pow(r,6)-60*pow(r,4)+30*pow(r,2)-4)*r*sin(theta);
    continue;
  d:corr+=val[i]*4*(35*pow(r,6)-60*pow(r,4)+30*pow(r,2)-4)*r*cos(theta);
    continue;
  e:corr+=val[i]*4*(21*pow(r,4)-30*pow(r,2)+10)*pow(r,3)*sin(3*theta);
    continue;
  f:corr+=val[i]*4*(21*pow(r,4)-30*pow(r,2)+10)*pow(r,3)*cos(3*theta);
    continue;
  g:corr+=val[i]*4*(7*pow(r,2)-6)*pow(r,5)*sin(5*theta);
    continue;
  h:corr+=val[i]*4*(7*pow(r,2)-6)*pow(r,5)*cos(5*theta);
    continue;
  i:corr+=val[i]*4*pow(r,7)*sin(7*theta);
    continue;
  j:corr+=val[i]*4*pow(r,7)*cos(7*theta);
    continue;
  k:corr+=val[i]*3*(70*pow(r,8)-140*pow(r,6)+90*pow(r,4)-20*pow(r,2)+1);}
  return corr;}
