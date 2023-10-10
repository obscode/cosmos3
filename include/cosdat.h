#ifndef  INCLUDE_COSDAT_H
#define  INCLUDE_COSDAT_H

typedef struct obsdef{
  char   instrument[10];
  int    gr_order;
  float  gr_angle;
  float  temp;
  float  alignrot;
  char   grating[20];
  char   dewoff[80];
  char   camera[6];
  char   mode[7];
  char   mask[80];
  int    nshuffle;
  float  ranod;
  float  decnod;
  char   dewar[15];
  char distor[80];} obsdef;

typedef struct fitsdef{
  int    naxis;
  long   naxes[3];
  int    bitpix;
  int    binning;
  int    ybinning;
  int    overscan;
  int    biaslins;
  int    subrstr;
  int    subx[2];
  int    suby[2];
  int    nshuffle;
  float  ranod;
  float  decnod;
  float  exptime;
  char   date[12];
  char   dewarori[3];} fitsdef;

typedef struct dewdat{
  int    nchip;
  int    xchip;
  int    ychip;
  int    xarsize;
  int    yarsize;
  int    hand;
  float  scale;
  float  theta;
  int    xz[9];
  int    yz[9];
  int    xloc[9];
  int    yloc[9];
  int    sx[9];
  int    sy[9];} dewdat;
 
  
#endif
