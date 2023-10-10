/* Utility programs for Spherical Astronomy - Header file */
/* Mask genertion edition */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
**                                                                   **
**                  Copyright (C) 2006 by                            **
**     Ken Clardy,  Carnegie Observatories,  Pasadena California.    **
**                                                                   **
**      This software is proprietary, and may be used or copied      **
**          only with the written consent of the author.             **
**              All copies must retain this notice.                  **
**                                                                   **
** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ---------------------------------------------------------------------- */

#ifndef	INCLUDE_KDCUTIL_H
#define	INCLUDE_KDCUTIL_H

/* INCLUDEs ------------------------------------------------------------- */

#include  <stdio.h>
#include  <math.h>
#ifdef  __sun
#include  <sunmath.h>  // sincos, nint (and others)
#else
/* This should really be ifdef LINUX; these items are
  not defined in the linux math.h file, but should be... */
#ifdef  HAVE_SC
extern void  sincos (double, double*, double*);
#else
void  sincos (double, double*, double*);
#endif
extern double  round (double);
/* An alternative is to write our own substitutes for these
  things and define them when needed */
#endif
#include  <time.h>     // get time

/* DEFINEs -------------------------------------------------------------- */

/* #define  pi  3.14159265358979323846 */
/* Size of alternate angular measures */
/* #define  DEGREE  (pi/180.0) */
/* #define  ARCSECOND (DEGREE/3600.0) */
/* #define  HOUR  (pi/12.0) */

#define  TwoPi      (M_PI*2.0)
#define  Degree     (M_PI/180.0)
#define  ArcMinute  (Degree/60.0)
#define  ArcSecond  (Degree/3600.0)
#define  Hour       (M_PI/12.0)
/* Note - M_PI is defined in <math.h> */
#define  Sidereal   1.002737909256
/* Ratio sidereal to solar interval */

#define  SysEpoch  2440587.5
/* Sysepoch is J.D. at epoch of Unix time, 0h GMT 1 Jan, 1970 */

/* TYPEDEFs ------------------------------------------------------------- */

struct  vector3  { double  x, y, z; };
typedef	struct	vector3	vect3;

struct	vector2  { double  x, y; };
typedef	struct	vector2	vect2;

struct	pvector3 { double  r, t, z; };
typedef struct  pvector3 pvec3;

struct  pvector2 { double  r, t; };
typedef struct  pvector2 pvec2;

struct	matrix2 { vect2  x, y; };
typedef	struct	matrix2  matrix2;

struct	matrix3 { vect3  x, y, z; };
typedef	struct	matrix3  matrix3;

/* GLOBALs -------------------------------------------------------------- */

/* FUNCTION PROTOTYPEs -------------------------------------------------- */

vect3	make3vect (double, double, double);
vect2	make2vect (double, double);
vect3	get3de2v (vect2, double z);
vect2	get2de3v (vect3);
vect3	sum3vect (vect3, vect3);
vect2	sum2vect (vect2, vect2);
vect3	sub3vect (vect3, vect3);
vect2	sub2vect (vect2, vect2);
vect3	mul3vect (vect3, double);
vect2	mul2vect (vect2, double);
double	dot3vect (vect3, vect3);
double	dot2vect (vect2, vect2);
double	vect3norm (vect3);
double	vect2norm (vect2);
vect3	norm3vect (vect3);
vect2	norm2vect (vect2);
vect3	cross3vect (vect3, vect3);
double	cross2vect (vect2, vect2);
vect2	flip2vect (vect2);
pvec2	polar2 (vect2);
vect2	rect2 (pvec2);
pvec3	polar3 (vect3);
vect3	rect3 (pvec3);
double	znorm (vect3);

double	dist3 (vect3, vect3);
double	dist2 (vect2, vect2);
double	nearp3 (vect3, vect3);
double	nearp2 (vect2, vect2);
double	dseg3 (vect3, vect3);
double	dseg2 (vect2, vect2);
double	pdis2 (vect2, vect2);
int	inside (double, double, double);
double	dparel (vect2, vect2,  vect2) ;

vect3	rotate (vect3, double, int);
vect2	rot2 (vect2, double);
vect3	sph2vec (double, double);
vect3	had2vec (double, double);
void	vec2sph (double*, double*, vect3);
void	vec2had (double*, double*, vect3);
vect2	pol2vec (double, double);
void	vec2pol (double*, double*, vect2);
void	disp3vect (char*, vect3, char*);
void	disp2vect (char*, vect2, char*);
vect3	tang3vec (vect3, double, double);
vect2	tang2vec (vect3);
vect3	vec2tang (vect2);
vect3	vec3tang (vect3, double, double);
void	smatrix (vect3*, vect3);
void	mrotate (vect3*, double, int);

void	dntri ( double* az, double* zd, double* pa, 
		double* da, double* dz, double* dp,
		double  h,  double  d,  double  g);

vect2	difref (vect2 fp, vect2 zp);
double	ndxref (double wvlen, double p, double t);
vect3	refracted (vect3 v, double cr, int z);

vect3	prec (double, vect3, double);
void	pprec  (double*, double*, double, double, double, double);
void	hdprec (double*, double*, double, double, double, double);

vect3	mmult3 (vect3*, vect3);
vect3	tmult3 (vect3*, vect3);
vect2	mmult2 (vect2*, vect2);
void	trpose3 (vect3*);
void	trpose2 (vect2*);
void	mcpy3 (vect3*, vect3*);
void	mcpy2 (vect2*, vect2*);
void	pmatrix (vect3*, double, double);
void	r2matrix (vect2*, double);
void	mxm3 (vect3*, vect3*, vect3*);

vect2	mmul2v (matrix2, vect2);
matrix2	mtr2   (matrix2);
matrix2	minv2  (matrix2);
matrix2	mul2m  (matrix2, matrix2);
matrix2	rot2m  (double);
vect3	mmul3v (matrix3, vect3);
matrix3	mtr3   (matrix3);
matrix3	minv3  (matrix3);
// matrix3	minv3g (matrix3);
matrix3	mul3m  (matrix3, matrix3);
matrix3	rot3m  (double, int);
matrix3	precmat (double, double);
matrix3	precmat3 (double, double);

double	fdifmod (double, double, double);
double  fmodsum (double, double, double);
double	fabsdif (double, double, double);

int	Nint (double a);

double  jdnow (void);
void	jd2date (double jd, int* year, int* month, int* day);
double	jdate ( int year, int month, int day);
double	Jyear (double  jd);
double	jdaye (double  year);
double	gst0 (double  jd);

char	exchar (int);	/* Formerly static */
char*	timestamp (time_t* tt);
time_t	unstamp (char* ts);

double	tranfw (double x, double zp, double sf);
double	tranrv (double x, double zp, double sf);

/* ---------------------------------------------------------------------- */
/* End of kdcutil.h package */

#endif	/* INCLUDE_KDCUTIL_H */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
