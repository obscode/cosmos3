/* Utility programs for Spherical Astronomy */
/* Mask Generation suplemential edition */

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

/* We need to be able to do these things:
 1) Precess catalog position to given epoch (coords or vector)
 2) Given a center, and catalog position, find projected vector
 3) Inverse, from projected point, find ra/dec
.. */

/* Other uses of these programs:
 1) Do l.s. solution for change x, y, rotation; then convert x, y
  changes into ra/dec changes.
 2) rotation of 2-d x,y values about a given center
  Uses 2-d dot product, vector add, subtract.
.. */

/* Radical idea -- rather than subscripting, we define a 3-vector
as a *structure* with x,y,z components.  Then we can transfer it as
a value, use it as an argument, and define a matrix as an array of
3 vectors.  */


/* ---------------------------------------------------------------------- */
/*  Defines  */
#ifdef    NEED_SINCOS
#undef    HAVE_SC
#else
#define   HAVE_SC
#endif
/* HAVE_SC is defined when the sincos function is available, and can
  save up to 40% execution time when computing both sin and cos of
  the same angle.  If this function is not available, undef here.  */
/* NEED_SINCOS is defined in those systems which need the sincos
  function, having no such thing in their math library.  A lame version
  using sin and cos is implemented here, to keep things compatable.  It
  in no way improves the computational efficiency as a real implementation
  of sincos would naturally do.  */

/* ---------------------------------------------------------------------- */

#include	"kdcutil.h"

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* Compatability code */

#ifdef   NEED_SINCOS
void	sincos (double angle, double* s, double* c) {
    *s = sin (angle);
    *c = cos (angle);
}
#endif

/* ---------------------------------------------------------------------- */

/* Some basic vector things... */

vect3	make3vect (double x, double y, double z) {
vect3	result;
    result.x = x;
    result.y = y;
    result.z = z;
    return result;
}

vect2	make2vect (double x, double y) {
vect2	result;
    result.x = x;
    result.y = y;
    return result;
}

vect3	get3de2v (vect2 v, double z) {
vect3	result;
    result.x = v.x;
    result.y = v.y;
    result.z = z;
    return result;
}

vect2	get2de3v (vect3 v) {
vect2	result;
    result.x = v.x;
    result.y = v.y;
    return  result;
}

vect3	sum3vect (vect3 a, vect3 b) {
vect3	result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return result;
}

vect2	sum2vect (vect2 a, vect2 b) {
vect2	result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    return result;
}

vect3	sub3vect (vect3 a, vect3 b) {
vect3	result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    return result;
}

vect2	sub2vect (vect2 a, vect2 b) {
vect2	result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    return result;
}

vect3	mul3vect (vect3 a, double b) {
vect3	result;
    result.x = a.x * b;
    result.y = a.y * b;
    result.z = a.z * b;
    return result;
}

vect2	mul2vect (vect2 a, double b) {
vect2	result;
    result.x = a.x * b;
    result.y = a.y * b;
    return result;
}

double	dot3vect (vect3 a, vect3 b) {
    return (a.x*b.x + a.y*b.y + a.z*b.z);
}

double	dot2vect (vect2 a, vect2 b) {
    return (a.x*b.x + a.y*b.y);
}

double	vect3norm (vect3 a) {
    return (sqrt(a.x*a.x + a.y*a.y + a.z*a.z));
}

double	vect2norm (vect2 a) {
/*     return (sqrt(a.x*a.x + a.y*a.y)); */
    return  (hypot (a.x, a.y));
}

vect3	norm3vect (vect3 a) {
vect3	result;
double	r;
    r = vect3norm (a);
    if (r == 0.0) return a;
    result.x = a.x/r;
    result.y = a.y/r;
    result.z = a.z/r;
    return  result;
}

vect2	norm2vect (vect2 a) {
vect2	result;
double	r;
/*  r = vect2norm (a);  */
    r = hypot (a.x, a.y);
    if (r == 0.0) return a;
    result.x = a.x/r;
    result.y = a.y/r;
    return  result;
}

vect3	cross3vect (vect3 a, vect3 b) {
vect3	result;
    result.x = a.y*b.z - b.y*a.z;
    result.y = a.z*b.x - b.z*a.x;
    result.z = a.x*b.y - b.x*a.y;
    return result;
}

double	cross2vect (vect2 a, vect2 b) {
    return (a.x*b.y - a.y*b.x);
}

vect2	flip2vect (vect2 a) {
/* Useful for supplemental angle in vec2pol calls, for example */
vect2	result;
    result.x = -a.x;
    result.y =  a.y;
    return  result;
}

pvec2	polar2 (vect2 a) {	/* Convert to polar */
pvec2	result;
/*  retult.r = vect2norm (a);  */
    result.r = hypot (a.y, a.x);
    result.t = atan2 (a.y, a.x);
    return  result;
}

vect2	rect2 (pvec2 a) {	/* Convert to rectangular */
vect2	result;
#ifdef  HAVE_SC
double	s, c;
    sincos (a.t, &s, &c);
    result.x = a.r * c;
    result.y = a.r * s;
#else
    result.x = a.r * cos(a.t);
    result.y = a.r * sin(a.t);
#endif
    return  result;
}

pvec3	polar3 (vect3 a)  {	/* Convert to polar3 */
pvec3	result;
    result.r = hypot (a.y, a.x);
    result.t = atan2 (a.y, a.x);
    result.z = a.z;
    return  result;
}

vect3	rect3 (pvec3 a)  {	/* Convert to rectangular */
vect3	result;
#ifdef  HAVE_SC
double	s, c;
    sincos (a.t, &s, &c);
    result.x = a.r * c;
    result.y = a.r * s;
#else
    result.x = a.r * cos(a.t);
    result.y = a.r * sin(a.t);
#endif
    result.z = a.z;
    return  result;
}

double	znorm (vect3 a)  {	/* Ray tracing special norm */
double	N;

    N = 1.0 - a.x*a.x - a.y*a.y;
    return ( N < 0.0 ? -a.z : sqrt(N) );
}
/* znorm is used in ray tracing; the "cannonical momentum" of a ray
  is a unit vector multiplied by the index of refraction.  When crossing
  a boundry, (z is normal to the boundry) the x and y components remain
  the same but z is changed to normalize to the new index.  Going from
  index i to index r, multiply the normalized unit vector by i/r and
  then apply znorm to find z; then re-normalize the unit vector.  If 
  the resulting z would be immaginary, that is total internal reflection,
  and we return -z instead of the refracted normal z.  */



/* ---------------------------------------------------------------------- */

/*  ---  Geometry Utilities  ---  */

double	dist3 (vect3 a, vect3 b) {
    return  vect3norm (sub3vect (a, b) );
}

double	dist2 (vect2 a, vect2 b) {
    return  vect2norm (sub2vect (a, b) );
}

/* Distance along vector a to nearest point to p */
double	nearp3 (vect3 a, vect3 p) {
double	result;
    result = vect3norm (a);
    if (result > 0.0) result = dot3vect (a, p) / (result*result);
    return  result;
}

double	nearp2 (vect2 a, vect2 p) {
double	result;
    result = vect2norm (a);
    if (result > 0.0) result = dot2vect (a, p) / (result*result);
    return  result;
}

/* Distance from p to segment 0,a; point on line limited by ends */
double	dseg3 (vect3 a, vect3 p) {
double	f;
    f = nearp3 (a, p);
    if (f < 0.0) f = 0.0;
    if (f > 1.0) f = 1.0;
    return  dist3 (mul3vect(a,f), p);
}

double	dseg2 (vect2 a, vect2 p) {
double	f;
    f = nearp2 (a, p);
    if (f < 0.0) f = 0.0;
    if (f > 1.0) f = 1.0;
    return  dist2 (mul2vect(a,f), p);
}


/* Signed perpendicular distance from point p to vector a */
double	pdis2 (vect2 a, vect2 p) {
double	A;
    A = vect2norm (a);
    if (A == 0.0) return  vect2norm (p);
    return  cross2vect (a, p) / A;
/* or.. return  cross2vect (norm2vect(a), p);  ?? */
}

/*  Interval comparison -- integer decision function;
  Return -1, 0 or 1 for x relative to directed interval a,b:
  Result is -1 (a), 0 (inside), +1 (b) side of interval... */
/* Compare is inclusive of the boundry */
int	inside (double a, double x, double b) {
    if (a < b) {
	if (x < a) return -1;
	return x > b ? 1 : 0;
    } else {
	if (x > a) return -1;
	return x < b ? 1 : 0;
    }
}

/*  Minimum distance to parallelogram:  
  A parallelogram is defined by side vectors a and b, from a corner
  at origin (0,0).  Find least distance from point p to this parallelogram,
  with the distance being negative inside the figure.  */
double	dparel (vect2 a, vect2 b,  vect2 p)  {
double	result;
double	d;
int	in;

/* If either side is zero, use distance to segment */
    if (vect2norm (a) == 0.0) return dseg2 (b, p);
    if (vect2norm (b) == 0.0) return dseg2 (a, p);

/* Find if point e is interior to parallelogram - in = 0 if so */
    in = inside (0.0, cross2vect (a, p), cross2vect (a, b) ) |
	 inside (0.0, cross2vect (b, p), cross2vect (b, a) );

/* Find least distance to any of 4 segments bounding parallelogram */
    result = dseg2 (a, p);
    d = dseg2 (b, p);
    if (d < result) result = d;
    d = dseg2 (a, sub2vect (p, b));
    if (d < result) result = d;
    d = dseg2 (b, sub2vect (p, a));
    if (d < result) result = d;
    if (!in) result = -result;
    return  result;
}




/* ---------------------------------------------------------------------- */

/* Spherical geometry vector things */

vect3	rotate (vect3 v, double angle, int axis) {
vect3	result;
double	s, c;
    if (angle == 0.0) return  v;  /* null case */
#ifdef  HAVE_SC
    sincos (angle, &s, &c);
#else
    s = sin (angle);
    c = cos (angle);
#endif
    switch (axis) {
	case 1:
	    result.x = v.x;
	    result.y = v.y*c - v.z*s;
	    result.z = v.z*c + v.y*s;
	    break;
	case 2:
	    result.x = v.x*c + v.z*s;
	    result.y = v.y;
	    result.z = v.z*c - v.x*s;
	    break;
	case 3:
	    result.x = v.x*c - v.y*s;
	    result.y = v.y*c + v.x*s;
	    result.z = v.z;
	    break;
	default:
	    result = v;
    }
    return result;
}

vect2	rot2 (vect2 v, double angle) {
vect2	result;
double	s, c;
#ifdef  HAVE_SC
    sincos (angle, &s, &c);
#else
    s = sin (angle);
    c = cos (angle);
#endif
    result.x = v.x*c - v.y*s;
    result.y = v.y*c + v.x*s;
    return result;
}

vect3	sph2vec (double lng, double lat) {
vect3	result;
#ifdef  HAVE_SC
double	sa, ca;
double	sd, cd;
/*    result = make3vect (cos(lat), 0.0, sin(lat)); */
/*    result = rotate (result, lng, 3); */
    sincos (lng, &sa, &ca);
    sincos (lat, &sd, &cd);
    result.x = cd * ca;
    result.y = cd * sa;
    result.z = sd;
#else
double	cd;
    cd = cos (lat);
    result.x = cd * cos (lng);
    result.y = cd * sin (lng);
    result.z = sin (lat);
#endif
    return result;
}

vect3	had2vec (double lng, double lat) {
/* Similar to sph2vec, using hours and degrees arguments */
vect3	result;
#ifdef  HAVE_SC
double	sa, ca;
double	sd, cd;
/*    result = make3vect (cos(lat), 0.0, sin(lat)); */
/*    result = rotate (result, lng, 3); */
    sincos (lng*Hour,   &sa, &ca);
    sincos (lat*Degree, &sd, &cd);
    result.x = cd * ca;
    result.y = cd * sa;
    result.z = sd;
#else
double	cd;
double	h, d;
    h = lng * Hour;
    d = lat * Degree;
    cd = cos (d);
    result.x = cd * cos (h);
    result.y = cd * sin (h);
    result.z = sin (d);
#endif
    return result;
}


void	vec2sph (double *lng, double *lat, vect3 v) {
    *lng = atan2 (v.y, v.x);
    if (*lng < 0.0) *lng += TwoPi;
    *lat = atan2 (v.z, hypot (v.y, v.x) );
}

void	vec2had (double *lng, double *lat, vect3 v) {
/* Similar to vec2sph, returning hours and degrees arguments. */
    *lng = atan2 (v.y, v.x) / Hour;
    if (*lng < 0.0) *lng += 24.0;
    *lat = atan2 (v.z, hypot (v.y, v.x) ) / Degree;
}

vect2	pol2vec (double r, double a) {
vect2	result;
#ifdef  HAVE_SC
double	s, c;
    sincos (a, &s, &c);
    result.x = r * c;
    result.y = r * s;
#else
    result.x = r * cos (a);
    result.y = r * sin (a);
#endif
    return  result;
}

void	vec2pol (double *r, double *a, vect2 v) {
/*  *r = vect2norm (v);  */
    *r = hypot (v.y, v.x);
    *a = atan2 (v.y, v.x);
    if (*a < 0.0) *a += TwoPi;
}

void	disp3vect (char *left, vect3 v, char *right) {
    printf ("%s%12.9f %12.9f %12.9f%s",
	left, v.x, v.y, v.z, right );
}

void	disp2vect (char *left, vect2 v, char *right) {
    printf ("%s%12.9f %12.9f%s",
	left, v.x, v.y, right );
}

/* 
  The tangent plane coordinate system here has X along the reference
  vector direction, Y a tangent pointing east, and Z a tangent pointing
  north.  We can transform both to the tangent (tang3vec) and from it
  (vec3tang).  This is used to project points to a tangent plane with
  a center in the reference direction.  For multiple points, a projection
  matrix may be defined, in the same way as a precession matrix below.

  For a projected vector in tangent plane coordinates, one may plot
  the x,y (right, up) pair as (-v.y/v.x, v.z/v.x) in radian units;
  scale appropriately to effective radius.
*/

vect3	tang3vec (vect3 v, double lng, double lat) {
/* From vector to tangent plane */
vect3	result;
    result = rotate (v, -lng, 3);
    result = rotate (result, lat, 2);
    return result;
}

vect2	tang2vec (vect3 v) {
/* From tangent plane vector, find 2-vector in plane */
vect2	result;
    result.x = -v.y/v.x;
    result.y =  v.z/v.x;
    return  result;
}

vect3	vec2tang (vect2 v) {
/* Inverse of above (normalization is tricky; the 2-vector
  MUST be in radian units here) */
vect3	result;
    result.y = -v.x;
    result.z =  v.y;
    result.x = 1.0;
    return (norm3vect(result));
}

vect3	vec3tang (vect3 v, double lng, double lat) {
/* From tangent plane to vector; transform X to ref. vector */
vect3	result;
    result = rotate (v, -lat, 2);
    result = rotate (result, lng, 3);
    return result;
}

/* NEED -- DERIVE TANGENT PROJECTION MATRIX FROM EITHER A REFERENCE
 DIRECTION (A,D) OR A REFERENCE VECTOR. */

void	smatrix (vect3 *a, vect3 v) {
/* Obtain projection matrix a from reference vector v */
    a[0] = norm3vect(v);	/* Along reference vector */
    a[1] = norm3vect (make3vect(-v.y, v.x, 0.0));	/* East */
/*  if (vect3norm(a[1]) == 0.0) a[1] = make3vect (0.0, 1.0, 0.0); */
    if (vect3norm(a[1]) < 0.001) a[1] = make3vect (0.0, 1.0, 0.0);
    a[2] = norm3vect (cross3vect(a[0], a[1]));	/* North */
}

void	mrotate (vect3 *a, double angle, int axis) {
/* Rotate the transform in matrix a by the angle about axis... */
vect3	v, w;
double	s, c;
    if (angle == 0.0) return;  /* null case */
#ifdef  HAVE_SC
    sincos (angle, &s, &c);
#else
    s = sin (angle);
    c = cos (angle);
#endif
    switch (axis) {
	case 1:
	    v = a[1];
	    w = a[2];
	    a[1] = sub3vect ( mul3vect(v,c), mul3vect(w,s) );
	    a[2] = sum3vect ( mul3vect(w,c), mul3vect(v,s) );
/* note -- we may need to store each mul3vect result independently
since they are made up in local storage in that program, and could
conflict during execution... */
	    break;
	case 2:
	    v = a[0];
	    w = a[2];
	    a[0] = sum3vect ( mul3vect(v,c), mul3vect(w,s) );
	    a[2] = sub3vect ( mul3vect(w,c), mul3vect(v,s) );
	    break;
	case 3:
	    v = a[0];
	    w = a[1];
	    a[0] = sub3vect ( mul3vect(v,c), mul3vect(w,s) );
	    a[1] = sum3vect ( mul3vect(w,c), mul3vect(v,s) );
	    break;
	default:
	    ;
    }
}

/* ---------------------------------------------------------------------- */

/*  Geometry of Alt-Az Telescope  */
/* For the supprot of alt-az telescopes, and instruments mounted
  at the Naysmith focus, we add some basic telescope geometry.
  The Parallactic Angle is defined as the angle at an observed
  object from the direction to the North Celestial Pole to the
  current Zenith, measured counterclockwize from North through East.
  The Zenith Distance is the angle between the current Zenith and
  the object being observed.
  The Azimuth is the angle measured from North on the horizon,
  through East to the object being observed.

  When the object's position vector in celestial coordinates is rotated
  about the polar axis by siderial time, we get a position vector
  based on current hour angle.  This rotated by latitude gives a
  position vector in altitude and azimuth.

  The current zenith can be rendered as a celestial position, its
  right ascension is the sideral time, and its declination is the
  latitude.

  To find both zeinth distance and parallactic angle, find the current
  zenith vector in a coordinate system containing unit vectors along
  the object direction, the tangent vector pointing north at the object,
  and the tangent vector pointing east.  The zenith vector's position
  in a coordinate system centered on the object may be found by starting
  with the negative hour angle and latitude location and rotating by
  the object's declination.
.. */

/* Solve entire p-z-s navigation triangle at one step */
void	dntri ( double* az, double* zd, double* pa, 
		double* da, double* dz, double* dp,
		double  h,  double  d,  double  g)
/* Finds az=azimuth, zd=zeinth_dist, pa=parallactic and
   the time (hour angle) derivatives da, dz, dp from
   h=hour_angle, d=declination, g=latitude.
   All angles are in radians!                  */
{
double  st, ct;
double  sp, cp;
double  sd, cd;
double  x, y, z;
double	e, n;
double	dx, dy;
double	de, dn;
double	s;
double	r, rr;

/* Look up all trug functions first. */
#ifdef  HAVE_SC
    sincos (h, &st, &ct);	/* Hour Angle  */
    sincos (d, &sd, &cd);	/* Declination */
    sincos (g, &sp, &cp);	/* Latitude    */
#else
    st = sin (h);
    ct = cos (h);
    sd = sin (d);
    cd = cos (d);
    sp = sin (g);
    cp = cos (g);
#endif

/* Find vector at star position pointing to zenith */
    x = cp*st;			/* East tangent  */
    y = sp*cd - cp*sd*ct;	/* North tangent */
    z = sp*sd + cp*cd*ct;	/* Toward star   */

/* Find star vector in horizon system */
    e = -cd*st;			/* east  */
    n = cp*sd - sp*cd*ct;	/* north */
/* z is the same in both systems, star <dot> vertical */

/* Find component perpendicular to zenith, star vector */
    s = hypot (x, y);	/* Or, (e,n) as well */

/* Compute angles from vector components */
    *az = atan2 (e, n);
    *zd = atan2 (s, z);
    *pa = atan2 (x, y);

/* Now, compute the derivatives */
/* Limit size of derivatives by limiting size of 1/s.  At a
  distance of 20.6265 arc seconds, s is 1.0e-4 */
#define  MINIMX 1.0e-4
    r = 1.0 / ( (s < MINIMX) ? MINIMX : s );
    rr = r * r;

    de = -cd*ct;
    dn = sp*cd*st;
    dx = cp*ct;
    dy = cp*sd*st;

    *da = (n*de - e*dn) * rr;
    *dz = cp*cd*st * r;
    *dp = (y*dx - x*dy) * rr;

}


/* ---------------------------------------------------------------------- */

/* Refraction things */

vect2	difref (vect2 fp, vect2 zp)
{
/* Given zenith point, having magnitude of zenith distance in radians;
  and a field point, having magnitude of tangent plane coordinate in
  radians; will return the differential refraction displacement of the
  field point in tangent plane coordinates.  This result is to be
  multiplied by the actual refractive index.
  NOTE that this result is independent of tangent plane coordinate
  system, whether it is right handed or left handed, and of the 
  actual refraction.  The parallactic angle is in the direction of
  the zenith point here.
  This subroutine is intended to be a general solution to the problem
  of differential refraction in many applications.               ... */

vect2	dr;		/* The result  */
double	z1, z2;		/* zenith distances in radians  */
double	p1, p2;		/* paralactic angles at center, fp  */
double	r, g;		/* polar coordinates of fp  */
double	B, cz;		/* auxiliary angles in triangle  */
double	sr, cr;		/* sin,cos of r  */
double	sz1, cz1;	/* sin,cos of z1  */
double	sA, cA;		/* sin,cos of auxiliary angle A  */
vect2	r1, r2;		/* refraction displacement at center, fp  */

    vec2pol (&z1, &p1, zp);	/* z and parallactic at center  */
    r1 = pol2vec (tan(z1), p1);	/* refraction of center  */
/* We could substitute, both here and in the computation of r2 below,
  a better approximation than tan(z1) should it be desired... */
    vec2pol (&r, &g, fp);	/* Plane polar coords of fp  */

/* We have a spherical triangle; its corners are:
	A : the center of the field
	B : the field point
	C : the zenith
and we have the equations, in which the lower case letters are the
oposite sides of the triangle:
	y = sin a sin B = sin b sin A
	x = sin a cos B = cos b sin c - sin b cos c cos A
	z = cos a       = cos b cos c + sin b sin c cos A
and we substitute the needed values as:
	a : z2,  b : z1,  c : r.
	A : angle at center = p1 - g;
	B : angle at field point: p2 + B = g + Pi;
	thus p2 = g + Pi - B.
Using vec2sph we derive B and cz from the x,y,z values above,
and note that cz is the complement of z2.		... */

#ifdef  HAVE_SC
    sincos (z1, &sz1, &cz1);	/* side b  */
    sincos ( r,  &sr,  &cr);	/* side c  */
    sincos (p1-g, &sA, &cA);	/* Angle A  */
#else
    sz1 = sin (z1);
    cz1 = cos (z1);
    sr  = sin (r);
    cr  = cos (r);
    sA = sin (p1-g);
    cA = cos (p1-g);
#endif
    vec2sph ( &B, &cz,  make3vect (
	cz1 * sr - sz1 * cr * cA,
	sz1 * sA,
	cz1 * cr + sz1 * sr * cA)  );
    p2 = g + M_PI - B;		/* Parallactic of fp  */
    z2 = M_PI/2.0 - cz;		/* Actual z2  */
    r2 = pol2vec (tan(z2), p2);	/* Refraction at fp  */
    dr = sub2vect (r2, r1);	/* Differential Refraction  */
    return  dr;
}

double	ndxref (double wvlen, double p, double t)
/* Approximate index of refraction of atmosphere */
/* Wavelength in angstroms, pressure in mm-hg, temp. in Celsius */
/* A good visual value of wvlen is 5500.0, and a good approximation
  for p at Las Campanas is 575.0  */
/* Returns index less 1.0 */
{
double	rwsq, q, c;

/* Refraction from Allen, Astrophysical Quantities, 3rd. edition */
    rwsq = (1.0e+8) / (wvlen*wvlen);
    q = 64.328 + 29498.1/(146.0-rwsq) + 255.4/(41.0-rwsq);
    c = p*(1.0 + (1.049-0.0157*t)*p*1.0e-6) / (720.88261*(1.0+0.003661*t));
    return (q*c*1.0e-6);
}
/* Using wvlen 5500, rwsq = 3.3, q = 277.8, and t = 0, p = 575.0,
  c = .798, refn = 221.737e-6;  .01/refn = .022 radian. */


vect3	refracted (vect3 v, double cr, int z)  {
/* Find refracted (horizon) vector fron non-refracted */
/* v is a 3 vector in horizon coordinates */
/* cr is coefecient of refraction, probably determined by ndxref above */
/* z is subscript of vertical component of v; 0=x, 1=y, 2=z */
#define  MINZ  0.01
#define  MAXP  (cr/MINZ)
vect3	r;
    r = norm3vect (v);
/*    if (z < 0 || z > 2) return r;  */
/*    if (r[z] < 0.01) r[z] += 0.1;  */
/*    else  r[z] += (cr / r[z]);  */
    switch  (z) {
	case 0:
	    if (r.x < MINZ)  r.x += MAXP;
	    else  r.x += (cr / r.x);
	    break;
	case 1:
	    if (r.y < MINZ)  r.y += MAXP;
	    else  r.y += (cr / r.y);
	    break;
	case 2:
	    if (r.z < MINZ)  r.z += MAXP;
	    else  r.z += (cr / r.z);
	    break;
	default:
	    return  r;
    }
    return  norm3vect (r);
}
/* Note refracted is not same refraction approximation used in
  the difref program above; yet. */


/* ---------------------------------------------------------------------- */

/* Precession things */

vect3	prec (double e, vect3 pf, double ef)
/* Result is pf at epoch ef precessed to epoch e. */
{
/* Compute vector position <p> at epoch <e> from position
   <pf> at epoch <ef>.
   Rigorous method -- precessional motion only.
   Uses IAU (1976) Precession constants
        (Astronomy and Astrophysics vol. 73, 282-284 (1979) )
   Epochs are in Julian years.  JE = 2000.0 + (JED - 2451545.0)/365.25
   Method is E.S.A.E. pages 30-31.   */
vect3	p;
double	t, tn;
static	double	zeta = 0.0;
static	double	z    = 0.0;
static	double	thta = 0.0;
static	double	es   = 0.0;
static	double	efs  = 0.0;

/* If same epochs as last time, skip the angle finding... */
    if ( (e != es) || (ef != efs) ) {
	t = (e - ef) * 1.0e-2;
	tn = (ef - 2000.0) * 1.0e-2;
     
	zeta = t * (11.1808609e-3 + tn*(6.77071e-6 - tn*674.e-12)
		    + t * (1.463556e-6 - tn*1.668e-9 + t*87.257e-9));
	z = zeta + t*t*(3.843603e-6 + tn*1.988e-9 + t*994.e-12);
	thta = t * (9.7171735e-3 - tn*(4.13692e-6 + tn*1.052e-9)
		    - t * (2.06846e-6 + tn*1.052e-9 + t*202.812e-9));
	es  = e;
	efs = ef;
    }

    p = rotate (pf, zeta, 3);
    p = rotate (p, -thta, 2);
    p = rotate (p, z, 3);
    return  p;
}

void	pprec (double *a, double *d, double e,
	double a1, double d1, double e1)  {
/* Precess from position a1, d1 at e1 to *a, *d at e.  Positions are
  in radians, and epoch in Julian years */
vect3	p;
    p = prec (e, sph2vec (a1, d1), e1);
    vec2sph (a, d, p);
}

void	hdprec (double *a, double *d, double e,
	double a1, double d1, double e1)  {
/* Precess from position a1, d1 at e1 to *a, *d at e.  Positions are
  in hours and degrees, and epoch in Julian years */
vect3	p;
    p = prec (e, had2vec (a1, d1), e1);
    vec2had (a, d, p);
}


/* ---------------------------------------------------------------------- */

/* (Simulated) Matrix things */

/* A matrix is an array of 3 3-vectors, or 2 2-vectors.  It is passed  by
  a pointer reference, and can be transposed, multiplied, copied, etc.  */

vect3	mmult3 (vect3 *a, vect3 v) {
vect3	result;
    result.x = dot3vect (a[0], v);
    result.y = dot3vect (a[1], v);
    result.z = dot3vect (a[2], v);
    return  result;
}

vect3	tmult3 (vect3* a, vect3 v) {
/* Multiply transpose of a by vector v */
vect3	result;
    result.x = a[0].x * v.x + a[1].x * v.y + a[2].x * v.z;
    result.y = a[0].y * v.x + a[1].y * v.y + a[2].y * v.z;
    result.z = a[0].z * v.x + a[1].z * v.y + a[2].z * v.z;
    return  result;
}


vect2	mmult2 (vect2 *a, vect2 v) {
vect2	result;
    result.x = dot2vect (a[0], v);
    result.y = dot2vect (a[1], v);
    return  result;
}

void	trpose3 (vect3 *a) {
/* Transpose a 3-matrix in place; this will invert a 3-d
  rotation matrix. */
double	t;
    t = a[0].y;   a[0].y = a[1].x;   a[1].x = t;
    t = a[0].z;   a[0].z = a[2].x;   a[2].x = t;
    t = a[1].z;   a[1].z = a[2].y;   a[2].y = t;
}

void	trpose2 (vect2 *a) {
/* Transpose a 2-matrix in place; this will invert a 2-d
  rotation matrix. */
double	t;
    t = a[0].y;   a[0].y = a[1].x;   a[1].x = t;
}

void	mcpy3 (vect3 *a, vect3 *b) {
/* Copy matrix b to location a */
int	i;
    for (i=0; i<3; i++) a[i] = b[i];
}

void	mcpy2 (vect2 *a, vect2 *b) {
/* Copy matrix b to location a */
int	i;
    for (i=0; i<2; i++) a[i] = b[i];
}

void	pmatrix (vect3 *p, double e, double ef) {
/* Make a 3-matrix p which will cause precession transformation to
  equinox e from equinox ef when multiplied by an equinox ef position. */
    p[0] = prec (e, make3vect (1.0, 0.0, 0.0), ef);
    p[1] = prec (e, make3vect (0.0, 1.0, 0.0), ef);
    p[2] = prec (e, make3vect (0.0, 0.0, 1.0), ef);
/* Transform the matrix to obtain correct multiplication.  */
    trpose3 (p);
/* To transform vector v1 at epoch e1 to vector v2 at epoch e2, do:
    vect3  pm[3];
    pmatrix (pm, e2, e1);
    v2 = mmult3 (pm, v1);       */
}

void	r2matrix (vect2 *a, double angle) {
/* Form a 2-d rotation matrix */
#ifdef  HAVE_SC
    sincos (angle, &(a[1].x), &(a[0].x));
#else
    a[0].x = cos (angle);
    a[1].x = sin (angle);
#endif
    a[0].y = -a[1].x;
    a[1].y =  a[0].x;
}

void	mxm3a (vect3* a, vect3* b, vect3* c) {
/* a = b x c, matrix wise... */
    a[0].x = b[0].x * c[0].x + b[0].y * c[1].x + b[0].z * c[2].x;
    a[0].y = b[0].x * c[0].y + b[0].y * c[1].y + b[0].z * c[2].y;
    a[0].z = b[0].x * c[0].z + b[0].y * c[1].z + b[0].z * c[2].z;
    a[1].x = b[1].x * c[0].x + b[1].y * c[1].x + b[1].z * c[2].x;
    a[1].y = b[1].x * c[0].y + b[1].y * c[1].y + b[1].z * c[2].y;
    a[1].z = b[1].x * c[0].z + b[1].y * c[1].z + b[1].z * c[2].z;
    a[2].x = b[2].x * c[0].x + b[2].y * c[1].x + b[2].z * c[2].x;
    a[2].y = b[2].x * c[0].y + b[2].y * c[1].y + b[2].z * c[2].y;
    a[2].z = b[2].x * c[0].z + b[2].y * c[1].z + b[2].z * c[2].z;
}

void	mxm3 (vect3* a, vect3* b, vect3* c) {
/* a = b x c, matrix wise... */
int	i;
    for (i=0; i<3; i++) {
	a[i].x = b[i].x * c[0].x + b[i].y * c[1].x + b[i].z * c[2].x;
	a[i].y = b[i].x * c[0].y + b[i].y * c[1].y + b[i].z * c[2].y;
	a[i].z = b[i].x * c[0].z + b[i].y * c[1].z + b[i].z * c[2].z;
    }
}

/* To both precess and project, we make a precession matrix b and
  a projection matrix a, and then the combined matrix c:
    smatrix (a, v)	-- projection, done second
    pmatrix (b, e, ecat)  -- precession, done first
    mxm3 (c, a, b)	-- combined precession and projection
.. */

/* ---------------------------------------------------------------------- */

/* Real Matrix things */

/* A (real) matrix is a structure, of 2 vect2 or 3 vect3 elements,
  labeled appropriately.  Functions will return a matrix value,
  with all the usual (needed) operations.  */


/* -- 2-dimensional Matrix things -- */

vect2	mmul2v (matrix2  mx,  vect2  v) {
/* Return the vector product of a 2-matrix and a 2-vector  */
vect2	result;
    result.x = dot2vect (mx.x, v);
    result.y = dot2vect (mx.y, v);
    return  result;
}

matrix2	mtr2 (matrix2 a) {
/* Return the transpose of a 2-matrix */
matrix2	result;
    result.x = make2vect (a.x.x, a.y.x);
    result.y = make2vect (a.x.y, a.y.y);
    return  result;
}

matrix2	minv2 (matrix2 a) {
/* Return the inverse of a 2-matrix */
matrix2	result;
double	d;
    d = a.x.x * a.y.y - a.x.y * a.y.x;
/* or d = cross2vect (a.x, a.y) */
    if (d == 0.0) {	/* Singular, return identity */
	result.x.x = result.y.y = 1.0;
	result.x.y = result.y.x = 0.0;
    } else {
	result.x.x =  a.y.y / d;
	result.y.y =  a.x.x / d;
	result.x.y = -a.x.y / d;
	result.y.x = -a.y.x / d;
/*	result.x.y = -a.y.x / d;  */
/*	result.y.x = -a.x.y / d;  */
    }
    return  result;
}

matrix2	mul2m (matrix2 a, matrix2 b) {
/* Standard matrix multiplication; direct. */
matrix2	result;
    result.x.x = a.x.x * b.x.x + a.x.y * b.y.x;
    result.x.y = a.x.x * b.x.y + a.x.y * b.y.y;
    result.y.x = a.y.x * b.x.x + a.y.y * b.y.x;
    result.y.y = a.y.x * b.x.y + a.y.y * b.y.y;
    return  result;
}

matrix2	rot2m (double a) {
/* Create a 2-d rotation matrix for angle a */
matrix2	result;
#ifdef  HAVE_SC
    sincos (a, &(result.y.x), &(result.x.x));
#else
    result.x.x = cos (a);
    result.y.x = sin (a);
#endif
    result.x.y = -result.y.x;
    result.y.y =  result.x.x;
    return  result;
}


/* -- 3-dimensional Matrix things -- */

vect3	mmul3v (matrix3  mx,  vect3  v) {
/* Return the vector product of a 3-matrix and a 3-vector  */
vect3	result;
    result.x = dot3vect (mx.x, v);
    result.y = dot3vect (mx.y, v);
    result.z = dot3vect (mx.z, v);
    return  result;
}

matrix3	mtr3 (matrix3 a) {
/* Return the transpose of a 3-matrix */
matrix3	result;
    result.x = make3vect (a.x.x, a.y.x, a.z.x);
    result.y = make3vect (a.x.y, a.y.y, a.z.y);
    result.z = make3vect (a.x.z, a.y.z, a.z.z);
    return  result;
}

/* PLEASE NOTE that for an orthonormal (rotation) matrix, the inverse
  and the transform are theoretically identical.  */

matrix3	minv3 (matrix3 a) {
/* Return the inverse of a 3-matrix, direct computation. */
matrix3	result;
matrix3 c;	/* The conjugate matrix */
double	d;	/* Determinant */

/* Find the conjugate matrix using cross products */
    c.x = cross3vect (a.y, a.z);
    c.y = cross3vect (a.z, a.x);
    c.z = cross3vect (a.x, a.y);
    d = dot3vect (a.x, c.x);	/* Determinant; y and z rows work too */
    if (d == 0.0) return rot3m (0.0, 0);	/* Return identity */
    d = (d == 0.0) ? 1.0 : 1.0 / d;	/* Scale by 1/d */
    result.x = mul3vect (c.x, d);
    result.y = mul3vect (c.y, d);
    result.z = mul3vect (c.z, d);
    return  mtr3 (result);	/* Actual result is transpose here */
}

matrix3	mul3m (matrix3 a, matrix3 b) {
/* Standard matrix multiplication, implemented by vectors */
matrix3	result;
matrix3	t;
    t = mtr3(b);
    result.x = mmul3v (a, t.x);
    result.y = mmul3v (a, t.y);
    result.z = mmul3v (a, t.z);
    return  mtr3(result);
}

matrix3	rot3m (double angle, int axis) {
/* Form rotation matrix which when multiplied by a vector will have
  the identical effect as the rotate subroutine */
/* If axis is not 1,2,3 the identity matrix is returned. */
matrix3	result;
double	s, c;
#ifdef  HAVE_SC
    sincos (angle, &s, &c);
#else
    s = sin (angle);
    c = cos (angle);
#endif
    switch (axis) {
	case 1:
	    result.x = make3vect (1.0, 0.0, 0.0);
	    result.y = make3vect (0.0,   c,  -s);
	    result.z = make3vect (0.0,   s,   c);
	    break;
	case 2:
	    result.x = make3vect (  c, 0.0,   s);
	    result.y = make3vect (0.0, 1.0, 0.0);
	    result.z = make3vect ( -s, 0.0,   c);
	    break;
	case 3:
	    result.x = make3vect (  c,  -s, 0.0);
	    result.y = make3vect (  s,   c, 0.0);
	    result.z = make3vect (0.0, 0.0, 1.0);
	    break;
	default:
	    result.x = make3vect (1.0, 0.0, 0.0);
	    result.y = make3vect (0.0, 1.0, 0.0);
	    result.z = make3vect (0.0, 0.0, 1.0);
    }
    return  result;
}

matrix3	precmat (double e, double ef) {
/* Return the 3-matrix p which will cause precession transformation to
  equinox e from equinox ef when multiplied by an equinox ef position. */
/* Update of pmatrix routine above. */
matrix3	result;
    result.x = prec (e, make3vect (1.0, 0.0, 0.0), ef);
    result.y = prec (e, make3vect (0.0, 1.0, 0.0), ef);
    result.z = prec (e, make3vect (0.0, 0.0, 1.0), ef);
/* This actually gives the transform of what is desired, so... */
    return  mtr3(result);
/* To transform vector v1 at epoch e1 to vector v2 at epoch e2, do:
    v2 = mmul3v (precmat(e2, e1), v1);       */
}

matrix3	precmat3 (double e, double ef) {
/* Return a precession matrix directly formed from precession constants. */
matrix3	result;
double	t, tn;
static	double	zeta = 0.0;
static	double	z    = 0.0;
static	double	thta = 0.0;
static	double	es   = 0.0;
static	double	efs  = 0.0;

/* Angles found exactly as in the prec subroutine */
/* If same epochs as last time, skip the angle finding... */
    if ( (e != es) || (ef != efs) ) {
	t = (e - ef) * 1.0e-2;
	tn = (ef - 2000.0) * 1.0e-2;
     
	zeta = t * (11.1808609e-3 + tn*(6.77071e-6 - tn*674.e-12)
		    + t * (1.463556e-6 - tn*1.668e-9 + t*87.257e-9));
	z = zeta + t*t*(3.843603e-6 + tn*1.988e-9 + t*994.e-12);
	thta = t * (9.7171735e-3 - tn*(4.13692e-6 + tn*1.052e-9)
		    - t * (2.06846e-6 + tn*1.052e-9 + t*202.812e-9));
	es  = e;
	efs = ef;
    }

/* Method 1...
    result = mul3m (rot3m (z, 3), rot3m (-thta, 2));
    result = mul3m (result, rot3m (zeta, 3));
.. */

/* Method 2 */
    result = mul3m (rot3m (-thta, 2), rot3m (zeta, 3));
    result = mul3m (rot3m (z, 3), result);

/* Method 3 ...
    result = rot3m (zeta, 3);
    result = mul3m (rot3m (-thta, 2), result);
    result = mul3m (rot3m (z, 3), result);
..*/

/* All 3 methods give the same result. */

    return  result;
}


/* ---------------------------------------------------------------------- */

/* R.A. overlap utilities */

double	fdifmod (double a, double b, double c) {
/* Difference limited to +/- c/2 */
double	result;
double	h;
    result = a - b;
    h = fabs (c / 2.0);
    while (result >  h) result -= c;
    while (result < -h) result += c;
    return  result;
}

double  fmodsum (double a, double b, double c) {
/* Sum, put into range 0<= x <c */
double	result;
double	h;
    result = a + b;
    h = fabs (c);
    while (result >=  h) result -= c;
    while (result < 0.0) result += c;
    return  result;
}

double	fabsdif (double a, double b, double c) {
/* Absolute value of difference in closest direction */
    return  fabs ( fdifmod(a, b, c));
}


/* ---------------------------------------------------------------------- */

/*  --  Arithmetic  --  */

int	Nint (double a)
/* Return nearest integer */
{
#ifdef  __sun
    return  nint(a);
#else
    return  (int)round(a);
#endif
}

/* ---------------------------------------------------------------------- */

/* -- Simple Time Programs -- */

double  jdnow (void)
/* Return current time as J.D. */
{
static long int        time_value;
/* static double	sysepoch = 2440587.5; */
/* 2440587.5 is J.D. at epoch of system time, 0h GMT on 1 Jan., 1970 */
/* SysEpoch is now defined in the kdcutil.h file. */
double          d;

    time (&time_value);		/* Seconds since the SysEpoch */
    d = (double)time_value / 86400.0;	/* Convert time to days */
    return (d + SysEpoch);	/* Convert clock time directly to J.D. */
}

/* Note that the returned JD and its associated GMT will be correct;
  however, local time would need the current time zone. */


void	jd2date (double jd, int* year, int* month, int* day)
/* Convert jd to Gregorian year, month and day (GMT date) */
{
long int        j, y, d, m;

    j = (long)(jd-0.5) - 1721118;	/* De-bias JD */
    y = (4*j-1)/146097;		/* Century */
    j = 4*j - 1 - 146097*y;	/* Quarter days in century */
    d = j/4;			/* Days in century */
    j = (4*d+3)/1461;		/* Year */
    d = 4*d + 3 - 1461*j;	/* Quarter days in year */
    d = (d+4)/4;		/* Days */
    m = (5*d-3)/153;		/* Month */
    d = 5*d - 3 - 153*m;	/* 1/5 days in month */
    d = (d+5)/5;		/* Day in month */
    y = 100*y + j;		/* Combine year, century */

    *year = y + (m+2)/12;	/* Correct year */
    *month = 1 + ((m+2) % 12); 	/* Correct month */
    *day = d;			/* Day of month */
    return;
}

double	jdate ( int year, int month, int day)
/* Convert Gregorian y-m-d to an actual J.D. at 0h GMT */
{

long int  j, c;
long int  jy, jm;

    month += 9;
    jy = year + ((month/12) - 1);
    jm = month % 12;

    c = jy/100;
    jy -= 100*c;
    j = (146097*c)/4 + (1461*jy)/4 + (153*jm+2)/5 + day + 1721118;
    return (double)j + 0.5;
}

/*  ===  Definitions for Julian Day conversions  ===  */
/* #define  JULIAN_YEAR  365.25  */
/* #define  JULIAN_EPOCH 1721045.0  */
/* #define  JYEAR(jd)  ((jd - JULIAN_EPOCH)/JULIAN_YEAR)  */
	/* JYEAR converts a julian date to standard epoch, years */
/* #define  JDEPOCH(epoch) (JULIAN_EPOCH + epoch*JULIAN_YEAR)  */
	/* JDEPOCH converts a standard epoch to julian date */
/* These are converted to subroutines below... */

double	Jyear (double  jd)
/* Convert a J.D. to a Julian Year epoch */
{
double	y;
    y = (jd - 1721045.0) / 365.25;
    return  y;
}

double	jdaye (double  year)
/* Convert a Julian Year epoch to a J.D. */
{
double	y;
    y = 1721045.0 + year * 365.25;
    return  y;
}


double	gst0 (double  jd)
/* Return the GST at 0 hours for the current J.D.  Just add the
  local solar time to get sidereal time.  Result in days. */
/* Also, lst = gst + east longitude */
{
double  tc;
double  gst;

    tc = (jd - 2451545.0) / 36525.0;
/* tc is time in Julian Centuries from the standard epoch 2000.0 */
    gst = 24110.54841
	+ tc*(8640184.812866
	+ tc*(0.093104
	- tc*6.2e-6));	/* GMST - UT1, seconds */
    gst /= 86400.0;	/* ST dif, days */
    gst -= (int)gst;
    if (gst < 0.0) gst += 1.0;	/* gst in days now */
    return  gst;
}


/* ---------------------------------------------------------------------- */

/* -- Time Stamp support -- */

// static	char	exchar (int c)
char	exchar (int c)
{
/* -- Returns the extended character coresponding to integer c;
  Should be in the ASCII sort compare sequence.  Normal return is
  0-9, A-Z, a-z.  A negative integer returns "."; there are
  62 values numbered 0-61; integers over 61 return "~" which is
  larger than any of the standard characters.  */

    if (c <  0) return '.';
    if (c < 10) return '0'+c;
    c -= 10;
    if (c < 26) return 'A'+c;
    c -= 26;
    if (c < 26) return 'a'+c;
    return '~';
}

/* --- Inverse extended character function... */
static	int	charex (char c)
{
/* -- Return an int which would result in c from exchar;
  The '.' returns -2, and '~' returns 63;
  Character values less than '.' give -3, between '.' and '0' give -1.
  Other unknown characters return 62.  '0'-'9', 'A'-'Z', 'a'-'z' give
  their proper extended values.  */
    if      (c <  '.') return -3;
    else if (c == '.') return -2;
    else if (c <  '0') return -1;
    else if (c <= '9') return c - '0';
    else if (c <  'A') return 62;
    else if (c <= 'Z') return 10 + c - 'A';
    else if (c <  'a') return 62;
    else if (c <= 'z') return 36 + c - 'a';
    else if (c == '~') return 63;
    else  return 62;
}

/* -- Time stamp formatting, from a unix system time... */
static	char*	tmstamp (time_t tx)
{
static	char	buffer[12];	/* Returns the value  */
struct	tm	*ut;

/* Find the ut time in a structure  */
    ut = gmtime (&tx);

/* Format the time into buffer  */
    sprintf (buffer, "%02d%02d%02d.%c%c%c",
	ut->tm_year % 100,
	ut->tm_mon + 1,
	ut->tm_mday,
	exchar(ut->tm_hour),
	exchar(ut->tm_min),
	exchar(ut->tm_sec)  );

    return  buffer;
}

/* -- A timestamp function, returning buffer with current time... */
char*	timestamp (time_t* tt)
{
time_t	now;

/* Find the current time, and return a time string for it, using
  the more general tmstamp routine above.  */
    time (&now);
    if (tt != NULL) *tt = now;
    return  tmstamp (now);
}

time_t	unstamp (char* ts)
{
/* Translate a time stamp, and provide the unix time; inverse of tmstamp */
int	y, m, d;
int	k;
char	ch, cm, cs;
int	days;
time_t	rt;		/* Result  */

    k = sscanf (ts, "%02d%02d%02d.%c%c%c",
	&y, &m, &d, &ch, &cm, &cs);
    if (k != 6) return -1;

/* Since there is no direct unix inverse to gmtime (mktime has the local
  time zone folded in, and is dodgy around daylight time changes) we find
  the unix time from our existing JD subroutines... */
    if (y < 70) y += 100;
    y += 1900;
    days = jdate (y, m, d) - SysEpoch;	/* Since System Epoch  */
    rt = days * 24 + charex (ch);	/* hours  */
    rt =  rt  * 60 + charex (cm);	/* minutes  */
    rt =  rt  * 60 + charex (cs);	/* seconds  */
    return  rt;
}


/* ---------------------------------------------------------------------- */

/* -- Linear Transforms -- */
/* Forward:  v = x*s + z;
   Inverse:  x = (v-z)/s;  */

double	tranfw (double x, double zp, double sf)
/* Forward transform */
{
    return (zp + x * sf);
}

double	tranrv (double x, double zp, double sf)
/* Inverse transform */
{
    if (sf == 0.0) return 0.0;
    return ( (x-zp) / sf);
}



/* ---------------------------------------------------------------------- */


/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

