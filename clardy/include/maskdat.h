/*  maskdat.h -- Definitions of mask data structures */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
**                                                                   **
**                   Copyright (C) 2006 by                           **
**     Ken Clardy, Carnegie Observatories, Pasadena California.      **
**                                                                   **
**      This software is proprietary, and may be used or copied      **
**          only with the written consent of the author.             **
**              All copies must retain this notice.                  **
**                                                                   **
** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef  INCLUDE_MASKDAT_H
#define  INCLUDE_MASKDAT_H

/*  ====  Global Macros  ====  */
/* Support for the data types defined here... */

#define  C_loop(q,a)  for(q=(a); q!=NULL; q=(q->next!=(a))?q->next:NULL)
/* C_loop (indexpointer, headpointer)
  Forms a "for" loop over all elements of the (possibly circularly)
  linked list, using the "next" element pointer.  Loop terminates with
  the index pointer NULL if a NULL pointer is encountered, or the
  head pointer is encountered in a "next" link.  Use of "continue"
  or "break" will work normally.  The head pointer may be modified,
  so long as the "next" pointer referencing it is also modified.      */


/* --NOTE-- problem with "passband" structure.  Defined in and needed
  in the optics definitions header, it is also used here in this
  header.  We probably should include the optics header to get it; or
  put it in a more general header used by both them and us */
/* -- optics.h header information is now here. */

/* May want to put vector definitions from kdcutils here... */

#include	"ioutils.h"
#include	"kdcutil.h"
#include	"mgfeats.h"

/*  ===  Optics Definitions  ===  */

#define  Angstrom  1.0e-7
#define  Micron  1.0e-3
#define  Inch	25.4
#define CAM_ANGLE	( (180.0 + 45.0) * Degree )


/* Other fundamental structures */

typedef	struct	intval {	/* ordered interval */
	double	lo;
	double	hi;
	} intval;

typedef	struct	bbox {		/* bounding box */
	intval	x;
	intval	y;
	} bbox;

/*  === Optics data structures from optics.h === */
/* optics.h is now obsolete, its definitions are here. */

/* passband applies to filters, detector sensitivity, spect. limits */
typedef  struct  pb {
	double	blue;	/* Shortest wavelength */
	double	red;	/* Longest wavelength */
	} passband;

//typedef  struct  filt {
//	struct  filt  *next, *last;
//	passband	pass;
                                    //	int		ident;		/* Internal identification */
//	char		name[32];
//	} filter;


/* detector must be defined in the detector plane, with
 pixels/mm and borders, areas between live elements, etc. */
typedef  struct  det {
//	passband	sens;
	bbox		bb;	/* Edges */
	vect2		pcount;	/* Pixel count each direction */
	vect2		psize;	/* Pixel sizes each direction */
	double		vrad;	/* Vignetting Radius */
/* This may change to a more detailed map of pixel space in the
  linear coordinates */
	bbox		bbs;	/* Bounds of nod&shuffle detector */
	bbox		bbn;	/* Bounds of normal detector */
	bbox		bbd;	/* Bounds of default detector */
	} detector;

/* Optics definitions for "optutils.h" here */

/*  ====  Optics data structures  ====  */

/* Generalized polynomial */

typedef struct poly {
        int     order;
        double  *coef;
        } poly;
/* coef will contain order+1 values */

typedef  struct  datlist {
	struct datlist	*next;
	char	*data;
	struct element *eg;
	} datlist;
// NOTE that datlist is now specific to element
// instrument queues.


// element IS SIMILAR TO xlist STRUCTURE IN optix.h
/* Structure(s) for the gmap ini data stuff... */
typedef  struct  element {
	struct element  *next, *last;
	char	*name;
	int	kw;	// Type of element. (a Keyword value)
	int	flag;		/* The default flag, maybe others */
	void   *data;
	struct  datlist  *head;
// CHANGE THIS TO POINTER TO (VOID) DATA STRUCTURE ?
	} element;

/* Bit values for "flag" word: 
	1 = Default element (not used here)
	2 = Alternate (x) axis modifier
..*/

typedef	struct	stgrism {
	double	lpmm;
	double	angle;
	char*	glass;
	element*  eg;
	} stgrism;
/* The glass string gets matched to a glass entry containing the
  index polynomial when this is used. */
/* Optionally, we could add a pointer to the polynomial, and fill it
  in after reading. */

typedef	struct	egrism {
	double	inangle;
	char	*inglass;
	element	*ig;
	double	gangle;
	double	lpmm;
	char	*eglass;
	element	*eg;
	double	eangle;
	} egrism;
/* The enhanced grism has a lead angle and glass, then a
  grating at some other angle, and an exit glass, ending
  at some exit angle.  Either glass may be null.  Zero angle
  is normal to coordinate system, usual input is negative
  angle and exit is positive angle.  */

typedef	struct	stechel {
	vect3	pris1;
	vect3	pris2;
	vect3	grnor;
	vect3	grdis;
	double	lpmm;
	char*	glass;
	element*  eg;
#ifdef  GTILTMOE
	double	gtilt;	// debug - tilt angle for grating.
	vect3	gtnor;	// debug tilted normal
	vect3	gtdis;	// debug tilted dispers
#endif
	} stechel;
/* Done similarly to grism with vector components
  rather than defining angles.  */

/*  (( Might need the focdat structure ?  ))  */

typedef	struct	gapdat {
	struct	gapdat	*next, *last;
	char	*name;
	char	*instr;
	int	ns;
	double	x1u, x1b;
	double	x2u, x2b;
	} gapdat;
/* Since the gaps need not be orthogonal to coordinates, they can't
  be stored as a bounding box.  The x-u values are at the upper edge
  of the detector, and the x-b at the bottom edge.  */

#define  GAPHEAD  "GapHead"
// GAPHEAD is name of the special gap header element


/* ... End of optutils data structures */

/*  =====  Mask space data  =====  */

typedef struct  avoids  {
        struct avoids  *next, *last;
//      vect2   c;
        double  x, y;
        double  r;
        }  avoids;

// ( This is from cututil.h )

typedef struct  {
        double  curve;  // diopters
        double  beam;   // diameter, mm.
        double  slew;   // rate mm/min
        double  cut;
        double  fine;
        double  afoc;   // autofocus range
        double  dx, dy, dz;
	double	flange;
	double	aangle;
	double	aclear;
	double	LongSlit;
	double	SlitGap;
	double	tool;
	double  cte;
	double	tcut;
        avoids  *av;    // avoidance queue
	int	mma;
        }  cdata;
/* cdata contains all the cutting parameters for any given
  instrument.  Typically used in an array by instrument, or
  passed as a (pointer) subroutine argument when used.  */

/* The object data structure -- data from object list */
typedef  struct  objdat {
	char	*name;		/* Name, dynamic string */
	double	ra, dec;	/* Celestial coordinates */
	double	priority;	/* For optimization */
	namlist	*precom;	/* Pre-comments */
	namlist	*postcom;	/* Post-comments */
	} objdat;

enum	slitshape {
	CIRCLE,
	SQUARE,
	RECTANGLE,
/* Add special, complex shapes here... */
	SPACED_CROSS
	};

typedef  struct  slit {
	int	shape;		/* See the shape enum */
	double	width;		/* Width or only dimension */
	double	alen;		/* Length to left */
	double	blen;		/* length to right */
	double	angle;		/* Position angle */
//	double	curve;		/* Curvature in diopters */
	} slit;

// NOTE -- These conflict types are yet unused
enum	conflict_types {
	NO_CONFLICT,
	OD_Y,
	OD_X,
	INT_SPECT,
	INT_IMAGE,
	DUPLICATE,
	BAD_SPOT
	};

/* Conflict node */
typedef  struct  conf_tag {
	struct  conf_tag  *next, *last;
	int	type;
	struct	conf_tag  *bak, *fwd;	/* 2-way per-object list */
	struct  objq  *obj;	/* Primary object back-pointer */
	struct  objq  *other;	/* Other object */
	} cfl;

/* Spectrum edge descriptor */
typedef  struct  edge_tag {
	vect2	p;	/* Position of center */
	vect2	r;	/* Red end */
	vect2	u;	/* blue end */
	double	a, b;	/* Coeficients of curve */
	} sp_edge;

/* Spectrum on Detector descriptor */
/* Used for detail conflict resolution computation */
typedef  struct  sod_tag {
	struct  sod_tag  *next;	/* Allow pop-list for images */
	int	order;		/* The order of this spectrum */
/* The spectrum edges, stored as above */
	sp_edge	e1, e2;
/* The bounding box limits for full compares */
	bbox	bb;
/* Second bounding box for on-detector test */
	bbox	bd;
	vect2	center;    /* DEBUG object at cw location */
	int	on_det;	/* see if any part is on the detector */
	} spect;

/* Object values -- data for the object's spectrum */
typedef  struct  objval {
/* Put here the projected spectra positions on detector as vect2 values */
	/* Pointer to main spectrum traced boundry */
	/* Queue to ghost image spectra */
	int		ncf;	/* Conflict count */
	cfl		*cfq;	/* Conflict queue pointer */
	spect		*sod;	/* Primary spectrum */
	spect		*img;	/* Other order image list */
	bbox		dimg;	/* Direct image limits */
	int		stat[2];	/* Extension status values */
	int		flag;	/* DEBUG flag, normally 0 */
	} objval;


enum  object_types {
	OBJ_UNKNOWN,
	OBJ_OBJECT,
	OBJ_REFERENCE,
	OBJ_GAP
	};

/*  Bits in objq.flag */
#define  OBJECT_FLAG    1
#define  OBJECT_ACTIVE  4
#define  OBJECT_SPECIAL 8
#define  OBJECT_VIRTUAL 16
#define  OBJECT_SERVICE (1<<5)

#ifdef  OBJCHECK
/* Bits in cflag check flag */
#define  OC_NEGUSE (1<<0)
#define  OC_PRIX   (1<<1)
#define  OC_NSHORT (1<<2)
#define  OC_PREQ   (1<<3)
#define  OC_WIDTH  (1<<4)
#define  OC_WIDEX  (1<<5)
#define  OC_SHAPE  (1<<6)
#define  OC_ALEN   (1<<7)
#define  OC_BLEN   (1<<8)
#define  OC_LENGTH (1<<9)
#define  OC_ANGLE  (1<<10)
#define  OC_UNCUT  (1<<11)
#endif

/* A queue of objects read in, and having slits assigned */
typedef  struct  objq {
	struct  objq	*next, *last;
	objdat		dat;
	vect2		smpos;	/* Position on the mask (was vect3) */
	slit		slit;
	objval		*val;	/* Conflict data */
	int		type;
	int		flag;
	int		use;	/* Use count */
	int		order;	/* Ordering to restore list */
	vect2		smp1, smp2;	/* For long slit ordering */
#ifdef  NNA
	struct	objq	*nn;	/* Neighbor pointer */
	struct	objq	*nx;	/* Neighbor pointer */
#endif
#ifdef  OBJCHECK
	int		cflag;	/* Check bits flag */
#endif
	} objq;

/* Slit edge for extension algorithm: */
typedef	struct	sleg_tag {
	struct	sleg_tag   *next, *last;
	objq	*objec; 	/* The object */
	objq	*confl;		/* The conflict */
	int	side;		/* A or B side, 0 or 1 index */
	double	free;		/* Free space (for sorting) */
	} sleg;

/* Status bits for stat[side] */
#define	Stat_Free  1
#define	Stat_Prov  2
#define	Stat_Mark  4

/*  ===  Observation Data  ===  */

//#ifdef  GSDATA
typedef	struct	gpos_tag {
	double	ra;	/* Right Ascension, in hours */
	double	dec;	/* Declination, in degrees */
	double	ep;	/* Epoch, years */
	int	cn;	/* Catalog Number */
	}  gpos;
//#endif

/* Observation flags (bits in Obs.Flag) */
#define  OFxrd   1
#define  OFpb2   2
#define  OFifu   4
#define  OFsxt   8
#define  OFmoe  16
#define  OFral  (1 << 11)
/* xrd = extra orders;  pb2 = passband 2 active;
  ifu = IFU mode on;  sxt = slit extend feature;
  moe = Multiple Object Echellette in use;
  ral = Remove ALignments for detector/gap;   */

/* Detector status flags (bits in Obs.Dstat) */
#define  DS_ns     1
#define  DS_gaps   2
#define  DS_orient 4

/* Observation structure -- the observing parameters */
typedef  struct  obs {

/* Identification */
	char	oname[16];	/* Observer name (required) */
	char	fname[8];	/* File name (required) */
	char		*title;		/* Observation title */

/* Celestial position and instrument orientation */
	double		ra, dec;	/* Field center */
	double		equinox;	/* Coordinate system */
	double		epoch;	/* Observation date in years and decimal */
	double		angle;	/* Instrument position angle */
	double		temp;	/* temperature expected */
	double		ha;	/* Hour angle expected */
//	vect3		mpos[3];	/* Position transform matrix */
	matrix3		mxpos;		/* Position transform matrix */

/* Guiders */
//#ifdef  GSDATA
	gpos		gp1;
	gpos		gp2;
//#endif

/* Wavelength space items */
	passband	pb;	/* Wavelength limits */
	passband	pb2;	/* Limits for on-detector */
	double		cw;	/* Center wavelength */
	double		sw;	/* Spectrum center (mean) wavelength */

/* Objects */
	namlist		*comments;	/* Global comments */
	namlist		*objfiles;	/* Object list files input */
	objq		*ob;	/* Object queue header */
	int		reuseob;	/* Re-use of objects */
	int		reuserf;	/* Re-use of reference */
	int		refsel;		/* Reference selector */
	int		reflimit;	/* Reference count limit */
	double		Pmusthave;	/* Must Have priority */
	double		Pdecide;	/* Decision parameter priority */

/* Mask space */
	element		*Tel_scope;
	slit		default_slit;
	slit		reference_slit;

/* Instrument */
	element		*Instr_mt;
	element		*Direct_inst;

/* Disperser */
	element		*Disperser;
	int		order;		/* Disperser order */
	double		D_angle;	/* Disperser (grating) angle */
	passband	ddir;	/* Dispersion direction */
/* definition - multiply ddir.blue by slit.width and add to 
  object position for the blue edge of spectrum; same for red;
  value stored is +0.5 or -0.5 depending on dispersion direction
  relative to slit mask coordinate direction. */
/* actually, ddir is not really used yet... */

/* Detector */
	detector	*Detect;
	int		Dstat;		/* Detector and gap status */
//	filter		*filter;

/*  New optical references... */
	element		*elist;		// Element list pointer

/* Features */
	int		ex_order;	/* Extra order conflicts flag */
	int		dlevel;		/* Debug message level */
	double		minsep;		/* Minimum separation spec. */
	double		dca, dcb;	/* Do not cut l/r in arcsec */
	int		pb2act;		/* Detector limit flag */
	int		IFUmode;	/* The IFU feature */
	double		IFUra, IFUdec;	/* Offset Center position */
	int		Flag;		/* MOE, later slext, exorder, etc. */
	int		slext;		/* Slit extension flag */
	namlist		*pcom;		/* Rotator flip, etc. comments */
	double		cref;		/* Refraction coeficient */
	vect2		zenith;		/* Zenith (radians) in mask xy */
	int		date;		/* Date as mjd for intgui */
	} Obs;


/* We need the grating used, passband limit, overlap used, filter,
and perhaps default slit parameters to be recorded for observing
time and reference */

/* A slit mask data structure for output (and input) containing
 all the data in the SMDF */
typedef  struct  slitmask {
	char	observer[16];	/* Standard identification */
	char	maskfile[8];
	Obs	obs;		/* Observing parameters */
/* Obs structure contains field center, eq, ep, pa values for observ'n */
	namlist	*comments;
/* type or structure containing size limits, curvature, etc. */
	objq	*slitlist;	/* List of slit (object) data blocks */
/* When writing, write only active objects */
/* When reading, the objval stuff is NULL pointer */

	} slitmask;
// above is unused in programs (yet?)


/*  ====  */
#endif     /*  INCLUDE_MASKDAT_H  */
/*  ====  */
