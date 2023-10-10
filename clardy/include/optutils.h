/*  optutils.h  ==  Generalized Optical Transform definitions  */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
**                                                                   **
**                  Copyright (C) 2005 by                            **
**     Ken Clardy,  Carnegie Observatories,  Pasadena California.    **
**                                                                   **
**      This software is proprietary, and may be used or copied      **
**          only with the written consent of the author.             **
**              All copies must retain this notice.                  **
**                                                                   **
** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ---------------------------------------------------------------------- */

#ifndef	INCLUDE_OPTUTILS_H
#define	INCLUDE_OPTUTILS_H

/* ---------------------------------------------------------------------- */

// #ifdef IS_OPT_C
// #define EXT_OPT  
// #else
// #define EXT_OPT extern
// #endif

/* INCLUDEs ------------------------------------------------------------- */

#include	"maskdat.h"
/* maskdat.h includes ioutils.h, kdcutil.h */
/* both are needed here */


/* DEFINEs -------------------------------------------------------------- */

/* This is defined elsewhere, should be incluced... (mgutils.c) */
/* Standard shorthand for an allocation */
#define  xalloc(b)  ( (b*) malloc (sizeof(b)) )

/* All 3 below are defined in maskdat.h by including optics.h */
/* #define  Angstrom  1.0e-7 */
/* #define  Micron  1.0e-3 */
/* #define  Inch  25.4 */


/* TYPEDEFs ------------------------------------------------------------- */


// PLEASE NOTE, THAT SOME OF THESE MAY DUPLICATE DEFINITIONS
// WHICH ARE AVAILABLE ELSEWHERE.  PERHAPS IN IOUTILS, WHICH
// WILL BE NEEDED HERE AND SHOULD BE INCLUDED...

/* Structure definitions */

/* THESE STRUCTURES MUST BE RECONCILED WITH THOSE IN optix.h
  AND optix.c FILES; THE GENERALIZED OPTICS USE CODE.  */

/*  NOTE -- Several structures have been moved to "maskdat.h" for
  common access to data types from other programs, which may not
  be using the optutils.c subroutines.  */

typedef	enum	Etype {
// SHOWS TYPE OF DATA BLOCK USED
	Void,
	Focuser,
	Grism,
//	DGrism,
	EGrism,
	Echelle,
	Instrument,
	Glass,
	FloatValue,
	Secondary,		// For level 2 keywords - change name
	Detector,
	Gap,
	Other
	} Etype;


// SIMILAR TO THE vmod STRUCTURE IN optix.h
typedef  struct  focdat {
	double	scale;
	poly	*rmod;
	poly	*tmod;
	poly	*wmod;
	poly	*curv;
	} focdat;


// typedef  struct  kwlist {
// 	struct kwlist	*next;
// 	char	*key;
// 	char	*data;
// 	} kwlist;



// The datlist pointer is used for instrument list, and for
// any other additional text data.
// The "flag" can contain information on the axis if needed, with
// a 0 bit for standard (y-axis) and 1 for non-standard (x-axis).

/* The element pointer will be used as argument to get-keyword value
  subroutines.  It will initialize the datlist pointer from its head
  value, and the datlist will be searched for keywords in the same
  way that files are read in the "older" version of the program. */


typedef enum  Keyword {
        VOID,
        END,
        Telescope,
        FOCuser,
        Colimator,
        Reimager,
        GRAting,
        GRIsm,
//	DGRIsm,
	EGRIsm,
	ECHelle,
        Reflection,
        ANgle,
        Filter,
        DETector,
        INSTrument,
        FOCLEN,
	SCale,
        FOCCURV,
        Field,
	Pixels,
        AXIS,
/*        ROTATION,  */
        Lines,
	GLASS,
	RMOD,
	TMOD,
	WMOD,
	ALIGNROT,
	PRIS1,
	PRIS2,
	GNORM,
	GDISP,
	GTILT,	/* debug value */
	GAP
        } Keyword;

/* GLOBALs -------------------------------------------------------------- */

extern	int	mdebug;	// global debug stuff for optutils.c

extern	int	dodebug;	// global for compatability.
extern	double*	GRangle;		// Grating angle[s]
extern	int*	GRorder;		// Grating or Grism order[s]
extern	element	*Cur_grating;
extern	element *Cur_grism;

#ifdef  OTDEBUG
extern	double	Cam_angle;	// special debug value...
extern	int	llong;		// Another debug thing
extern	int	direct;		// Another debug thing
extern	double	cur_den;	// Also debug for ntest
extern	double	grindex;	// Also debug for ntest
#endif

// Initialization of all globals appears in optutils.c

/* FUNCTION PROTOTYPEs -------------------------------------------------- */

// DATA TYPE UTILITIES
Keyword		Decode_Keyword (char* key);
//char*		Keyname (Keyword k);
char*		typename (Etype k);
Etype		type_of (Keyword k);
element*	newelement (char* name);
focdat*		newfdat (void);
datlist*	newdat (char* line);

// DESTRUCTOR-LIKE ROUTINES
datlist*	freedat (datlist* p);
void		freedatlist (datlist* p);
void		freepoly (poly* p);
void		freestgrism (stgrism* p);
//void		freesdgrism (sdgrism* p);
void		freexgrism (egrism* p);
void		freestechel (stechel* p);
void		freefdat (focdat* p);
element*	pop_element (element* p);
element*	push_element (element* b, element* a);
element*	free_element (element* p);


// THESE ARE ACTUALLY IN IOUTILS
// int     intvalue (char*  field);
// double  floatvalue (char*  field);

double*	newfv (char* c);

// THE PARSING MIGHT BE PUT INTO IOUTILS
// char*	skipover (char* line, char* delim);
// char*	parse (char* line, char* token, char* delim);
vect3	readv3 (char* text);

// THE POLYNOMIAL THINGS MIGHT BE PUT ELSEWHERE
poly*	read_poly (char* text);
void	showpoly (poly* p, FILE* f);
double  eval (poly*  p, double arg);

void*	optstore (char* line, void* head);

// ELEMENT QUEUE UTILITIES
element*	find_element (element* head, char* name);
//int		count_elements (element* head);
int		count_types (element* head, Keyword kw);

// ERROR HANDLING FOR OPTICAL DATA
// void	error_title (int te);
int	find_errors (element* head);

// Optical computations
vect3	snell (vect3 a, vect3 s, double ndx);
vect3	diffract (vect3 a, vect3 s, vect3 g, double v);

// The general programs...
double  Gr_Angle (double lambda, double density, double Camangle, int order);

vect3	Op_transform (vect3 A, element* p, double wavl, double temp);

vect3   Inv_transform (vect3 B, element* p, double wavl, double temp);

// Disperser Utilities
double	grwavc (element* te, double cur_temp);

/* ---------------------------------------------------------------------- */
/* End of optutils.h package */

#endif	/* INCLUDE_OPTUTILS_H */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
