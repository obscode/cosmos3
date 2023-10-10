/*  mgutils.c == Mask Generation Utility Subroutines */

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

/* These utility routines are entered into the IMACS library for
  general use.  */

/* %*%*% -- Use to find genlist.c include section */

/*  ===  Standard Includes  ===  */
#include        <stdio.h>
#include        <string.h>
#include        <stdlib.h>
#include	<math.h>
#include	<ctype.h>	/* For "toupper" */
#include	<time.h>	/* For clock - debug use */

#undef  debug

#undef  debug_slex

// //#define  RLSFULL
// #undef  RLSFULL
/* Symbol RLSFULL allows the full computation of reordering using
  actual end points of each slit.  If not on, actual start point is
  used for both start and end of the slit. */
// Moved to mgfeats.h for future coordination with cututils.c ...

#define  Tname  "Magellan"
// "Magellan" has the ADC, "MagNoADC" does not.

#define  MODULE_NAME  "MGutils"
const	char	MGUNAME[] = "Mask Gen Utils";
const	char	MGHIVD[]  = "1";
const	char	MGUVERS[] = ".46";


/*  ===  Local Includes  ===  */

// #include	"mgfeats.h"
// Features now included by maskdat.h

//#include	"mgutils.h"
#include	"mgutils.h"

// #include	"mobtem.h"
// Include the monthly Tav array

/* -- see below -- #include	"genlist.c"  */

/* Note:  mgutils.h includes all of ioutils, kdcutils, optutils, so
  they don't need to be included here. */

/*  -- Timing utility and time poiners... */
#include	"timepak.c"
#define  TX_size  30

// Assignment of TX subscripts to easily change them
// At least the early ones are done here, move useful ones
// here when stable...
#define  TS_whole   1	// whole thing
#define  TS_sorts   2	// sorting in slit exten
#define  TS_freet   3	// free space computations
#define  TS_freeb   4	// Free space binary search
#define  TS_conxs   5	// conflict-x searches

#define  TS_initc   6	// Initial conf-x		**unused
#define  TS_augms   7	// Augment side in conf-x	**unused
#define  TS_geteg   8	// Get edge subroutine

#define  TS_confx   9	// conf-x full search

/* Currently unused symbols */
#define  TS_autst  10	// augment side, total
#define  TS_augxo  11	// augment side extra order
#define  TS_augsb  12	// augment side boundaries
#define  TS_copsd  13	// copy sod
#define  TS_inifc  14	// initial find conflict


//#define  TS_dotime	// Timer output on/off
#undef   TS_dotime	// Timer output on/off


/*  ===  Definitions  ===  */

#define max(a, b)	((a) < (b) ? (b) : (a))
#define min(a, b)	((a) > (b) ? (b) : (a))
#define abs(x)		((x) >= 0 ? (x) : -(x))
/* Above should be in general.h ? from <macros.h> */
/* Perhaps they could be in mgfeats.h maybe? */

#define	true  1
#define false 0
/* --NOTE-- Find where to define above things !! */

#define  ever  ;;
/* for the for (ever) loop structure... */

#define  XordLim 20
/* XordLim is how many extra orders we allow as a limit */

#define  Latitude  (-29.01418)
#define  Pressure  575.0
/* Latitude and atmospheric pressure in mm Hg are defined
  statically here for Las Campanas, Magellan telescopes.
  When needed, these quantities are per-telescope things which
  should be indexed by telescope name.   */

// #define  DefTemp  20.00390625
/* Default temperature is 20 + 1/256 degrees, an exact value */
// #define  DefTemp  11.99609375
#define  DefTemp  12.0390625
/* Better suggestion is 12 - 1/256 degrees, also exact */
/* Actual average is 12.28437 (weighted) or 12.03818 -- use nearest
  1/256 fraction to make it exact.  Use 12 + 73/256 or 12 + 10/256
  i.e. 12.28515625 and 12.0390625 */
// MOVE THIS TO MASKDAT, then use it in maskcut to turn off the
// CTE use if it exactly matches.  Also, initialize in defobs
// and in the gui reading of the obs data.
// Best would be to use the maskcut flag to turn CTE on/off.


/*  ===  Global and static values  === */
static int	refcount = 0;
double	D_maskrad = D_maxradi;	// Global here...

FILE	*dbfile;		// Actually present here.

int	Nfaroff = 0;	/* Far off mask messages */
int	Noutdet = 0;	/* Outside detector messages */
int	Nduppos = 0;	/* Duplicate mask positions */
int	Nconflx = 0;	/* Conflicts reported */


static objq	*Tobj = NULL;	// Test-only object (for add_default)
/* The test-only object pointer is used in slit extension, where it is
  NOT desired to create a conflict queue entry for the test object, and
  this is to be communicated to the conflict queue subroutines. */

/*  ===  Debug function prototypes === */

/* These debug things are used from msdebug.c, which is
  included in the main program.  */
// void	dump_conf (objq *ob);
// void	check_cfcount (objq* ob, char* title);
// Above now included in this program for linking ease.


	/* - - - - - - - - - - *\
	|			|
	|   Generalized  List	|
	|			|
	\* - - - - - - - - - - */


/*  ====  Generalized List Handling  ====  */
/* The following is "genlist.c", a logically separate file, which is
  reproduced here rather than included as a temporary measure to make
  this package complete for source distribution.  */

/*  %*%*%*%*%*%*%  Include of genlist.c  %*%*%*%*%*%*%  */
/*  genlist.c == General list handling code  */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
**                                                                   **
**                   Copyright (C) 2005 by                           **
**     Ken Clardy, Carnegie Observatories, Pasadena California.      **
**                                                                   **
**      This software is proprietary, and may be used or copied      **
**          only with the written consent of the author.             **
**              All copies must retain this notice.                  **
**                                                                   **
** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// --- %$%$% Start genlist.h %$%$%
// This section is code for genlist.h

/* #include	<stdlib.h>	// For structure, NULL */
/* From /usr/include/stdlib.h -- */
#ifndef NULL
#define NULL    0
#endif

/*  ====  General Data Structure  ====  */

typedef	struct	g_link {
	struct	g_link	*next, *last;
	} g_link;
/* Please note that the order of last, next must be the same in
this and all covered structures, and that these pointers must
come first in the structure.  */
/* It is recommended that the next pointer come first for
compatability with structures which are singly linked. */

/*  ====  Function Prototypes  ====  */

/* See %#%#% function prototype section in mgutils.h */

// --- %$%$% End genlist.h %$%$%


/*  ====  General List Handling  ====  */

void*	G_pop_c (void* x)
/* Remove structure from circular queue, return pointer to next */
{
g_link	*a = (g_link*)x;
g_link	*p, *q;
    if (a == NULL) return a;
    p = a->next;
    if (p == a) p = NULL;
    q = a->last;
    if (q == a) q = NULL;
    if (p == NULL) p = q;
    if (p != NULL) p->last = q;
    if (q != NULL) q->next = p;
    a->last = a->next = a;
    return p;
}

void*	G_pop_l (void* x)
/* Remove structure from linear queue, return pointer to next */
{
g_link	*a = (g_link*)x;
g_link	*p, *q;
    if (a == NULL) return a;
    p = a->next;
    if (p == a) p = NULL;
    q = a->last;
    if (q == a) q = NULL;
    if (p != NULL) p->last = q;
    if (q != NULL) q->next = p;
    a->last = a->next = NULL;
    if (p == NULL) p = q;
    return p;
}

void*	G_pop (void* x, int c)
/* Remove structure from queue, return pointer to next */
/* c is true if circularly linked */
{
g_link	*a = (g_link*)x;
g_link	*p, *q;
    if (a == NULL) return a;
    p = a->next;
    if (p == a) p = NULL;
    q = a->last;
    if (q == a) q = NULL;
    if (p == NULL && c) p = q;
    if (p != NULL) p->last = q;
    if (q != NULL) q->next = p;

    if (c)
	a->last = a->next = a;
    else
	a->last = a->next = NULL;

    if (p == NULL) p = q;
    return p;
}



void*	G_head (void* a)
/* Find head of queue, circular or linear */
{
g_link	*p;
    if (a == NULL) return a;
    for (p=a; p->last != NULL; p = p->last) {
	if (p->last == a) return a;
// If circular, return a as the head of the list.
    }
    return  p;
}

void*	G_tail (void* a)
/* Find tail of queue, circular or linear */
{
g_link	*p;
// if circular, use a->last; if linear, follow path
    if (a == NULL) return a;
    for (p=a; p->next != NULL; p = p->next) {
	if (p->next == a) break;
// break if found circular, return a is not exactly right.
// If circular, return the value before a as the tail.
    }
    return  p;
}


// push into queue, possibly 2 ways
// Scan list to check for proper links (circular and linear)
// and optionally repair?
// circularize a linear list


void	G_circ (void* a)
/* Check for linear list a, and make circular if so */
{
g_link	*p, *q;
    if (a == NULL) return;
    p = G_head (a);
    q = G_tail (a);
    if (p->last != NULL && q->next != NULL) return;	// redundant?
    if (p->last == NULL) p->last = q;
    if (q->next == NULL) q->next = p;
/* Should (recursively?) loop until both head and tail are found
  to be without NULL pointers, although logically that can happen
  only once?  */
}


void*	G_push_l (void* a, void* b)
/* Put linear list a before linear list b and return head */
{
g_link	*p, *r;
    if (b == NULL) return a;
    if (a == NULL) return b;
/* Scan to tail of a (make that p) a..p - 0 */
    p = G_tail (a);
/* Scan to head of b (make that r) 0 - r..b */
    r = G_head (b);
/* Hook up ...p - r... */
    p->next = r;
    r->last = p;
/* Find the new head and return it... */
    return  G_head (a);
}

void*	G_push_cl (void* y, void* x)
/* Put queue or element b (y) into queue a (x), so as to preceede the
 element a, and return the new head (a (x)) pointer.  */
/* b is head of queue inserted before a */
{
g_link	*a = (g_link*)x;
g_link	*b = (g_link*)y;
g_link	*p, *r;

/* Find queue ends such that a...p and b...r are merged, into
  something like  ...p - b...r - a... ; this replaces the links
  ...p - a... and ...r - b... with the p - b and r - a links. */
/* If a has null predicessor, we link ...0 - b...r - a...; if
  b has null predicessor, find the next end and use that. */
/* Take care of null cases */
    if (a == NULL) return b;
    if (b == NULL) return a;
/* Link it up */
    p = a->last;		/* ...p - a... */
    if (p == NULL) p = G_tail (a);
    r = b->last;		/* ...r - b... */
    if (r == NULL) {	/* Scan to end */
	for (r = b; r->next != NULL; r = r->next) {
	    if (r->next == b) break;
	}
    }
    p->next = b;
    b->last = p;
    r->next = a;
    a->last = r;		/* ...p - b...r - a... */
    return a;
}


void*	G_push_c (void* y, void* x)
/* Put element b (y) before element a (x), and other circularly linked
  items get linked together.  Return new head (a) pointer. */
/* b is tail of list inserted before a */
{
g_link	*a = (g_link*)x;
g_link	*b = (g_link*)y;
g_link	*p, *r;

/* Find queue ends such that a...p and r...b are merged, into
  something like  ...p - r...b - a... ; this replaces the links
  ...p - a... and ...b - r... with the p - r and b - a links. */
/* If a has null predicessor, or b a null follower, circularize
  that list and proceed as before.  */
/* Take care of null cases */
    if (a == NULL) return b;
    if (b == NULL) return a;
/* Link it up */
    if (a->last == NULL) G_circ (a);
    p = a->last;		/* ...p - a... */
    if (b->next == NULL) G_circ (b);
    r = b->next;		/* ...b - r... */
    p->next = r;
    r->last = p;
    b->next = a;
    a->last = b;		/* ...p - r...b - a... */
    return a;
}

int	listcount (void* x)
/* Return number of list nodes in a general list -- count nodes */
{
int	n=0;
g_link	*a;
    C_loop (a, x) n++;
    return  n;
}

/*  %*%*%*%*%*%*%  End Include of genlist.c  %*%*%*%*%*%*%  */

/*  ===  Memory Allocation utilities  ===  */

#define  xalloc(b)  ( (b*) malloc (sizeof(b)) )

//	b = xalloc (objq);
//	memset (b, 0, sizeof(objq));	/* CLEAR TO ZERO */

#define  zalloc(b) ( (b*) calloc (1, sizeof(b)) )
// allocate and return structure set to zero.

/* Another way to do this would be this:  */
// #define  balloc(b)  ( (b*) mxclr (sizeof(b)) )
// void*	mxclr (size_t b)  {
// void*	p;
//     p = malloc (b);
//     memset (p, 0, b);
//     return  p;
// }


	/* - - - - - - - - - - - - - - *\
	|				|
	|  Interval and Bounding Box	|
	|				|
	\* - - - - - - - - - - - - - - */

/*  ===  Interval and Bounding Box support  ===  */

// This may be relocated into a more general utility
// code module and included above...

intval	interval (double a, double b) {
/* Return interval with the two values */
intval	r;

    if (a < b) {
	r.lo = a;
	r.hi = b;
    } else {
	r.lo = b;
	r.hi = a;
    }
    return r;
}

int	inintv (double a, intval v) {
/* Compare a to interval, result is 0 if a is outside low,
  1 if inside (inclusive of end points) and 2 if outside high */
/* Change to "within" when the old one is obsoleted */
    if (a < v.lo) return 0;
    if (a > v.hi) return 2;
    return 1;
}

int	ovr_lap (intval a, intval b) {
/* Find overlap state of intervals; 0 for no overlap,
  1 if overlap without inclusion,
  2 if a within b, 4 if b within a.  */
/* Change to "overlap" when the old one is obsoleted */

const int	rval[3][3] = { {0,1,4}, {1,2,1}, {4,1,0} };
int	i, j;

    i = inintv (a.lo, b);
    j = inintv (a.hi, b);
    return  rval[i][j];
}

int	bbover (bbox a, bbox b) {
/* Compare bounding boxes for overlap, inclusion.  Returns
  0 if no overlap, 1 if intersected,
  2 if a is included within b and 4 if b is included within a */
int	i, j;
    i = ovr_lap (a.x, b.x);
    if (i == 0) return 0;
    j = ovr_lap (a.y, b.y);
    if (j == 0) return 0;
    if (i == j) return i;
    return 1;
}

bbox	boundbox (intval x, intval y) {
/* Form a bounding box from the x and y extent intervals */
bbox	r;
    r.x = x;
    r.y = y;
    return r;
}

bbox	boundbv2 (vect2 a, vect2 b) {
/* Form bounding box from coordinates of 2 points */
bbox	r;
    r.x = interval (a.x, b.x);
    r.y = interval (a.y, b.y);
    return  r;
}

intval	intex (intval v, double p) {
/* Extend interval v to include point p */
intval	r;
    r = v;
    if (p < r.lo) r.lo = p;
    if (p > r.hi) r.hi = p;
    return  r;
}

bbox	boundex (bbox b, vect2 p) {
/* Extend bounding box b to include point p */
bbox	r;
    r.x = intex (b.x, p.x);
    r.y = intex (b.y, p.y);
    return  r;
}

bbox	boundor (bbox a, bbox b) {
/* Form inclusive bounding box combining a and b */
bbox	r;
    r = a;
    if (b.x.lo < r.x.lo) r.x.lo = b.x.lo;
    if (b.x.hi > r.x.hi) r.x.hi = b.x.hi;
    if (b.y.lo < r.y.lo) r.y.lo = b.y.lo;
    if (b.y.hi > r.y.hi) r.y.hi = b.y.hi;
    return  r;
}

bbox	bbedge (sp_edge e) {
/* Form a bounding box from a spectrum edge */
bbox	r;
    r = boundbv2 (e.r, e.u);
    r = boundex (r, e.p);
    return  r;
}

static	bbox	bbspec (spect* s) {
/* Form a bounding box from a spectrum */
bbox	r;
    r = boundor (bbedge(s->e1), bbedge(s->e2));
    return  r;
}

double	intlen (intval x) {
/* Return length of interval */
    return  fabs (x.hi - x.lo);
}


/*  ===  Maskgen specific I/O  ===  */

// THIS SHOULD BE MOVED TO IOUTILS AS IT IS AN IO UTILITY
// TO SELECT A FILE IN GLOBAL DATA AREAS...
// NO, NO, NO.  THE FINDFILE STUFF IS SPECIFIC TO MASKGEN...

void	data_file (char** fspec, char* filename, char* ev, int mf)
/* Used by mask generation software to locate a filespec from
  a given filename and a given environment name.  The filespec is
  returned in the provided buffer, and may be null if nothing
  was found.  */
{

/* First, look in the local directory, and parents */
    findfile (fspec, "r", filename, NULL);
/* Then, look at recommended environment variable. */
    findfile (fspec, "r", filename, ev);
/* If still not found, look in the default directories. */
    findfile (fspec, "r", filename, "/usr/local/etc/maskcut");
    findfile (fspec, "r", filename, "/usr/local/etc/maskgen");
    findfile (fspec, "r", filename, "/usr/local/etc");
/* All out of ideas.  If not found by now, it is lost. */

/* Give a utility message if file is not found, and mf is on. */
    if (fspec == NULL  &&  mf > 0) {
	printf (" ** Could not find \"%s\" anywhere.\n", filename);
	printf ("  Please define %s as the proper directory.\n", ev);
    }
}

	/* - - - - - - - - - - *\
	|			|
	|   Get Optical Data	|
	|			|
	\* - - - - - - - - - - */

/*  ===  Subroutine to read optical data  ===  */

element*  get_optic_data (char* filename, int dl)
{
char	*opfspec = NULL;
element	*gs = NULL;	/* Optic details queue pointer */
// ==getdata==


/* For debug, when feature is on, try for the test version first. */
#ifdef  GAPVO
    findfile (&opfspec, "r", "optictest.dat", NULL);
#endif

    data_file (&opfspec, filename, "OPTICDEF", 1);

//    if (opfspec == NULL) {
//	printf (" ** Could not find \"%s\" anywhere.\n", filename);
//	printf ("  Please define OPTICDEF as the proper directory.\n");
//    } else {

    if (opfspec != NULL) {
/* Here's where we read the optical definitions... */
	if (dl > 0)
	    printf (" Optic data from \"%s\".\n", opfspec);	// debug

/* Get the optics data from the optics data file into an element queue. */
	mdebug = 0;		// Turn off optutils debug messages.
	gs = (element*) get_data (
	    opfspec,		// File to be read
	    *optstore,	// User subroutine supplied here
	    WHITESPACE,	// Delimiter character set
	    "!#",		// Initial character for comment
	    "!",		// Introduces comment anywhere
	    "#;",		// Comment if preceeded by delimiter
	    "",		// Comment if followed by delimiter
	    "" 		// Comment if pre & follow by delimiter
	    );

	if (gs == NULL && dl >= 0) {
	    printf (" ** Could not read \"%s\".\n", opfspec);
	}

/* Add globals required.  Will need to put in real values later */

// The global pointers are assigned when obs structure is loaded,
// otherwise they will be NULL
//	if (GRangle == NULL) GRangle = xalloc (double);
//	*GRangle = 45.0 * Degree;
//	if (GRorder == NULL) GRorder = xalloc (int);
//	*GRorder = 2;

	free (opfspec);
    }
    return  gs;
}


	/* - - - - - - - - - - *\
	|			|
	|   Instrument  Name	|
	|			|
	\* - - - - - - - - - - */

/*  ===  Globalization  ===  */

Indx	InstIndx (char* name)
/* Return an instrument index from Indx enum for a given string. */
{
int	l;
Indx	k = -1;
    if (name == NULL) return k;
    l = strlen (name);
    if (l < 1) return k;
    if      (!strncasecmp(name, "IMACS",	min( 5,l))) k = I_imacs;
    else if (!strncasecmp(name, "Centerfield",	min(11,l))) k = I_cfg;
    else if (!strncasecmp(name, "LDSS-2",	min( 6,l))) k = I_lds2;
    else if (!strncasecmp(name, "LDSS2",	min( 5,l))) k = I_lds2;
    else if (!strncasecmp(name, "LDSS3",	min( 5,l))) k = I_lds3;
    else if (!strncasecmp(name, "LDSS-3",	min( 6,l))) k = I_lds3;
    else if (!strncasecmp(name, "LDSSNEW",	min( 7,l))) k = I_lds3;
    else if (!strncasecmp(name, "Other",	min( 5,l))) k = I_other;
    else  k = -1;
    return  k;
}
/* Above is the only valid use of the alternate ldss-3 names */

/* The extended index, mostly differing in distinguishing long/short
  imacs cameras, is not used.  See instindx.otl  */



char*	InstName (Indx i)
/* Return pointer to literal name of the index value. */
{
static	char	*c;
    switch (i) {
	case I_imacs:	c = "IMACS";		break;
	case I_cfg:	c = "Centerfield";	break;
	case I_lds2:	c = "LDSS-2";		break;
	case I_lds3:	c = "LDSS-3";		break;
	case I_other:	c = "Other";		break;
	default:	c = "Unknown";		break;
    }
    return  c;
}

	/* - - - - - - - - - - *\
	|			|
	|   Queue  utilities	|
	|			|
	\* - - - - - - - - - - */


/*  ===  Object Queue manipulation utilities  ===  */
/* Push and pop queue routines */

objq*	pop_obj (objq* a)
/* Remove structure from queue, return pointer to next */
{
objq	*p, *q;
    if (a == NULL) return a;
    p = a->next;
    if (p == a) p = NULL;
    q = a->last;
    if (q == a) q = NULL;
    if (p == NULL) p = q;
    if (p != NULL) p->last = q;
    if (q != NULL) q->next = p;
    a->last = a->next = a;
    return p;
}

spect	*kill_spect (spect *s)
{
spect	*n;
    if (s == NULL) return s;
    n = s->next;
    free (s);
    return  n;
}

objval	*kill_objval (objval *v)
/* Kill the extra storage for a permanently removed object */
{
spect	*s;
cfl	*c;	// had..  *d;
    if (v == NULL) return v;
    for (s = v->sod; s != NULL; ) s = kill_spect (s);
    for (s = v->img; s != NULL; ) s = kill_spect (s);

//    for (c=v->cfq; c != NULL; ) {
//	c = pop_cfl (d=c);
//	free (d);
//    }

/* The above is wrong.  If for some reason the cfq pointer is
  not null, we have a problem.  Before this is called, the object
  which contains v should have been set inactive.  If that is not
  done, pointers which point to v (usually the per-object chain
  from another cfl) could still exist.  If we had a back pointer
  for the secondary chain, we could fix this.  We would pop only
  the current cfl, fix its secondary chain if any, and be done.
  This program is called by kill_obj, and by fill_object (twice).
  The latter 2 are done with inactive object only. */

    c = v->cfq;
    if (c != NULL) {
	pop_cfl (c);	/* Return value of pop_cfl is not used ?? */
/* Ignore the return, which is probably another conflict indexed
  in another node.  We just free this one.  If used properly, there
  will be no conflict anyway, as the node should have been set to
  inactive, which removes the conflict queue item correctly. */
	free (c);
	v->cfq = NULL;
    }
/* The above takes care of any cfl which might have existed;
  it removes it from both queues cleanly, and frees it.  */

    free (v);
    return  NULL;
}

objq	*kill_obj (objq *a)
/* Remove the given object, return the "next" one */
{
objq	*b;
objq	*r;
objval	*v;
int	nc;

    if (a == NULL) return  a;
    b = a;	/* Looks ugly to do this to a, use a copy */
    r = pop_obj (b);	/* The next victim is here */

/* Free each of the dynamic elements in struct b, then b itself */
    if (b->dat.name != NULL) free (b->dat.name);
    kill_name (b->dat.precom);
    kill_name (b->dat.postcom);
/* --NOTE-- now need to free all objval, including conflict queue... */

/* We should insure that kill_objval is only done on an inactive
  node.  Thus, we need to set a inactive above... */
/* But, we will set it inactive here... */
    v = b->val;
    nc = set_inactive (b);
    if (v != NULL) {
	v->ncf -= nc;
	b->val = kill_objval (v);
    }
    free (b);
    return  r;
}

objq	*kill_objq (objq *a)
/* Recursively remove all object storage, ideally, then returns NULL */
{
objq	*b;
    b = a;
    while (b != NULL)  b = kill_obj (b);
    return  b;
}

objq	*push_obj (objq *b, objq *a)
/* Put the queue or element b into queue a, so that it
  preceeds a.  Return the new head (a) pointer.  */
{
objq *p;
objq *r;

/* Take care of null cases */
    if (b == NULL) return a;
/* Find queue ends such that a...p and b...r are merged, into
  something like  ...p - b...r - a... ; this replaces the links
  ...p - a... and ...r - b... with the p - b and r - a links. */
/* If a has null predicessor, we link ...0 - b...r - a...; if
  b has null predicessor, find the next end and use that. */
/* Link it up */
    if (a == NULL) return b;
    if (b == NULL) return a;
    p = a->last;		/* ...p - a... */
    if (p == NULL) {	/* Scan to end */
	for (p=a; p->next != NULL; p = p->next) {
	    if (p->next == a) break;
// BAD CONSTRUCT BREAK
	}
    }
//    if (p == NULL) p = a;	/* maybe a pointers got nulled ? */
    r = b->last;		/* ...r - b... */
    if (r == NULL) {	/* Scan to end */
	for (r = b; r->next != NULL; r = r->next) {
	    if (r->next == b) break;
// BAD CONSTRUCT BREAK
	}
    }
//    if (r == NULL) r = b;	/* Maybe b pointers got nulled ? */
//    if (p != NULL) p->next = b;
    p->next = b;
    b->last = p;
    r->next = a;
    a->last = r;		/* ...p - b...r - a... */
    return a;
}


/*  ===  Other dynamic storage  ===  */

void	kill_ptr (void** p)
/* Kill pointer storage if not NULL, and clean pointer */
{
    if (p == NULL) return;
    if (*p == NULL) return;
    free (*p);
    *p = NULL;
}

void	clean_obs (Obs* a)
/* Remove contents from an Obs structure */
{
/* Check for trivial case first */
    if (a == NULL) return;

/* Free all dynamic storage in Obs structure */
    kill_ptr ((void*)&(a->title));
    a->elist = NULL;
    a->Tel_scope = NULL;
    a->Disperser = NULL;
// Those are in the element queue, handled separately!
    a->Instr_mt    = NULL;
    a->Direct_inst = NULL;
    kill_ptr ((void*)&(a->Detect));
//    kill_ptr ((void*)&(a->filter));

/* Kill all lists linked to the Obs structure */
    a->comments = kill_name (a->comments);
    a->ob       = kill_objq (a->ob);
    a->objfiles = kill_name (a->objfiles);
    a->pcom     = kill_name (a->pcom);

/* Drop GRangle and GRorder if they point here */
    if (GRangle == &(a->D_angle)) GRangle = NULL;
    if (GRorder == &(a->order)) GRorder = NULL;
}


Obs*	kill_obs (Obs* a)
/* Kill a dynamcially allocated Obs structure; return NULL */
{
/* Check for trivial case first */
    if (a == NULL) return a;

/* Clean the innards of the structure */
    clean_obs (a);

/* Free the obs structure itself */
    free (a);

    return  NULL;
}



/* ====  Debug routines (temporary)  ==== */

void	dbroq (char* title, Obs* ob)
{
objq	*p, *q;
objval	*v;
int	n, a, s, c, t, b;

// See if we do a report...
    if (dbfile == NULL) return;
    if (ob == NULL) return;
// dbfile comes from "DEBUGFILE name" in .obs file
    if (ob->dlevel < 2) return;	// debug level required

// Looks like we do the report.
    n = a = s = c = t = b = 0;
    p = ob->ob;
    C_loop (q, p) {
	n++;				// Count all nodes in queue
	if (q->flag & OBJECT_ACTIVE) a++;	// Active nodes
	v = q->val;
	if (v != NULL) {
	    if (v->sod != NULL) s++;	// sod present
	    if (v->cfq != NULL || v->ncf > 0) {	// some conflicts
		c++;	// conflicts
		if (!(q->flag & OBJECT_ACTIVE)) b++;
	    }
	    t += v->ncf;		// conflict counts
	}
    }		// End of loop through object queue

    fprintf (dbfile, " #@%18s, ob=%d, a=%d, s=%d, c=%d",
	    title, n, a, s, c);
    fprintf (dbfile, " tc=%d, -ac=%d\n", t, b);
    fflush (dbfile);
}


/*  ===  Subroutines for Mask Generation  ===  */


/* Some basic geometry routines */

static	double	pbmean (passband  p)
/* Return mean wavelength of a passband */
{
    return  (0.5 * (p.blue + p.red) );
}


/*  ====  Generalized Storage Support  ====  */

gsto	*new_gsto (void)
/* Create storage, clean and return */
{
gsto	*r;
    r = xalloc (gsto);
    r->next = r->last = NULL;
    r->indx = r->type = 0;
    r->data = NULL;
    return  r;
}

gsto	*kill_gsto (gsto* a)
/* Kill the node a and all its subsidiary storage, and return
  the pointer to where it was removed from its queue.  If needed,
  operate recursively... */
{
gsto	*r;	// What we return
gsto	*b;	// scratch
void	*v, *w;

/* Deal with any null case calls */
    if (a == NULL) return NULL;

/* Rip a from any queue it is in, and remember where for return */
    r = G_pop_l (a);

/* Deal with the data, depending on type */
    switch (a->type) {
	case GS_gsto:
/* For another gsto node, we recursively call ourself */
	    for (b=a->data; b != NULL; ) {
		b = kill_gsto (b);
	    }
	    break;
	case GS_int:
	case GS_double:
/* For scaler items, we just free the storage itself */
	    if (a->data != NULL) free (a->data);
	    break;
	case GS_strx:
/* For a standard queue with no dynamic elements, we make a loop
  to pop and free all the queued elements. */
	    for (v=a->data; v != NULL; ) {
		w = G_pop_l (v);
		free (v);
		v = w;
	    }
	    break;
    }

/* Free the node itself */
    free (a);

/* Return pointer to what is left. */
    return  r;

}

int	get_CDtype (int item)
/* CD specific function, returns the GS type for the CD item */
{

    switch (item) {

/* All the special service items */
// None at the present time...
//	    return  GS_gsto;

/* The double precision scalers */
	case  CD_curve:
	case  CD_beam:
	case  CD_slew:
	case  CD_cut:
	case  CD_fine:
	case  CD_afoc:
	case  CD_dx:
	case  CD_dy:
	case  CD_dz:
	case  CD_flange:
	case  CD_aangle:
	case  CD_aclear:
	case  CD_LongSlit:
	case  CD_SlitGap:
	case  CD_tool:
	case  CD_cte:
	case  CD_tcut:
	case  CD_dmx:
	    return  GS_double;

/* The integer values */
	case  CD_refr:		// We store int subscript there
	case  CD_mma:
	    return  GS_int;

/* Special queue elements */
	case  CD_av:
	    return  GS_strx;

/* Other things... */
    }

    return  -1;
}

static	gsto	*gsfind (gsto* head, int item)
/* Service routine -- look through the gsto queue for the first
  item matching item, and return its pointer */
{
gsto	*g;	// scratch pointer
    if (head != NULL) {
	C_loop (g, head) {
	    if (g->indx == item) return g;
	}
    }
    return  NULL;
}

static	void	gsmakeptr (gsto* a)
/* Service routine -- make a pointer of the appropriate (type) and
  put it in the data field of a.  */
{
// create data - gsto/int/double/void
    if (a == NULL) return;
    if (a->data != NULL) return;
    switch (a->type) {
	case GS_gsto:
/* Actually, we don't want to do this.  We will get the entry into
  the proper place by a push into the queue, modifying the head. */
//	    a->data = malloc (sizeof( gsto* ));
	    break;
	case GS_int:
	    a->data = malloc (sizeof( int ));
	    break;
	case GS_double:
	    a->data = malloc (sizeof( double ));
	    break;
	case GS_strx:
	    a->data = malloc (sizeof( void** ));
	    *( (void**)a->data ) = NULL;
	    break;
/* Otherwise, we would need these:
// void	**pp;	// a dynamically generated pointer pointer
//	    pp = malloc (sizeof( void** ));
//	    *pp = NULL;
//	    a->data = pp;
.. */
    }

}


void*	GS_lookup (gsto* *head, int index, int item, int flag)
/* Returns a pointer, the type of which is dependent on the type
of the item, to the item storage.  If flag is set, storage for the
appropriate type is created, otherwise if flag is not set, and no
storage exists, a NULL is returned.  For scaler items, the pointer
returned points to the data storage; for gsto items it points to
the gsto list head, for queue items it points to a pointer location,
and that pointer is set by the user.  May create *head. */
{
//void*	r=NULL;	// Return value
gsto	*g;	// primary queue element heading secondary queue
gsto	*b;	// secondary queue pointer
gsto	*s;	// item in secondary queue locating data
int	*rx, ref;
int	dtyp;
static	void*  mynull = NULL;

/* If item is not CD_refr, we get reference recursively; if it
 is, we will naturally look it up, using no reference */
    if (item != CD_refr) {
	rx = (int*) GS_lookup(head, index, CD_refr, 0);
	if (rx == NULL) ref = -1;
	else  ref = *rx;
    } else  ref = -1;

/* Find the data type we are using from the item */
    dtyp = get_CDtype (item);

/* Find the secondary queue for this item */
    g = gsfind (*head, item);

/* If no secondary queue, and we are inserting, make a secondary queue
  and insert it in the primary queue.  If not inserting, return NULL. */
    if (g == NULL) {	// No secondary queue found
	if (flag) {	// Insertion mode -- make secondary queue
	    g = new_gsto();
	    g->indx = item;
	    g->type = GS_gsto;	// not dtyp;
	    *head = G_push_l (*head, g);
	} else {	// Not inserting, no secondary, bummer...
	    return NULL;
	}
    }
    b = (gsto*) g->data;
/* b is now the pointer to the secondary queue */

/* See if we have a direct storage in the secondary queue, and
  return it if so. */
    s = gsfind (b, index);	// Secondary queue element, or null

/* If we have it, return the data pointer. */
    if (s != NULL) {		// Direct secondary storage exists
	if (flag)  gsmakeptr (s);	// make it if needed
	return  s->data;
    }
/* After this, we have no direct storage element. */

/* If inserting, make a direct storage element, put it into the
  secondary queue, and return its location. */
    if (flag) {
	s = new_gsto();
	s->indx = index;
	s->type = dtyp;
	gsmakeptr (s);
	b = (gsto*)G_push_l ((gsto*)(g->data), s);	// Insert in queue
	g->data = b;
	return  s->data;
    }
/* After this, we have no direct storage, and are not inserting. */

/* If request is for CD_refr, we have NULL. (no refer or default) */
    if (item == CD_refr) return NULL;

/* Look for reference storage in secondary queue, and if we do, return
  that storage  (reference index feature) */
    s = gsfind (b, ref);
    if (s != NULL) return s->data;

/* As a special thing, return null for avoid request, and for
  anything which should not be global here.   */
    if (item == CD_av) return  &(mynull);

/* Finally, see if we have any storage in the secondary queue, and if
  we do, return the first one found.  (default global feature) */
    if (b != NULL) return b->data;

/* If no storage has been found, return NULL */
    return  NULL;
}



/* === Cutting parameters and avoidance support === */

/* The newavoid should become static when all uses have been resolved */

/* Obtain new avoidance section for given data. */
static	avoids*	newavoid (double x, double y, double r)
{
avoids	*p;
    p = (avoids*) malloc (sizeof(avoids));
    p->next = p->last = NULL;
    p->x = x;
    p->y = y;
    p->r = r;
    return  p;
}


/* A static (local) copy of the traditional routine, where head
  points to the expanded cfig array, an array of cdata structures.
  Data decoded from the file lines is stored into the array.  */
/* If head is NULL, get a pointer to the array and clean the array.
  That pointer is used for the return value. */
/* Keep istrument index as a local static variable. */
/* Use configstore from maskcut.c as a model */

static	void*   cconfsto (char* line, void* head)
/* Data storage for configuration file (former configstore) */
/* Parameters:  line is the file line read, de-commented;
  head is the data pointer.  Returns head as modified.  When called
  with a NULL head, a new data array is obtained.  When called with
  a NULL line, that is a signal for end of file processing.  */
{
#define  KEYSIZE  32
char	keyword[KEYSIZE];
char	name[KEYSIZE];
char	*c;
int	i;
double	x, y, r;
avoids	*ax;
avoids	*px;
int	*mp;
avoids	**ppx;
gsto	*CDhead;
static	int	init=0;
static	int	inst=0;

/* For each case, store data in proper cfig structure. */

/* Check for final call here */
    if (line == NULL) {
double	dmax;
/* Compute the maximum distance to move for cututil and others,
  which replaces the large comutation */
double	*cv;
double	*af;
double	*dm;
	CDhead = head;
	for (i=0; i <= I_other; i++) {
	    cv = (double*)GS_lookup (&CDhead,i,CD_curve,0);
	    af = (double*)GS_lookup (&CDhead,i,CD_afoc,0);
	    dm = (double*)GS_lookup (&CDhead, i, CD_dmx, 0);

//	    *(double*)GS_lookup (&CDhead, i, CD_dmx, 1) =
	    dmax = (cv == NULL || af == NULL || *cv == 0.0) ? 1.0e6 :
		0.75 * sqrt ( fabs (8.0e3 * (*af) / *cv ));

	    if (dm == NULL || dmax != *dm) {
		dm = (double*)GS_lookup (&CDhead, i, CD_dmx, 1);
		*dm = dmax;
	    }

	}
//	D = (cfig[inst].curve == 0.0) ? 1.0e6 :
//	    0.75 * sqrt (fabs(8.0e3 * cfig[inst].afoc / cfig[inst].curve));
	head = CDhead;
	inst = I_other;		// Set instrument index to unused one
	return head;
    }

/* Make sure we have storage defined */
    CDhead = head;

    c = parse (line, keyword, WHITESPACE);

/* Decide basic keywords in a cascaded if statement here */
    if        (!strcasecmp (keyword, "CASE")) {
// find case name, hence the instrument index...

	c = parse (c, name, WHITESPACE);

/* Try to find instrument name... */
	inst = InstIndx (name);

	if (init == 0) {	// Initialize some (global) values
	    *(double*)GS_lookup(&CDhead, inst, CD_curve, 1) =   0.0;
	    *(double*)GS_lookup(&CDhead, inst, CD_beam, 1)  =   0.05;
	    *(double*)GS_lookup(&CDhead, inst, CD_slew, 1)  =2032.0; // 80 ipm
	    *(double*)GS_lookup(&CDhead, inst, CD_cut, 1)   = 254.0; // 10 ipm
	    *(double*)GS_lookup(&CDhead, inst, CD_fine, 1)  = 127.0; //  5 ipm
	    *(double*)GS_lookup(&CDhead, inst, CD_afoc, 1)  =   7.5;
	    *(double*)GS_lookup(&CDhead, inst, CD_dx, 1)    =   0.0;
	    *(double*)GS_lookup(&CDhead, inst, CD_dy, 1)    =   0.0;
	    *(double*)GS_lookup(&CDhead, inst, CD_dz, 1)    =   1.270;
	    *(double*)GS_lookup(&CDhead, inst, CD_flange, 1) =313.05;
	    *(double*)GS_lookup(&CDhead, inst, CD_aangle, 1) =  0.0;
	    *(double*)GS_lookup(&CDhead, inst, CD_aclear, 1) =  0.020;
	    *(double*)GS_lookup(&CDhead, inst, CD_LongSlit,1)= 40.0;
	    *(double*)GS_lookup(&CDhead, inst, CD_SlitGap, 1) = 0.400;
	    *(double*)GS_lookup(&CDhead, inst, CD_tool, 1)   = 15.0;
	    *(double*)GS_lookup(&CDhead, inst, CD_cte, 1)    = 17.0e-6;
	    *(double*)GS_lookup(&CDhead, inst, CD_tcut, 1)   = 20.0;
	    *(int*)GS_lookup(&CDhead, inst, CD_mma, 1)   = 0;
/* Also, set all avoid pointers to null... */

	    init++;
	}


	if (inst < 0) {	/* Unknown case message */
	    printf (" Case \"%s\" not recognized.\n", name);
	    printf ("  ( Valid cases include: ");
	    for (i=0; i < Instruments; i++) {
		if (i > 0) printf (", ");
		printf ("\"%s\"", InstName(i));
		fflush (stdout);
	    }
	    printf (" ).\n");
	}

    } else if (!strcasecmp (keyword, "reference")) {
	*(int*)GS_lookup(&CDhead, inst, CD_refr, 1) = InstIndx (c);

// Representative storage of double keywords...
    } else if (!strcasecmp (keyword, "curve")) {
	*(double*)GS_lookup(&CDhead, inst, CD_curve, 1) = floatvalue (c);
    } else if (!strcasecmp (keyword, "beam")) {
	*(double*)GS_lookup(&CDhead, inst, CD_beam, 1)     = floatvalue (c);
    } else if (!strcasecmp (keyword, "slew")) {
	*(double*)GS_lookup(&CDhead, inst, CD_slew, 1)     = floatvalue (c);
    } else if (!strcasecmp (keyword, "cut")) {
	*(double*)GS_lookup(&CDhead, inst, CD_cut, 1)      = floatvalue (c);
    } else if (!strcasecmp (keyword, "fine")) {
	*(double*)GS_lookup(&CDhead, inst, CD_fine, 1)     = floatvalue (c);
    } else if (!strcasecmp (keyword, "afoc")) {
	*(double*)GS_lookup(&CDhead, inst, CD_afoc, 1)     = floatvalue (c);
    } else if (!strcasecmp (keyword, "flange")) {
	*(double*)GS_lookup(&CDhead, inst, CD_flange, 1)   = floatvalue (c);
    } else if (!strcasecmp (keyword, "Aangle")) {
	*(double*)GS_lookup(&CDhead, inst, CD_aangle, 1)   = floatvalue (c);
    } else if (!strcasecmp (keyword, "Aclear")) {
	*(double*)GS_lookup(&CDhead, inst, CD_aclear, 1)   = floatvalue (c);
    } else if (!strcasecmp (keyword, "LongSlit")) {
	*(double*)GS_lookup(&CDhead, inst, CD_LongSlit, 1) = floatvalue (c);
    } else if (!strcasecmp (keyword, "SlitGap")) {
	*(double*)GS_lookup(&CDhead, inst, CD_SlitGap, 1)  = floatvalue (c);
    } else if (!strcasecmp (keyword, "RTOOL")) {
	*(double*)GS_lookup(&CDhead, inst, CD_tool, 1)     = floatvalue (c);
    } else if (!strcasecmp (keyword, "CTE")) {
	*(double*)GS_lookup(&CDhead, inst, CD_cte, 1)      = floatvalue (c);
    } else if (!strcasecmp (keyword, "Tcut")) {
	*(double*)GS_lookup(&CDhead, inst, CD_tcut, 1)     = floatvalue (c);

    } else if (!strcasecmp (keyword, "PIN") ||
               !strcasecmp (keyword, "LATCH")) {
/* Avoidance item to be put into queue... Read x,y,r */
	c = parse (c, name, WHITESPACE);
	x = floatvalue (name);
	c = parse (c, name, WHITESPACE);
	y = floatvalue (name);
	c = parse (c, name, WHITESPACE);
	r = floatvalue (name);
	ax = newavoid (x, y, r);	// New element to enqueue
	ppx = (avoids**)GS_lookup(&CDhead, inst, CD_av, 1);
//ppx is the pointer to the list address
	px = G_tail (*ppx);
	if (px == NULL) {
	    *ppx = ax;	// New list
	} else {	// Add it (ax) to the end of the list
	    px->next = ax;
	    ax->last = px;
	}
    } else if (!strcasecmp (keyword, "OFFSET")) {
/* An offset for the tool positioning... dx, dy, dz */
	c = parse (c, name, WHITESPACE);
	*(double*)GS_lookup(&CDhead, inst, CD_dx, 1) = floatvalue (name);
	c = parse (c, name, WHITESPACE);
	*(double*)GS_lookup(&CDhead, inst, CD_dy, 1) = floatvalue (name);
	c = parse (c, name, WHITESPACE);
	*(double*)GS_lookup(&CDhead, inst, CD_dz, 1) = floatvalue (name);
    } else if (!strcasecmp (keyword, "Mmode")) {
char	target[] = "Allowed";
	c = parse (c, name, WHITESPACE);
	mp = (int*)GS_lookup(&CDhead, inst, CD_mma, 1);
	*mp = strcmp (name, upcase(target)) ? 0 : 1;	// metric mode allowed
    } else {  // Unknown keyword, complain about it
	printf (" Unknown keyword \"%s\" in cut configuration.\n", keyword);
    }
    return  CDhead;
}


	/* - - - - - - - - - - - - - - *\
	|				|
	|    Get Mask Cutting Data	|
	|				|
	\* - - - - - - - - - - - - - - */


gsto*	get_cut_data (char* filename, int dl)
{
gsto	*ret;
// ==getdata==

char	*fspec = NULL;

/* Program to take as input a file spec of a known existant cutting
  configuration file.  It will return a pointer to a dynamic array
  which holds the decoded data.  It calls a local copy of configstore
  from maskcut.c and others, which stores the data into the array.  */

/* NOTE:  We should take a file name rather than filespec, and do the
  search for the file here, similarly to the way opticdef does it.
  However, also call cconfsto with null value to get the defaults
  put into a basic structure so that some things can proceede such
  as the ncplot program.  May want to return NULL instead to maskcut
  since the parameters are more needed there... */

/* Perhaps a second parameter indicating whether a default set would be
  acceptable; if so, call with dummy case imacs line, and dummy null
  line to obtain a default structure.  Otherwise, give null.  Then,
  maskgen and ncplot would allow defaults, and maskcut would not.  */

    data_file (&fspec, filename, "CUTDEF", 1);

    if (fspec != NULL) {
/* Here's where we read the mask cutting definitions... */
	if (dl > 0)
	    printf (" Mask cutting data from \"%s\".\n", fspec);   // debug

/* Call get_data on the file, and save the pointer returned. */
//	printf (" Reading: \"%s\".\n", fspec);
	ret =  get_data (
	    fspec,		// File to be read
	    *cconfsto,	// User subroutine supplied here
	    WHITESPACE,	// Delimiter character set
	    "!#",		// Initial character for comment
	    "!",		// Introduces comment anywhere
	    "#;",		// Comment if preceeded by delimiter
	    "",		// Comment if followed by delimiter
	    ""		// Comment if pre & follow by delimiter
	    );
	if (ret == NULL && dl >= 0) {
	    printf (" ** Could not read \"%s\".\n", fspec);
	}
	free (fspec);
    }
/* Return that pointer, the file was closed by get_data. */

    return  ret;	// The dynamic array
}


/* -- Geometry of slit:  Find if slit intersects a given circle */

// /* The Incircle is not used.  We will replace it with Xcircle */
// int	Incircle (objq* q, vect2 p, double r)
// /* Return true (1) if slit for object q intersects the circle at
//   point p of radius r.  Otherwise, return 0.
//   Used for intersection decision with
//   avoidance areas around edge for cutting machine.  */
// /* However, we need Rmax to insure slit is entirely contained within
//   a IMACS mask radius. */

/* A subroutine finds intersection of rectangle, described by an end
  point, width and length, with a circle described by a center and
  radius.  Returns -1 if inside, 0 if intersecting, 1 if outside circle. */

static	int	Xcrecd (vect2 s, vect2 a, double w, vect2 p, double r)
/* Find intersection of rectangle with circle.
  s is the location of the center of one side of the rectangle
  a is the vector distance from s to the center of oposite side
  w is the width of the rectangle perpendicular to a
  p is the center of the circle
  r is the radius of the circle
 Return value is as follows:
 -1  Rectangle is entirely contained within circle
  0  Rectangle intersects the circle, or circle is within rectangle
  1  Rectangle is entirely outside the circle
Used to implement Xcircle for rectangular shapes    */
{
vect2	b;	// width vector
vect2	e;	// corner location
vect2	f;	// another corner
vect2	c;	// center displacement
double	d;	// distance

/* Generate the width vector perpendicular to a, length w, call it b */
    b = rot2 (a, M_PI_2);
    d = vect2norm (b);
    b = mul2vect (b, w/d);

/* Find circle center displacement from parallelogram corner, call it c */
    e = sub2vect (s, mul2vect(b, 0.5));
    c = sub2vect (p, e);

/* Find parallelogram distance of circle center, call that d */
    d = dparel (a, b, c);

/* Translate that into intersection value */
    if (d > r) return 1;	// circle outside rectangle
    if (d < 0.0) {	// Center inside rectangle
	d = -d;
	if (r < d) return 0;	// circle within rectangle
    }

// Now, circle is not outside or inside rectangle.  Can intersect, or
// rectangle can be within circle.  If any corner of the rectangle is
// outside the circle we have an intersection, otherwise the rectangle
// is entirely within the circle.  Test all 4 corners in turn.
    if (dist2 (p, e) > r) return 0;
    if (dist2 (p, sum2vect(e,b)) > r) return 0;
    f = sum2vect (e, a);
    if (dist2 (p, f) > r) return 0;
    if (dist2 (p, sum2vect(f,b)) > r) return 0;

    return  -1;	// Must be contained if we havent left yet.
}

int	Xcircle (objq* q, vect2 p, double r)
/* Find intersection of object q with circle of radius r about point p.
  Return values as follows:
  -1  q is entirely contained within the circle
   0  q intersects the circle, or q is null
   1  q is entirely outside the circle
 This routine is used to implement Rmax and in_avoid; also any other uses
 comparing an object with a circular border.     */
{
vect2	h;		// the object position
vect2	v, u;		// unit vectors
//vect2	y;		// side vector
vect2	f;		// position or center side
double	s, c;		// rotation sine, cosine for unit vector
double	d;		// A distance

//double	br = 0.0;	// The biggest radius; what we return.
double	w;		// fundamental width/2
double	a, b;		// alength and blength
int	i;		// Intersection indicator
int	k;		// Intersection sub-indicator
int	j;		// utility counter

/*  How to proceede:
  Circle intersection exists if sum of radii exceeds distance between centers.
  For Rectangle, Square find least distance from point to paralleogram, and
  intersects if that is less than circle radius.
  Spaced-cross is modeled by circle and 4 rectangular slits.    */

/* First, take care of the null and trivial cases */
    if (q == NULL) return 0;

    sincos (q->slit.angle, &s, &c);	// angle is in radians
    h = q->smpos;			// Position of object
    w = q->slit.width / 2.0;		// Fundamental half-width


    switch (q->slit.shape)  {
    case  CIRCLE:
/* Simple, compare sum of radii to center distance */
	d = dist2 (p, h);	// distance between centers
	if (d > (r+w)) return 1;	// outside
	if (d < (r-w)) return -1;	// inside
	return  0;	// Neither, must intersect.
//	break;		// redundant, actually
    case  RECTANGLE:
	a = q->slit.alen;
	b = q->slit.blen;
	break;
    case  SQUARE:
	a = b = w;
	break;
    case  SPACED_CROSS:
	a = -(w + q->slit.blen);
	b = q->slit.alen - a;
/* Slit of alen, spaced blen past circle of radius w. */
	break;		// redundant, actually
    }	// end of case

/* Rectangle is located at (center-a, -w) with sides
  2w and (b+a).  Rotate 4 times and add circle for spaced-cross */
/* The a,b,w vectors are derived from s,c and a,b,w lengths. */
/* Derive these vectors: */

    v = make2vect (c, s);	// Unit vector along slit
    u = sub2vect (h, mul2vect(v,a));	// Location of end
    v = mul2vect (v, b+a);	// Slit vector
    i = Xcrecd (u, v, 2.0*w, p, r);
    if (q->slit.shape == RECTANGLE || q->slit.shape == SQUARE) return i;

/* We have a spaced cross, do the rotations and etc. */
    if (q->slit.shape == SPACED_CROSS) {
	if (i == 0) return  i;
	for (j=0; j<3; j++) {	// Check other 3 slits...
	    f = rot2 (sub2vect (u, h), M_PI_2);
	    u = sum2vect (f, h);
	    v = rot2 (v, M_PI_2);
	    k = Xcrecd (u, v, 2.0*w, p, r);
	    if (k == 0) return  0;	// intersected
	    if (k != i) return  0;	// in to out or out to in, intersect
	}
/* If not returned yet, check the center circle too... */
	d = dist2 (p, h);
	if (d > (r+w) && i > 0) return 1;	// outside
	if (d < (r-w) && i < 0) return -1;	// inside
    }
// We are also returning 0 if shape is unknown or strange...
    return  0;		// intersection is all that's left
}

//#define DEBUGAVOID
#undef  DEBUGAVOID
/* DEBUGAVOID is used for debug things below */

int	avoidanc (objq* ob, avoids* av, double tool)
/* Return 1 if object is within an avoidance area */
/* LATER, CHANGE CALL SEQUENCE to use a pointer to a cdat structure which
  will contain the cdat->av avoidance pointer, and cdat->tool tool value;
  this will reduce the calling sequence by one element. */
{
avoids	*px;
vect2	c;
//vect2	pc;
#ifdef  DEBUGAVOID
static	int	kk=0;		// debug only
int	k=0;		// debug only
int	special;	// debug only
#endif

#ifdef  DEBUGAVOID
    special = !strcmp(ob->dat.name, "xavd");	// debug only
    if (special) {	// debug only section
	printf (" %s is at %7.2f %7.2f\n",
		ob->dat.name, ob->smpos.x, ob->smpos.y);
    }
#endif

    C_loop (px, av) {
	c = make2vect (-px->x, px->y);
/* Above minus sign:  The avoidance coordinates are in mask cutting
  coordinates, which are mirror imaged with respect to the actual mask
  locations contained in the object structure.  This will put both in
  the same coordinate system for the check.  */

#ifdef  DEBUGAVOID
/* debug section */
	if (special) {	// debug only section
double r;
	    r = dist2 (c, ob->smpos);
	    printf (" %s is %7.2f from %7.2f %7.2f of radius %6.2f\n",
		ob->dat.name, r,
		c.x, c.y, px->r + tool);
	}
/* end debug section */
#endif

	if (Xcircle (ob, c, px->r + tool) < 1) return 1;
#ifdef  DEBUGAVOID
	k++;	// debug only
#endif
    }

#ifdef  DEBUGAVOID
    if (kk < 3) {	// debug only section
	printf (" Searched %d avoidances.\n", k);
	kk++;
    }
#endif

    return  0;
}



/* -- Obsolete grating/grism angle/wavelength programs... */

// static  double  Gr_Angle_2 ( double lambda, double density, double Camangle, int order)
// /* Compute the grating angle, radians */
// Replaced by Gr_Angle in optutils.c

// static  double  Gangle (int Mode, double lambda, grating* g)
// /* Compute the grating angle, radians */
// Obsolete, unused.  Use Gr_Angle now.

// double  Gwavl (int Mode, grating* g, double lambda)
// /* Find the center wavelength of a grism.  If none possible,
//   use the given lambda instead.  */
//  Replaced by "grwavc" in optutils.c

/*  ---  Optical utility (local to mgutils)  --  */

/* foclen is used only within mgutils, so is declared static; we really
  don't want outside users to do this bad method anyway.  */

static	double	foclen (Obs* ob)
/* Find a focal length by some means */
{
focdat  *fd;
double	FocalLength;

// We prefer a focal length from the optics data; failing that, we
// substitute an approximate value.

/* First, take care of any NULL things... */
    FocalLength = 71089.123;	// default value (Magellan)
    if (ob == NULL)  return  FocalLength;

/* Look for the good stuff first */
    if (ob->Tel_scope != NULL) {	// An actual element
	fd = (focdat*) ob->Tel_scope->data;
	if (fd == NULL) {
	    FocalLength += 0.25;	// debug ident.
	    if (ob->dlevel > 0)
		printf ("  No data, focal length is %.3f\n", FocalLength);
	} else {
	    FocalLength = fd->scale;
	}
/* Last resort is to use the default value established */
    } else {
	if (ob->dlevel > 0)
	    printf ("  Using old Focal Length of %.3f\n", FocalLength);
    }
    return  FocalLength;
}


/*  ====  Checking the object list warning flags  ====  */

#ifdef  OBJCHECK
static	int	objc_total (objq* ob, int mask)
/* Return number of active objects with masked error bits set */
{
int	r = 0;
objq	*p;
    C_loop (p, ob) {
	if ((p->flag & OBJECT_ACTIVE) && (p->cflag & mask)) r++;
    }
    return  r;
}

static	int	objc_test (objq* ob, int mask, int px, char* msg)
/* Utility to compute error percent and give message, return status */
{
int	x = 0;
int	b;

    b = (objc_total (ob, mask) * 100) / active_objs (ob);
    if (b >= px) {
	printf ("  * %3d%% have %s.\n", b, msg);
	x++;
    }
    return  x;
}


int	object_check (objq* ob)
/* Return >0 if we would like to halt due to errors. */
{
int	r = 0;
int	n;
int	b;

/* Find total objects, see if we have very few */
    n = active_objs (ob);
    if (n < 7) return r;	// Trivial

/* Find number of possible errors in total; see if we have large number */
    b = (objc_total (ob, -1) * 100) / n;
//    if (b < 1) return r;	// Few errors (debug only)
    if (b < 10) return r;	// Few errors

/* Errors found, report that they exist. */
    printf (
	"\n ** Warning: %3d%% of active objects have unusual parameters!\n",
	b);
    r++;

/* Do individual error reports on each bit in cflag field here */

    r += objc_test (ob, OC_PREQ,  10, "pre-comments (these are unusual)");
    r += objc_test (ob, OC_NSHORT,10, "names less than 2 characters long");
    r += objc_test (ob, OC_NEGUSE, 2, "negative use counts");
    r += objc_test (ob, OC_PRIX,  10, "priorities <-10 or >30");
    r += objc_test (ob, OC_WIDTH,  5, "slit width under 0.3 arcseconds");
    r += objc_test (ob, OC_WIDEX, 10, "slit width over 15 arcseconds");
    r += objc_test (ob, OC_SHAPE,  3, "shape codes out of range");
    r += objc_test (ob, OC_ALEN | OC_BLEN, 10, "slit length out of range");
    r += objc_test (ob, OC_LENGTH, 5, "length less than width");
    r += objc_test (ob, OC_ANGLE, 10, "tilt angles over 20 degrees");

// This is older way replaced by objc_test subroutine...
    b = (objc_total (ob, OC_UNCUT) * 100) / n;
    if (b >= 5) {
	printf ("  * %3d%% have non-default lengths in N&S mode\n", b);
	r++;
    }
    printf ("\n");
    return  r;
}

#endif  // on OBJCHECK

/*  =====  Filling in some structures  =====  */

	/* - - - - - - - - - - *\
	|			|
	|   Get object struc	|
	|			|
	\* - - - - - - - - - - */

objq*  get_object (char* data, int type, namlist* preq, Obs* ob)
/* Obtain an object queue element, and fill with data from line,
  put in preq pointer, leaving postq pointer null for now */
{
objq*	b;
// double	scale;
char	*c;
char	token[64];

/* Get an object data structure and clean it */
    b = xalloc (objq);
    memset (b, 0, sizeof(objq));	/* CLEAR TO ZERO */
    b->next = b->last = b;

/* Fill in some default and required values */
//    b->dat.priority = 12.0;		// default value
    if (type == OBJ_REFERENCE) b->dat.priority = ob->Pmusthave - 0.5;
    else  b->dat.priority = ob->Pdecide + 5.0;

//    b->dat.wlimit.red = b->dat.wlimit.blue = 0.0;
    b->dat.precom  = preq;
    b->dat.postcom = NULL;

    if (type == OBJ_REFERENCE)	b->slit = ob->reference_slit;
    else			b->slit = ob->default_slit;
#ifdef  OBJCHECK
    b->cflag = 0;	// Redundant
/* Check the default slit parameters here */
    if (b->slit.width <  0.3) b->cflag |= OC_WIDTH;
    if (b->slit.width > 15.0) b->cflag |= OC_WIDEX;
    if (b->slit.shape < 0 || b->slit.shape > SPACED_CROSS)
		b->cflag |= OC_SHAPE;
    if (b->slit.alen < 0.0) b->cflag |= OC_ALEN;
    if (b->slit.blen < 0.0) b->cflag |= OC_BLEN;
    if ((b->slit.alen + b->slit.blen) < (0.99 * b->slit.width))
		b->cflag |= OC_LENGTH;
    if (fabs(b->slit.angle) > (20.0 * Degree)) b->cflag |= OC_ANGLE;
    if (preq != NULL) b->cflag |= OC_PREQ;
#endif  // on OBJCHECK
    b->val = NULL;
    b->type = type;
    b->flag = 0;	// Redundant
    b->use = 0;

/* Read what can be found in the object record... (first 3 required) */
//    c = parse (data, token, " ");
/* Many instances of " " are replaced by WHITESPACE here 05/3-15 */
    c = parse (data, token, WHITESPACE);
    b->dat.name = dynamstr (token);
    c = parse (c, token, WHITESPACE);
    b->dat.ra  = posvalue (token) * Hour;
    c = parse (c, token, WHITESPACE);
    b->dat.dec = posvalue (token) * Degree;
/* NOTE scaling of R.A. to other than hours is done in the calling
  program if an appropriate "&RADEGREE" record was found */
#ifdef  OBJCHECK
    if (strlen(b->dat.name) < 2) b->cflag |= OC_NSHORT;
#endif

/* The rest of the input record is actually optional */
    c = parse (c, token, WHITESPACE);  if (*token == AN) return b;
    b->dat.priority = atof (token);
#ifdef  OBJCHECK
    if (b->dat.priority < -10.0) b->cflag |= OC_PRIX;
    if (b->dat.priority >  30.0) b->cflag |= OC_PRIX;
#endif
    c = parse (c, token, WHITESPACE);  if (*token == AN) return b;
/* Add the use count decoding.  If we have an integer, it is a use
  count; if floating point it is the slit width and use count is 0.
  Currently, decide that on basis of existing decimal point in field.  */
    if (strchr(token, '.') == NULL) {	// integer, find use
	b->use = atoi (token);
#ifdef  OBJCHECK
	if (b->use < 0) b->cflag |= OC_NEGUSE;
#endif
	c = parse (c, token, WHITESPACE);  if (*token == AN) return b;
    }
    b->slit.width = atof (token);
#ifdef  OBJCHECK
    b->cflag &= ~(OC_WIDTH | OC_WIDEX);
    if (b->slit.width <  0.3) b->cflag |= OC_WIDTH;
    if (b->slit.width > 15.0) b->cflag |= OC_WIDEX;
#endif
    c = parse (c, token, WHITESPACE);  if (*token == AN) return b;
    b->slit.shape = atoi (token);
#ifdef  OBJCHECK
    b->cflag &= ~OC_SHAPE;
    if (b->slit.shape < 0 || b->slit.shape > SPACED_CROSS)
		b->cflag |= OC_SHAPE;
#endif
    c = parse (c, token, WHITESPACE);  if (*token == AN) return b;
    b->slit.alen  = atof (token);
#ifdef  OBJCHECK
    b->cflag &= ~(OC_ALEN | OC_LENGTH);
    if (b->slit.alen < 0.0) b->cflag |= OC_ALEN;
    if ((b->slit.alen + b->slit.blen) < (0.99 * b->slit.width))
		b->cflag |= OC_LENGTH;
/* ADD CHECK OF UNCUT LENGTH, FLAG ANY ALEN DIFFERENT FROM DEFAULT,
  AND ANY UNCUT GREATER THAN DEFAULT */
    if ( (b->slit.alen != ob->default_slit.alen)
	&& (ob->dca != 0.0) ) b->cflag |= OC_UNCUT;

#endif
    c = parse (c, token, WHITESPACE);  if (*token == AN) return b;
    b->slit.blen  = atof (token);
#ifdef  OBJCHECK
    b->cflag &= ~(OC_BLEN | OC_LENGTH);
    if (b->slit.blen < 0.0) b->cflag |= OC_BLEN;
    if ((b->slit.alen + b->slit.blen) < (0.99 * b->slit.width))
		b->cflag |= OC_LENGTH;
    if ( (b->slit.blen != ob->default_slit.blen)
	&& (ob->dcb != 0.0) ) b->cflag |= OC_UNCUT;
// NOTE THAT WE DON'T CHECK THE SIDE WHICH IS NOT UNCUT,
// BUT THAT MAY NOT BE THE RIGHT THING TO DO...
#endif
    c = parse (c, token, WHITESPACE);  if (*token == AN) return b;
    b->slit.angle = atof (token) * Degree;
#ifdef  OBJCHECK
    b->cflag &= ~OC_ANGLE;
    if (fabs(b->slit.angle) > (20.0 * Degree)) b->cflag |= OC_ANGLE;
#endif

/* Add any more optional things here, decode as above */
    return  b;

/* Slit parameters were to be scaled in fillobject program, when the
  scale factor (focal length) is known. */

}

/* The get_grating thing is totally obsolete */
//grating*	get_grating (char* name)
///* Fake set-up to simulate a grating. */

// in dropping get_grating, we are missing the default order of
// some gratings.  That should only affect intgui in practice.


	/* - - - - - - - - - - *\
	|			|
	|   Get Detector struc	|
	|			|
	\* - - - - - - - - - - */


//detector*	get_detector (char* name, element* gs, int dl)
static	detector*	get_detector (element* inst, element* gs, int dl)
/* Set up a detector structure from an instrument, for real. */
{
detector	*d;	// New structure
element		*e = NULL;
datlist		*p;
double	dsize;
double	dxtent;
double	psiz;
//double	vrad;
double	pixls;
intval	de;

/* Get a relatively clean detector pointer */
    d = xalloc (detector);

/* Attempt to obtain the optical data for a detector part of the
  supplied instrument.  If not available, use default IMACS stuff. */

/* Find an element pointer to a detector element (or NULL) */
    if (inst != NULL) {
	if (type_of(inst->kw) == Detector) e = inst;
	else if (type_of(inst->kw) == Instrument) {
	    C_loop (p, inst->head) {	// loop instrument members
		if (p->eg == NULL) p->eg = find_element (gs, p->data);
		if (p->eg == NULL) continue;
		if (type_of(p->eg->kw) == Detector) {
		    e = p->eg;
		    break;
		}
	    }

	}
    }

// debug section - show we found detector
//    if (e != NULL && dodebug > 0) {
//	printf (" --Detector %s found in instrument %s.\n",
//		e->name, inst->name);
//    }
// end debug only stuff

/* If we did not find the detector, fill in default values and go */
    if (e == NULL) {		// Fill with defaults

	dsize  = 125.0;
	dxtent = dsize / 2.0;
	pixls  = 8192.0;
	psiz   = dsize / pixls;
	de = interval (-dxtent, dxtent);

	d->bb     = boundbox (de, de);
	d->pcount = make2vect (pixls, pixls);
	d->psize  = make2vect (psiz, psiz);
	d->vrad   = 100.0;
	d->bbs = d->bbn = d->bbd = d->bb;

	if (inst != NULL && dl > 0) {
	    printf (" ** Default detector information used for %s.\n",
		inst->name );
	}

    }  else {			// Fill from element
	*d = *((detector*) e->data);	// get the data

/* Compute the pixel size, if needed... */
	dsize = intlen (d->bb.x);
	pixls = d->pcount.x;
	if (pixls > 0.0) d->psize.x = dsize / pixls;

	dsize = intlen (d->bb.y);
	pixls = d->pcount.y;
	if (pixls > 0.0) d->psize.y = dsize / pixls;

/* If no vigneting radius, set it high for compatability */
	if (d->vrad == 0.0) d->vrad = hypot (intlen(d->bb.x), intlen(d->bb.y));

	if (dl > 0) printf (" --Detector %s found in instrument %s.\n",
		e->name, inst->name);

    }

    return  d;
}


static	void	get_type (char* t, objq* ob, Obs* obs)  {
/* Put type into the title location */
    if        (ob->type == OBJ_REFERENCE) {
//	strcpy (t, "Reference");
	strcpy (t, "Alignment");
    } else if (ob->type == OBJ_GAP)       {
	strcpy (t, "Gap");
    } else if (ob->type == OBJ_OBJECT &&
	       ob->dat.priority < obs->Pmusthave)  {
	strcpy (t, "Must-Have");
    } else if (ob->type == OBJ_OBJECT)    {
	strcpy (t, "Object");
    } else {
	strcpy (t, "Unknown");
    }
}


/*  ===  Reading files  ===  */

	/* - - - - - - - - - - *\
	|			|
	|  Default Obs struct.	|
	|			|
	\* - - - - - - - - - - */


/* Subroutine to create and return a default obs structure */
Obs*	defobs (element* gs)
{
Obs	*a;

/* Get storage for the structure(s) */
    a = xalloc (Obs);

/* Establish default values in new structure */
    strncpy (a->oname, "Annonymous", 16);
    strncpy (a->fname, "Unnamed", 8);
    a->title = NULL;
/*    a->object = NULL; */
//    a->instrument = NULL;
    a->ra = 24.0 * Hour;	/* Indicate no data */
    a->dec = 90.0 * Degree;	/* Indicate no data */
    a->equinox = 2000.0;
    a->epoch = Jyear (jdnow() + 7.0);
    a->angle = 0.0;  /* Seems resonable */
//    a->temp = 15.0;  /* -NOTE- improve to seasonal default */
//    a->temp = 5.0 + 5.0 * ( (int) (2.0 * acos(TwoPi *
//			(a->epoch - (int)a->epoch - 0.15) )) );

//    a->temp = 5.0 + 5.0 * ( floor (0.5 + 2.0 * acos(TwoPi *
//			(a->epoch - floor(a->epoch) - 0.15) )) );
//    a->temp = 5.0;	// See what's wrong with linux?
    a->temp = DefTemp;	// 20 + 1/256 is default no-temp value
/* Produces an approximate seasonal temperature default */
#ifdef  DIFREF
//    a->ha = 0.0;	// Default is meridian; could also use -25.0 ?
    a->ha = -25.0;	// Default if not specified is "off" when reading!
#else
    a->ha = -25.0;	// Set condition for no d.r. when not in use
#endif
/* ?? mxpos is not initialized... */

/* Initialize guiders to all null */
//#ifdef GSDATA
    a->gp1.ra = a->gp1.dec = a->gp1.ep = 0.0;
    a->gp1.cn = 0;
    a->gp2.ra = a->gp2.dec = a->gp2.ep = 0.0;
    a->gp2.cn = 0;
//#endif

    a->pb.red = 7900.0625;
    a->pb.blue = 4000.0625;	/* Default values flagged with 1/16 */
/* Keep less than factor of 2 for 2nd order self conflict */
    a->pb2 = a->pb;
    a->pb2act = 0;
    a->cw = 5800.0625;  /* Just a guess */
//    a->sw = pbmean (a->pb);
    a->ob = NULL;

/* Set reuse and select counts */
    a->reuseob = 999;
    a->reuserf = 999;
    a->refsel  = 1;
    a->reflimit = 99;
    a->order = 0;
/* Order is 0 until it is read.  If not read, it is defaulted to 1
  in the normalize_obs routine later. */

/* Paremeterized decision priorities */
    a->Pmusthave = -2.0;
    a->Pdecide = 0.0;

    a->D_angle = 0.0;

/* Set the global pointers to this obs structure... */
    GRangle = &(a->D_angle);
    GRorder = &(a->order);
/* This is done here, as this program is intended to be called only in
  read_obsfile and read_smdf to create the master obs structure; if
  other obs structures are designed, these initializations above should
  be moved to the read_obsfile and read_smdf programs which call defobs.  */

// //
/* Set the current element list into the structure... */
    a->elist = gs;
    a->Tel_scope   = NULL;
    a->Disperser   = NULL;
    a->Instr_mt    = NULL;
    a->Direct_inst = NULL;

//    a->Detect = get_detector ("default");
    a->Detect = get_detector (NULL, gs, 0);
/* -NOTE- Add some IMACS values as defaults here */
    a->ddir.blue = 0.5;
    a->ddir.red = -0.5;
/* -NOTE- above is a guess; compute it below */
    a->Dstat = 0;
//    a->filter = NULL;
    a->default_slit.shape = RECTANGLE;
    a->default_slit.width = 1.0;
    a->default_slit.alen  = 6.0;
    a->default_slit.blen  = a->default_slit.alen;
    a->default_slit.angle = 0.0;
//    a->default_slit.curve = 0.0;
/* --NOTE-- Find actual units for the above... */
    a->reference_slit.shape = SQUARE;
    a->reference_slit.width = D_refsiz;		// 2.0;
    a->reference_slit.alen  = D_refsiz/2.0;	// 2.0;
    a->reference_slit.blen  = D_refsiz/2.0;	// 2.0;
// or, use alen for blen...
    a->reference_slit.angle = 0.0;
//    a->Det_mode = Grating | LongCam;
    a->ex_order = 0;	// Do extra order conflicts
    a->dlevel = 0;		// Debug off.
    a->minsep = D_minsep;	// Min. (-overlap), pixels
    a->dca = a->dcb = 0.0;	// Do not cut distances
    a->objfiles = NULL;
/* --NOTE-- the objfiles will be required upon reading; if none
  are given, a prompt will be done  (i.e. DO the prompt!) */
    a->comments = NULL;
    a->IFUmode = 0;	// Off as standard.
    a->IFUra   = a->ra;
    a->IFUdec  = a->dec;
    a->Flag  = 0;	// All flags start as off...
    a->slext = 0;	// no extension as default
    a->pcom  = NULL;	// no pass-through comments until read
    a->cref  = 0.0;	// Start with nothing...
    a->zenith = make2vect (0.0, 0.0);
    a->date = 0;	// default 0 date.
    return  a;
}

void	get_mpos (Obs* a)
/* Use ra, dec, angle to form the mpos, invp matrices in a */
/* Use ra, dec, angle to form the mpos matrix in a */
/* Use ra, dec, angle to form mxpos matrix in a */
{
vect3	v;
vect3	t[3];
    v = sph2vec (a->ra, a->dec);
    smatrix (t, v);
    mrotate (t, a->angle, 1);
    a->mxpos.x = t[0];
    a->mxpos.y = t[1];
    a->mxpos.z = t[2];
// old way:
//    smatrix (a->mpos, v);
//    mrotate (a->mpos, a->angle, 1);
//    mcpy3 (a->invp, a->mpos);
//    trpose3 (a->invp);
}


	/* - - - - - - - - - - *\
	|			|
	|	Normalize	|
	|			|
	\* - - - - - - - - - - */

static	element*	valid_epointer (element* e)
/* Return element pointer e if valid, NULL if not */
{
/* Drop any NULL pointer */
    if (e == NULL) return NULL;
/* Any self-referent is a single name holder */
    if (e->next == e) return NULL;
    if (e->last == e) return NULL;
/* If not linked to a queue, it is not valid */
    if (e->next == NULL) return NULL;
    if (e->last == NULL) return NULL;
/* Whatever is left must be a valid pointer */
    return  e;
}

void	normalize_obs (Obs* a)
/* When optic element names and data are present, set up the older
optic-element data in the Obs structure to be compatable as much as
possible with the new things.  All this will become obsolete when
the old elements are removed and all optic computations are done by
the generalized transform. */
/* This is called by read_obsfile after reading all data in the
  file; and by interface gui after modifications to elements and
  before writing the Obs data out to the file.  */
{
element	*gs = NULL;
int	imx = I_imacs;	// default 1;
char	*iname;

// NOTE -- we really should generate an instrument index here, and
// store it in the obs structure in addition to the names...
// it would be useful in the intgui, and elsewhere...

    iname = (a->Instr_mt == NULL) ? NULL : a->Instr_mt->name;
//    imx = ExtIndx (iname);
    imx = InstIndx (iname);

//    if (iname != NULL) {
//	if (!strcmp(iname, "IMACS_sc")) imx = 2;
//	if (!strcmp(iname, "IMACS_lc")) imx = 1;
//	if (!strncmp(iname, "LDSS", 4)) imx = 0;	// or I_imacs ?
//    } else  imx = 1;	// default

/* Find an element queue item, if any */
    if (gs == NULL) gs = a->elist;
// if (gs == NULL) printf (" ** Normalize_obs, no elist pointer!\n");  //debug
    if (gs == NULL) gs = valid_epointer (a->Tel_scope);
    if (gs == NULL) gs = valid_epointer (a->Instr_mt);
    if (gs == NULL) gs = valid_epointer (a->Disperser);
/* The pointers above might have been created with a fake element
  structure to hold the element name only when reading SMDF file.  */


// Need to get detector dimensions until we do a better job of
// this by examining the instrument elements...
    kill_ptr ((void*)&(a->Detect));
//    a->Detect = get_detector(iname, gs);
    a->Detect = get_detector(a->Instr_mt, gs, a->dlevel);

/* Apply the detector defaults based on Dstat word */
// WANT TO RESTRICT ORIENTATION FEATURE TO IMACS INSTRUMENTS ONLY
    if ((a->Dstat & DS_orient) && (InstIndx(a->Instr_mt->name) == I_imacs)) {
	a->Detect->bb = (a->Dstat & DS_ns) ? a->Detect->bbn : a->Detect->bbs;
    }

/* Use default telescope if necessary */
    if (a->Tel_scope == NULL && gs != NULL) {
//	if (imx == IX_ldss) a->Tel_scope = find_element (gs, "Magellan2");
	if (imx == I_lds3) a->Tel_scope = find_element (gs, "Magellan2");
	else		   a->Tel_scope = find_element (gs, Tname);
    }
// Should look to noadc or instrument to default the telescope!

/* The GRangle and GRorder pointers were set when the Obs structure
  was initialized in defobs.  They point to D_angle and order items.  */

    if (a->order == 0) a->order = 1;	// actual default.

    if (a->Disperser != NULL) {
element	*e;
//stgrism	*sg;
double	den;
	e = a->Disperser;
	if (e->kw == GRAting) {		// If grating:
// Density is a element property.
// Set angle using Gangle subroutine
	    den = *((double*)e->data);
	    a->D_angle = Gr_Angle (a->cw, den, CAM_ANGLE, a->order);
	    Cur_grating = e;
	} else if (e->kw == GRIsm || e->kw == EGRIsm) {	// If grism:
// Angle is a grism property in element, also density
	    Cur_grism = e;
// Set grism angle and order globals correctly
//	    if (e->kw == DGRIsm) {
//		a->D_angle = ((sdgrism*)e->data)->grism.angle * Degree;
	    if (e->kw == EGRIsm) {
// change definition to max (gangle, eangle) - inangle here:
// or simply the max (gangle, eangle) would work too...
// This is entirely obsolete, anyway...
// Defined as grating angle -- use gangle only
		a->D_angle = ((egrism*)e->data)->gangle  * Degree;
	    } else {		// otherwise assume a regular grism...
		a->D_angle = ((stgrism*)e->data)->angle * Degree;
	    }
	    a->cw = grwavc (e, a->temp);
	} else if (e->kw == ECHelle) {
	    a->D_angle = 0.0;
	    a->order = 0;
	    Cur_grating = e;
	}	// Grating/Grism
    }

/* Find sw.  If non-default passband values, it is the mean.
  Otherwise, it is the cw value. */
    if (
	( (a->pb.blue - floor(a->pb.blue) ) != 0.0625 ) &&
	( (a->pb.red  - floor(a->pb.red ) ) != 0.0625 )  ) {
	a->sw = pbmean (a->pb);
// Default passband values have 0.0625 angstroms as the fractional
// part to identify them.  If both have been changed, use mean
    } else {
// Otherwise, with a default value at red or blue end, use cw
	a->sw = a->cw;
    }


/* If MOE is being used, turn off incompatable features */
    if (a->Flag & OFmoe) {
	a->Flag  &= ~(OFsxt | OFxrd);
	a->ex_order = 0;
	a->slext    = 0;
	a->Dstat &= ~DS_gaps;
    }

/* Set the D_maskrad value based on current instrument */
//    D_maskrad = (imx > 0) ? D_maxradi : D_maxradl;
//    D_maskrad = (imx == IX_ldss) ? D_maxradl : D_maxradi;
    D_maskrad = (imx == I_lds3 || imx == I_lds2) ? D_maxradl : D_maxradi;

/* Set hard limits on temperature field */
    if (a->temp < -40.0) a->temp = -40.0;
    if (a->temp >  40.0) a->temp =  40.0;
/* Note that since write_obs calls this normalize routine before
  writing, this also limits the temp values which can be output
  by the interface gui program. */

/* Not making the section below depend on DIFREF, since it should and
  can be computed in all cases. */
/* Get the differential refraction things computed once here... */
    if (a->ha < -24.0) {	// Turn dif_ref stuff "off"
	a->cref = 0.0;
	/* a->zenith does not matter now */
    } else {		// Real dif_ref computations...
/* Local variables within these { } here */
double	az, zd, pa;
double	da, dz, dp;
double	p, s, c;
double	w;

/* Find the refraction value from known constants. */
	w = pbmean (a->pb);
	if (w < 100.0) w = 5500.0;	// Safety statememt.
	a->cref = ndxref (w, Pressure, a->temp);

/* Find the zenith position vector,
  from declination, latitude, hour angle */
	dntri ( &az, &zd, &pa,  &da, &dz, &dp,
		a->ha * Hour,	// ha is kept in hours, need radians
		a->dec,		// dec is kept in radians
		Latitude * Degree  );

/* Error to have too large a zenith distance... */
	if ((zd/Degree) > 80.0) {	// Bad zenith distance
	    a->cref = 0.0;
	    a->zenith = make2vect (0.0, 0.0);
	    printf (" ! Zenith distance of %.1f exceeds 80 Deg.\n",
		zd / Degree );
	    printf (" ! No differential refraction computed!\n");
	} else {		// Good zenith distance
	    p = pa - a->angle;		// Hope this is right!
/* The pa is an angle on the sky from North through East.  The a->angle
  is the orientation of the local y axis on the sky.  Our result is
  the position angle in local x,y system.  Think -- as angle increases,
  our y axis displaces reducing the Zenith vector pa.  */
	    sincos (p, &s, &c);
	    a->zenith = mul2vect ( make2vect (-s, c), zd);
/* The -s comes from the fact that the East direction is in the -x
  direction in a coordinate system where y is North, and the position
  angle increases from North to East, or ccw from the y axis. */
	    if (a->dlevel > 0)
		printf ("<DR> %11.4g pa=%.2f p=%.2f z=(%.2g, %.2g)\n",
		    a->cref, pa/Degree, p/Degree,
		    a->zenith.x, a->zenith.y  );
	} // End zenith distance resonableness conditional
    }	// End dif_ref conditional
    if (a->slext && a->dlevel > 0) printf (" --Slit extension requested.\n");
/* End of normalize_obs */
}


	/* - - - - - - - - - - - - - - *\
	|				|
	|   Read  Observation  File	|
	|				|
	\* - - - - - - - - - - - - - - */

Obs*	read_obsfile (FILE* input, element* gs)
/* Read observation file (see outline and data definitions) */
{
Obs	*a;

int	length;
int	l;
char	buffer[OBSLIST_ENTRY+1];
int	k, j;	// scratch
int	dbw = 0;
int	gpw = 0;
// char	cval;
char	bf;
char	*key;
char	*rest;
char	*cc;	// debug only
char	*st;
char	*cp;
char	token[64];
char	buf[64];
double	val;
double	bv;
// double	scale;
element	*ep;

namlist*	comq = NULL;	/* Remove others above... */
namlist*	pcomq = NULL;
namlist*	bn;

//if (DEBUG_UTL) {
// printf ("Entering read_obsfile.\n"); fflush (stdout); //debug
//}

#ifdef  TS_dotime
	TX_con (TX_size);	// Debug of getedge will be started...
#endif

/* Get a default obs structure to fill with data... */
    a = defobs (gs);

/* Loop reading the file */
/* How done:
  Read data records, identify first blank separated field
  If comment, see if .-= flag sends it to comment queue
  Use switch statement fo field record types
..*/

//if (DEBUG_UTL) {
// printf ("In read_obsfile, start reading.\n"); fflush (stdout); //debug
//}

/* Use the example read, and fix it up here... */

    while ( fgets (buffer, OBSLIST_ENTRY, input) != NULL) {
/* Clean out terminating newline from buffer */
	trimend (buffer, NEWLINES);
	length = strlen (buffer);
	if (length < 1) continue;	/* Drop null lines */
/* Test for and deal with comment records */
	if (*buffer == '!' || *buffer == '#') {	/* comment */
	    bf = buffer[1];
	    if (bf == '.' || bf == '-' || bf == '=') {  /* keeper */
		comq = pushtxt (buffer+2, comq);
/* Pass through comments starting #! or !# -- output will be !#  --
  that is done in write_comments which always writes ! comments. */
	    } else if (*buffer == '#' && bf == '!') {	/* pass-thru */
		pcomq = pushtxt (buffer+2, pcomq);
	    } else if (*buffer == '!' && bf == '#') {	/* pass-thru */
		pcomq = pushtxt (buffer+2, pcomq);
	    }
	    continue;		/* Otherwise ignore the comment */
	}

/* Deal with individual records by first character */
/* Obtain the keyword and separate from rest of buffer */
	key = strtok (buffer, WHITESPACE);	/* keyword */
	rest = strtok (NULL, "");
	if (key == NULL) continue;	/* ? may want to debug print? */
	if (rest == NULL) {
	    l = strlen (key);
	    rest = key + l;
	}

/* Trim whitespace from start and end of rest */
	trims (rest, WHITESPACE);

/* Debug print of commend... */
//	if (DEBUG_UTL) {
//	    printf ("Ob.key: \"%s\", rest: \"%s\".\n", key, rest);
//	    fflush (stdout);
//	}

/* We uppercase the "key" string at this time,
  so the compares will be non-case-sensitive */
	for (st=key; *st != '\0'; st++) *st = toupper(*st);
/* OR -- we could use "strcasecmp" instead of "strcmp"  */

/* Decision tree on keyword... */
	if        (strcmp (key, "OBSERVER") == 0) {
/* OBSERVER  put 2-16 characters into observer name  */
	    while ((st=strchr(rest, ' '))) *st = '_';
// Replace any embeded blanks with underscores
	    strncpy (a->oname, rest, 16);
/* -NOTE- need to enforce 2 char. minimum! */
	    while (strlen (a->oname) < 2) strcat (a->oname, "-");
	} else if (strcmp (key, "FILENAME") == 0) {
/* FILENAME  put 1-8 characters into file name  */
	    strncpy (a->fname, rest, 8);
/* -NOTE- need to enforce 1 char. minimum! */
	    if (strlen(a->fname) < 1) strcpy (a->fname, ".");
	} else if (strcmp (key, "TITLE") == 0) {
/* TITLE  put rest into dynamic string into title  */
	    free (a->title);
	    a->title = dynamstr (rest);
	} else if (strcmp (key, "HEADER") == 0) {
/* HEADER  Add to the commentary text with prefix of "  "  */
	    l = strlen (rest);
	    st = (char *) malloc (l + 3);
	    strcpy (st, "  ");
	    strcat (st, rest);
	    bn = get_nam (st);
	    a->comments = push_name (bn, a->comments);
	} else if (strcmp (key, "DATEOBS") == 0) {
/* DATEOBS  sets epoch of observation  */
/* -NOTE- need time program here; decode dateobs into year, month,
  day -- then to decimal years and store that in a->epoch */
/* Format of rest should be %d-%d-%d read as year-month-day; then year
  needs 2000 added if < 1800. */

	} else if (strcmp (key, "CENTER") == 0) {
/* CENTER  read ra and dec from line  */
	    a->ra  = posvalue (strtok (rest, WHITESPACE)) * Hour;
	    a->dec = posvalue (strtok (NULL, WHITESPACE)) * Degree;
	    get_mpos (a);	/* Get transform matrix */
	} else if (strcmp (key, "EQUINOX") == 0) {
/* EQUINOX  set equinox from line  */
	    sscanf (rest, " %lf", &(a->equinox));
	} else if (strcmp (key, "POSITION") == 0) {
/* POSITION  set angle, in degrees  */
	    val = 0.0;
	    sscanf (rest, " %lf", &val);
	    a->angle = val * Degree;
	    get_mpos (a);	/* Get transform matrix */
	} else if (strcmp (key, "GS1") == 0) {
/* GS1 is guide star 1 */
//#ifdef  GSDATA
		a->gp1.ra  = posvalue   (strtok (rest, WHITESPACE));
		a->gp1.dec = posvalue   (strtok (NULL, WHITESPACE));
		a->gp1.ep  = floatvalue (strtok (NULL, WHITESPACE));
//		a->gp1.cn  = intvalue   (strtok (NULL, WHITESPACE));
//#endif
	} else if (strcmp (key, "GS2") == 0) {
/* GS2 is guide star 2 */
//#ifdef  GSDATA
		a->gp2.ra  = posvalue   (strtok (rest, WHITESPACE));
		a->gp2.dec = posvalue   (strtok (NULL, WHITESPACE));
		a->gp2.ep  = floatvalue (strtok (NULL, WHITESPACE));
//		a->gp2.cn  = intvalue   (strtok (NULL, WHITESPACE));
//#endif
	} else if (strcmp (key, "SLITSIZE") == 0) {
/* SLITSIZE  set default slit parameters */
/* Allow the actual defaults to be replaced only with
  real data from the input line... */
//	    cp = parse (rest, token, " ");
/* Several instances of " " replaced by WHITESPACE 05/3-15 */
	    cp = parse (rest, token, WHITESPACE);
	    if (*token != AN) a->default_slit.width  = atof (token);
// Does not include a shape code at present time...
	    cp = parse (cp, token, WHITESPACE);
	    if (*token != AN) a->default_slit.alen   = atof (token);
	    cp = parse (cp, token, WHITESPACE);
	    if (*token != AN) a->default_slit.blen   = atof (token);
	    cp = parse (cp, token, WHITESPACE);
	    if (*token != AN) a->default_slit.angle  = atof (token) * Degree;
/* Scaling of slit size to millimeters is done when actual object
  data is read as input; defaults are kept in arc seconds. */
	} else if (!strcmp (key, "REFHOLE")) {
/* REFHOLE  set reference hole default parameters */

/* Debug - show the reference slit widths...  */
//    printf (" $$2.5 Ref. Slit w,a,b: %.4f, %.4f, %.4f\n",
//	a->reference_slit.width,
//	a->reference_slit.alen, a->reference_slit.blen);
//    printf (" For REFHOLE, rest = <%s>\n", rest);

/* Do reference as above... */
	    cp = parse (rest, token, WHITESPACE);
	    if (*token != AN) a->reference_slit.width  = atof (token);
	    cp = parse (cp, token, WHITESPACE);
	    if (*token != AN) a->reference_slit.shape  = atoi (token);
	    cp = parse (cp, token, WHITESPACE);
	    if (*token != AN) a->reference_slit.alen   = atof (token);
	    cp = parse (cp, token, WHITESPACE);
	    if (*token != AN) a->reference_slit.blen   = atof (token);
	    cp = parse (cp, token, WHITESPACE);
	    if (*token != AN) a->reference_slit.angle  = atof (token) * Degree;
/* Also, reference slit sizes are scaled only when objects are read. */
/* If shape is square or circle, set alen and blen to width/2; this will
  set the proper spectrum size... */
	    if (a->reference_slit.shape == SQUARE ||
		a->reference_slit.shape == CIRCLE )
		a->reference_slit.alen = a->reference_slit.blen =
		a->reference_slit.width / 2.0;
/* NOTE: Above should be done for all; AND, don't allow slits shaped
  as square or circle to escape when read in later, too! */

/* Debug - show the reference slit widths...  */
//    printf (" $$3 Ref. Slit w,a,b: %.4f, %.4f, %.4f\n",
//	a->reference_slit.width,
//	a->reference_slit.alen, a->reference_slit.blen);

	} else if (strcmp (key, "SLEXTEND") == 0) {
/* SLEXTEND  slit extension algorithm is allowed. */
	    sscanf (rest, " %d", &(a->slext));
//	    if (a->slext) printf ("  - - Slit extension requested.\n");
	} else if (strcmp (key, "OVERLAP") == 0) {
/* OVERLAP  read overlap value  */
//	    sscanf (rest, " %lf", &(a->minsep));
// minsep operates in inverse fashion from overlap as a parameter.
		a->minsep = -atof (rest);
	} else if (strcmp (key, "UNCUTLEFT") == 0) {
/* UNCUTLEFT read amount in arcseconds to not cut on left slit end */
	    sscanf (rest, " %lf", &(a->dca));
	} else if (strcmp (key, "UNCUTRIGHT") == 0) {
/* UNCUTRIGHT read amount in arcseconds to not cut on right slit end */
	    sscanf (rest, " %lf", &(a->dcb));
	} else if (strcmp (key, "MUSTHAVE") == 0) {
/* MUSTHAVE read priority of must have object */
	    sscanf (rest, " %lf", &(a->Pmusthave));
	} else if (strcmp (key, "PDECIDE") == 0) {
/* PDECIDE read priority where decision algorithm changes from number
  of conflicts (below it) to priority preference (above it) */
	    sscanf (rest, " %lf", &(a->Pdecide));
	} else if (strcmp (key, "INSTNAME") == 0) {	// OBSOLETE
//	    printf ("-- Instrument name currently unused.\n");
	    strcpy (buf, strtok (rest, WHITESPACE));
/* NOTE - ADD SECTION - see if instrument name is an element in the
  current element queue (gs), with different error message if not. */
	    ep = find_element (gs, buf);
	    if (ep == NULL) {
//		printf (" **Instrument name \"%s\" not recognized.\n", buf);
		printf (" **Obsolete INSTNAME name \"%s\" not recognized.\n",
			buf);
	    } else {
		a->Instr_mt = ep;
		printf (" *INSTNAME is obsolete.  Please change to INSTRUMENT");
		printf ("\n  -- Used \"%s\" as Instrument.\n", buf);
	    }

	} else if (strcmp (key, "INSTRUMENT") == 0) {
/* INSTRUMENT -- Set optical element based on name */
	    a->Instr_mt = find_element (gs, rest);
// NOTE that may fail if comment follows name??
/* ALSO NEED == read Telescope, grating/grism/disperser names here, and
  validate them with the element queue!! */
	} else if (strcmp (key, "TELESCOPE") == 0) {
/* TELESCOPE  read telsecope name, find element  */
	    cc = strtok (rest, WHITESPACE);				//debug
//	    printf (" -- Telescope name \"%s\".\n", cc);	//debug
	    strcpy (buf, cc);
// Fix up the telescope name, if Mag2NoADC, use Magellan2
	    if (!strcmp(cc, "Mag2NoADC")) strcpy (buf, "Magellan2");
	    a->Tel_scope = find_element (gs, buf);
	    if (a->Tel_scope == NULL)
		printf (" **Telescope %s not found.\n", buf);

	} else if (strcmp (key, "DISPERSER") == 0) {
/* DISPERSER -- Set optical element based on name */
	    a->Disperser = find_element (gs, rest);
/* Recognize the special disperser for echelle support */
	    if (strcasecmp (a->Disperser->name, "MOE") == 0) {
		printf (" * MOE recognized. *\n");	// debug only
		a->Flag |= (OFmoe | OFpb2);
		a->Flag &= ~(OFsxt | OFxrd);
// may need to set Det_mode properly here too...  but, that
// should have been done by long camera instrument line.
	    }
	} else if (strcmp (key, "GRATING") == 0) {
/* GRATING  read grating name (unused here)  */
////	    printf ("-- Grating name currently unused.\n");
	    printf (" **Obsolete Grating keyword for %s found - dropped.\n", rest);
	} else if (strcmp (key, "FILTER") == 0) {
/* FILTER  read filter name (unused here)  */
	    printf ("-- Filter name currently unused.\n");
	} else if (strcmp (key, "WAVELENGTH") == 0) {
/* WAVELENGTH  read center wavelength  */
	    sscanf (rest, " %lf", &(a->cw));
	} else if (strcmp (key, "WLIMIT") == 0) {
/* WLIMIT  read wavelength limits, set as red and blue  */
	    l = sscanf (rest, " %lf %lf", &val, &bv);
	    if (l == 2) {
		a->pb.blue = min (val, bv);
		a->pb.red  = max (val, bv);
/* If no secondary limit has been set, fill in from here */
		if (a->pb2act == 0) a->pb2 = a->pb;
//		a->sw = pbmean (a->pb);  // done in normalize
	    }
	} else if (strcmp (key, "DLIMIT") == 0) {
/* DLIMIT  read wavelength limits, set as red and blue  */
/* If found, this is a secondary limit set, and we flag that */
	    l = sscanf (rest, " %lf %lf", &val, &bv);
	    if (l == 2) {
		a->pb2.blue = min (val, bv);
		a->pb2.red  = max (val, bv);
		a->pb2act = 1;
	    }
	} else if (strcmp (key, "EXORDER") == 0) {
/* EXORDER  Read flag to de-conflict extra orders */
	    sscanf (rest, " %d", &(a->ex_order));
	    if (a->ex_order) a->Flag |= OFxrd;
	} else if (strcmp (key, "TEMPERATURE") == 0) {
/* TEMPERATURE  read expected temperature  */
//	    sscanf (rest, " %lf", &(a->temp));
	    sscanf (rest, " %lf", &val);
/* Limit temperature read to resonable values... */
	    if (fabs(val) <= 40.0) a->temp = val;
	} else if (strcmp (key, "HANGLE") == 0) {
/* HANGLE  read expected hour angle in hours */
	    if ( (sscanf (rest, " %lf", &(a->ha))) != 1) a->ha = -25.0;
	} else if (strcmp (key, "OBJFILE") == 0) {
/* OBJFILE  read an object filespec, add to list of filespecs.  */
	    a->objfiles = pushtxt (rest, a->objfiles);
// problem -- dlevel has not been defined at this point...
if (dodebug > 1) {
	printf (" +Object file name = %s\n", a->objfiles->name);
}
	} else if (strcmp (key, "REPOBJ") == 0) {
/* REPOBJ  read repeat object limit */
	    sscanf (rest, " %d", &(a->reuseob));
	} else if (strcmp (key, "REPREF") == 0) {
/* REPREF  read repeat reference object limit */
	    sscanf (rest, " %d", &(a->reuserf));
	} else if (strcmp (key, "REFSEL") == 0) {
/* REFSEL  read reference object selection value */
	    sscanf (rest, " %d", &(a->refsel));
	} else if (strcmp (key, "REFLIMIT") == 0) {
/* REFLIMIT  read reference object limit value */
	    sscanf (rest, " %d", &j);
	    if (j < 4 && j > 0) j = 4;
	    if (j < 99) a->reflimit = j;
	} else if (strcmp (key, "IFU") == 0) {
/*  IFU   read IFU mode value */
	    sscanf (rest, " %d", &(a->IFUmode));
	} else if (strcmp (key, "OFFCENTER") == 0) {
/* OFFCENTER  read IFU object ra and dec from line  */
	    a->IFUra  = posvalue (strtok (rest, WHITESPACE)) * Hour;
	    a->IFUdec = posvalue (strtok (NULL, WHITESPACE)) * Degree;
	} else if (strcmp (key, "DEBUGFILE") == 0) {
/* DEBUGFILE -- create a debug output file */
	    dbfile = iopen (rest, "w");
	} else if (strcmp (key, "ORDER") == 0) {
/* ORDER  read a requested disperser order */
	    sscanf (rest, " %d", &(a->order));
	} else if (strcmp (key, "GAPS") == 0) {
/* GAPS -- process the gaps feature value */
	    k = sscanf (rest, " %d", &j);
	    if (k == 1) gpw = dbw = j;
	} else if (strcmp (key, "DATE") == 0) {
	    a->date = intvalue (rest);
// 	} else if (strcmp (key, "DBOUND") == 0) {
// /* DBOUND -- process the detector boundries value */
// 	    k = sscanf (rest, " %d", &j);
// 	    if (k == 1) dbw = j;

//	} else if (strcmp (key, "XXX") == 0) {
	} else {	/* unrecognized */
	    printf ("** Keyword \"%s\" unrecognized in obs file.\n", key);
	}  /* End decisions on keyword */
    }	/* End of input loop */

//xxooxx

    a->comments = push_name (comq, a->comments);
    a->pcom = push_name (pcomq, a->pcom);

/* Rectify the bounds and gaps values... */
    j = 0;
    if (gpw > 0) j = gpw - 1;
//    if (dbw > 0) j = dbw - 1;
    j &= 1;
    a->Dstat = 0;
    k = (a->Instr_mt == NULL) ? -1 : InstIndx(a->Instr_mt->name);
    if (k == I_imacs) {
//    if (InstIndx(a->Instr_mt->name) == I_imacs) {
	if (j)   a->Dstat |= DS_ns;
	if (gpw) a->Dstat |= DS_gaps;
	if (dbw) a->Dstat |= DS_orient;
    }

/* Note that this program is used by various other applications to
  read an .obs file to a structure.  Some don't need the object
  queues to be filled in, so we leave that as a separate step.
  [ just call  read_objectqueue (*Obs) for that. ]  */

// normalize the obs structure; optics, telescope, disperser
    normalize_obs (a);

    return  a;
}

/* Program to check/set default temperature */
/* Called only by maskgen; put here to get DefTemp definition */
void	check_temp (Obs* obs, gsto* Chead)
{
int	myinst;
double	*ct;

/* In the case that the instrument is ldss (or others in future), and
  the temp in the obs struct is near DefTemp, we set the temp in obs
  structure to the cutting temperature found from Chead so that the
  temperature compensation is defaulted to null.  */

// Get rid of null pointers...
    if (obs == NULL) return;
    if (Chead == NULL) return;

/* Return early if it is not the LDSS-3 instrument */
    myinst = (obs->Instr_mt == NULL) ? -1 : InstIndx (obs->Instr_mt->name);
    if (myinst != I_lds3) return;

/* If the temperature currently in obs is not sufficiently near
  the default value, we return; it is probably being set on
  purpose to cause temperature compensation.  */
    if (fabs (obs->temp - DefTemp) > 0.02) return;

// Set the obs structure temperature to the "right" value.
    ct = (double*)GS_lookup (&Chead, myinst, CD_tcut, 0);
    if (ct != NULL) obs->temp = *ct;
    else            obs->temp = 20.0;

}


	/* - - - - - - - - - - *\
	|			|
	|   Read Object list	|
	|			|
	\* - - - - - - - - - - */

// This below -----
// is designed to be internal to the read_objectlist routine,
// where it is used.  See how to do that.  Pascal can have
// internal subroutines, why not C?
static	void	dumq (namlist* h, char* p)   /* Local Internal subroutine */
    {
    namlist	*r;
	if (h == NULL) return;
	C_loop (r, h) { fprintf (stderr, "%s%s\n", p, r->name); }
	h = kill_name (h);
    }  /* End local subroutine */
//  end -----


objq*  read_objectlist (FILE*  input, Obs* ob)
/* Read the object list file, return an object queue.  Comments are
  repeated on stderr messages, as is the title */
{

/* How we do it:
  Read data records, identify by first character.
  If comment, see if pre (.) or post (-) queue needs to be made;
  At first object record, all post queue comments are linked to
  the header comment queue.
  Hold a pre and post comment queue, move these to the object queue
  header when an object is detected.
  At object or reference star, decode information into a new
  objq entry.
  At end, put all post comments into last object.
..*/
objq*	a = NULL;
objq*	b = NULL;
int	length;
int	type;
int	radflag=0;
double	scale;
char	buffer[OBJECTLIST_ENTRY+1];
namlist*	obcomq = NULL;
namlist*	preq = NULL;
namlist*	postq = NULL;

/* Check for nulls;  */
    if (input == NULL) return NULL;

/* Get value for scaling slit dimensions... */
    scale = foclen (ob) * ArcSecond;
/* "scale" is size of one arc second in mm. */

    radflag = 0;
/* Loop to read data from object list file */
    while ( fgets (buffer, OBJECTLIST_ENTRY, input) != NULL) {
/* Clean out terminating newline from buffer */
	trimend (buffer, NEWLINES);
/* Trim any leading whitespace, and check for null line */
	trimbeg (buffer, WHITESPACE);
	length = strlen (buffer);
	if (length < 1) continue;
/* Deal with possible in-line comment.  If a ! or # character is
  found that is not the first non-blank character, set it null. */
	length = strcspn (buffer, "!#");
	if (length > 0) buffer[length] = AN;
/* Trim any trailing whitespace also, and check for null */
	trimend (buffer, WHITESPACE);
	length = strlen (buffer);
	if (length < 1) continue;

/* Deal with individual records by first character */
	type = OBJ_UNKNOWN;
	switch (*buffer) {
	case '$':		/* Title */
#ifdef  debug
		printf ("-Title: <%s>\n", buffer+1);
#endif
	    if (ob->dlevel > 0) fprintf (stderr,
		"Object file title: %s\n", buffer+1);
/* Title is not otherwise used or remembered */
	    break;
	case '!':		/* Comment */
/* Check for and save pre (.), post (-) and global (=) comment queues, all
  others are mercifully forgotten */
	    if (buffer[1] == '.') {
#ifdef  debug
		printf ("-Precomment: <%s>\n", buffer+1);
#endif
		preq = pushtxt (buffer+2, preq);
	    } else if (buffer[1] == '-') {
#ifdef  debug
		printf ("-Postcomment: <%s>\n", buffer+1);
#endif
		postq = pushtxt (buffer+2, postq);
	    } else if (buffer[1] == '=') {
#ifdef  debug
		printf ("-Globalcomment: <%s>\n", buffer+1);
#endif
		obcomq = pushtxt (buffer+2, obcomq);
	    } else {
/* Ignore me */
	    }
	    break;
	case '#':		/* Comment */
#ifdef  debug
		printf (" Alternate comment <%s>\n", buffer+1);
#endif
	    ;
	    break;
	case '&':		/* Keyword introducer */
	    if (!strncasecmp(buffer+1, "RADEGREE", 8)) radflag = 1;
	    break;
/* NOTE:  NEW FEATURE.
  Add a special entry line to flag the case in which right ascension will
  be entered in degrees.  This is sometimes done in a file using decimal
  degrees for both coordinates.
  Assign a special leading character to keyword entries, try "&"
  Follow that with a keyword, e.g. "RADEGREE"
  This flag will be turned on by the special entry; it is normally off.  */
	case 'B':		/* Object */
	case '@':		/* Object */
	    type = OBJ_OBJECT;	/* Intentional fall-through */
	case '*':		/* Reference Star */
	    if (type == OBJ_UNKNOWN) type = OBJ_REFERENCE;
#ifdef  debug
		printf ("-Object %d:  <%s>\n", type, buffer+1);
#endif
	/* Put existing postcomment queue in last or global location */
	    if (a == NULL) obcomq = push_name (postq, obcomq);
	    else      a->last->dat.postcom = postq;
	    postq = NULL;
	/* Fill in a new object entry in the header queue */

	    b = get_object (buffer+1, type, preq, ob);	/* New object */
	    if (radflag) b->dat.ra /= 15.0;	/* Convert scale */
// Scale the slit values for this object; slit values for objects
// are in millimeters; the defaults and input are in arc seconds
	    b->slit.width *= scale;
	    b->slit.alen  *= scale;
	    b->slit.blen  *= scale;
/* Enforce shape restriction on alen, blen; since these are used to
  actually compute the slit ends for spectrum placement... */
	    if (b->slit.shape == SQUARE ||
		b->slit.shape == CIRCLE )
		b->slit.alen = b->slit.blen =
		b->slit.width / 2.0;
	    a = push_obj (b, a);
	    preq = NULL;
	    break;
/* Note -- only difference between @ and * cases is OBJ_OBJECT vs
  OBJ_REFERENCE; thus they are combined with statements selecting
  a variable to the correct type... */
	default:		/* Unknown data record */
/* Did not find character - error message needed */
	    printf ("-Unidentified object error record: <%s>\n", buffer);
	}  /* End switch on first character */
    }	/* End of input loop */

/* Assign any remaining postq comments to proper place */
    if (a == NULL) obcomq = push_name (postq, obcomq);
    else      a->last->dat.postcom = postq;
    postq = NULL;

/* Send the obcomq to the stderr messages, and then delete
  the list pointer queue. */
/* Also same for preq */

/* Send comments remaining in global and precomment queues to the
  stderr file with a prefix; then kill the queued commentary. */
// Want local routine to be here -----
    dumq (preq,   " !!-");
    dumq (obcomq, " !!=");

/* Input is done, return the dynamic object queue */
    return  a;
}


objq*	read_objectqueue (Obs* ob)
{
namlist	*r;		// Scratch
objq	*q = NULL;	// result
objq	*f;		// scratch
int	j=0;
FILE	*INF;		// scratch

    if (ob == NULL) return NULL;
/* Loop over the namlist of object files, reading each file into
  an object queue.  Those queues are made into the actual object
  queue which is returned.      */
    C_loop (r, ob->objfiles) {
	INF = iopen (r->name, "r");
	if (INF != NULL) {
	    f = read_objectlist (INF, ob);
	    fclose (INF);
	    q = push_obj (f, q);	// Put found into queue
	}
    }

/* Fill in the order values in the result queue */
/* This sets the ordering for the object queue.  It should be globally
  correct, since the entire queue has been read in here, and no other
  object reading is scheduled at this time. */
    C_loop (f, q) { f->order = j++; }

/* Result is the object queue. */
    return  q;
}

void	fill_objects (Obs* ob)
/* Read all objects specified into the obs structure's object
  queue pointer from the obs structure's object file specifications. */
{

/* Take care of any trivial cases */
    if (ob == NULL) return;

/* Clean out any existing object queue */
    ob->ob = kill_objq (ob->ob);

/* Obtain the new object queue, put in structure. */
    ob->ob = read_objectqueue (ob);
}


	/* - - - - - - - - - - *\
	|			|
	|   Conflict geometry	|
	|			|
	\* - - - - - - - - - - */


/*  ===  Conflict geometry  ===  */


int	within (double a, double x, double b)
/* Return for x relative to ab interval will be
  0 (a), 1 (in) or 2 (b) side of the interval. */
/* x equal to a or b is condisered inside */
{
    if (a < b) {
	if (x < a) return 0;
	return x > b ? 2 : 1;
    } else {
	if (x > a) return 0;
	return x < b ? 2 : 1;
    }
}


bool	overlap (double x, double y, double a, double b)
/* True if xy interval overlaps ab interval */
{
/* The simple, yet inefficient way, is to use the above */
/*    if (within (a, x, b) == 1) return true; */
/*    if (within (a, y, b) == 1) return true; */
/*    return false; */

/* We can speed things by not comparing what we know: */
    if (a < b) {
	if (x < a) {		/* x outside a side */
	    return (y >= a);
	} else if (x > b) {	/* x outside b side */
	    return (y <= b);
	} else {		/* x inside */
	    return true;
	}
    } else {
	if (x > a) {
	    return (y <= a);
	} else if (x < b) {
	    return (y >= b);
	} else {
	    return true;
	}
    }
//    return  0;
// a return statement here can't be reached...
}


	/* - - - - - - - - - - - - - - *\
	|				|
	|   Differential  Refraction	|
	|				|
	\* - - - - - - - - - - - - - - */

/* Differential Refraction support utilities */
#ifdef  DIFREF
static	vect3	drefractn (vect3 pos, vect2 zp, double rc)
/* Return pos as augmented by differential refraction */
{
vect2	p;
vect2	d;
vect3	r;

/* Short work of trivial case... */
    if (rc == 0.0) return pos;

/* Find the differential */
    p = get2de3v (pos);
    d = difref (p, zp);

/* Add to pos, and normalize */
    p = sum2vect (p, mul2vect(d, rc));
    r = get3de2v (p, sqrt (1.0 - p.x*p.x - p.y*p.y) );

/* Return result */
    return  r;
}
#endif	// on DIFREF

#ifdef  DIFREF
static	vect3	urefractn (vect3 pos, vect2 zp, double rc)
/* Return vector which, when sent to drefractn, would return pos. */
{
vect2	p;
vect2	d;
vect2	s;
vect2	q;
vect3	r;

/* Short work of trivial case... */
    if (rc == 0.0) return pos;

/* Find the differential */
    q = p = get2de3v (pos);
    d = make2vect (rc, rc);
    s = make2vect (0.0, 0.0);

/* Find the un-refracted vector */
    while (vect2norm(sub2vect(d,s)) > 1.4e-8) {
	s = d;
	d = mul2vect (difref (q, zp), rc);
	q = sub2vect (p, d);
    }

    r = get3de2v (q, sqrt (1.0 - q.x*q.x - q.y*q.y) );

/* Return result */
    return  r;
}
#endif  // on DIFREF

	/* - - - - - - - - - - *\
	|			|
	|  Transform Geometry	|
	|			|
	\* - - - - - - - - - - */

#undef  MPDBUG
// Above debugs mask vector position comparison...

vect2	maskvect (Obs *obs, double ra, double dec)
/* Obtain the slit mask position vector for a celestial position */
{
vect3	a, v;
vect2	r;
const	vect2	bad = { 1.e5, 0.0 };
vect3	p;
// debug quantities here...
#ifdef  MPDBUG
static  int  k = 0;	//debug
double	f;	//debug
#endif

/* Get the celestial position vector */
    v = sph2vec (ra, dec);  /* equinox, 6H, NCP */

/* Use transform matrix... */
//    v = mmult3 ( obs->mpos, v );
    v = mmul3v ( obs->mxpos, v );

/* Telescope optics are needed HERE to change v, a collimated position
  vector, into r, a focal plane vector... */
/* It will also need an inverse for the bplot stuff? */

/* Now, x,y,z is really z,-x,y, and we need to multiply by proper focal
 length, and re-compute z based on curvature */
    a.x = -v.y;
    a.y =  v.z;
    a.z =  v.x;

// Really should fix this transform in the obs->mpos matrix
// to eliminate this step...
// That would involve the get_mpos routine...

/* Check for wildly wrong positions, over 45 degrees or so */
    if (a.z < 0.70) return bad;

/* Apply differential refraction to vector a here */
#ifdef  DIFREF
    a = drefractn (a, obs->zenith, obs->cref);
#endif

/* Compute in (one) of two ways depending on element type */
#ifndef  MPDBUG
    if (obs->Tel_scope == NULL) {	// missing element*
#endif
	v = mul3vect (a, foclen(obs));
// We allow focal_length here as new element is not present...
	r = get2de3v (v);
#ifdef  MPDBUG
    f = sqrt (a.x*a.x + a.y*a.y) / ArcSecond;	// debug
    if (obs->Tel_scope != NULL) {
#else
    } else {		// element* present
#endif

/* Do the accurate computation if we have a non-null element here... */
	p =  Op_transform (a, obs->Tel_scope, obs->sw, obs->temp);
// Now using correct OPTICS transform
#ifdef  MPDBUG
// REPORT OLD/NEW POSITIONS FOR FIRST 20 OR SO CASES
	if (++k < 10)  {
	    printf (" -- Old %8.3f %8.3f (%8.3f) = (%8.3f arcsec.)\n",
		    r.x, r.y, vect2norm (r), f );
	    printf ("    New %8.3f %8.3f (%8.3f)\n",
		    p.x, p.y, vect2norm (get2de3v(p)) );
	}
// THE DEBUG PRINT APPEARS TO BE CORRECT 03/3-4 KDC
#endif
	r = get2de3v (p);
    }
    return  r;
}

/* NOTE! -- The above will need to use the optic general transform
  for the telescope to put in distortion terms.  We also want to
  have an inverse, obtaining ra/dec from a mask position.  See the
  plotting (bplot.c) for an example of inverting it. */

/* NOTE -- we may need a routine to make a mask 3-vector by adding
  curvature data from telescope.  We may also need routines to make
  2 from 3 or 3 from 2 dimensional vectors. */

/* Following is the inverse of maskvect... */

void    unmaskvect (Obs* ob, vect2 p, double* ra, double* dec) {
/* Inverse of maskvect routine; find ra/dec from object position */
vect3	a, r, v;

// Inverse of the method by which ra/dec are obtained in maskvect
    r = get3de2v (p, 0.0);
    if (ob->Tel_scope == NULL) {
	a = mul3vect (r, 1.0/foclen(ob));
	a.z = sqrt (1.0 - a.x*a.x - a.y*a.y);
    } else {
	a = Inv_transform (r, ob->Tel_scope, ob->sw, ob->temp);
    }

/* Remove the differential refraction here */
#ifdef  DIFREF
    a = urefractn (a, ob->zenith, ob->cref);
#endif

    v.x =  a.z;
    v.y = -a.x;
    v.z =  a.y;
//    v = tmult3 (ob->mpos, v);
    v = mmul3v (mtr3 (ob->mxpos), v);
    vec2sph (ra, dec, v);
// Above obtained from wing_set; and could be used by it.
}

void	pbnormal (passband *pb)
/* Make sure red > blue, as in wavelength */
{
double	w;
    if (pb == NULL) return;
    if (pb->red >= pb->blue) return;
    w = pb->red;
    pb->red = pb->blue;
    pb->blue = w;
    return;
}


#define	DIRECT_ELEMENT_NAME  "*Direct*"
// ==direct==

static	element	*Dir_Inst (element *ex)
/* Recursive routine:  Returns direct version of the given
  instrument element.  If element is not an instrument,
  or is null, simply return the given pointer.
  When an instrument, make an element to hold the new
  instrument, construct the new name by appending DIRECT_ELEMENT_NAME
  and cycle through the data elements, copying all those which are
  not instrument elements.  If an instrument element is found, call
  this routine recursively to obtain the direct version.  Put all
  the things and the direct version of instruments into the list
  by name.  Then, push the new element after the given instrument
  and return its element pointer.		--  See drimage.otl  --
..*/
{
char	*c;	// Scratch
int	j, n;	// Scratch
char	dname[128];
element	*re=NULL;	// Return value
datlist	*d;	// for element list
datlist	*q;	// index pointer
datlist	*dle;	// For list construction
element	*ep;	// An element pointer
element	*rx;	// recursing pointer

/* Get rid of null cases */
    if (ex == NULL) return ex;		// No null pointer
    if (ex->kw != INSTrument) return ex;	// Not instrument
    n = strlen (DIRECT_ELEMENT_NAME);	// Check for special name,
    j = strlen (ex->name);		// at end of instrument name
    if (j > n) {
	c = ex->name + (j-n);
	if (!strcmp(c, "DIRECT_ELEMENT_NAME")) return ex;
    }

/* Check if this has already been done -- construct special name of
  direct element and see if such an element already exists, and
  return it if so. */
    sprintf (dname, "%s%s", ex->name, DIRECT_ELEMENT_NAME);
    re = find_element (ex, dname);
    if (re != NULL) return re;

/* Now, get serious.  We need to make this one, so we get an element
  block and put it into re; then set values appropriately. */
    re = xalloc (element);
    re->next = re->last = re;
    re->name = dynamstr (dname);
    re->kw   = INSTrument;
    re->flag = ex->flag & (~1);
    re->data = NULL;
    re->head = NULL;
    dle = NULL;

/* Construct the instrument by copying the element list from the
  source, with appropriate changes where found... */
    C_loop (q,ex->head) {
// Look for special names here -- grating/grism
	if (!strcasecmp(q->data, "GRATING")) {
	    d = newdat ("IMACS_direct");
	} else if (!strcasecmp(q->data, "GRISM")) {
	    d = newdat ("IMACS_direct_grism");
	} else {
/* We have an actual element pointer... */
	    ep = q->eg;
	    if (ep == NULL) {	// find element if any
		ep = find_element (ex, q->data);
		q->eg = ep;
	    }
	    if (ep == ex) {	// bad error - re-entrant
	printf (" ** Re-entrant directed instrument %s canceled.\n",
				ex->name);
		re = free_element (re);
		break;
	    }
	    if (ep == NULL) {	// Not a real element, fake it.
		d = newdat (q->data);
	    } else {		// Real element, see if instrument, disperser
		if (ep->kw == INSTrument) {	// is instrument
// Here's the recursive part!
		    rx = Dir_Inst (ep);
		    d = newdat (rx->name);
		    d->eg = rx;
		} else if (ep->kw == GRAting) {
		    d = newdat ("IMACS_direct");
		} else if (ep->kw == GRIsm || ep->kw == EGRIsm) {
		    d = newdat ("IMACS_direct_grism");
		} else {		// not instrument or special
		    d = newdat (ep->name);
		    d->eg = ep;
		}
	    }
	}	// Test for special names
	if (dle == NULL) re->head = d;
	else            dle->next = d;
	dle = d;
	if (d->eg == NULL) {	// Get element if possible
	    d->eg = find_element (ex, d->data);
	}
    }	// Loop constructing instrument list
    re->data = re->head;

/* Put the constructed direct element (re) into the element queue
  following the given element (ex).  */
//    rx = push_element (ex, re);
    rx = G_push_c (ex, re);	// rx becomes re here

/* Return the newly constructed direct instrument */
    return  re;
}

// ==direct==
static	element	*Dinstr (Obs *obs)
/* Return pointer to the direct instrument coresponding to
 the actual instrument indicated by the obs structure.  If
 necessary, manufacture such an instrument and put it into
 the optic element queue.  Otherwise, return the pointer.  */
{
element	*r=NULL;

/* Normally, the element pointer is stored in obs, we get that
  one if it is there. */
    if (obs->Direct_inst != NULL) return obs->Direct_inst;

/* If we don't have it (yet), we obtain the direct instrument
  from the current instrument in the obs structure.  This will
  construct an appropriate instrument if necessary. */

    r = Dir_Inst (obs->Instr_mt);
    obs->Direct_inst = r;

    return  r;
}


/* The "det_pos" program is where we go to find the mapping from a
  slit mask position (sm) to a detector position ( a 2-vector ) using
  the optical transform code and the instrument defined in the Obs
  structure (pointed to by obs->Instr_mt).  */
/* The other use of Op_transform is the "maskvect" program, above,
  which is used to find the slit mask position from celestial
  coordinates and the Obs structure.  */

// #undef	DPOLD

#undef	DPDBUG

vect2	det_pos (Obs *obs, double lambda, vect3 sm, int order)
/* Return is detector position, via optics code. */
// Correct to use OPTICS here
{
vect2	r;
vect3	p;
int	*nord;
// debug variables
#ifdef	DPDBUG
static	int	k=0;	// debug
#endif

    if (obs->Instr_mt == NULL) {
/* If instrument is null, we can't transform, return zero */
	r = make2vect (0.0, 0.0);
    } else {
/* Save, and later restore, the global order pointer. */
	nord = GRorder;
	GRorder = &(order);
	p = Op_transform (sm, obs->Instr_mt, lambda, obs->temp);

#ifdef  DPDBUG
	if (++k < 10) {
	    printf (" -- Old d.p. %8.3f %8.3f\n",
		r.x, r.y );
	    printf ("    New d.p. %8.3f %8.3f\n",
		p.x, p.y );
	}
#endif
	r = get2de3v (p);
/* Restore global order pointer before returning */
	GRorder = nord;
    }
    return  r;
}

/* Inverse of det_pos (pos_det, eh?) */

vect2	pos_det (Obs *obs, double lambda, vect2 dp, int order)
/* Return is slit mask position, via optics inversion. */
{
vect3	vp;
vect2	r;
vect3	p;
int	*nord;

    if (obs->Instr_mt == NULL) {
/* If instrument is null, we can't transform, return zero */
	r = make2vect (0.0, 0.0);
    } else {
/* Save, and later restore, the global order pointer. */
	vp = get3de2v (dp, 0.0);
	nord = GRorder;
	GRorder = &(order);
	p = Inv_transform (vp, obs->Instr_mt, lambda, obs->temp);

	r = get2de3v (p);
/* Restore global order pointer before returning */
	GRorder = nord;
    }
    return  r;
}

// ==direct==
static	vect2	dir_pos (Obs *obs, vect3 sm)
/* Return direct image position... */
{
vect2	r={0.0,0.0};
vect3	p;
int	*nord;
int	dord=0;
double	*nang;
double	dang=22.5*Degree;	// NOTE: LITERAL ASSUMED VALUE

/* Check null case */
    if (obs == NULL) return r;

/* Check that a direct instrument actually exists */
    if (obs->Direct_inst == NULL) obs->Direct_inst = Dinstr (obs);

// Much below was copied from det_pos, but with changes...

    if (obs->Direct_inst != NULL) {
/* Save, and later restore, the global order pointer. */
	nord = GRorder;
	GRorder = &(dord);
/* Save, and later restore, the global angle pointer */
	nang = GRangle;
	GRangle = &(dang);
	p = Op_transform (sm, obs->Direct_inst, obs->sw, obs->temp);
	r = get2de3v (p);
/* Restore global order, angle pointers before returning */
	GRorder = nord;
	GRangle = nang;
    }
    return  r;	// r is zero if instrument was null
}

	/* - - - - - - - - - - *\
	|			|
	|  Spectrum Geometry	|
	|			|
	\* - - - - - - - - - - */

sp_edge	getedge (Obs *obs, objq *ob, vect2 slend, double ovl,
		double xcent, int order, passband wl)
/* Compute slit edge structure from obs, ob and slit end vector; allow
  overlap by ovl, move edge -right- by that amount */
{
sp_edge	v;
vect2	s;
vect3	ss;
vect2	rd, bd;
double	y;
double	d;
double	cwavl;		/* Center wavelength */
/* Obtain the vector position of the slit end */
	TX_beg (TS_geteg);
    s = sum2vect (ob->smpos, slend);
    ss = get3de2v (s, 0.0);	/* For optics. */
    cwavl = pbmean (wl);	/* Center of passband used here. */
/* ==NOTE== Alternatively, find z from image surface curvature */
    y = s.y;
/* Compute the center, red and blue spectrum points */
    v.p = det_pos (obs, cwavl, ss, order);
/* Find the proper ovl sign.  Negate it if the center of the edge is
  less than the center of the spectrum; otherwise it is left as is.  */
    if (v.p.x < xcent) ovl = -ovl;

	v.p.x += ovl;
    s.y = y + obs->ddir.red * ob->slit.width;
/* NOTE s IS NOT USED HERE; ss IS; DO WE MEAN ss ABOVE???? */
    v.r = det_pos (obs, wl.red, ss, order);
	v.r.x += ovl;
    s.y = y + obs->ddir.blue * ob->slit.width;
    v.u = det_pos (obs, wl.blue, ss, order);
	v.u.x += ovl;
/* Find the polynomial coeficients */
    if (order == 0) {
/* We use a=0 and simplify computation here */
	v.a = 0.0;
	rd = sub2vect (v.r, v.u);
	v.b = (rd.y == 0.0) ? 0.0 : rd.x / rd.y;
    } else {
/* Poly is x-x0 = a(y-y0)^2 + b(y-y0) */
	rd = sub2vect (v.r, v.p);
	bd = sub2vect (v.u, v.p);
	d = rd.y*bd.y*(rd.y - bd.y);
	v.a = (rd.x*bd.y - bd.x*rd.y) / d;

//	v.b = (rd.x*bd.y*bd.y  - bd.x*rd.y*rd.y) / d;
/* DEBUG -- The above computation is wrong!  */
/* How to solve correctly:
   dx = a * dy * dy + b * dy; and we have two cases of dx and dy;
  these cases are the r and b cases, in vectors rd, bd.  The 2 equations are
	rdx = (a * rdy + b) * rdy
	bdx = (a * bdy + b) * bdy
  we want to solve for a, b here.
	rdx/rdy = a * rdy + b
	bdx/bdy = a * bdy + b
  subtracting these two, we get
	rdx/rdy - bdx/bdy = a * (rdy - bdy)
  which yields a = (rdx/rdy - bdx/bdy) / (rdy - bdy)
  and b may be found from wither equation then.
temporarily, put this solution in and see if it works: yes.
so, take d = rdy*bdy*(rdy-bdy) as done above,
  and a * d = rdx*bdy - bdx*rdy, so a was correct, b was not.
  Test that by dropping the v.a computation below as not needed: yes.
  So, we find b as either:
	b = rdx/rdy - a * rdy
	b = bdx/bdy - a * bdy
  Multiply both by d and substitute a to get:
	b*d = rdx*bdy(rdy-bdy) - rdy*(rdx*bdy - bdx*rdy)
	b*d = bdx*rdy(rdy-bdy) - bdy*(rdx*bdy - bdx*rdy)
  Multiply out both equations:
	b*d = rdx*bdy*rdy - rdx*bdy*bdy - rdy*rdx*bdy + rdy*bdx*rdy
	b*d = bdx*rdy*rdy - bdx*rdy*bdy - bdy*rdx*bdy + bdy*bdx*rdy
  Canceling terms yields:
	b*d = rdy*rdy*bdx - rds*bdy*bdy
	b*d = bdx*rdy*rdy - bdy*rdx*bdy
  Which are the same, giving a better expression for b, below:
  (the sign on b was inverted)
... */
//	v.a = (rd.x/rd.y - bd.x/bd.y) / (rd.y - bd.y);
//	v.b = rd.x/rd.y - v.a * rd.y;
/* Fix it up later if this is correct...  Fixed below: */

	v.b = (bd.x*rd.y*rd.y - rd.x*bd.y*bd.y) / d;

/* Compute poly as x-x0 = (y-y0)*(a*(y-y0) + b) */
    }
	TX_enm (TS_geteg, "getedge subroutine.");
    return  v;
}


vect2  lslit (slit s)
/* Return the left end of the slit */
{
vect2	p;
    p.x = -s.alen * cos(s.angle);
    p.y = -s.alen * sin(s.angle);
    return  p;
}

vect2  rslit (slit s)
/* Return the right end of the slit */
{
vect2	p;
    p.x = s.blen * cos(s.angle);
    p.y = s.blen * sin(s.angle);
    return  p;
}

// wslit used in smdfps, smplot and dbplot
vect2	wslit (slit s)
/* Return the (half) width vector of slit */
{
vect2	p;
    p.x = 0.5 * s.width * -sin(s.angle);
    p.y = 0.5 * s.width *  cos(s.angle);
    return  p;
}

// ==direct==
static	bbox	direct_loc (Obs *obs, objq *ob)
/* Return a bounding box in detector space for the direct image of
  the given object, using direct instrument */
{
bbox	rb;	// Result to be returned.
vect2	dp;	// Detector position
vect2	se;	// slit end
vect2	cv;	// corner vector
double	w;	// scratch
int	k;	// scratch

/* Initialize things.  We hope our pointers are not null... */
    dp = dir_pos (obs, get3de2v (ob->smpos, 0.0));	// like det_pos
    rb = boundbv2 (dp, dp);

/* Some prelim. */
    w = 0.5 * ob->slit.width;
//    cv.x = -w * sin(ob->slit.angle);
//    cv.y =  w * cos(ob->slit.angle);
//    cv = rot2 (pol2vec (w, ob->slit.angle), M_PI/2.0);
    cv = rot2 (pol2vec (w, ob->slit.angle), M_PI_2);

/* What is done depends on what kind of mask hole we have */
    switch (ob->slit.shape) {
	case CIRCLE:
	    dp = dir_pos (obs, get3de2v(sum2vect(ob->smpos, cv), 0.0));
	    rb = boundex (rb, dp);
	    dp = dir_pos (obs, get3de2v(sub2vect(ob->smpos, cv), 0.0));
	    rb = boundex (rb, dp);
	    se.x =  cv.y;		// simply rotating cv by
	    se.y = -cv.x;		// pi/2 here, quickly
	    dp = dir_pos (obs, get3de2v(sum2vect(ob->smpos, se), 0.0));
	    rb = boundex (rb, dp);
	    dp = dir_pos (obs, get3de2v(sub2vect(ob->smpos, se), 0.0));
	    rb = boundex (rb, dp);
	    break;
	case SQUARE:
// Fall through, since these things are the same anyway
	case RECTANGLE:
	    se = sum2vect (ob->smpos, lslit (ob->slit));
	    dp = dir_pos (obs, get3de2v(sum2vect(se, cv), 0.0));
	    rb = boundex (rb, dp);
	    dp = dir_pos (obs, get3de2v(sub2vect(se, cv), 0.0));
	    rb = boundex (rb, dp);
	    se = sum2vect (ob->smpos, rslit (ob->slit));
	    dp = dir_pos (obs, get3de2v(sum2vect(se, cv), 0.0));
	    rb = boundex (rb, dp);
	    dp = dir_pos (obs, get3de2v(sub2vect(se, cv), 0.0));
	    rb = boundex (rb, dp);
	    break;
	case SPACED_CROSS:
/* Consists of hole with diameter of slitwidth, and a slit of length alen
  starting blen outside the hole edge.  4 slits are around the hole */
	    w = 0.5 * ob->slit.width + ob->slit.blen + ob->slit.alen;
	    cv = pol2vec (w, ob->slit.angle);
	    for (k=0; k<4; k++) {
		dp = dir_pos (obs, get3de2v(sum2vect(ob->smpos, cv), 0.0));
		rb = boundex (rb, dp);
//		cv = rot2 (cv, M_PI/2.0);
		cv = rot2 (cv, M_PI_2);
	    }
    }	// switch on slit shape

    return  rb;
}

/* --NOTE-- check out that this is what is desired; may have to
  fix the angles, or add a dimension scale factor.  Also, put the
  prototypes of these functions into the header file */

spect	*getspect (Obs *obs, objq *ob, int order, passband wl)
/* Create and return spectrum description structure */
{
spect	*p;
vect3	ss;	/* DEBUG */
double	spoverlap;
double	cwavl;

    if (!(ob->flag & OBJECT_ACTIVE)) return NULL;
// NOTE -- APPARENTLY if we omit spectra for inactive objects, many
// conflicts are allowed onto the mask.  Some test for sod null is
// possibly short-circuiting a loop!  Look at all val->sod stuff...

/* If ob->type == OBU_GAP or if (ob->flag & OBJECT_VIRTUAL) the spectrum
  should already be computed, and we should not be called... */

/*    p = (spect*) malloc (sizeof(spect)); */
    p = xalloc (spect);
    p->next = NULL;
    p->order = order;
    cwavl = pbmean (wl);

/* DEBUG - get the object center position */
/* Actually, the center position is needed for edge computation */
    ss = get3de2v (ob->smpos, 0.0);  /* For optics. */
    p->center = det_pos (obs, cwavl, ss, order);
// Above uses Op_transform; as does getedge

    spoverlap = 0.5 * obs->Detect->psize.x * obs->minsep;
/* We allow spectra to overlap by minsep by decreasing the spectrum widths
  by spoverlap on each side.  While lslit is on the sky (mask), the value
  of spoverlap is on the detector. */
    p->e1 = getedge (obs, ob, lslit(ob->slit),
	spoverlap, p->center.x, order, wl);
    p->e2 = getedge (obs, ob, rslit(ob->slit),
	spoverlap, p->center.x, order, wl);
// getedge uses det_pos
/* The effect of a positive minsep is to INCREASE the width of the
  spectrum-on-detector */

//    p->bb = boundor (bbedge(p->e1), bbedge(p->e2));
// or could also use:  p->bb = bbspec (p);
    p->bb = bbspec (p);

/* See if bounding box is anywhere on detector */
//    p->on_det = overlap (xl, xh, obs->Detect->low.x, obs->Detect->high.x)
//	&& overlap (yl, yh, obs->Detect->low.y, obs->Detect->high.y);

    p->on_det = (bbover (p->bb, obs->Detect->bb) & 7);
/* Overlap is 2 if bb is entirely on detector, 1 if overlaps, and
  is 4 if detector is entirely within bb (really should not happen).
  So we try & 3 to see what happens... */

/* Special section for secondary bounding box */
    if (obs->pb2act) {	// Find second bounding box
sp_edge	e1, e2;
	e1 = getedge (obs, ob, lslit(ob->slit),
	    spoverlap, p->center.x, order, obs->pb2);
	e2 = getedge (obs, ob, rslit(ob->slit),
	    spoverlap, p->center.x, order, obs->pb2);
// Use similar code to find bounds...
	p->bd = boundor (bbedge(e1), bbedge(e2));
    } else {		// Copy primary to secondary bounds
	p->bd = p->bb;
    }

    return  p;
}


	/* - - - - - - - - - - *\
	|			|
	|    M.O.E.  Support	|
	|			|
	\* - - - - - - - - - - */

/* === Multiple Object Echellette support === */

/* #define  ECON  45575.0  -- Now obsolete, dynamically derived */
// Constant ECON yields center wavelength of any order

/* Replace Echelle CONstant with a subroutine...  These subroutines
  and constants are global, but only in this module; use is by the
  MOE support routines only.  */

static	double	yctran (double wav, double temp, element* es, element* an)
/* Return y coordinate of translated input ray to disperser es */
/* Also use transform by angle an if not null */
{
vect3	inv={0.0, 0.0, 1.0};
vect3	a, b;
    a = Op_transform (inv, es, wav, temp);
    b = (an == NULL) ? a : Op_transform (a, an, wav, temp);
    return  b.y;
}

static double  gcwavx (element *te, double temp, int ord)
/* Return center wavelength of disperser te at order ord */
/* This may be used to replace grwavc at a later time. */
{

//... copied from grwavc (optutils.c) here...

// vect3	inv={0.0, 0.0, 1.0};
double  ah, al;
double  yh, yl;
double  a, y;
double  b;
int     k;
int*	xord;
element*  an;

/* Check for current element (should be disperser) */
    if (te == NULL) {
	printf (" **No current disperser.\n");
	return  0.0;
    }

/* Check for proper element type */
//    if (te->kw != GRIsm && te->kw != DGRIsm) {
//	printf (" **Element %s is not a grism.\n", te->name);
//	return 0.0;
//    }
/* Set an to null unless element type is grating or echelle, then
  set it to imacs angle */
    if (te->kw == GRAting || te->kw == ECHelle) {
	an = find_element (te, "IMACS_angle");
    } else {
	an = NULL;
    }

/* Save and later restore GRorder pointer... */
    xord = GRorder;
    GRorder = &(ord);

/* Transforming from on-axis vector inv through the grism te, yields
  a y-value:   Op_transform (inv, te, cur_wavl, cur_temp).y
  which should be zero when the center wavelendth is current.
  We vary cur_wavl here to find the zero y value.  */


/* Roll through solution space by 500 Angstrom increment, searching for
  the y value to change sign; then do binary search about the interval */
    al = 300.0;		/* Start truly low */
//    yl = Op_transform (inv, te, al, cur_temp).y;
//static	double	yctran (double wav, double temp, element* es, element* an)
	yl = yctran (al, temp, te, an);
    if (yl == 0.0) {GRorder = xord; return al; }
    for (ever) {
        ah = al + 500.0;
//	yh = Op_transform (inv, te, ah, cur_temp).y;
	yh = yctran (ah, temp, te, an);

/* Debug report solution space */
/*      printf ("  ..Wavelength %f y = %f\n", ah, yh); */  /* debug */
        if (yh == 0.0) {GRorder = xord; return  ah; }
        if ( (yl < 0.0) && (yh > 0.0) ) break;
        if ( (yl > 0.0) && (yh < 0.0) ) break;
        if (ah > 15000.0) {
            printf ("  **No wavelength solution found for %s.\n", te->name);
	    GRorder = xord;
            return 0.0;
        }
        al = ah;
        yl = yh;
    }
/* If we get through, yl and yh are different signs, and the
  interval al to ah brackets the solution. */
    k = 0;
    while ( (ah - al) > 1.0e-6 ) {
/* To speed convergence, we combine linear interpolation with a
  binary search.  The proportion of the interpolate coming from the
  linear method increases with each iteration. */

        k++;            /* Iteration counter */
        if (fabs(yh) > fabs(yl)) {      /* Linear interpolation */
            a = al - yl * (ah - al) / (yh - yl);
        } else {
            a = ah - yh * (ah - al) / (yh - yl);
        }
        b = (ah + al) / 2.0;            /* Binary search */
//        a = (a*(double)(2*k) + b) / (double)(2*k+1);
        a = (a*(double)(4*k) + b) / (double)(4*k+1);
/* This usually results in about 12 iterations, rather than the 30 or
  so a strictly binary operation would give.  A strictly linear method
  can get "stuck" making very small changes in the interval. */

//	y = Op_transform (inv, te, a, cur_temp).y;
	y = yctran (a, temp, te, an);

/* Debug trace of solution */
/*      printf ("  ...Wavelength %f y = %f\n", a, y); */ /* debug */
        if (y == 0.0) {GRorder = xord; return  a; }
/* Check which side it is on */
        if ( ((yl < 0.0) && (y < 0.0)) ||
             ((yl > 0.0) && (y > 0.0)) ) {      /* Advance al */
            al = a;
            yl = y;
        } else {                                /* Lower ah */
            ah = a;
            yh = y;
        }
    }
/* Finish off with linear interpolation */
/* (a-al)/(0-yl) = (ah-a)/(yh-0) = (ah-al)/(yh-yl) */
    if (fabs(yh) > fabs(yl)) {
        a = al - yl * (ah - al) / (yh - yl);
    } else {
        a = ah - yh * (ah - al) / (yh - yl);
    }
// printf (" Gwavl Wavelength = %7.3f,", a);  /* debug */
// printf (" Interpolation %d iterations.\n", k);  /* debug */
    GRorder = xord;
    return  a;

}

static double  ewcon (element *te, double temp)
/* Return echelle wavelength constant for given element; use a stored
  computed value if called with same parameters. */
{
static	double	EWC = 0.0;
static	element*  EP = NULL;
static	double	ECT = 0.0;
double	zw;
double	sum;
int	k;
int	ord;

/* Characterize the calls for element and temperature */
    if (EWC != 0.0 && EP == te && ECT == temp) return  EWC;

/* If we fall through that, do a detailed computation,
  saving the parameters for characterization... */
    EP = te;
    ECT = temp;

    for (ord=4,k=0,sum=0.0; ord<11; ord++) {
	zw = gcwavx (te, temp, ord);
	sum += ord * zw;
	k++;
    }

    EWC = sum / (double)k;
    return  EWC;
}

/* NOW, GLOBALLY REPLACE "ECON" WITH ewcon (find_element(gs, "MOE"), temp)
  or something to that effect. */
// NEED TO SUPPLY ELEMENT ADDRESS AND TEMP WHEN GETTING THE CONSTANT.
// HOW TO DO IT:
// in calls from moespect specify element pointer and temp obtained
// at the start of that program
// in calls from the subroutines, use element pointer and temp passed.
// add element pointer and temp in moecw, moeorder, moepb calls.

// keep old stuff behind gtiltmoe symbol undefine; put new ones
// in the gtiltmoe symbol define section...
// same for ECON define, put behind symbol.
// Once made static, remove moecw, moeorder, moepb from mgutils.h !!!

static	double	moecw (int n, element* e, double t)
/* Return center wavelength of given order */
{
//    return ( ECON / (double)n );
    return ( ewcon(e,t) / (double)n );
}

static	int	moeorder (double  wavl, element* e, double t)
/* Return order for given wavelength */
{
int	n;
//    n = floor ( (ECON/wavl) + 0.5);
    n = floor ( (ewcon(e,t)/wavl) + 0.5);
/* We use here the "half order" method. */
/* Alternate computation would be: (the "half way" method)--
	n = floor (ewcon(e,t)/wavl);
	n = (wavl <  (0.5 * (moecw(n,e,t)+moecw(n+1,e,t)) ) ? n : n+1;
.. without using extra storage for 1/2 way point..  */
    return  n;
}

static	passband	moepb (int n, element* e, double t)
/* Return wavelength range for a given order */
{
passband	p;
//    p.red  = ECON / ((double)n - 0.5);
    p.red  = ewcon(e,t) / ((double)n - 0.5);
//    p.blue = ECON / ((double)n + 0.5);
    p.blue = ewcon(e,t) / ((double)n + 0.5);
/* We use here the "half order" method. */
/* Alternate computation would be: (the "half way" method)--
	p.red  = 0.5 * ( moecw(n,e,t) + moecw(n-1,e,t) );
	p.blue = 0.5 * ( moecw(n,e,t) + moecw(n+1,e,t) );
.. just as defined in moeorder above... */
    return  p;
}

spect	*moespect (Obs *obs, objq *ob)
/* Get a Multiple Object Echellette spectrum description for the
  given object.  Return as a spectrum structure... */
{
spect	*p;
vect3	ss;	/* DEBUG */
double	spoverlap = 0.0;
int	lowrd, hiord;
int	midord;
vect2	dpl, dpr;
vect2	cpos;
int	cordr;
int	rr;
double	d;
double	obcw;
static	int	kd=0;	// debug counter
double	t=obs->temp;
element	*e;

/* similar to getspect here */

/* Bypass null cases here...*/
    if (!(ob->flag & OBJECT_ACTIVE)) return NULL;
    p = xalloc (spect);
    p->next = NULL;
//    p->order = order;
// What do we put in for the order?  mid - order ? */
// MAYBE WE MAKE ORDER A FLAG FOR MOE - 0 OR EVEN -1 ?

    e = find_element (obs->elist, "MOE");

/* NOTE getedge uses obs->cw as center wavelength for the order;
  we need to fudge that to being center of this order!  */
    obcw = obs->cw;		// save for later

/* Find the order range... */
    lowrd = moeorder (obs->pb.red, e, t);
    hiord = moeorder (obs->pb.blue, e, t);
    midord = (lowrd+hiord+1)/2;
// If even number of orders, mid is on the blue side,
// Since the spectrum ends bend toward the blue cross-dispersion

/* DEBUG - get the object center position */
/* Actually, the center position is needed for edge computation */
    ss = get3de2v (ob->smpos, 0.0);  /* For optics. */
//    p->center = det_pos (obs, obs->cw, ss, order);
    p->center = det_pos (obs, moecw(midord, e, t), ss, midord);

/* Find the red spectrum edge.  It is the longer edge. */
// get the slit ends of red spectrum, and
// find the end farther from the center position.
// That defines left/right side...
    cordr = (midord > lowrd) ? midord : lowrd + 1;
    cpos = det_pos (obs, moecw(cordr,e,t), ss, cordr);
// That defines a center position at least one order more than lowrd
    dpl = det_pos (obs, moecw(lowrd,e,t),
	get3de2v (sum2vect(ob->smpos,lslit(ob->slit)), 0.0), lowrd);
    dpr = det_pos (obs, moecw(lowrd,e,t),
	get3de2v (sum2vect(ob->smpos,rslit(ob->slit)), 0.0), lowrd);
    rr = ( fabs(dpr.x - cpos.x) > fabs (dpl.x - cpos.x) );
// rr is true if red is right; i.e. if right end of redest order
// is to the right of left end of same order...
// Now, if rr, we get red end from right side order lowrd,
// and if not rr, get red end from left side of same order.

    if (kd++ < 2 && obs->dlevel > 0) {	// debug print of orders...
	printf (" Orders %d (red) to %d (blue); mid=%d; rr=%d\n",
		lowrd, hiord, midord, rr );
//	printf ("  -- cx: %.2f, lx: %.2f, rx %.2f\n",
//		cpos.x, dpl.x, dpr.x );
    }

    obs->cw = moecw(lowrd,e,t);
//    obs->sw = moecw(lowrd,e,t);
    if (rr) {
	p->e1 = getedge (obs, ob, rslit(ob->slit),
	    spoverlap, cpos.x, lowrd, moepb(lowrd,e,t));
    } else {
	p->e1 = getedge (obs, ob, lslit(ob->slit),
	    spoverlap, cpos.x, lowrd, moepb(lowrd,e,t));
    }

/* WE NEED A LITTLE MORE HERE.
  The total y extent, fabs(p->e1.u.y - p->e1.r.y) should be over 1/2 of
  our total detector size, obs->Detect->bb.lo - hi
  If it is not, we need to extend y to make it so.
Add a little to both ends to get to 0.50001 * extent.
Then extend the y and x as done below for side 2.
  This does not seem in practice to occur.
.. */


/* Now, get the blue spectrum edge.  It is shorter, and will be
  of the oposite sense to the red edge.  */

    obs->cw = moecw(hiord,e,t);
//    obs->sw = moecw(hiord,e,t);
    if (rr) {
	p->e2 = getedge (obs, ob, lslit(ob->slit),
	    spoverlap, p->center.x, hiord, moepb(hiord,e,t));
    } else {
	p->e2 = getedge (obs, ob, rslit(ob->slit),
	    spoverlap, p->center.x, hiord, moepb(hiord,e,t));
    }
/* NOTE -- in the above we may want to change p->center to cpos for the
  center position, or choose another cpos in a smaller order than
  the current hiord, especially if cordr >= hiord, which could be true */


/* The extent of the blue (2) edge is insufficient, but its polynomial
  coeficients are presumed correct.  So we extend it to match the
  height of the red (1) edge. */
/* The position of p->e2.r and p->e2.u need to be extended in y to
  equal p->e1.r.y and p->e1.u.y; their x values need to be extended
  by the polynomial x-x0 = (y-y0)*(a*(y-y0) + b) where the coefficients
  a and b are p->e2.a and p->e2.b.  */

    p->e2.r.y = p->e1.r.y;
    d = p->e2.r.y - p->e2.p.y;
    p->e2.r.x = p->e2.p.x + d * (d * p->e2.a + p->e2.b);
    p->e2.u.y = p->e1.u.y;
    d = p->e2.u.y - p->e2.p.y;
    p->e2.u.x = p->e2.p.x + d * (d * p->e2.a + p->e2.b);

/* ALL THE PROCESS FOR MOE SPECTRUM GOES HERE...
   check slit length and shorten if needed
   find order range
   find center right and left ends
   find top and bottom for low order side
   find top and bottom for high order side and extrapolate
   to the same top and bottom as low order side.
   need to duplicate all the getedge stuff here, too.
   find the bounding box -- use bbspec
.. */

    obs->cw = obcw;		// Restore it
//    obs->sw = obcw;		// Restore it
/* Now finish up with the bounding box things */
//    p->bb = boundor (bbedge(p->e1), bbedge(p->e2));
    p->bb = bbspec (p);
    p->on_det = (bbover (p->bb, obs->Detect->bb) & 3);
    p->bd = p->bb;

/* Set vertical extent of bd separately.  Use y values derived
  from the rr edge of slit and the secondary wavelength range.
  If pb2 is not active, we set a default size based on detector
  size; allowing for object placement over middle 2/3 of height
  of the detector, or limit y to +/- 1/6 detector total size. */

    if (obs->pb2act) {
/* Get positions for slit images at secondary wavelength range */
	if (rr) {
	    ss = get3de2v (sum2vect(ob->smpos,rslit(ob->slit)), 0.0);
	    dpr = det_pos (obs, obs->pb2.red, ss, moeorder(obs->pb2.red,e,t));
	    dpl = det_pos (obs, obs->pb2.blue,ss, moeorder(obs->pb2.blue,e,t));
	} else {
	    ss = get3de2v (sum2vect(ob->smpos,lslit(ob->slit)), 0.0);
	    dpr = det_pos (obs, obs->pb2.red, ss, moeorder(obs->pb2.red,e,t));
	    dpl = det_pos (obs, obs->pb2.blue,ss, moeorder(obs->pb2.blue,e,t));
	}

/* Set secondary bounding box y values based on those y values */
	if (dpr.y > dpl.y) {
	    p->bd.y.lo = dpl.y;
	    p->bd.y.hi = dpr.y;
	} else {
	    p->bd.y.lo = dpr.y;
	    p->bd.y.hi = dpl.y;
	}

    } else {
/* Restrict y values of bd to total length of 1/3 range of the
  obs->Detect.bb range.  If that range is exceeded, decrease each
  end equally to get within the range.  */
double	dy, by, y;
	dy = (obs->Detect->bb.y.hi - obs->Detect->bb.y.lo)/3.0;
	by = p->bd.y.hi - p->bd.y.lo;
	y = (by - dy) / 2.0;	// Amount to move each end if positive
	if (y > 0.0) {
	    p->bd.y.lo += y;
	    p->bd.y.hi -= y;
	}
    }

    return  p;
}

	/* - - - - - - - - - - *\
	|			|
	|  Conflict  Geometry	|
	|			|
	\* - - - - - - - - - - */


/*   ===   Conflict geometry   ===   */


static	int	spect_edge_intersect (vect2 e1, vect2 e2, double ed,
		sp_edge p,  double xd, int flag)
/* Compute score of intersection points for line e1-e2 with active
y side in the ed direction; to spectral edge p <was... parabola about
 point p with coefs. a and b> and active x side in the xd direction.  */
{
double	m, v, x, y;
int	kl, kp;	/* Scores for line and parabola */
int	kx, ky;

/* In a loop, iteratively find x and y of intersection, starting with
  y at line, and x on parabola at that y.  For each x, y find the
  interval location on the line for x and parabola for y.  Continue
  until two successive y and x values are the same.  The correct the
  scores for active side of opposite line.  Sum scores to return. */

    m = (e2.y - e1.y) / (e2.x - e1.x);
    for (kx=ky=-1, kl=kp=-2, y=e1.y;
	(kx != kl) || (ky != kp);  )  {
	v = y - p.p.y;
	x = p.p.x + v*(p.a*v + p.b);
	kx = kl;
	kl = within (e1.x, x, e2.x);
	y = e1.y + m*(x - e1.x);
	ky = kp;
/*	kp = within (p.r.y, y, p.u.y);  */
/* -NOTE- confused here -- p is really a position, not a slit edge;
 yet we have formulated a slit edge use for it.  We need the y extent
 and it is not supplied in the calling parameters...
 -- So, we are changing p,a,b in call sequence to sp_edge p! */
	kp = within (p.r.y, y, p.u.y);
    }

/* Correct the values for active side of other line... */
/* If going from r to b or e1 to e2 is not the same direction as
  going from the other line to the reference point, we invert */
    if ( ((e1.x - e2.x)  * (p.p.x - xd)) < 0.0) kl = 2 - kl;
    if ( ((p.r.y - p.u.y) * (e1.y - ed)) < 0.0) kp = 2 - kp;
/* ==NOTE== Check carefully the signs above */
/* Change score in 0-1-2 space to 1-2-0 space here */
    kl = (kl + 1) % 3;
    kp = (kp + 1) % 3;
//if (flag) {
//double c, h;
//	printf ("   From %.2f %.2f to %.2f %.2f active in %.1f\n",
//		e1.x, e1.y,  e2.x, e2.y,  ed);
//	printf ("  Other %.2f %.2f to %.2f %.2f active in %.1f\n",
//		p.r.x, p.r.y,  p.u.x, p.u.y, xd);
//	printf ("  Cross %.2f %.2f, center %.2f %.2f to %.2f\n",
//		x, y,  p.p.x, p.p.y, xd);
//	printf ("  Edge_intersect scores l %d, p %d.\n", kl, kp);
//// Compute x at r,p,u values and write with actual values... */
//	c = p.r.y - p.p.y;
//	h = p.p.x + c*(p.a*c + p.b);
//	printf ("  Red x comp: %.2f vs. %.2f;", h, p.r.x);
//	c = p.u.y - p.p.y;
//	h = p.p.x + c*(p.a*c + p.b);
//	printf (" Blue x comp: %.2f vs. %.2f\n", h, p.u.x);
//}
    if (flag) fflush (stdout);	// fake use to placate lint
    return  (kl + kp);
}


static	boolean	spect_intersect (vect2 e1, vect2 e2, vect2 en, spect *b,
		 int flag)
/* Score intersection of line segment e1-e2 with both sides of
  other spectrum b, return true if intersection detected. */
/* flag is added for DEBUG print only. */
{
int	k1, k2;

/* Score intersection with edge 1 of b */
    k1 = spect_edge_intersect (e1, e2, en.y,
	b->e1, b->e2.p.x, flag);
    if (k1 > 3) return 1;
//    if (flag && k1 > 3) {	// debug version;
//	printf ("  spect_intersect k1 = %d, intersection.\n", k1);
//	return  1;
//    }

/* Score intersection with edge 2 of b */
// was..    k1 = spect_edge_intersect (e1, e2, en.y,
// and...	b->e2, b->e1.p.x);
// suspect that this should be "k2" not "k1" here!
// also, the final 2 parameters are probably wrong...
    k2 = spect_edge_intersect (e1, e2, en.y,
	b->e2, b->e1.p.x, flag);
    if (k2 > 3) return 1;
//    if (flag && k2 > 3) {	// debug version;
//	printf ("  spect_intersect k2 = %d, intersection.\n", k2);
//	return  1;
//    }
/* DEBUG note -- we note that the only values k1 or k2 can take
  are 0 through 4... */
//    if (flag) printf (" spect_intersect k1 = %d, k2 = %d\n", k1, k2);

/* See if total indicates intersection */
    return ( (k1+k2) > 5 );
}


static	boolean	spect_overlap (spect *a, spect *b, int flag)
/* Determine if the two given spectral structures overlap on the
  detector space; true if conflicted.  */
/* flag is added for DEBUG only */
{

/* NOTE -- WE HAVE TO ADD OVERLAP CRITERIA.  FIND AMOUNT OF SPECTRAL
  OVERLAP PERPENDICULAR TO DISPERSION, AND SEE THAT IT IS LESS THAN
  THE OVERLAP CRITERIA SPECIFIED IN THE OBS DATA! */

/* Look for bounding-box overlaps.  If no overlap in either
  direction, we may safely return false.  */
    if (ovr_lap(a->bb.x, b->bb.x) == 0) return 0;
    if (ovr_lap(a->bb.y, b->bb.y) == 0) return 0;
// We don't use bbover above for increased speed...

// if (flag) printf (" spect_overlap past bound test (does overlap).\n");

/* Look in detail at the ends of the spectra which lie within
  the bounding box of the other spectrum in y.  For each such,
  test in detail the intersection.  If a conflict is found,
  immediately return true.  If not found, finally return false. */
    if (overlap (a->e1.r.y, a->e2.r.y, b->bb.y.lo, b->bb.y.hi)) {
// if (flag) printf (" spect_overlap pri red overlaps sec box.\n");
	if (spect_intersect (a->e1.r, a->e2.r, a->e1.u, b, flag))
	return 1;
    }

    if (overlap (a->e1.u.y, a->e2.u.y, b->bb.y.lo, b->bb.y.hi)) {
// if (flag) printf (" spect_overlap pri blue overlaps sec box.\n");
	if (spect_intersect (a->e1.u, a->e2.u, a->e1.r, b, flag))
	return 1;
    }

// if (flag) printf (" spect_overlap location 2.\n");

    if (overlap (b->e1.r.y, b->e2.r.y, a->bb.y.lo, a->bb.y.hi)) {
// if (flag) printf (" spect_overlap sec red overlaps pri box.\n");
	if (spect_intersect (b->e1.r, b->e2.r, b->e1.u, a, flag))
	return 1;
    }

    if (overlap (b->e1.u.y, b->e2.u.y, a->bb.y.lo, a->bb.y.hi)) {
// if (flag) printf (" spect_overlap sec blue overlaps pri box.\n");
	if (spect_intersect (b->e1.u, b->e2.u, b->e1.r, a, flag))
	return 1;
    }

// if (flag) printf (" spect_overlap location 4.\n");

/* Passed all the required tests, so no conflict... */
    return  0;
}


	/* - - - - - - - - - - *\
	|			|
	|  Object data filling	|
	|			|
	\* - - - - - - - - - - */

/*  ===  Fill in object data from optics data  ===  */

static	int	fillobject (objq* ob, Obs* obs)
/* Fill data for object ob.  Get an objval structure.  Compute
  the mask position.  Note that object reading has already
  filled in the slit values, objdat structure, type and flag;
  and put ob into the object queue.  This program is called
  in the first queue traversal to find all spectra. */
/* NOTE - could return 0 if not active, 1 if active */
{
int	ord;
spect	*sp;
// double	scale;
int	ul;
int	x;
char	typnam[16];

/* Note - ob is not queue here, but individual object */
    if (ob == NULL) return 0;

/* Set the object active status here, active unless inhibited by
  the re-use count or by reference selection criteria. */

/* See if this object should be active -- test use count. */
    if (ob->type == OBJ_OBJECT) ul = obs->reuseob;
    else if (ob->type == OBJ_REFERENCE) ul = obs->reuserf;
    else  ul = 999;
    if (ob->use <= ul) ob->flag |= OBJECT_ACTIVE;

/* If a reference object and selection is active, de-select those
objects which have been selected and should not be... */
    if ( (ob->type == OBJ_REFERENCE) && (obs->refsel > 1) &&
	    (ob->flag & OBJECT_ACTIVE) ) {
	if ((refcount++ % obs->refsel) != 0) {
//	    ob->flag &= !OBJECT_ACTIVE;
	    ob->flag &= ~OBJECT_ACTIVE;
	}
    }

/* Obtain a type name for this object to use in the messages */
    if (ob->type == OBJ_OBJECT) {
	sprintf (typnam, "Object");
    } else if (ob->type == OBJ_REFERENCE) {
	sprintf (typnam, "Align.");
    } else {
	sprintf (typnam, "Unknown");
    }


/* Compute smpos as a vector from telescope data, and obs */
/* To do this computation, we derive a position matrix for the
  field center position.  Then, using celestial position for the
  object, convert that into a vector.  Project that vector onto
  the tangent plane using the position matrix.  Finally, determine
  the z-value using focal plane curvature. */

    ob->smpos = maskvect (obs, ob->dat.ra, ob->dat.dec);
// Above computes mask position using OPTICS
// Also, sets special position if far from center position.

/* ==NOTE== assuming that ra, dec have been precessed to the equinox
  of the obs structure by the object reading code! */
/*  NOTE -- we (yet) do not do any precession, assuming that the input
  star lists are all in the equinox specified in the .obs file.  This
  is not completely realistic, although a utility to precess a star list
  file to a given equinox is easily constructed. */
// The maskvect subroutine uses transform matrix and opt transform code.
// NOTE -- check that instrument position angle is included in the
// transform matrix which is used there...

/* Check for object FAR off the mask or detector.  If so, set it to
  non-active status, and skip further images as it will not be used. */
    if (vect2norm(ob->smpos) > 500.0) {		// NOTE ARBITRARY LENGTH
	if (Nfaroff < MESCOUNT && obs->dlevel > 0)
		printf (" -%s %s FAR off mask.\n", typnam, ob->dat.name);
//	ob->flag &= !OBJECT_ACTIVE;
	ob->flag &= ~OBJECT_ACTIVE;
	Nfaroff++;
    }

/* MOVE TO THIS POINT CODE TO DETECT:
	OBJECT FAR OFF MASK, AND SET IT TO INACTIVE
	OBJECT OFF MASK AT ALL, OR IN AVOID REGION, AND SET INACTIVE
	OBJECT EXCEEDING USE COUNT, AND SET INACTIVE
  THEN, ALL OBJECTS SET INACTIVE SHOULD NOT GET OBJVAL OR SOD STORAGE.
WE MAY HAVE TO WORRY ABOUT A NULL OBJVAL OR SOD POINTER BEING REFERENCED.
 .. */

/* Temporarily, null out any object at all over 15 arc minutes from
  the center; with no ADC this is 309.44 mm, with ADC use 310.18 mm */
//    if (vect2norm(ob->smpos) > 330.0) {  // was.. 309.44) {
// DUE TO CONFLICT WITH AVOIDANCES IN MASKCUT, WE TEMPORARILY CLEAN
// ALL OUTSIDE 302 MM.  LATER, WE WILL ALLOW A LARGER RADIUS (ABOVE),
// AND SEARCH OUT ANY IN OR NEAR AVOIDANCE ZONES...

/* The radius allowed depends on instrument... */
    if (strcmp(obs->Instr_mt->name, "LDSS3") == 0) {
	if (vect2norm(ob->smpos) > D_maxradl) {
	    if (obs->dlevel > 0)
		printf (" -%s %s is off mask.\n", typnam, ob->dat.name);
	    ob->flag &= ~OBJECT_ACTIVE;
	}
    } else {
	if (vect2norm(ob->smpos) > D_maxradi) {  // was.. 309.44) {
//	if (vect2norm(ob->smpos) > 302.0) {  // was.. 309.44) {
// NOTE: FIX ABOVE TO D_maxradi WHEN AVOIDANCES ARE RATIONALIZED
	    if (obs->dlevel > 0)
		printf (" -%s %s is past 15'.\n", typnam, ob->dat.name);
//	    ob->flag &= !OBJECT_ACTIVE;  /* NEVER USE ! TO INVERT BITS!! */
	    ob->flag &= ~OBJECT_ACTIVE;
	}
    }

/* Temporarily, null out any alignment object which is over 12 arc minutes
  from the center of the mask...   That translates to about 251 mm. */
/* Do this only if the telescope is MagNoADC ... */
    if ( (ob->type == OBJ_REFERENCE) && (vect2norm(ob->smpos) > 251.0)
	&& ( !strcmp(obs->Tel_scope->name, "MagNoADC") )
// Removed Mag2NoADC in the above line
	) {
	if (obs->dlevel > 0)
	    printf (" -%s %s is past 12'.\n", typnam, ob->dat.name);
//	ob->flag &= !OBJECT_ACTIVE;
	ob->flag &= ~OBJECT_ACTIVE;
    }
// THE ABOVE SHOULD BE MOVED TO mask_space PROGRAM IN MASKGEN IF IT
// IS STILL NEEDED.

/* At this point, if we have made the object inactive due to obs count,
  or off mask, we kill any objval stuff and return as there is nothing
  further to do to these unworthy ones... */
    if (!(ob->flag & OBJECT_ACTIVE)) {
	ob->val = kill_objval (ob->val);
	return 0;
    }
// Following this, objects are active, and get spectra computed.

/* Check for existing objval; if none, get dynamic storage
  and create one.   Initialize conflict count and queue pointers
  to null.  Same for spectrum pointers. */
    if (ob->val == NULL) {
/*	ob->val = (objval *) malloc (sizeof(objval)); */
	ob->val = xalloc (objval);
/* Note that this is only done in maskgen, not when reading
  the SMDF file objects... */
	ob->val->ncf = 0;
	ob->val->cfq = NULL;
	ob->val->sod = NULL;
	ob->val->img = NULL;
/* DEBUG section - turn on flag for interesting nodes  */
//	if      (!strcmp(ob->dat.name, "13073_25.4")) ob->val->flag = 1;
//	if      (!strcmp(ob->dat.name, "obj1174")) ob->val->flag = 1;
//	else if (!strcmp(ob->dat.name, "obj1184")) ob->val->flag = 1;
//	else if (!strcmp(ob->dat.name, "obj447")) ob->val->flag = 1;
//	else if (!strcmp(ob->dat.name, "obj448")) ob->val->flag = 1;
//	else
//	else if (!strcmp(ob->dat.name, "13566_20.7")) ob->val->flag = 1;
//	else if (!strcmp(ob->dat.name, "12578_24.1")) ob->val->flag = 1;
//	else

	ob->val->stat[0] = 0;
	ob->val->stat[1] = 0;

	ob->val->flag = 0;
    }


/* If not null, obtain the primary spectrum.  Then, working both
  ways, find those secondary spectra which are on detector.  Push
  these onto the list in ob. */
/* If not found yet, obtain the primary spectrum */
    if (ob->val->sod == NULL) {
	if (obs->Flag & OFmoe) {
	    ob->val->sod = moespect (obs, ob);
	} else {
	    ob->val->sod = getspect (obs, ob, obs->order, obs->pb);
// Computes spectrum positions using OPTICS -- uses getedge and det_pos
// Was setting object active here; now its far up...
	}
    }

    if (ob->val->sod != NULL) {
/* This extra test is probably not needed any longer... */
	if (vect2norm(ob->val->sod->center) > 200.0) {
	    if (obs->dlevel > 0)
		printf (" -%s %s FAR off detector.\n", typnam, ob->dat.name);
//	    ob->flag &= !OBJECT_ACTIVE;
	    ob->flag &= ~OBJECT_ACTIVE;
/* Remove the sod and objval too, and return as nothing more remains. */
	    ob->val = kill_objval (ob->val);
	    return 0;
	}
// NOTE - not quite right for alignment objects!
// Alignment objects do not need to be on the detector...
// See computation of direct bound box for objects below...
    }

/* If not found yet, push all possible images onto queue */
/* Skip this section if ob->ex_order is not on; this is a selectable
  feature set in the observation file to de-conflict other orders. */
/* NOTE that if disperser is "direct" and not dispersing, we should NOT
  do this loop.  Detect that and bypass these by returning here. */
    if (ob->val == NULL) return (ob->flag & OBJECT_ACTIVE) ? 1 : 0;

/* Planned enhancement:
  When implemented by intgui, we will have 2 flags in ex_order, 1
  for the OBJ_REFERENCE and 2 for others, 3 means both.  Check the
  correct bit with the object type to see if we have enabled extra order
  images for that object type.  */
#ifdef  EXORDX
    if (ob->type == OBJ_OBJECT) x = 2;
    else if (ob->type == OBJ_REFERENCE) x = 1;
    else  x = 0;
#else
    x = 3;	// Compatable with old methods
#endif
/* x is now the reference/object bbit matching ex_order mask... */
    if (ob->val->img == NULL && (ob->flag & OBJECT_ACTIVE)
	&& (obs->ex_order & x) && !(obs->Flag & OFmoe) ) {
/* NOTE that we are NOT getting alignment extra-order spectra... */
/* Do this each way; first going up... */
	for (ord=obs->order + 1; ord < XordLim; ord++) {
	    sp = getspect (obs, ob, ord, obs->pb);
	    if (sp->on_det) {
		sp->next = ob->val->img;
		ob->val->img = sp;
	    } else {
		sp = kill_spect (sp);
		break;
	    }
	}
/* Now, going down... */
	for (ord=obs->order - 1; ord > -XordLim; ord--) {
	    sp = getspect (obs, ob, ord, obs->pb);
	    if (sp->on_det) {
		sp->next = ob->val->img;
		ob->val->img = sp;
	    } else {
		sp = kill_spect (sp);
		break;
	    }
	}
    }	/* End of image finding operation */

/* Compute here the direct image bounding box, placing it into the
  dimg element of the objval structure.  This will be a direct image
  for the alignment object (or the slit) and will be compared to the
  detector bound for alignments in find single conflicts later. */
    ob->val->dimg = direct_loc (obs, ob);
// ==direct==

    return (ob->flag & OBJECT_ACTIVE) ? 1 : 0;
}	// end of fillobject routine


int	object_setup (Obs* b)
/* Fill in data on all objects in the queue, and return a count
  of the objects so filled. */
{
objq	*a;	/* Head of object queue */
int	n = 0;
objq	*r;	/* scratch */

// ===debug section
// /* checl grangle... */		// debug debug
//     if (GRangle == NULL) {
// 	printf (" NULL GRangle pointer found. \n");
//     }
//     printf (" GRangle value = %.5f\n", *GRangle/Degree);
// ===debug section

/* Check for null parameter */
    if (b == NULL) {
	printf (" ** Observation data missing.\n");
	n++;
    }
    a = b->ob;

/* Check for a null object queue */
    if (a == NULL) {
	printf (" ** Object queue missing.\n");
	n++;
    }

/* If any errors were detected, return a zero. */
    if (n > 0) return  0;

#if (DEBUG_UTL > 0)
    printf (" ..Setting up objects for %.8s.\n", b->fname);
#endif

/* Check for any missing components... */
    if (b->Detect == NULL) {
	printf (" ** Observation data missing Detector data.\n");
	n++;
    }

/* If any errors were detected, return a zero. */
    if (n > 0) return  0;

// Need to make direct version of current instrument...
//    b->Direct_inst = Dinstr (b);
// above is not needed; use of dir_loc should reference
// Dinstr if the object pointer is null the first time, and
// use the pointer on subsequent calls.
// ==direct==

    C_loop (r, a) n += fillobject (r, b);
/* NOTE could (should?) have fillobject return 1 only if a good, active
  object is being filled.  Then n += fillobject() would count the number
  of active objects created here. */
    return  n;
}


	/* - - - - - - - - - - *\
	|			|
	|  GAP virtual objects	|
	|			|
	\* - - - - - - - - - - */

/*  ====  Generating virtual objects for gaps  ====  */

// /* A debug subroutine to trace gap existance */
// #ifdef  GAPVO
// void	Gap_count (int k, objq* head) {
// /* Count and report number of gap objects in queue */
// objq	*p;
// int	n=0;
// int	a=0;
//     C_loop (p, head) {
// 	if (p->type == OBJ_GAP) {
// 	    n++;
// 	    if (p->flag & OBJECT_ACTIVE) a++;
// 	}
//     }
//
// /* Report number, and number active */
//     printf (" At %d, there are %d Gap objects, %d are active.\n", k, n, a);
// }
// #endif

#ifdef  GAPVO
static	objq*	gap_queue (Obs* ob)  {
objq	*head = NULL;	// head queue to return
objq	*b;		// object pointer
gapdat	*p;
element	*hg;
sp_edge	edg;
spect	*sd;
objval	*v;
double	toright = 0.0;
double	tored = 0.0;
double	xmid;
double	x1, x2;
//double	top, bot;	// detector bounds
//intval	xint;
intval	yint;
vect2	de;

/* Need the gap support in Dstat, and need gap data from opticdef too */
    if (!(ob->Dstat & DS_gaps)) return NULL;
    hg = find_element (ob->elist, GAPHEAD);
    if (hg == NULL) return NULL;		// Gap header element

//    printf (" > Looking for Gap objects for %s.\n",
//	ob->Instr_mt->name);					// debug only

/* Save the detector Y-interval for bounds and edge positions */
    yint = ob->Detect->bb.y;

/* Find the directions of the object slits, both left/right and
  in the blue/red direction.  Look through objects until a real
  object is found, and then find positions in its spectrum.  */
    C_loop (b, ob->ob) {
	if (b->type != OBJ_OBJECT) continue;	// Find an object
	v = b->val;
	if (v == NULL) continue;		// with objval structure
	sd = v->sod;
	if (sd == NULL) continue;		// and spectrum
	toright = sd->e2.p.x - sd->e1.p.x;	// Get x difference
	tored = sd->e1.r.y - sd->e1.u.y;	// And y difference
	break;					// Split when found.
    }	// looking for an example object

/* Scan the gap queue for each useful gap. */
    C_loop (p, hg->data) {
//	printf (" > Obj. %s, Inst. %s ?\n",
//		p->name, p->instr);				// debug only
//	if (strcmp(p->instr, ob->Instr_mt->name)) continue;
	if (strncasecmp (p->instr, ob->Instr_mt->name,
	    min (strlen(p->instr), strlen(ob->Instr_mt->name)) ) ) continue;
//	printf (" > Dstat %d, object->ns %d.\n",
//		ob->Dstat & DS_ns, p->ns);			// debug only
	if ((ob->Dstat & DS_ns) != p->ns) continue;
/* If we pass those tests, we have instrument and ns orientation matched. */

/* Create an object (b) for gap structure (p) here. */
	b = zalloc (objq);	// Allocate and clear an objq structure

/* Fill in the elements of the objq structure... */
	b->next = b->last = b;
	b->dat.name = dynamstr (p->name);
//	printf (" > Constructing gap %s.\n", b->dat.name);	// debug only
	b->dat.priority = -9.99e99;
/* Will fill in ra, dec a little later, also smpos */

	b->slit.shape = RECTANGLE;
/* Leave rest of the slit null, and construct an objval structure. */
	v = zalloc (objval);
	b->val = v;

/* Construct an objval structure, mostly the sod element */
	sd = zalloc (spect);
	sd->order = ob->order;
	v->sod = sd;

/* Get the elements of the sd structure from the gap position on
  the detector, and the signs for direction. */
	x1 = 0.5 * (p->x1u + p->x1b);
	x2 = 0.5 * (p->x2u + p->x2b);
	xmid = 0.5 * ( x1 + x2 );
	sd->center = make2vect (xmid, 0.0);

/* Make an edge for gap side 1 */
	edg.p = make2vect (x1, 0.0);
	if (tored < 0.0) {
	    edg.r = make2vect (p->x1b, yint.lo);
	    edg.u = make2vect (p->x1u, yint.hi);
	} else {
	    edg.r = make2vect (p->x1u, yint.hi);
	    edg.u = make2vect (p->x1b, yint.lo);
	}

/* Come up with the edge curve coeficients, this for a line. */
	edg.a = 0.0;	// Linear curve...
	de = sub2vect (edg.r, edg.u);
	edg.b = de.x / de.y;
//	sd->bb = bbedge (edg);


/* See which sod edge to assign, depending on whether the sign of
 gap edges 2-1 is the same as spectra edges 2-1 */
	toright *= (x2 - x1);	// only the sign counts
	if (toright < 0.0) {
	    sd->e2 = edg;
	} else {
	    sd->e1 = edg;
	}

/* Now, make and assign the edge of gap side 2 */
	edg.p = make2vect (x2, 0.0);
	if (tored < 0.0) {
	    edg.r = make2vect (p->x2b, yint.lo);
	    edg.u = make2vect (p->x2u, yint.hi);
	} else {
	    edg.r = make2vect (p->x2u, yint.hi);
	    edg.u = make2vect (p->x2b, yint.lo);
	}
	if (toright < 0.0) {
	    sd->e1 = edg;
	} else {
	    sd->e2 = edg;
	}
	edg.a = 0.0;	// Linear curve...
	de = sub2vect (edg.r, edg.u);
	edg.b = de.x / de.y;
//	sd->bb = boundor (sd->bb, edg);

// /* Make up the bounding boxes for the gap.  First, get the x interval
//   which includes all the gap x positions, and then combine that with
//   the y interval we have been using for detector extent.  */
// 	xint = interval (p.x1u, p.x2u);
// 	xint = intex (xint, p.x1b);
// 	xint = intex (xint, p.x2b);
// 	sd.bb = boundbox (xint, yint);
// // use bbedge and boundedge and boundor as in make spect...
// // use boundor and bbedge as in getspect
// // or use bbspec
// //    p->bb = boundor (bbedge(p->e1), bbedge(p->e2));

/* Make up bounding boxes for the gap, use bbspec */
	sd->bb = bbspec (sd);
	sd->bd = sd->bb;

/* That completes the v->sod element, and others are null.  So val is
  done, and we continue with filling the object structure. */
	b->type = OBJ_GAP;
	b->flag = OBJECT_ACTIVE | OBJECT_VIRTUAL;

/* Obtain the position for this fake object using the center position
  which was computed above in the sd element.  Back project this from
  the detector to the mask, and then from the mask to the sky.  */
	b->smpos = pos_det (ob, ob->sw, sd->center, ob->order);
	unmaskvect (ob, b->smpos, &(b->dat.ra), &(b->dat.dec) );

/* Finally, put this object into the local queue which we will return. */
	head = push_obj (b, head);
//	printf (" > Finished %d gap %s.\n",
//	    listcount(b), b->dat.name);   			// debug only
    }	// loop over all possible gaps
//    printf (" ..Returning %d gap objects\n", listcount(head));	// debug only
    return  head;
}

void	add_gaps (Obs* ob)  {
/* If the gap feature is on, add gaps to the object queue */
objq	*g;

/* Leave quietly if there is no gap request */
    if (!(ob->Dstat & DS_gaps)) return;

/* This should be the only call to the program above.   Get the
  gap object queue, and put it at the end of the regular queue. */
    g = gap_queue (ob);
//    printf (" ..Adding %d gap objects\n", listcount(g));	// debug only
    ob->ob = push_obj (g, ob->ob);
}
#endif	// on GAPVO


	/* - - - - - - - - - - *\
	|			|
	|    Sort  Objects	|
	|			|
	\* - - - - - - - - - - */

/*  ====  Sorting the object queue  ====  */


/* Sort the object queue by position on detector across dispersion */

int	sort_decide (objq *a, objq *b)
/* Return -1, 0 or 1 depending on a>b, a=b, a<b */
{
int	n;
int	ua, ub;				// Object is used

/* If either is NULL, return 0 */
    if (a == NULL) return 0;
    if (b == NULL) return 0;

/* Put all unused (inactive) objects after all used (active) ones */
/* This is needed to avoid some bad interactions in conflict finding
  and resolution which can miss conflicts and leave them on mask. */
// This may be omitted to locate these interactions...
// /*
    ua = a->flag & OBJECT_ACTIVE;
    ub = b->flag & OBJECT_ACTIVE;
    if (ua > ub) return  1;
    if (ub > ua) return -1;
// */

/* Sort by Bx if active, but by position if not active */
    n = a->val == NULL || b->val == NULL;
    if (n == 0) n = a->val->sod == NULL || b->val->sod == NULL;
    if (n) {
	if (a->smpos.x < b->smpos.x) return  1;
	if (a->smpos.x > b->smpos.x) return -1;
// Note we test by interchanging sense of smpos due to image inversion
    } else {
	if (a->val->sod->center.x < b->val->sod->center.x) return  1;
	if (a->val->sod->center.x > b->val->sod->center.x) return -1;
    }

    return 0;
}

int	sort_check (objq* a, int(decision)(objq*,objq*))
{
/* Check for proper sort order among active nodes, and report the
 number of errors found.  */
objq	*p, *q;
int	k;
    for (k=0,p=a; p != NULL; ) {		// loop on p
	if (p->flag & OBJECT_ACTIVE) {		// need p active
	    for (q=p->next; q != NULL; q=q->next) {	// loop on q
		if (q == a) break;			// stop if wraparound
// BAD CONSTRUCT BREAK (no)
		if (q->flag & OBJECT_ACTIVE) {		// need q active
		    if (decision (p, q) < 0) k++;	// the meat
		}		// end q decide
	    }			// end q loop
	    p = q;		// p is next active node (or a)
	} else p = p->next;	// p not active, keep scanning
	if (p == a) break;	// test at end of loop
    }				// end loop on p
    return  k;
}

objq*	sort_merge (objq *a, objq *b, int(decision)(objq*,objq*))
/* Return pointer to head of merged queue with all nodes
in sort order. */
/* ==NOTE== the decision function is implicitly global here */
{
objq	*p;
objq	*q;

/* Binary merge operation:
  Start a new queue at a null head.
  Loop while both a, b are not null:
	Test element at head of a and b queues
	If a comes first, pop that element from a and
	push it onto the result queue.
	Otherwise, b comes first, pop that element and
	push it onto the result.
  The loop exited when a or b became null.  Push each of these
  onto the result queue, a null being pushed is a null operation.

.. */
    p = NULL;
    while ( (a != NULL) && (b != NULL) ) {
/* Binary merge, decide which one is pulled */
	if (decision (a, b) < 0)	b = pop_obj (q=b);
	else				a = pop_obj (q=a);
	p = push_obj (q, p);
    }
    p = push_obj (a, p);
    p = push_obj (b, p);
    return  p;
}


/* Push and pop queue routines */
sortq*	push_sq (sortq* a, sortq* b)
/* Push a before b, and return new b (head) which may change */
{
sortq	*p, *r;
    if (b == NULL) return a;
    if (a == NULL) return b;
    p = b->last;
    if (p == NULL) p = b;	/* ...p - b... */
    r = a->last;
    if (r == NULL) r = a;	/* ...r - a... */
    a->last = p;
    p->next = a;
    r->next = b;
    b->last = r;	/* ...p - a...r - b... */
    return  b;
}


sortq*	pop_sq (sortq* a)
/* Remove structure from queue, return pointer to next */
{
sortq	*p, *q;
    if (a == NULL) return a;
    p = a->next;
    if (p == a) p = NULL;
    q = a->last;
    if (q == a) q = NULL;
    if (p == NULL) p = q;
    if (p != NULL) p->last = q;
    if (q != NULL) q->next = p;
    a->last = a->next = a;
    return p;
}

sortq*	kill_sq (sortq* a)
/* Remove structure from queue, free it, return pointer to next */
{
sortq	*p;
    p = pop_sq (a);
    if (a != NULL) free (a);
    return  p;
}



objq*	sort_obj (objq *head, int(decision)(objq*,objq*))
/* Return pointer to head of sorted queue */
{
sortq	*sqh, *sqa, *sqb;
objq	*cur;

/* Allocate meta-queue heads in a list, each pointing to a
sort-order queue.  Split the queues between these, probably
with NULL split ends. */

/* How to start:
  If the object queue is null, return null.

  Create outer loop - initialize by obtaining a current object by
  popping it from the object queue, and set sort queue head null.
  Outer loop (while current object not null):
    Obtain a new sort queue element
    Do the inner loop as follows--
    Inner loop (forever):
	Put the current object into the sort object queue.
	Obtain next object from object queue; if null, break.
	Test current object w.r.t. last item in sort object q;
	if correct order, push into queue, if not, break.
    Push the sort queue element into the sort queue.
  -- end of outer loop.
  Result is the sort queue.

.. */

    sqh = NULL;
    head = pop_obj (cur=head);
    while (cur != NULL) {
/*	sqa = (sortq*) malloc (sizeof(sortq)); */
	sqa = xalloc (sortq);
	sqa->q = NULL;
	sqa->next = sqa->last = NULL;
//	while (1) {		// Should use for (ever), #deine ever ;;
	for (ever) {
	    sqa->q = push_obj (cur, (objq*)(sqa->q));
// NOTE (objq*)sqa->q also works, at least in solaris.
	    head = pop_obj(cur=head);
	    if (cur == NULL) break;
//	    if (sort_decide (((objq*)(sqa->q))->last, cur) < 0) break;
	    if (decision (((objq*)(sqa->q))->last, cur) < 0) break;
	}
	sqh = push_sq (sqa, sqh);
    }

/* The meta-queue list is circular.  Traverse it, popping off
two items at a time, merging them, and pushing the result back onto
the meta-queue.  When second pop gives NULL, the first one
is the result, which is returned. */
    if (sqh == NULL) return NULL;
    while (sqh != NULL) {
	sqh = pop_sq (sqa=sqh);
	if (sqh == NULL) break;
	sqh = pop_sq (sqb=sqh);
	sqa->q = sort_merge ((objq*)(sqa->q), (objq*)(sqb->q), decision);
	sqb = kill_sq (sqb);
	sqh = push_sq (sqa, sqh);
    }
    cur = (objq*)(sqa->q);
// NOTE  (objq*)sqa->q; also works, at least in solaris
    kill_sq (sqa);

    return  cur;
//    return  sqa->q;
}


	/* - - - - - - - - - - - - - - *\
	|				|
	|   Conflict Queue Management	|
	|				|
	\* - - - - - - - - - - - - - - */


/*  ===  Conflict queue management  === */


/* Push and pop queue routines */
cfl*	push_cfl (cfl* a, cfl* b)
/* Push a before b, and return new b (head) which may change */
{
cfl	*p, *r;
    if (b == NULL) return a;
    if (a == NULL) return b;
    p = b->last;
    if (p == NULL) p = b;	/* ...p - b... */
    r = a->last;
    if (r == NULL) r = a;	/* ...r - a... */
    a->last = p;
    p->next = a;
    r->next = b;
    b->last = r;	/* ...p - a...r - b... */
    return  b;
}

/* ==NOTE== may have to keep cfl list head as a global, and check for
it in both push and pop operations so a header is always available */

cfl*	pop_cfl (cfl* a)
/* Remove structure from queue, return pointer to next */
{
cfl	*p, *q;
cfl	*r, *s;
cfl	**pri;
    if (a == NULL) return a;
    p = a->next;
    if (p == a) p = NULL;
    q = a->last;
    if (q == a) q = NULL;	/* ...q - a - p... */
    if (p == NULL) p = q;
    if (p != NULL) p->last = q;
    if (q != NULL) q->next = p;
    a->last = a->next = a;
//    a->last = a->next = NULL;	// debug try -- no change
/* Do the same for the secondary per-node links... */
    r = a->fwd;
    if (r == a) r = NULL;
    s = a->bak;
    if (s == a) s = NULL;	/* ...s - a - r... */
    if (r != NULL) r->bak = s;
    if (s != NULL) s->fwd = r;
    a->bak = a->fwd = NULL;	/* non-circular */
/* Find our primary per-node pointer, and maintain it... */
    if (r == NULL) r = s;	/* Get a better pointer */
    pri = &(a->obj->val->cfq);
    if (*pri == a) *pri = r;
    return p;
}

cfl*	spop_cfl (cfl* a)
/* Remove structure from queue, returning next per-node pointer */
{
cfl	*p, *q;
cfl	*r, *s;
cfl	**pri;
    if (a == NULL) return a;
    p = a->next;
    if (p == a) p = NULL;
    q = a->last;
    if (q == a) q = NULL;	/* ...q - a - p... */
    if (p == NULL) p = q;
    if (p != NULL) p->last = q;
    if (q != NULL) q->next = p;
    a->last = a->next = a;
//    a->last = a->next = NULL;	// debug try -- no change
/* Do the same for the secondary per-node links... */
    r = a->fwd;
    if (r == a) r = NULL;
    s = a->bak;
    if (s == a) s = NULL;	/* ...s - a - r... */
    if (r != NULL) r->bak = s;
    if (s != NULL) s->fwd = r;
    a->bak = a->fwd = NULL;	/* non-circular */
/* Find our primary per-node pointer, and maintain it... */
    if (r == NULL) r = s;	/* Get a better pointer */
    pri = &(a->obj->val->cfq);
    if (*pri == a) *pri = r;
/* This is identical to pop_cfl except for the return. */
    return r;
}

cfl*	find_cfl (objq* ob)
/* Find a cfl pointer from the object queue. */
{
objq	*p;
cfl	*r=NULL;
    if (ob == NULL) return NULL;
    C_loop (p, ob) {
	if (p->val != NULL) {
	    r = p->val->cfq;
	    if (r != NULL) break;
	}
    }
    return  r;
}


/*  ===  Add and drop conflicts  ===  */

/*  Add and drop conflicts with objects */

cfl*	add_conflict (objq *cur, objq *other, int type, bool inhibit)
/* Add conflict described to conflict list */
/* Return conflict pointer; also add to queue at head */
/* When called, cur and other have been checked for active, and
  the existance of the objval pointer. */
{
cfl	*p, *q;
cfl	*head;
char	*c1, *c2;

    if (cur == NULL) return NULL;   /* is this needed? (it is now) */
    q = cur->val->cfq;
    if (cur == Tobj) return q;	// Suppress adding for test object
    head = find_cfl (cur);

/* Start by getting the structure up */
    p = xalloc (cfl);
    p->last = p->next = NULL;
    p->type = type;

/* Put into current object's push down list */
//    q = cur->val->cfq;	// done earlier
    if (q != NULL) q->bak = p;
    p->fwd = q;
    cur->val->cfq = p;
    p->bak = NULL;

/* Set object and other pointers */
    p->obj = cur;
    p->other = other;

/* Push into total queue in standard way */
    head = push_cfl (p, head);

/* If not inhibited, call self for secondary */
    if (! inhibit)  {
	q = add_conflict (other, cur, type, 1);
	if (dodebug > 1) {
	    c1 = (other == NULL) ? "NULL" : other->dat.name;
	    c2 = (cur == NULL) ? "NULL" : cur->dat.name;
	    printf (" -- Secondary add_conflict");
	    printf (" %s - %s type %d\n", c1, c2, type);
	}
    }

/* Increment conflict counter in object */
    cur->val->ncf++;
/* DEBUG section - report any flagged nodes */
//     if (cur->val->flag) {
// 	if (other != NULL && other->val->flag) {
// 	    c1 = (other == NULL) ? "NULL" : other->dat.name;
// 	    c2 = (cur == NULL) ? "NULL" : cur->dat.name;
// // fudge the type for conflicts between special nodes
// 	    type += 50;
// 	    p->type = type;
// 	    printf (" ++ %s - %s type %d inhib. %d; count = %d\n",
// 		c2, c1, type, inhibit,  cur->val->ncf);
// 	}
//     }	// end DEBUG
    return  p;
/* ==NOTE== how do we return new head value???  Use **head ?? */
}


cfl*	drop_conflict (cfl *a, bool inhibit)
/* Remove the conflict, return the next per-object conflict */
{
cfl	*p, *q;
objq	*cur;
int	j;
//char	*pri, *sec;	/* Debug */
//char	nn[] = "NULL";	/* Debug */

    if (a == NULL) return NULL;
    cur = a->obj;

/* DEBUG section - check on secondary, and if flagged report
  on the conflict being reduced... */
// Could also look for a->type > 40 here!
//     if (cur->val->flag) {
// objq	*sec;
// char	*c1, *c2;
// 	sec = a->other;
// 	if (sec != NULL && sec->val->flag) {
// 	    c1 = (cur == NULL) ? "NULL" : cur->dat.name;
// 	    c2 = (sec == NULL) ? "NULL" : sec->dat.name;
// 	    printf (" -- %s - %s type %d inhib. %d; count = %d\n",
// 		c2, c1, a->type, inhibit,  cur->val->ncf);
// 	}
//     }	// end DEBUG section

/* If not inhibited, we go remove this conflict from secondary */
    if (! inhibit && a->other != NULL) {
/* ==NOTE== may need nested if to test *a->other->val not null */
/* We search the secondary object's per-object conflict queue to
  find one which contains the current object as its secondary.  If
  found, call ourselves with inhibit to pluck it out. */
	if (a->other->val != NULL) {
	    for (j=0,p=(a->other->val->cfq); p != NULL; j++,p=(p->fwd)) {
		if (p->other == cur) {
/* Found it!  Change the pointer to here to result, thus removing
  this from the queue; and inhibit secondary search */
		    p = drop_conflict (p, 1);
/* If p were equal to a->other->val->cfq, we should also reset
  the value of a->other->val->cfq to the new p.
  This is now done in spop_cfl used by drop_conflict directly, below. */
/* We do it here for good form... */
		    if (j == 0) a->other->val->cfq = p;
		    a->other = NULL;
		    break;
		}
	    }  /* Loop to search secondary */
	}  // Test for non-null a->other->val for safety
    }  /* Test for inhibited or null */

/* Pop from the list */
    q = spop_cfl (a);	/* Allow to pop head of queue */
//    if (*head == a) *head = q;	/* Re-establish the head */

/* Decrement conflict count of current object */
    if (cur->val != NULL) cur->val->ncf--;

/* return next by-object pointer */
//    q = a->q;
//  done by using spop above

/* NOTE - may need debug print of what was dropped? */
//    if (dodebug > 0) {
//	printf (" -drop_conflict %s/%s -%d- ", pri, sec, inhibit);
//	printf (" [%08x]", a);
//	printf ("\n");
//    }

/* also, clean storage of the freed node */
    free (a);
    return  q;
}


	/* - - - - - - - - - - *\
	|			|
	|    Find Conflicts	|
	|			|
	\* - - - - - - - - - - */


/*  ===  Find conflicts  ===  */

double	DeltaX = 2.5;	// A global here...

// The "delta" defined here is the width to which we search across the
// dispersion direction for possible conflicts.  It can depend on the
// camera focal length (see instrument ident), the default slit widths
// and specific slit widths.  Currently it is a global parameter; yet
// it really should be computed rather than defined!

// #define  deltaX  2.5
/* At 125mm/8192 pixels, 9 mm is 590 pixels.  NOTE - we really need to
compute this limit using possible curvature and slit width parameters
projected onto the detector space. */
/* The bounding box vs straight spectrum distance is roughly 0.12 mm
 in IMACS, or about 8 pixels.  The slit size is defined in the obs file,
 and some scale can be defined.  Here in IMACS, the slit is 12 arcsec and
 images to about 1.62 mm.  So, we need a delta of 0.135 * slit size +
 0.12 mm and take a factor of 1.5 additional. */
/* We can compute this later, for now, use 2.5 mm. */
/* This value could be put in the parameter list of find_conflict */

// NOTE -- we should determine deltaX from the width of a
// bounding box on the actual detector, and multiply that by
// about 1.2 to 1.5 for safety; then add 2x slit width.

// ALSO NOTE that if sod->center is removed from debug use, a similar
// location must be put into the objval structure.

//    DeltaX = 1.1 * bigest_bounding_box (head);

/* Program to find widest bounding box on the detector */
static	double	bigest_bounding_box (objq* head)
{
double	b = 0.0;
double	e;
objq	*p;
objval	*v;
spect	*s;
    b = 0.0;	// Do this actively as well as above...
    C_loop (p, head) {
	if (p->flag & OBJECT_ACTIVE) {	// Only for active ones
	    v = p->val;
	    if (v != NULL) {
		s = v->sod;
		if (s != NULL) {
		    e = fabs (s->bb.x.hi - s->bb.x.lo);
//		    if (e > b) b = e;
// debug: reject any overly gross bb values here...
		    if (e > b && e < 1.e3) b = e;
if (e > 1.e3) printf (" *Node %s bbx = %.4g %.4g\n",
	p->dat.name, s->bb.x.lo, s->bb.x.hi );	// debug only
		}		// valid sod
	    }			// valid objval
	}			// active object test
    }
    return  b;
}


//#define NOCONFLRAL
#undef NOCONFLRAL
// Symbol NOCONFLRAL is used to include the code to directly remove
// a reference node which is off the array rather than letting it form
// a conflict which will be warned later.


int	find_conflict (objq *pri, objq *sec, Obs *obs)
/* Find conflict status of given objects, store if found.  Return is
 -1 if Bx out of range, 0 if no conflict, 1 if conflict found. */
{
//bool	on_x, on_y;
int	k;	/* result */
double	w;	/* Slit width */
int	ks;	/* Emergency stop - debug */
cfl	*cf;	/* Conflict added; put in queue, not actually used here */
int	cc = 0;
double	xs, ys, rv;
double	x1, x2, y1, y2;
int	on;
int	good;
//static  int	kmsg=0;
/* ==NOTE== We should use cf to put into chead queue here, rather
  than in the add_conflict code, or else not return cf from add_conflict */
spect	*p;
int	special;	/* DEBUG status */
// char	*c1, *c2;	// debug pointers

    if (pri == NULL) return 0;
    if (sec == NULL) return 0;
    k = 0;

// Omit any finding of a conflict for inactive objects
    if (!(pri->flag & OBJECT_ACTIVE)) return k;

/* If pri and sec are same, we only check for spectrum off detector
  and for ghost images conflicting with primary images.  Actually,
  the latter (ghost self-conflict) should have been checked, as it
  will be a global problem if it exists at all. */
    if (pri->val == NULL) {
	printf (" *!(FC) Pri %s active without objval!\n", pri->dat.name);
	fflush (stdout);
	return  k;
    }
    if (pri->val->sod == NULL) {
	printf (" *!(FC) Pri %s active without sod!\n", pri->dat.name);
	fflush (stdout);
	return  k;
    }

    if (pri == sec) {
/* Check bounding box and detector */

/* Revised to use bounding box code; and secondary bound box */
/* Only use of secondary bounding box here */
	if (pri->type == OBJ_REFERENCE) {
/* Use direct image bb for the reference objects */
	    on = bbover (pri->val->dimg, obs->Detect->bb);
	} else {
	    on = bbover (pri->val->sod->bd, obs->Detect->bb);
	}
// result is 2 if spectrum is included within detector, other values
// are 0 if non-intersect, 1 if intersect, 4 if detector within spectrum
// only value of 2 is acceptable here
//	good = (on == 2);
	good = (on == 2 || (pri->type == OBJ_GAP));
/* We must let gap objects stay good, although overlap may 1 or 4 there. */

/* Check also for vignetting radius (for IMACS short camera) */
/* But do this only for actual objects, not reference */
	if (pri->type == OBJ_OBJECT) {
	    x1 = pri->val->sod->bd.x.lo;
	    x2 = pri->val->sod->bd.x.hi;
	    y1 = pri->val->sod->bd.y.lo;
	    y2 = pri->val->sod->bd.y.hi;

	    xs = max (x1*x1, x2*x2);
	    ys = max (y1*y1, y2*y2);
	    rv = obs->Detect->vrad;
	    good = good && ((rv*rv) > (xs+ys));
	}

/* Also check for slit being within allowable IMACS radius -- this is
  IMACS specific, and we should really test with a mask limit routine... */
/* CHECK THAT MASKRAD HAS BEEN SET PROPERLY */
	good = good && (Xcircle(pri, make2vect(0.0,0.0), D_maskrad) < 0);
// Above was added for slit extension algorithm compatability
/* Changed to use global D_maskrad which is set in normalize routine;
  should be used for all instruments... */

/* Keep gap objects as good... */
	good = good || pri->type == OBJ_GAP;

	if ( (!good) && (pri->type == OBJ_OBJECT) ) {
/* NOTE -- We may want to keep priority < -2.0 objects??? */
/* See if we want objects < Pmusthave to be retained ?? */
/* Only check detector vs. sod for objects, not alignment stars... */
/* Appears spectrum bounding box is outside detector boundary */
	    cf = add_conflict (pri, NULL, 1, 1);
	    if (cf != NULL) cc++;  // nonsense for lint
	    if (obs->dlevel > 0 && Tobj == NULL && Noutdet < MESCOUNT) {
// Add Tobj test to suppress when testing extensions
		printf (" Spect. %s outside (detector/mask) boundary.\n",
		    pri->dat.name);
//		kmsg++;
	    }
	    Noutdet++;
/* 1 = conflict type; 1 = inhibit here */
	    k = 1;
	} else if ( (!good) && (pri->type == OBJ_REFERENCE) ) {
// Currently warn for this, later add it as a conflict, or inactivate.
#ifdef  NOCONFLRAL
	    if (obs->Flag & OFral) {
#endif
/* Removing alignments wrt detector... */
		cf = add_conflict (pri, NULL, 1, 1);
		if (obs->dlevel > 0 && Tobj == NULL) {
		    printf (" Reference %s outside detector/mask.\n",
			pri->dat.name);
/* Reference removal not limited by message count (yet). */
		}
		Noutdet++;
/* In this section, always allow off-detector refernce*/
#ifdef  NOCONFLRAL
	    } else {
		printf (" **WARNING: Reference %s falls off detector!\n",
		    pri->dat.name);
//printf ("-bbx: X %.2f %.2f; Y %.2f %.2f\n",
//	pri->val->dimg.x.lo, pri->val->dimg.x.hi,
//	pri->val->dimg.y.lo, pri->val->dimg.y.hi );	//debug ==direct==
	    }
#endif
/* In above section, we now always allow off-detector reference to have
  a conflict created.  We only remove the item for this conflict if
  automatic removal is allowed (later). */

	}	// Conditional for not on detector
	return  k;
    }	// end pri, sec same case.

// Since sec is different from primary here, check it also
// to see if it is active too.  If not, we don't check further...
    if (!(sec->flag & OBJECT_ACTIVE)) return k;
// COULD be placed after check for Bx range below...
// NO, that's a bad idea, check Bx range only for actives...

    if (sec->val == NULL) {
	printf (" *!(FC) Sec %s active without objval!\n", sec->dat.name);
	fflush (stdout);
	return  k;
    }
    if (sec->val->sod == NULL) {
	printf (" *!(FC) Sec %s active without sod!\n", sec->dat.name);
	fflush (stdout);
	return  k;
    }

/* DEBUG -- set special status */
    special = 0;
//    special = pri->val->flag && sec->val->flag;
//    if (special) {
//	c1 = (pri == NULL) ? "NULL" : pri->dat.name;
//	c2 = (sec == NULL) ? "NULL" : sec->dat.name;
//    }

/* Check for Bx out of range */
/* ==NOTE== need to compute b1 and b2 for each, or a standard
  delta b allowed (globally) to detect this.  We take delta B allowed
  as simply a global value here */
    if (fabs(sec->val->sod->center.x - pri->val->sod->center.x) > DeltaX) {
//	if (special) printf (" xx %s - %s far apart.\n", c1, c2);
	return  -1;
    }


/* Check mask space conflict - finds duplicates too; type 5 */
    w = pri->slit.width;
    if (w < sec->slit.width) w = sec->slit.width;
/* NOTE assuming width is positive for all slits! */
    w /= 2.0;	// Check for really close, under 1/2 slit width here.
/* If separated by under 1/2 largest slit width in both x and y,
  the two slits would intersect on the mask; we don't wish to
  allow that, so a conflict exists on the mask.  */
/* This allows slits which overlap in length; it is not a fully
  rogorous mask intersection computation. */

/* Omit test for duplicates if either is a gap object... */
    if (pri->type != OBJ_GAP && sec->type != OBJ_GAP) {
	if (fabs(sec->smpos.x - pri->smpos.x) < w) {
	    if (fabs(sec->smpos.y - pri->smpos.y) < w) {
		cf = add_conflict (pri, sec, 5, 0);
		if (cf != NULL) cc++;  // nonsense for lint
		if (obs->dlevel > 0 && Nduppos < MESCOUNT) {
#ifdef  PRIORDER
		    if (pri->order < sec->order) {
#else
		    if (pri->order > sec->order) {
#endif
			printf (" Duplicate mask position for %s and %s.\n",
			    pri->dat.name, sec->dat.name);
		    } else {
			printf (" Duplicate mask position for %s and %s.\n",
			    sec->dat.name, pri->dat.name);
		    }
		}
		Nduppos++;
	    }  // test for close in y
	}  // test for close in x
    }  // drop gap objects dup test

/* Check types; if both are reference objects, then we really don't
  care about conflicting spectra, so we just return a zero. */
    if (pri->type == OBJ_REFERENCE && sec->type == OBJ_REFERENCE) {
//	if (special) printf (" xx %s - %s both reference objs.\n", c1, c2);
//	return  0;
/* No, actually we do care!  Check overlap of direct images */
	if (bbover (pri->val->dimg, sec->val->dimg) ) {
// 0 if no intersection, non-zero indicates conflict
	    cf = add_conflict (pri, sec, 5, 0);
	    if (obs->dlevel > 1) {
		printf (" Image conflict between %s and %s.\n",
		pri->dat.name, sec->dat.name);
	    }
	    k++;
	}
	return  k;
    }

/* Check types, if both are gap objects, we don't care */
    if (pri->type == OBJ_GAP && sec->type == OBJ_GAP) return 0;
// Actually, the above should never happen, but we test anyway!

/* Check types; if reference or gap at both pri and sec, we need to
  compare the direct image bound box of the refernce and the spectrum
  of the gap.  Otherwise, fall through to check the spectra.  For
  now, we just do this directly... */
    if (pri->type == OBJ_REFERENCE) {
/* Case of secondary being reference was handled above */
	if (sec->type == OBJ_GAP) {
	    if (bbover (pri->val->dimg, sec->val->sod->bb) ) {
		cf = add_conflict (pri, sec, 3, 0);
		if (obs->dlevel > 1)
		    printf (" Alignment %s falls in gap %s.\n",
			pri->dat.name, sec->dat.name);
		k++;
	    }
	    return  k;
	}
    } else if (pri->type == OBJ_GAP) {
/* Case of secondary being gap was handled above */
	if (sec->type == OBJ_REFERENCE) {
	    if (bbover (pri->val->sod->bb, sec->val->dimg) ) {
		cf = add_conflict (pri, sec, 3, 0);
		if (obs->dlevel > 1)
		    printf (" Alignment %s falls in gap %s.\n",
			sec->dat.name, pri->dat.name);
		k++;
	    }
	    return  k;
	}
    }


/* Check primary spectrum against secondary */
    if (spect_overlap (pri->val->sod, sec->val->sod, special) ) {
/* NOTE -- NEED TO ADD OVERLAP CONSIDERATION TO SPECT-OVERLAP STUFF */
	cf = add_conflict (pri, sec, 3, 0);
	if (cf != NULL) cc++;  // nonsense for lint
/* 3 = conflict type, pri spectra, 0 = inhibit off. */
	if (obs->dlevel > 1 || special) {
	    printf (" Spectrum overlap between %s and %s.\n",
	    pri->dat.name, sec->dat.name);
	}
	k++;
    }
// /* DEBUG section */
// 	else if (special) {
// 	printf (" xx %s - %s no spectral overlap.\n", c1, c2);
//     }  // end DEBUG section */

/* Now, check image spectra.  If we are not checking extra orders, the
  image spectra have not been computed and the ->img pointers are null.
  The overhead added here is equivalent to 3 "if" statements. */
/* Check primary images against secondary spectrum */
/* Skip this if secondary is reference hole, or a gap */
    if (sec->type != OBJ_REFERENCE && sec->type != OBJ_GAP) {
	for (ks=0,p=pri->val->img; p != NULL; p = p->next) {
	    if (spect_overlap (sec->val->sod, p, special)) {
		cf = add_conflict (pri, sec, 4, 0);
		if (cf != NULL) cc++;  // nonsense for lint
/* 4 is image conflict */
    if (obs->dlevel > 1) {
	    printf (" Spectrum overlap between image of %s and %s.\n",
	    pri->dat.name, sec->dat.name);
    }
		k++;
	    }
	    if (p->next == pri->val->img) break;
	    if (++ks > 10) {
		printf (" *** Debug exit 1!!\n");
if (ks > 50) exit (0);		// debug loop killer
		break;	// debug
	    }
/* NOTE -- The above error exit, and one below are effective when
  there is a large error between the object coordinates and the
  center coordinate.  This produces many objects off mask and off
  detector, and many error messages.  This exit cuts down the
  excessive number of errors.  A better fix is desired. */
	}
    }

/* Check secondary images against primary spectrum */
/* Skip this if primary is reference hole, or a gap */
    if (pri->type != OBJ_REFERENCE && pri->type != OBJ_GAP) {
	for (ks=0,p=sec->val->img; p != NULL; p = p->next) {
	    if (spect_overlap (pri->val->sod, p, special)) {
		cf = add_conflict (pri, sec, 4, 0);
		if (cf != NULL) cc++;  // nonsense for lint
    /* 4 is image conflict */
	    if (obs->dlevel > 1) {
		printf (" Spectrum overlap between %s and image of %s.\n",
		pri->dat.name, sec->dat.name);
	    }
		k++;
	    }
	    if (p->next == sec->val->img) break;
	    if (++ks > 10) {
		printf (" *** Debug exit 2!!\n");
		break;	// debug
	    }
	}
    }
    if (obs->dlevel > 3 && cc > 2)  // nonsense for lint
	printf (" %d conflicts between %s and %s.\n", cc,
		pri->dat.name, sec->dat.name);
    return  k;
}


int	seek_conflict (objq *pri, bool single, Obs *obs)
/* Scan nearby objects to the primary to find conflicts.  Scan only one
  direction (forward) if single is true.  Enter found into conflict
  queue and object structures.  Return number of conflicts entered. */
{
objq	*p;
int	k;	/* Result */
int	i;

    Tobj = NULL;
    k = 0;
/* Check zero-th for inactive object - not there => does not conflict. */
    if (!(pri->flag & OBJECT_ACTIVE)) return k;
    if (pri->val == NULL) return k;		// no spectrum

/* Check first for being within the mask space.  This is IMACS specific,
  and NOTE - it should be replaced by a mask descirption check.  */

/* Check for being on the detector?  If fails, set it to inactive,
  and forget about it. */
/* This check is (still) done in the general find_conflict thing,
  and stays there except for gross errors found already in fillobject */

/* Start with a loop to scan from this object forward until we get a
  status of Bx out of range. */
    C_loop (p, pri) {
	if ((i = find_conflict(pri, p, obs)) < 0) break;
// NOTE -- we should also break if Bx decreases
	k += i;
    }

/* If not single sided, do a second loop to scan from the last object
  backward until Bx out of range. */
    if (single) return k;

    for (p=pri->last; p != NULL; p = (p->last == pri) ? NULL : p->last) {
	if ((i = find_conflict(pri, p, obs)) < 0) break;
	k += i;
//	if (p->last == pri) break;
// BAD CONSTRUCT BREAK
    }
    pri->val->ncf += k;	/* Enter count in node searched */
    return  k;
}


int	scan_conflicts (Obs *obs)
/* Scan the entire object queue for conflicts.  Do this by scanning
the queue and calling seek_conflict for each object in the queue.
  Return is total number of conflicts found. */
{
objq	*p;
int	k=0;	/* return value */
int	n;	/* debug method */

    DeltaX = 1.1 * bigest_bounding_box (obs->ob);
    if (obs->dlevel > 1)
	printf ("  Scan delta = %.2f vs. old 2.5 mm.\n", DeltaX);  //debug
    C_loop (p, obs->ob) {
	if (obs->dlevel > 1) {	/* Debug method */
	    n = seek_conflict (p, 1, obs);
	    if (n > 0)
		printf (" Got %d conflicts for \"%s\".\n", n, p->dat.name);
	    k += n;
	} else {	/* Normal, non-debug method */
	    k += seek_conflict (p, 1, obs);
	}
    }
    return  k;
}

	/* - - - - - - - - - - *\
	|			|
	|  Conflict Management	|
	|			|
	\* - - - - - - - - - - */

/*  ===  Debug only routine for maskgen  ===  */
/* This program is called in maskgen and in 2 other places here too.
  It was in msdebug.c which is included by maskgen, but was not
  included when these utility programs were included from library
  in stand-alone test routines.  To prevent an undefined symbol in
  linking such routines, it is placed here, and referenced in the
  header file for the benefit of maskgen.   */

void	check_cfcount (objq* ob, char* title)
/* Check and report any conflict count discrepencies... */
{
objq	*p;
objval	*v;
cfl	*ch;
cfl	*cp;
int	csq;
int	tsq=0;
int	tcf=0;
int	tco=0;

/* Loop over all the objects */
    C_loop (p, ob) {
/* Find conflict count for node p in two ways, the p->val->ncf and
  a direct count of the per-object queue itself.  Report by node. */
	v = p->val;
	if (v != NULL) {
// Count the secondary queue from cfq here
	    for (csq=0, cp = v->cfq; cp != NULL; cp = cp->fwd) csq++;
// Total the counts
	    tsq += csq;
	    tcf += v->ncf;
// Compare the counts
	    if (csq != v->ncf) {
		printf (" %s Node %s type %d, ncf=%d, count=%d\n", title,
		p->dat.name, v->cfq->type, v->ncf, csq);
	    }
	}
    }
/* May want to keep a count of conflicts found from objects, and
  compare that to the total size of the conflict queue itself... */
	ch = find_cfl (ob);
// Count all the cfl nodes in main queue
    C_loop (cp, ch) tco++;
/* Replace with listcount function, tco = listcount(ch) */
/* Report those numbers */
    if (tco == tcf) {
	printf (" %s Total conflicts = %d", title, tco);
    } else {
	printf (" %s Total cfl=%d, total ncf=%d", title, tco, tcf);
    }
    if (tcf == tsq) {
	printf (".\n");
    } else {
	printf (", total csq=%d.\n", tsq);
    }
//    printf (" %s Total cfl=%d, total ncf=%d, total csq=%d.\n",
//	title, tco, tcf, tsq);
}

/* -- End debug only routines for maskgen -- */


/*  ===  Manage conflicts - drop, count, resolve  ===  */


double	hcpriority (objq* ob, double pv)
/* Return highest (numerically smallest) priority of a conflicting
  object w.r.t. the given object.  If no conflicting object, return
  the default value pv. */
/* Used only by conflictedest, drop_single_confs, and
    other test routines below */
{
double	rv = pv;	// Return value.
double	sp;		// Secondary priority
objq	*pri, *sec;	// Conflict object pointers
objval	*v;	// object's value pointer
cfl	*c;	// object's conflict queue pointer
cfl	*cp;	// scratch queue pointer

// Take care of null cases
    if (ob == NULL) return rv;
    v = ob->val;
    if (v == NULL) return rv;
    c = v->cfq;

// Look at the conflict queue
    for (cp=c; cp!=NULL; cp=cp->fwd) {
	pri = cp->obj;
	if (pri != ob) continue;
	sec = cp->other;
	if (sec == NULL) continue;
	sp = sec->dat.priority;
	if (cp == c || sp < rv) rv = sp;
// If secondary is type OBJ_REFERENCE, should return musthave ?
// Take care of that by setting priority of reference objects.
    }
    return  rv;
}


#undef  CHKMHP
/* Use symbol CHKMHP to check for conflicts with musthave nodes during
  the checking of single conflicts.  This is only strictly correct if
  the queue is ordered in decreasing priority, which would also have
  the side effect of making all duplicates remove the "other" node, with
  a slight performance hit.  Keep it undefined to do the prioritized
  conflict resolution in the main conflict resolution code, which is
  normally executed later.  */
/* NOTE -- the reason we do not do the removal of conflicts with
  musthave nodes in the first pass is that it is possible for
  those musthave nodes to be removed for single conflicts such as
  off detector during that first pass.  Since they could be removed
  later in the pass than a conflicting node, that conflicting node
  would have been removed for a cause which later ceases to exist,
  and it should not have been removed.  All the single conflicts
  need to be removed on the first pass, and then the conflicts with
  musthave nodes (which are dual, not single conflicts) then get
  resolved in a second pass.  */

int	drop_single_confs (Obs *obs)
/* Drop all objects in queue with single conflicts -- i.e. those that
conflict with themselves, mostly being off the detector.
  Return is number of objects found and removed. */
/* Note -- no calls to this routine yet */
/* NOTE - CALLED FROM maskgen IN MAIN-LINE PROGRAM */
{
int	k=0;	/* result value */
int	ns=0;	/* Singles */
int	nd=0;	/* Duplicates */
int	np=0;	/* Priority drops */
objq	*p;
cfl	*cf;
//#ifdef  CHKMHP
double	musthave;
double	dp;
//#endif

    if (obs == NULL) return 0;
//#ifdef  CHKMHP
    musthave = obs->Pmusthave;
//#endif
//	printf (" -- Begpoint in drop_single_confs\n");
//	fflush (stdout);
    C_loop (p, obs->ob) {
	if ((p->flag & OBJECT_ACTIVE) && (p->val != NULL)) {
/* Only done if active... */
/* Check the conflict queue of this object for any single
  conflicts -- null other node, off array, or conflicting with
  a must have while not a must have ourselves */
/* The label here is referenced within the loop to restart the loop,
  since we can potentially modify the queue being traversed. */
retry:	    for (cf = p->val->cfq; cf != NULL; cf = cf->fwd) {
//		if (  (cf->other == NULL)	// Single
		if (  (cf->other == NULL)	// Single
		    || (cf->type == 1)		// off array
		    || (cf->other == p)		// same node
#ifdef  CHKMHP
		    || (cf->other->dat.priority <= musthave &&
		                p->dat.priority > musthave)
#endif
		   ) {
		    if ( (p->type != OBJ_REFERENCE) || (obs->Flag & OFral) ) {
/* That says either it is not a reference, or, if it is, the flag
  must indicate automatic removal. */
			set_inactive (p);
			k++;  ns++;
		    }
/* Allowing off-array conflict with reference to stay */
		    break;
		}

/* Check the conflict queue element for a duplicate.  If found, we
  drop the lower (numerically higher) priority node.  That is either
  node p, or cf->other.  If we drop p, we break the loop, otherwise
  we do not break, but restart the loop here .  The set_inactive
  affects the cfq. */
/* If priority of duplicates is equal, decide by ordering value */
		if (cf->type == 5) {
objq	*r = NULL;
double	dpr;
		    dpr = cf->other->dat.priority - p->dat.priority;
		    if      (dpr > 0.0)	r = cf->other;
		    else if (dpr < 0.0)	r = p;
		    else {
#ifdef  PRIORDER
			if (cf->other->order > p->order) r = cf->other;
#else
			if (cf->other->order < p->order) r = cf->other;
#endif
			else r = p;
		    }

/* Temporary fix for possible problem.   ..KADEBUG.. */
// This section should not be needed or used.
// Without it, the program MAY hang in an infinite loop here...
		if (r == NULL) continue;
// Remove line above for testing; should be OK logically

		    set_inactive (r);
		    k++;  nd++;
		    if (r == p) break;   else  goto retry;

//debug		    if (r == p) break;  // else  goto retry;
// NOTE -- had to change from simple else above to test below  ..KADEBUG
// without the test it will hang...
//debug		    else if (r != NULL) goto retry;

		}

/* Note that we are not checking for circular cf->fwd queue here. */
	    }
	}
    }
/* That completes all the single node conflicts. */

    if (ns > 0 && obs->dlevel > 0)
	printf (" %4d objects removed for single conflicts.\n", ns);
    if (nd > 0 && obs->dlevel > 0)
	printf (" %4d objects removed as duplicates.\n", nd);

//	printf (" -- Midpoint in drop_single_confs\n");  //KADEBUG
//	fflush (stdout);

/* As a separate step, remove any nodes which are active yet below
  the must-have priority, which conflict with any node which has the
  must-have or higher priority.  Such nodes will always be eliminated
  by later stages in conflict resolution, and it is best that they
  not be present to possibly confuse resolution of otherwise possible
  nodes or potentially cause removal of an otherwise available node. */

/* NOTE THAT the terms high-low here refer to priority, while terms
  greater-less refer to the value of the priority number.  Thus,
  less is high and greater is low.  Try to be consistent on this. */

/* Two ways to do this.
    1) For each active node lower than Pdecide, find hcpriority,
	and if that priority is higher than Pdecide, drop the node.
    2) For each active node higher than Pdecide, scan the
	conflict queue, dropping all nodes in that queue.
It would appear that (1) requires fewer conflict sub-list restarts,
and so is potentially simpler.  We need to scan all nodes once to find
either removal candidates or all high priority nodes.  Eliminating the
scan of the conflict list too simplifies things.
   */

    dp = musthave + 5.0;	// Actually anything greater is good
    C_loop (p, obs->ob) {
// NOTE above is a new method of ending the object loop...  KADEBUG
/* Skip nodes which are not removal candidates... */
/*  ..Check this condition, and when found  --skip--    ..Leaving.. */
	if (!(p->flag & OBJECT_ACTIVE))    continue;	// Active
	if (p->val == NULL)		   continue;	// Has conflicts
	if (p->flag & OBJECT_VIRTUAL)      continue;	// non-gaps
	if (p->dat.priority <= musthave)   continue;	// Low priority
	if (hcpriority (p, dp) > musthave) continue;	// High conflict
/* Now, we conflict with a musthave, are not one ourselves, so the
  current node is dead meat... */
	k++;  np++;
	set_inactive (p);
    }
/* That completes all the simple conflicts, either single or by any
  removable node with a must-have node.  Return the object count. */

    if (np > 0 && obs->dlevel > 0)
	printf (" %4d objects removed for conflict with must-have.\n", np);

//	printf (" -- Endpoint in drop_single_confs\n");   // KADEBUG
//	fflush (stdout);
    return  k;
}


static	int	count_conobj (objq *head)
/* Return total conflict count from object queue */
{
objq	*p;
cfl	*q;
int	k=0;	/* Result */
int	n=0;	/* Alternate result */

    C_loop (p, head) {
	if (p->val != NULL) {
	    k += p->val->ncf;
/* ==NOTE== could also run the val->cfq list, counting it separately,
  then compare that count to k.  Report at end if differs... */
	    for (q=p->val->cfq; q != NULL; q = q->fwd) n++;
	}
    }
    if (n != k) printf (" ** %d conf. counts, %d nodes.\n", k, n);
    return  k;
}

/* This is where "conflictedest" was, long ago... */


//#define  DECONFERRP
#undef  DECONFERRP
// Symbol DECONFERRP compiles a error report section into the
// deconflict program.
// This section is obsoleted by program repconf
// Leave DECONFERRP undefined to use repconf properly.

/* Special edition of C_loop for the conflict per-node queue; this may
  become useful elsewhere here... */
#define  Cf_loop(q,a) for(q=(a); q!=NULL; q=(q->fwd!=(a))?q->fwd:NULL)

int	repconf (Obs *obs)
/* Report any important outstanding conflicts, return number found.  */
{
int	nc=0;	// conflicts found
objq	*p;
objval	*v;
cfl	*c;
cfl	*cp;
cfl	*head=NULL;
cfl	*q;
int	d;
#ifndef  DECONFERRP
objq	*pri, *sec;
char	prilab[24];
char	seclab[24];
#endif

// Check for null case
    if (obs == NULL) return nc;

/* Loop over all the objects in list to find conflicts */
    C_loop (p, obs->ob) {
	v = p->val;
	if (v == NULL) continue;
	c = v->cfq;
	if (c == NULL) continue;

/* Now we have a per-object conflict queue, loop through it... */
	Cf_loop (cp, c) {

/* A conflict has been found; we need to:
  1) remove duplicates pri/sec vs sec/pri
  2) enter a copy on reportable queue  .. */

/* Look in local queue if any for duplicate if any */
	    if (head != NULL && cp->other != NULL) {
		d = 0;
		C_loop (q, head) {
		    if (cp->obj == q->other && cp->other == q->obj) {
			d++;
			break;
		    }
		}
		if (d > 0) break;	// duplicate is ignored.
	    }

/* No head or not duplicate, make a copy into local queue */
	    q        = xalloc (cfl);
	    q->next  = q->last = q->bak = q->fwd = NULL;
	    q->type  = cp->type;
	    q->obj   = cp->obj;
	    q->other = cp->other;
	    if (q->other == NULL) q->other = q->obj;	// self-conflicts
	    head = push_cfl (q, head);
	    nc++;
	}	// End loop over per-object queue
    }	// End loop over all objects.

/* Loop over local reportable queue to write reports */
/* Remove all entries in local queue at same time */
    while (head != NULL) {

/* Pop the reportable conflict from local queue. */
	c = head;
	head = pop_cfl (head);

/* Report on the conflict in local conflict c.  */
#ifndef  DECONFERRP
	pri = c->obj;
	get_type (prilab, pri, obs);
	sec = c->other;
	if (pri == sec || sec == NULL) {
	    printf (" ** CONFLICT -- %s \"%s\" is off detector.\n",
		    prilab, pri->dat.name );
	} else {	// NOTE that here sec cannot be NULL
	    get_type (seclab, sec, obs);
	    printf (" ** CONFLICT between %s \"%s\" <and> %s \"%s\".\n",
		    prilab, pri->dat.name,
		    seclab, sec->dat.name  );
	}
#endif
	free (c);
    }	// End loop removing local queue entries
    return  nc;
}


/*  ===  Object queue status -- set active/inactive  === */

// NOTE -- set active is not used anywhere!!  (yet)
int	set_active (objq *cur, Obs *obs)
/* Set an object active; scan for conflicts and put any on.  Return
number of conflicts found.  */
{
int	ul;

/* Check if already active; if so, quick exit -- be sure to set
the return depending on whether this node has conflicts */

/* We need a subroutine to count conflicts? */
    if ((cur->flag & OBJECT_ACTIVE) && (cur->val != NULL)) {
	return  cur->val->ncf;
    }

/* Check the re-use count to see if use of this object is allowed... */
    if (cur->type == OBJ_OBJECT) ul = obs->reuseob;
    else if (cur->type == OBJ_REFERENCE) ul = obs->reuserf;
    else  ul = 999;

/* Scan the node for conflicts, return status found */
    if (cur->use <= ul) cur->flag |= OBJECT_ACTIVE;
    if (cur->val == NULL) fillobject (cur, obs);
    return  seek_conflict (cur, 0, obs);
}
/*  NOTE THAT set_active SHOULD NOT BE USED.  AN OBJECT MAY BE SET
  INACTIVE IN THE INITIAL FILL_OBJECT ROUTINE IF IT IS TOO FAR OFF
  THE MASK OR DETECTOR, AND IMAGES OF THAT OBJECT ARE THEN NOT
  COMPUTED.  SO IT IS INCOMPLETE.  */


int	set_inactive (objq *cur)
/* Set an object inactive.  Remove all its conflicts from the queues,
possibly causing conflict queue to become null...
  Return is number of conflicts resolved.  */
{
int	k;	/* Return value */
int	nc=0;

/* If not currently active, return zero... */
    if (!(cur->flag & OBJECT_ACTIVE)) return  0;

/* OK, we have active node.  If no conflict count, we can return. */
//    cur->flag ^= OBJECT_ACTIVE;  /* Toggle bit off */
    cur->flag &= ~OBJECT_ACTIVE;  /* Set bit off */
    k = cur->val != NULL ? cur->val->ncf : 0;
    if (k == 0) return 0;	// cur->val is now not null...

/* We need to remove some conflicts. */

// debug-3-
//    if (!(cur->type == OBJ_OBJECT))
//	printf (" !! set_inactive on REFERENCE star\n");

    while (cur->val->cfq != NULL) {
	cur->val->cfq = drop_conflict (cur->val->cfq, 0);
	nc++;
//	cur->val->ncf--;
// Above already done in drop_conflict
    }
/* NOTE== As a DEBUG check, verify that ncf is now zero, squak if not */
    if (cur->val->ncf != 0) {
	printf (" *Setting %s inactive, conflict count remains %d!",
		cur->dat.name, cur->val->ncf);
	printf (" (was %d, count %d)\n", k, nc);
    }
    return  k;
}


	/* - - - - - - - - - - *\
	|			|
	|  Conflict Resolution	|
	|			|
	\* - - - - - - - - - - */


/*  ===  Mark-II Conflict Resolution Algorithm  ===  */

objq*	max_conflict (objq** loq, objq** hiq)
{
// Works similarly to conflictedest
// Scan the loq queue
// If any node has no conflicts, push it onto the hiq queue
// Find most conflicted node
// If any tie, take one with lowest priority
// Return most conflicted node
// If none found, null return
objq	*q=NULL;	// Candidate to return
objq	*p;		// pointer for scan
objq	*r;		// Aux. pointer.
// objq	*s;		// Start point
objval	*v;		// Extra pointer
int	k=0;		// most conflicts found
int	j;		// Loop control
int	nc;		// Conflict count for node

    if (loq == NULL) return q;
    if (*loq == NULL) return q;

    for (j=0,p=*loq;
	(p != NULL) && !((p == *loq) && j); ) {

	v = p->val;
	nc = (v == NULL) ? 0 : v->ncf;	// Conflict count of p

/* This compound if decides whether to drop a non-conflicting node,
  or to save a new candidate for maximum conflicts. */
	if (nc == 0 && hiq != NULL)	// not a conflict
	    {			// Drop node
	    r = pop_obj (p);	// Move from loq to hiq
	    *hiq = push_obj (p, *hiq);
	    if (*loq == p) *loq = r;
	    p = r;
	    continue;
	    }
	else if (nc == 0)		// not a conflict, no hiq
	    ;			// Do nothing, pass the node
	else if (nc > k || q == NULL)	// New candidate
	    {			// Save node and count
	    q = p;
	    k = nc;
	    }
	else if (nc == k && p->dat.priority > q->dat.priority)
	    q = p;
	else if (nc == k && p->dat.priority == q->dat.priority &&
#ifdef  PRIORDER
		p->order > q->order)
#else
		p->order < q->order)
#endif
	    q = p;

/* End of loop, advance p to next node */
	p = p->next;
	j++;	// do not increment if continued
    }
    return  q;
}


int	order_decide (objq* a, objq* b)
/* Return -1, 0 or 1 depending on a>b, a=b, a<b */
/* Check order of nodes. */
{
int	d;

/* If either is NULL, return 0 */
    if (a == NULL) return 0;
    if (b == NULL) return 0;

/* Get order diff */
    d = b->order - a->order;

/* Check the d value */
    if      (d > 0) return  1;
    else if (d < 0) return -1;
    else            return  0;
}


int	priority_decide (objq *a, objq *b)
/* Return -1, 0 or 1 depending on a>b, a=b, a<b */
/* Use numeric priority, lower comes first. */
/* Ties are ordered depending on order field */
{
double	d;

/* If either is NULL, return 0 */
    if (a == NULL) return 0;
    if (b == NULL) return 0;

/* Get priority difference */
    d = b->dat.priority - a->dat.priority;

/* Check the d value */
//    if (d == 0.0)	return 0;
//    if (d > 0.0)	return 1;
//    else		return -1;

/* Change this to return 1 or -1, but if d is 0, return
  the result of order_decide on the nodes. */

    if (d > 0.0) return  1;
    if (d < 0.0) return -1;
#ifdef  PRIORDER
    return  order_decide (a, b);
#else
    return  (- order_decide (a, b));
#endif
//return  0;  //debug   KADEBUG
}


int	deconflict (Obs *obs)
/* Remove conflicts.  Split out all objects below decision priority,
and loop to set those inactive.  Sort all objects on priority, and
loop to remove their conflicts.  Restore sort order. */
/* Return the number of objects rendered inactive */
{
objq	*head;
objq	*lowq=NULL;
int	k=0;	/* Count removed */
int	k1;
int	n;	/* Conflict count */
// int	nc;	/* intermediate */
int	j;	/* Conflict count at start */
int	ks;	/* scratch for ploop2 */
// int	dun;
objq	*p;	/* point to object */
objq	*r;	/* scratch */
double	Pd;
objval	*v;
cfl	*cq;
cfl	*qq;
#ifdef  DECONFERRP
objq	*pp = NULL;
objq	*po = NULL;
#endif
// int	nl;		// debug
int	tar;	// Total active references
int	erc;	// Excess references
int	dc;	// count delta
int	kd;	// skip count

/* Timing things for test/debug */
clock_t	begn;
clock_t	start, stopt;
clock_t	cftot, rmtot;
clock_t	begtt;
clock_t	sortt;
double	scale;
double	msecs;

/* Initialize */
    head = obs->ob;
//    debug_show_conflicts (head, "Start deconflict");
    j = count_conobj (head);	// PROBABLY NOT NEEDED

	/* - - - - - - - - - - - - - - *\
	|				|
	|  Removing Reference Objects	|
	|				|
	\* - - - - - - - - - - - - - - */

/* New Feature.  Ability to remove conflicting alignment objects.
  If the (obs->Flag & OFral) flag bit is on, we remove reference objects
  which have their direct image bounding boxes conflicting with any
  gap object, or with a higher priority reference object's direct
  image bounding box.  These conflicts are (now) found by the find_conflict
  program when reference objects and gaps are detected there.  These
  removals must be done first, so other objects are allowed to exist
  which do not conflict with other things but which may have conflicted
  with a reference object which would be removed.  */
/* Also, excess reference objects, exceeding reflimit value are removed
  in this section after dropping any off-detector and conflicts. */
    if (obs->Flag & OFral) {
	C_loop (p, head) {	// over objects
/* Locate all active reference objects, which have conflicts with
  gap objects or other reference objects. */
	    if (p->type != OBJ_REFERENCE) continue;
	    if (!(p->flag & OBJECT_ACTIVE)) continue;
/* Now, p is an active reference object, look for conflicts... */
	    v = p->val;
	    if (v == NULL) continue;	// No objval
	    if (v->ncf == 0) continue;	// No conflicts
/* Doing this in a similar way to the priority removal section below */
/* Loop over conflicting nodes; set p inactive if the conflicting
  node is a gap or a higher priority reference.  Only.  */
	    for (cq = v->cfq; cq != NULL; cq = qq) {	// over conflicts
		r = cq->other;
#ifndef  NOCONFLRAL
		if (r == NULL) r = cq->obj;
#endif
		qq = cq->fwd;
		if ( (r->type == OBJ_GAP) ||
		     ( (r->type == OBJ_REFERENCE) &&
		       (r->order < p->order) )
// may need to decide on ordering test above?
		   ) {
		    k += set_inactive (p);
		    break;	// exit conflict loop
		}
	    }	// loop over conflicts
	}	// loop over objects

	/* - - - - - - - - - - *\
	|			|
	|   Reference  Limit	|
	|			|
	\* - - - - - - - - - - */

/* Limit the reference objects using the reflimit parameter in the
  obs structure.  Count number of active references here, and if that
  exceeds the limit, remove them one at a time until the limit is met. */

/* Count number of active reference objects as tar */
	tar = 0;		// Total Active Reference objects
	C_loop (p, head) {
	    if (p->type != OBJ_REFERENCE) continue;
	    if (!(p->flag & OBJECT_ACTIVE)) continue;
	    tar++;		// count the active ones
	}

/* Find number of excess reference objects as erc */
	erc = tar - obs->reflimit;	// Excess Reference Count
	if (obs->reflimit >= 99 && erc >= 0) erc = -1;	// omit if default limit

/* Set a count delta as dc, either method here: */
//	dc = (erc > 0) ? ceil ((double)tar/(double)erc) : 1;
	dc = (erc > 0) ? (tar+erc-1)/erc : 1;

/* Initialize the loop pointer */
	p = head;

/* Put out a message if reference objects are being removed... */
	if (erc > 0) printf (" *- Removing %d excess alignment objects.\n",
			erc);

/* Loop while the erc is over zero */
	while (erc > 0) {
	    for (kd=dc; ; p = p->next) {
/* Skip around the non-reference, non-active objects */
		if (p->type != OBJ_REFERENCE) continue;
		if (!(p->flag & OBJECT_ACTIVE)) continue;
/* Skip while counting to find the dc'th one */
		if (--kd > 0) continue;
/* At this point we have counted dc active objects and have one... */
		k += set_inactive (p);	// remove an extra reference
		erc--;		// track the excess for loop control
		break;	// Get to the outer loop on erc > 0 here
	    }	// counting the next dc active references
	}	// while excess reference count is over zero
/* Now, the excess reference count is no longer over zero */

    }	// Removal of conflicting reference objects



	/* - - - - - - - - - - - - - - *\
	|				|
	|   Splitting Object queue	|
	|				|
	\* - - - - - - - - - - - - - - */


/* Split the object queue into two at Pdecide.  Loop through, pulling
  out nodes below Pdecide and putting them in the low queue. */
    Pd = obs->Pdecide;
    if (obs->Pmusthave > obs->Pdecide) Pd = obs->Pmusthave;

#define PLOOP2
// #undef  PLOOP2
#ifdef  PLOOP2

    for (ks=0,p=head;
	(p != NULL) && !((p == head) && ks); ) {

//	k++;  // wrong placement, only count if p gets incremented
	if (    (p->dat.priority > Pd)  &&
		(p->type == OBJ_OBJECT) ) {	// Pull it
// debug-2- make change if still bad -- yes it was!
// Don't put in low queue if object is a reference!
	    r = pop_obj (p);
	    if (p == head) head = r;
	    lowq = push_obj (p, lowq);
	    p = r;
	    continue;
	}
	p = p->next;
	ks++;
// BAD STRUCTURE LOOP - FIX UP?  See newloop.txt for hints.
    }
#else
// This one works, one above does not...
    p = head;
    while (head != NULL && p != NULL) {
	if (    (p->dat.priority > Pd)  &&
		(p->type == OBJ_OBJECT) ) {	// Pull it
	    r = pop_obj (p);
	    if (p == head) head = r;
	    lowq = push_obj (p, lowq);
	    p = r;
	    continue;
	}
	p = p->next;
	if (p == head) break;
// BAD CONSTRUCT BREAK
    }
#endif
// {		// debug only section
// int	nl, nh;
// 	nl = listcount (lowq);
// 	nh = listcount (head);
// 	printf (" * Pulled %d objects leaving %d in high queue.\n",
// 		nl, nh);
// 	printf (" * The decision level is %.2f.\n", Pd);
// 	printf (" # There are %d conflicting objects\n", j);
// }

	/* - - - - - - - - - - - - - - *\
	|				|
	|   Conflict count algorithm	|
	|				|
	\* - - - - - - - - - - - - - - */

/* Run the conflict count algorithm on the low queue, using the
  max_conflict decision support routine */

/* This is the conflict count algorithm. */
//    n = k = 0;
    n = k;
    cftot = rmtot = 0;				// debug timing
// We should time the following loop
    start = clock();				// debug timing
//    while (1) {	/* A forever loop */
    for (ever) {
// Time the following subroutine call
	begn = clock();				// debug timing
	p = max_conflict (&lowq, &head);
	cftot += clock() - begn;		// debug timing
	if (p == NULL) break;
	if (n < 100*(obs->dlevel-1)) {	// Added arbitrary count limit
	    printf ("  --Remove node %s\n", p->dat.name);
	}
// Time the following subroutine call
	begn = clock();				// debug timing
	n += set_inactive (p);
	rmtot += clock() - begn;		// debug timing
	k++;
    }
    stopt = clock() - start;			// debug timing

// Set scales for time reports
    scale = 1.0 / CLOCKS_PER_SEC;
    msecs = 1.0e3 * scale;

    if (k > 0 && obs->dlevel > 0) {	// debug, should be 0
	printf ("  Conflict count removed %d objects in %.3f seconds.\n",
	    k, scale*stopt);
	if (k > 1) printf ("  Logic used %.3f sec., removal %.3f millisec.\n",
	    scale*cftot, msecs*rmtot);
    }
    k1 = k;

	/* - - - - - - - - - - - - - - *\
	|				|
	|   Priority object algorithm	|
	|				|
	\* - - - - - - - - - - - - - - */

/* Push any remaining nodes back into the main queue, and then sort
  the main queue on priority alone */
    head = push_obj (lowq, head);
//    begn = clock();
    head = sort_obj (head, priority_decide);
//    sortt = clock() - begn;
//    printf ("  Sorted by priority in %.3f millisec.\n", scale*(double)sortt);
    obs->ob = head;	// Put in new pointer.

// Debug - check the conflict counts...
    if (obs->dlevel > 1) check_cfcount (head, "dcf:1");

//    rmtot = 0.0;
    begtt = clock();

/* The revised scan loop to inactivate all conflicting nodes unless
  they are of sufficient priority.  Scan once by priority. */
//    for (j=0,p=head;
//	    (p != NULL) && !(p == head && j > 0);
//		p = p->next) {
    j = 0;	// PROBABLY NOT NEEDED
    C_loop (p, head) {
//int	interesting = 0;	// debug
//	if (p->val != NULL) interesting = p->val->flag;
	j++;
	v = p->val;
	if (v == NULL) continue;	// No objval
	if (v->ncf == 0) continue;	// No conflicts

//	if (interesting) printf (" Conflicts for %s (%d).\n",
//		p->dat.name, v->ncf); // debug

/* Loop over all conflicting nodes and remove any which are not either
  musthave or reference objects. */
	for (cq = v->cfq; cq != NULL; cq = qq) {
	    r = cq->other;
//	if (interesting) printf ("  -- %s type %d\n",
//		r->dat.name, cq->type );	// debug
	    qq = cq->fwd;	// NOTE how we handle loop here!
// The cq->fwd and cq itself goes away during the set_inactive call;
// We save the next node here to avoid the problem.

#ifndef  NOCONFLRAL
	    if (r == NULL) continue;	// allow single conflict
#endif

	    if (r->type != OBJ_OBJECT) continue;	// debug try -
	    if (r->dat.priority > obs->Pmusthave) {
//		begn = clock();			// debug timing
		n += set_inactive (r);
//		rmtot += clock() - begn;		// debug timing
		k++;
		qq = v->cfq;  // We might have to start over here.
/* The above line appears to be needed in some cases, depending on
  the ordering of the conflicts and how they were scanned.  Since it
  only executes when a conflicting node is set inactive, there should
  not be an infinite loop problem because of it.

  Actually, when we allow "extra orders" feature, some conflicting
  nodes can have 2 conflicts of 2 different types, with the spectrum
  and with an extra order of that spectrum.  Removing the node also
  removes the extra orders, and the second conflict.  Since that conflict
  was the next in the queue, the qq pointer is killed, terminating
  the loop before it can remove the other nodes.  The only safe way
  around this is to start this loop over after each removal.   */
	    }
	}	// Loop over conflicted nodes
    }	// Loop scanning the total object list
    sortt = clock() - begtt;

	/* - - - - - - - - - - - - - - *\
	|				|
	|   Conflict queue debugging	|
	|				|
	\* - - - - - - - - - - - - - - */

//    printf (" Scanned %d objects in %.3f millisec.\n", j, scale*(double)rmtot);
//    printf (" Removed %d conflicts.\n", n);
//    printf ("  Removed by priority in %.3f millisec.\n", scale*(double)sortt);

    if (k > k1 && obs->dlevel > 1) {
	printf (" Removed %d objects in %.3f millisec.\n",
		k-k1, msecs*sortt);
    }

// Debug - check the conflict counts...
    if (obs->dlevel > 1) check_cfcount (head, "dcf:2");

#if (DEBUG_UTL > 0)
/* Debug section - check the conflict queue itself */
if (obs->dlevel > 1) {			// debug section
cfl	*ch;
cfl	*cp;
int	nc=0;
    ch = find_cfl (head);
    C_loop (cp, ch) {
	nc++;
	if (cp->next == ch) {
	    printf (" **Conflict queue is circular. ");
	    break;
	}
    }
    printf ("  Found %d nodes in conflict queue.\n", nc);


/* If a small number of conflicts (say < 10) we should dump the conflict
  queue here... */
//    if (nc < 10) dump_conf (head);

}	// end debug section
#endif

#ifdef  DECONFERRP
/* -- Error message section -- */
/* Scan the objects and conflicts.  If any must-have object conflicts
  with something (other must-have, reference, gap) we need to issue
  a message to the user about that... */
    C_loop (p, head) {
/* Sorted by priority, leave if we are past the musthave break point */
	if (p->dat.priority > obs->Pmusthave) break;
/* If no conflict, slide on by */
	v = p->val;
	if (v == NULL) continue;
/* Now, report only conflicts with higher (numerically less) priority
  objects, to avoid duplicate reports. */
	C_loop (cq, v->cfq) {
char	t1[32], t2[32];
#ifndef  NOCONFLRAL
	    if (cq->other == NULL) continue;
#endif
	    if (cq->obj->dat.priority < cq->other->dat.priority) continue;
	    if (pp == cq->obj && po == cq->other) continue;
	    get_type (t1, cq->obj, obs);	// object type for message
	    get_type (t2, cq->other, obs);
	    printf ("** CONFLICT between %s \"%s\" and %s \"%s\"\n",
		t1, cq->obj->dat.name, t2, cq->other->dat.name);
	    pp = cq->obj;	// To eliminate duplicate
	    po = cq->other;	// conflict messages.
	}
    }
#endif


/* Restore the original sort order for objects */
//    begn = clock();
    head = sort_obj (head, order_decide);
//    sortt = clock() - begn;
//    printf ("  Sorted by order in %.3f millisec.\n", msecs*(double)sortt);
    obs->ob = head;	// Put in new pointer.

/* Return the total count of nodes removed */
    return  k;
}	// end of deconflict


/*  ---  Remove inactive objects into separate queue  --  */

/* Upgraded version of remove inactive objects */
objq*	inact_obj (Obs *obs)
/* Return a queue of inactive objects, removed from the obs->ob
  queue in the observation structure. */
{
objq	*head;
objq	*p;
objq	*q;
objq	*b = NULL;
    for (p = head = obs->ob;  p != NULL;  )  {
//	if (!(p->flag & OBJECT_ACTIVE) ||
//		(p->flag & OBJECT_VIRTUAL) ) {  /* Pop/push p here */
/* Leave virtual objects here if active, to enable plotting of gap
  objects in final debug plot of the object queue. */
	if (!(p->flag & OBJECT_ACTIVE) ) {  /* Pop/push p here */
	    p = pop_obj (q = p);   /* Save p in q, pop and set p to next */
	    if (p == head) p = NULL;	/* To terminate */
	    if (head == q) head = p;	/* Head popped, replace it. */
	    b = push_obj (q, b);	/* push q into b list */
	} else {		/* No action, just follow link in p */
	    p = (p->next == head) ? NULL : p->next;
	}
    }
    obs->ob = head;	/* Set (possibly new) active list head */
    return  b;
}



/*  ===  Writing SMDF file from obs structure  ===  */
/* Write the Slit Mask Definition File from the mask data structure */

static	void	write_comment (FILE* file, namlist* comments, char c)
/* Write !c followed by comments onto file */
{
namlist	*r;
    if (comments == NULL) return;
    C_loop (r, comments) {
	fprintf (file, "!%c%s\n", c, r->name);
    }
}


	/* - - - - - - - - - - *\
	|			|
	|   Write  SMDF file	|
	|			|
	\* - - - - - - - - - - */

#ifdef  GSDATA
static	void	gpout (gpos g, FILE* f1, FILE* f2, int wl)
{
char	rbuf[16];
char	dbuf[16];
double	e;

    sexigw (rbuf, g.ra,   3, 1);
    sexigw (dbuf, g.dec, -3, 0);	/* Special feature */
    e = (g.ep == 0.0) ? 2000.0 : g.ep;

    fprintf (f1, " %s %s %.1f", rbuf, dbuf, e);
    if (wl) fprintf (f2, " %s %s %.1f", rbuf, dbuf, e);
}
#endif  // on GSDATA


void	write_smdf (Obs* obs)
/* Write the SMDF file; all data needed is in obs structure */
{
FILE	*f;
double	scale;
double	uca, ucb;
double	alen, blen;
double	roff, angl;
objq	*r;
char	name[16];
char	*cq;
char	*nv;
//namlist	*np;
char	line[128];
char	buf[128];
int	yr, mo, dy;
char	rabuf[32];
char	decbuf[32];
int	wl;
// spect	*sp;	/* NOTE debug only */
extern  char PROGNAME[];
extern	char VERSION[];
/* The extern reference must use [] rather than char* ... */

    wl = (obs->dlevel > 1);	// or debug level when replaced

/* Open the file */
    strncpy (name, obs->fname, 8);
    name[8] = '\0';
    strcat (name, ".SMF");
    f = fopen (name, "w");
    if (f == NULL) {
	strcpy (line, "Opening output file \"");
	strcat (line, name);
	strcat (line, "\"");
	perror (line);	/* to stderr */
	return;
    }

    if (obs->dlevel > 0) printf (" Writing SMDF file to %s\n", name);

    scale = foclen (obs) * ArcSecond;
    uca = scale * obs->dca;
    ucb = scale * obs->dcb;

    if (wl) fprintf (stderr, "\n - - - SMDF file - - -\n\n");

/* Write the basic data */
    fprintf (f, "NAME %.8s\n", obs->fname);
    if (wl) fprintf (stderr, "NAME %.8s\n", obs->fname);

// Note -- These strings should print left justified.

    fprintf (f, "OBSERVER %.16s\n", obs->oname);
    if (wl) fprintf (stderr, "OBSERVER %.16s\n", obs->oname);

    if (obs->title == NULL) {
	fprintf (f, "# Untitled\n");
	if (wl) fprintf (stderr, "# Untitled\n");
    } else {
	fprintf (f, "TITLE %s\n", obs->title);
	if (wl) fprintf (stderr, "TITLE %s\n", obs->title);
    }

/* note -- come up with date in the line field here... */
    jd2date (jdnow()-(8.0/24.0), &yr, &mo, &dy);
/* Changed to produce the actual date that the program was run, in
  the Pacific time zone... */
    sprintf (line, "%d-%02d-%02d%c", yr, mo, dy, AN);

// Add a standard time stamp
    strcpy (buf, timestamp (NULL));

    fprintf (f, "MADE %s By %s version %s%s = %s\n",
			line, PROGNAME, VERSION, MGUVERS, buf);
    if (wl) fprintf (stderr, "MADE %s By %s version %s%s = %s\n",
			line, PROGNAME, VERSION, MGUVERS, buf);

//  printf (" Location 2\n"); fflush (stdout);  //debug

/* Previous writing of telescope data was totally bogus; we now do it
  the right way, using the obs->Tel_scope pointer to the telescope
  element data structure, and some slight assumptions... */

    fprintf (f, "TELESCOPE %7.1f %7.5f %d %s\n",
	foclen (obs),
	eval ( ((focdat*)obs->Tel_scope->data)->curv, 0.0),
	3,	/* Reflections -- NOTE -- get that right in optics! */
	obs->Tel_scope->name);
    if (wl) fprintf (stderr, "TELESCOPE %7.1f %7.5f %d %s\n",
	foclen (obs),
	eval ( ((focdat*)obs->Tel_scope->data)->curv, 0.0),
	3,	/* Reflections -- NOTE -- get that right in optics! */
	obs->Tel_scope->name);


//  printf (" Location 3\n"); fflush (stdout);  //debug

    if (obs->Instr_mt == NULL) {
	cq = "Unknown";		// Error stop
    } else {
	cq = obs->Instr_mt->name;
    }

/* NOTE -- need to get proper names into structure */
    fprintf (f, "INSTRUMENT %s \n", cq);
    if (wl) fprintf (stderr, "INSTRUMENT %s \n",
		cq);

// printf (" Location 4\n"); fflush (stdout);  //debug

/* The GRATING output item is changed to DISPERSER, with a changed
  information content of <name, order, angle> with grating density
  obtained from the optic element of that name.  The GRATING output
  is then obsolete, and will not be done. */

/* Write the DISPERSER information -- name, order, and possibly the
  angle used.  */

    if (obs->Disperser == NULL) {
	printf (" ** NULL disperser field in obs!\n");
	fflush (stdout);
    } else {
	fprintf (f, "DISPERSER  %s  %d",
	    obs->Disperser->name,
	    obs->order);
	if (GRangle != NULL) fprintf (f, "  %7.5f", *GRangle/Degree);
	fprintf (f, "\n");
    }

/* Add a comment on grating angle.  If we really have a grating rather
  than a grism, the angle desired is angle - 22.5 degrees, the angle
  which the grating is tilted past the zero-dispersion angle.  */
// GratingAngle
    if (obs->Disperser != NULL) {
	if (obs->Disperser->kw == GRAting) {
	    fprintf (f, "# Grating Angle = %.4f Degrees.\n",
		*GRangle/Degree - 22.5);
// Could also use obs->D_angle which is pointed to by GRangle
	}
    }

    sprintf (line, "%s", sexig(obs->ra/Hour, 3, 3));
    sprintf (buf, "%s", sexig(obs->dec/Degree, 3, 2));
// CHANGE TO USE sexigw HERE!!

// printf (" Location 5\n"); fflush (stdout);  //debug


    fprintf (f, "POSITION 0 %s %s %6.1f %4.1f 0.0 0.0\n",
	line, buf, obs->equinox, obs->angle/Degree);
    if (wl) fprintf (stderr, "POSITION 0 %s %s %6.1f %4.1f 0.0 0.0\n",
	line, buf, obs->equinox, obs->angle/Degree);
/* NOTE -- need that sub-mask subscript, and field center x,y !! */

/* Add comments on rotator flip passed from intgui, if any */
//    for (np=obs->pcom; np != NULL; np = np->next) {
//	fprintf (f, "#!%s\n", np->name);
//	if (np->next == obs->pcom) break;
//    }
    write_comment (f, obs->pcom, '#');
    if (wl) write_comment (stderr, obs->pcom, '#');

	/*-----------------------------*\
	|				|
	|    Observing Line comment	|
	|				|
	\*-----------------------------*/

/* Add comment replicating some observing line parameters */
    roff = (InstIndx(obs->Instr_mt->name) == I_lds3) ? LDSSROA : CHUECO;
    angl = fmod (obs->angle/Degree + 90.0 + roff, 360.0) - 180.0;
    sexigw (line, obs->ra /Hour,   3, 3);
    sexigw (buf, obs->dec/Degree, -3, 2);	// Special feature!
    fprintf (f, "!.OC %.8s %s %s %6.1f 0.0 0.0 %7.2f OFF ",
	obs->fname, line, buf, obs->equinox, angl);
    if (wl) fprintf (stderr, "!.OC %.8s %s %s %6.1f 0.0 0.0 %7.2f OFF ",
	obs->fname, line, buf, obs->equinox, angl);
/* That puts observing catalog parameters 2-9 in the comment following
  the lead in of !.OC for a script to make catalog entries from all
  .SMF files present or named. */
/* Add observing parameters 10-15 to comment */
//#define  FULLOBC
// We want to define that at the next release...
// change FULLOBC to GSDATA
//#ifdef  FULLOBC
#ifdef  GSDATA
// These 6 things are ra/dec/epoch of guider probe 1, and same for gp2
// with ra,dec in hh:mm:ss and +dd:mm:ss, epoch in years.
// Could be filled in if intgui put guider positions into the .obs
// file, and we read them in here...
// That could also require or encourage guider star selection in intgui...
//    fprintf (f, "0.0 0.0 0.0 0.0 0.0 0.0 \n");
//    if (wl) fprintf (stderr, "0.0 0.0 0.0 0.0 0.0 0.0 \n");
    gpout (obs->gp1, f, stderr, wl);
    gpout (obs->gp2, f, stderr, wl);
//#else
#endif
/* or not... */
    fprintf (f, " \n");
    if (wl) fprintf (stderr, " \n");


    fprintf (f, "WAVELENGTH %6.1f\n", obs->cw);
    if (wl) fprintf (stderr, "WAVELENGTH %6.1f\n", obs->cw);

    fprintf (f, "TEMPERATURE %4.1f\n", obs->temp);
    if (wl) fprintf (stderr, "TEMPERATURE %5.2f\n", obs->temp);

    fprintf (f, "EPOCH %10.5f\n", obs->epoch);
    if (wl) fprintf (stderr, "EPOCH %10.5f\n", obs->epoch);

    fprintf (f, "WLIMIT %6.1f %6.1f\n", obs->pb.blue, obs->pb.red);
    if (wl) fprintf (stderr, "WLIMIT %6.1f %6.1f\n",
				obs->pb.blue, obs->pb.red);

/* Check if pb2 has been set and set pb2act if so... */
//    if ( (obs->pb2.blue != obs->pb.blue) ||
//	 (obs->pb2.red  != obs->pb.red ) )  obs->pb2act = 1;

    if (obs->pb2act) {
	fprintf (f, "DLIMIT %6.1f %6.1f\n", obs->pb2.blue, obs->pb2.red);
	if (wl) fprintf (stderr, "DLIMIT %6.1f %6.1f\n",
				obs->pb2.blue, obs->pb2.red);
    }

/* NOTE -- we also need obs->epoch computed, obs->pb somewhere, and
  several other things... */

#ifdef  DIFREF
//    if (!(obs->ha < -24.0)) {
	fprintf (f, "HANGLE %8.3f\n", obs->ha);
	if (wl) fprintf (stderr, "HANGLE %8.3f\n", obs->ha);
//    }
// Always output ha
#endif

/* Write the global comments */
    write_comment (f, obs->comments, '=');
    if (wl) write_comment (stderr, obs->comments, '=');

/* Write a comment if slits have been extended */
    if ( (obs->slext & 3) == 3) {
	fprintf (f, "# Slits have been extended.\n");
	if (wl) fprintf (stderr, "# Slits have been extended.\n");
    }

/* Loop to write objects */
    C_loop (r, obs->ob) {
	if (r->flag & OBJECT_ACTIVE) {	// only do active ones
	    if (r->type == OBJ_OBJECT) nv = "SLIT";
	    else if (r->type == OBJ_REFERENCE) nv = "HOLE";
	    else continue;		// will drop any gap objects
	    write_comment (f, r->dat.precom, '.');
	    if (wl > 1) write_comment (stderr, r->dat.precom, '.');

/* .. Old way of object output ...
	    fprintf (f, "%s %s %9.7f %9.6f %7.3f %7.3f %5.3f",
		    nv, r->dat.name, r->dat.ra/Hour, r->dat.dec/Degree,
		    r->smpos.x, r->smpos.y, r->slit.width);
	    if (wl > 1) fprintf (stderr, "%s %s %9.7f %9.6f %7.3f %7.3f %5.3f",
		    nv, r->dat.name, r->dat.ra/Hour, r->dat.dec/Degree,
		    r->smpos.x, r->smpos.y, r->slit.width);
... */

	    sprintf (rabuf, "%s%c", sexig(r->dat.ra/Hour, 3, 3), AN);
	    sprintf (decbuf, "%s%c", sexig(r->dat.dec/Degree, 3, 2), AN);
	    fprintf (f, "%s %-8s %s", nv, r->dat.name, rabuf);
	    if (wl > 1) fprintf (stderr, "%s %-8s %s",
			nv, r->dat.name, rabuf);
	    fprintf (f, " %12s", decbuf);
	    if (wl > 1) fprintf (stderr, " %12s", decbuf);
	    fprintf (f, " %8.3f %8.3f %5.3f",
		r->smpos.x, r->smpos.y, r->slit.width);
	    if (wl > 1) fprintf (stderr, " %8.3f %8.3f %5.3f",
		r->smpos.x, r->smpos.y, r->slit.width);

/* NOTE - NEW FEATURE -- SCALE THE SLIT DIMENSIONS OF WIDTH AND ALEN,
  BLEN BY SCALE = FOCAL LENGTH * ARCSECOND WHEN WRITING THEM OUT HERE!! */

    /* NOTE -- potential conflict between OBJ_REFERENCE and the hole
      shape code in r->slit.shape -- fix that */
/* NOTE - we may want SMDF positions to be sexigisimal, could handle
  like the POSITION keyword is done above */

    /* IF SHAPE IS RECTANGLE, ALSO WRITE LENGTHS AND ANGLE */
	    alen = r->slit.alen;
	    blen = r->slit.blen;
	    if (r->type == OBJ_REFERENCE) {
		fprintf (f, " %d", r->slit.shape);
		if (wl > 1) fprintf (stderr, " %d", r->slit.shape);
	    } else if (r->type == OBJ_OBJECT) {
// Reduce slit lengths being written by the uncut length parameter
// Intended primarily for nod/shuffle applications.
		alen -= uca;
		blen -= ucb;
// NOTE we may want to modify to -= uca/cos(r->slit.angle) if the angle
// is under some amount like 80 degrees...
		if (alen < 0.0) alen = 0.0;
		if (blen < 0.0) blen = 0.0;
	    }
/* Always write alen, blen, angle for reference as well as slits */
	    fprintf (f, " %5.3f %5.3f %5.1f",
		    alen, blen, r->slit.angle/Degree);
	    if (wl > 1) fprintf (stderr, " %5.3f %5.3f %5.1f",
		    alen, blen, r->slit.angle/Degree);

	    putc ('\n', f);
	    if (wl > 1) putc ('\n', stderr);

	    write_comment (f, r->dat.postcom, '-');
	    if (wl > 1) write_comment (stderr, r->dat.postcom, '-');

	}	// active object conditional
    }	// End of object writing loop

/* Close the file and return */
    fclose (f);
}


	/* - - - - - - - - - - *\
	|			|
	|   Write  .obf  file	|
	|			|
	\* - - - - - - - - - - */


void	write_obf (Obs* obs)
/* Write an object file from the object queue in the obs file; use
  the observation name with .obw extension, and increment the use
  count for all active objects.  */
{
objq	*r;
char	filename[32];
char	rafld[16];	// 13 needed
char	defld[16];	// 13 needed
char	tc;		// for type
int	usec;
FILE	*f;
double	scale;
double	alen, blen;	// to compare defaults properly
char	sbuf[64];	// slit parameter buffer
char	dbuf[64];	// slit default parameters
char	rbuf[64];	// reference default parameters
char	*ds;		// default string pointer
slit	*dfs;		// default slit pointer
int	ld;		// last difference index.
const	char	spform[] = " %.3f %d %.3f %.3f %.2f";

// Go away if null pointer
    if (obs == NULL) return;

// Come up with file name from obs name and .obw
    strncpy (filename, obs->fname, 8);
    filename[8] = AN;	// Just in case
    strcat (filename, ".obw");
// Open that file for writing (no worry about overwrite yet
    if (obs->dlevel >= 0)
	printf (" Writing object file with use counts to %s.\n", filename);
    f = iopen (filename, "w");
    if (f == NULL) return;	// Error already written

/* Get value for scaling slit dimensions... */
    scale = foclen (obs) * ArcSecond;
/* "scale" is size of one arc second in mm. */

// Come up with the default slit parameter strings (arcseconds)
    sprintf (dbuf, spform,
		obs->default_slit.width, obs->default_slit.shape,
		obs->default_slit.alen, obs->default_slit.blen,
		obs->default_slit.angle/Degree );
    sprintf (rbuf, spform,
		obs->reference_slit.width, obs->reference_slit.shape,
		obs->reference_slit.alen, obs->reference_slit.blen,
		obs->reference_slit.angle/Degree );

// Write global comments (were read from object files)
    write_comment (f, obs->comments, '=');

// Make a loop over all the object queue in standard way

    C_loop (r, obs->ob) {
// Find type symbol
	if      (r->type == OBJ_OBJECT)    tc = '@';
	else if (r->type == OBJ_REFERENCE) tc = '*';
	else if (r->type == OBJ_GAP) continue;	// drop virtual ones
	else  tc = '#';   // comment if unknown
// Find ra, dec fields
	sexigw (rafld, r->dat.ra /Hour,   3, 3);
	sexigw (defld, r->dat.dec/Degree, 3, 2);
/* NOTE -- have had problems with defld coming up blank, in a random
  way.  If blank, issue error message and re-do it... */
    if ((defld[0] < ' ') || (defld[1] < '-')) {	// possible blank declination - PANIC
	printf (" *** Dec. output error: %s: \"%s\" !!\n", r->dat.name, defld);
    }
// Find the use count
	usec = r->use;
	if (r->flag & OBJECT_ACTIVE) usec++;

// ADD WRITING OF PRE-COMMENTS HERE
	write_comment (f, r->dat.precom, '.');	// Any pre-comments

// Write the object line
	fprintf (f, "%c%-10s %s  %s  %6.3f  %d",
		tc, r->dat.name, rafld, defld,
		r->dat.priority, usec);

// Write, if needed, any non-default slit parameters for this object;
// fill a buffer with the current slit parameters to compare them...
/* Slit parameters are stored individually in mm, need to divide by
  scale (mm/arcsec) to get arc seconds comparable to defaults. */
	dfs = (r->type == OBJ_REFERENCE) ? &(obs->reference_slit) :
		&(obs->default_slit);
/* Objects with square or circle shapes have their alen, blen set
  to half width upon being read; if this is so, write the default
  values to avoid having the alen, blen put into the output file. */
	alen = r->slit.alen/scale;
	blen = r->slit.blen/scale;
	if (r->slit.shape == SQUARE || r->slit.shape == CIRCLE) {
	    if (alen == 0.5*r->slit.width/scale) alen = dfs->alen;
	    if (blen == 0.5*r->slit.width/scale) blen = dfs->blen;
	}
	sprintf (sbuf, spform,
		r->slit.width/scale, r->slit.shape,
		alen, blen,
		r->slit.angle/Degree );
// for reference objects, need to compare to correct default value alen/blen!!
/* For reference objects -- alen,blen are forced to width/2.  So, if
  shape is square or circle, we need to write for alen and blen the
  appropriate default values directly.
  Therefore, generate alen, blen and put them into above write, right.
*/
// Locate the default string to compare to.
	ds = (r->type == OBJ_REFERENCE) ? rbuf : dbuf;

// Find last difference between sbuf and ds strings...
// This could be a subroutine returning index ld...
	{
int	j;
	    ld = strlen (sbuf);
	    j  = strlen (ds);
// would need error checking on pointers and lengths in subroutine
	    while (ld >= 0) {
		if (j < 0) break;
		if (sbuf[ld] != ds[j]) break;
		ld--; j--;
	    }
// basically returning ld here...
	}

// Replace first space following sbuf[ld] with null, and write sbuf...
	if (ld >= 0) {
	    while (sbuf[ld] != '\0') {
		if (sbuf[ld] == ' ') sbuf[ld] = '\0';
		else ld++;
	    }
	    fprintf (f, "%s", sbuf);
	}
/* 	  else  fprintf (f, "%s!", sbuf);	// debug display...  */

// End the line
	fprintf (f, "\n");

// ADD WRITING OF POST-COMMENTS HERE
	write_comment (f, r->dat.postcom, '-');

    }	// End of object loop.

// At end of loop, close file.
    fclose (f);
}


	/* - - - - - - - - - - - - - - *\
	|				|
	|    Object Queue Utilities	|
	|				|
	\* - - - - - - - - - - - - - - */

objq*	find_object (objq* a, char* name)
/* Find and return named object in queue; NULL if none */
{
objq	*p;	// Pointer to queue

/* First, clear null cases */
    if (a == NULL) return NULL;
    if (name == NULL) return NULL;

/* Standard loop to find things */
    C_loop (p, a) {
/* If we have right one, return it immediately */
//	if (!strcasecmp(p->dat.name, name)) return p;
	if (!strcmp(p->dat.name, name)) return p;
    }

/* Did not find anything, return a NULL */
    return  NULL;
}


/* Former  int kobjs (objq* ob) is totally obsolete; it returned a
  count of objects in list; use  listcount (ob) instead.   */

int	active_objs (objq* a)
/* Return count of active nodes in object queue */
{
objq	*p;
int	n=0;
    C_loop (p, a) {
	if (p->flag & OBJECT_ACTIVE) n++;
    }
    return  n;
}


/* Report on numbers of objects active, slits, alignment;
  and also inactive slits and alignment holes. */
/* May add a unit to report on here, prepare for it... */
void	repobj (objq* ob)
{
FILE*	f = stdout;
int	ks, ka, ko;
int	ns, na, no;
int	kg, ng;
objq	*p;
    ks = ka = ns = na = 0;
    ko = no = 0;
    kg = ng = 0;
    if (ob == NULL) return;
    C_loop (p, ob) {
	if (p->flag & OBJECT_ACTIVE) {	// do not use == here
	    if (p->type == OBJ_OBJECT) ks++;
	    else if (p->type == OBJ_REFERENCE) ka++;
	    else if (p->type == OBJ_GAP) kg++;
	    else ko++;
	} else {
	    if (p->type == OBJ_OBJECT) ns++;
	    else if (p->type == OBJ_REFERENCE) na++;
	    else if (p->type == OBJ_GAP) ng++;
	    else no++;
	}
    }
// Not yet reporting ng and kg values for gaps; would add
// a column to the table if we really wanted to do this.
/* Gathered numbers, now put them in a nice table... */
    fprintf (f, " Objects - Slits  Alignment Other\n");
    fprintf (f, " On Mask:  %4d    %4d    %4d\n", ks, ka, ko+kg);
    fprintf (f, " Removed:  %4d    %4d    %4d\n", ns, na, no+ng);
    fprintf (f, "  Totals:  %4d    %4d    %4d\n", ks+ns, ka+na, ko+no+kg+ng);
}



	/* - - - - - - - - - - - - - - - - - - *\
	|					|
	|   Reordering  (traveling salesman) 	|
	|					|
	\* - - - - - - - - - - - - - - - - - - */

/*  ===  Object queue ordering  ===  */

#define  LongSlit  40.0
#define  SlitGap  0.500
// Above follow parameters in maskcut...
// they appear in cututil.c ...
// There, they are variables, however.

/* Symbol NNA definition moved to mgfeats.h file */

static	void	crect_ends (vect2 p, vect2 a, double w, vect2* s, vect2* e)
/* Find start s and end e of cutting for a rectangle with starting
  end at p and centerline a, with width w perpendicular to a. */
{
double	x, y;	// Lengths to move
vect2	ux, uy;	// Unit vectors of rectangle
vect2	r;	// Reference location

/* Compute the corner and length and width vectors, using the beam
  width; also the turn-on point.  */
    x = vect2norm (a); // - bwid;
    if (x < 0.0) x = 0.0;		// Movement length


/* Test to see if the slit is "long", i.e. exceeds a pre-determined
  maximum length specified in parameter LongSlit defined above.  If it
  is, we recursively call ourselves to cut two slits, separated by
  a gap of length SlitGap defined above.  */
    if (x > (LongSlit) ) {
/* The recursive slit splitting algorithm */
/* NOTE for end point things --
  We call with a NULL pointer for start when recursively cutting, so
  don't store start if null.  Store start point as e also.  */
double	gap;
double	hx;
vect2	half;
vect2	p2;
// Find new x length of each slit part, and compute new
// slit length (just under a/2) and second slit start location.
	gap = SlitGap; // * SF;		// Actual length of gap
	hx = (x - gap) / 2.0;		// Length of new slit part
	half = mul2vect (a, hx/x);	// Centerline for each part
//	z = crect (p, half, w, hp, fr, u);	// Cut first half
	crect_ends (p, half, w, s, e);
	p2 = sum2vect (p, mul2vect (a, (x-hx)/x) );	// Location of part 2.
//	z = crect (p2, half, w, z, fr, u);	// Cut second half
	crect_ends (p2, half, w, NULL, e);
//	return  z;
	return;
    }
/* NOTE that the above fix for excessively long slits produces a "slit"
  which is not precisely defined by the SMDF file, in that it has gap(s)
  along its length...  */

    y = w;  // - bwid;
    if (y < 0.0) y = 0.0;		// Movement width
    ux = norm2vect (a);
//    uy = make2vect (-ux.y, ux.x);	// Rotated left 90 deg.
    uy = make2vect (ux.y, -ux.x);	// Rotated right 90 deg.
// rotation in reverse due to reflection in cutting...
// Reference point is lower left corner; it is displaced from
// start p by bwid in +ux, by y/2 in -uy.
    r = p;  // sum2vect (p, mul2vect (ux, bwid*0.5));
    r = sub2vect (r, mul2vect (uy, y*0.5));

// r is our start and end point...
    if (s != NULL) *s = r;
    *e = r;

}

static	void	endpoints (objq* a)
/* Compute and store the slit end points for cutting in the
  object queue structure. */
{
vect2	c;	// object location
vect2	start;	// location of start corner, relative to c
vect2	slth;	// length
//vect2	wide;	// width
double	rs, rc;	// rotation sine, cosine
double	la, lb, w;	// corrected length, width
double	ls;	// length of slit

// Follows the methods of the maskcut program.
	c = a->smpos;	// object

	sincos (a->slit.angle, &rs, &rc);

// w is width/2 less tool, not less than zero
// la is length a, less tool, not less than zero
// lb is like la.
// start is -la, -w rotated by rc,rs
// slth is la+lb,0 rotated from x by rc,rs
// wide is 0,2w, rotated from y by rc,rs
	w = a->slit.width;
	if (w < 0.0) w = 0.0;

// Find the slit length unit vector
	slth = make2vect (rc,rs);

// We do not consider the tool size in this program (yet?)

	switch (a->slit.shape)  {
	case  RECTANGLE:
	    la = a->slit.alen;  // - tool;
	    lb = a->slit.blen;  // - tool;
/* la is movement left from center, lb is movement right */
	    ls = la + lb;
	    if ( ls <= 0.0) {
		if (la < 0.0) la = 0.0;
		if (lb < 0.0) lb = 0.0;
		ls = la + lb;
	    }
// Find the slit start position, including beam size
	    start = sub2vect (c, mul2vect (slth, la));  // la+tool));
// Get right length for the slit length vector including beam
	    slth = mul2vect (slth, ls);  // ls+bwid);
/* Plot the rectangle using these parameters */
//	    zlast = crect (start, slth, w, zlast, slow, f);
	    crect_ends (start, slth, w, &(a->smp1), &(a->smp2) );

	    break;
	case  CIRCLE:
// //	    ls = w/2.0 - tool;
// 	    ls = w;  // - bwid;
// //	    zlast = chole (c, ls, zlast, vslow, f);
// NOTE - REALLY SHOULD USE STARTING POINT NEAR ONE SIDE?
	    a->smp1 = a->smp2 = c;
	    break;
	default:
	case  SQUARE:
	    la = lb = w/2.0; // - tool;
	    ls = la + lb;
// Find the slit start position, including beam size
	    start = sub2vect (c, mul2vect (slth, la));  // +tool));
// Get right length for the slit length vector including beam
	    slth = mul2vect (slth, ls);  // +bwid);
/* Plot the rectangle using these parameters */
//	    zlast = crect (start, slth, w, zlast, slow, f);
	    crect_ends (start, slth, w, &(a->smp1), &(a->smp2) );
	    break;
// ADD HERE OTHER SPECIAL SHAPES...
// E.G. THE CIRCLE AND 4 SLITS CASE...
//	case  SPACED_CROSS:
//	    zlast = chole (c, w-bwid, zlast, vslow, f);
//// wide is c to slit displacement
//// slth is slit length
//	    wide = mul2vect (slth, w*0.5 + SF*r->slit.blen);
//	    slth = mul2vect (slth, SF*r->slit.alen);
//	    for (j=0; j<4; j++) {
//		start = sum2vect (c, wide);
//		zlast = crect (start, slth, w, zlast, slow, f);
//// Rotate wide and slth left by 90 degrees...
//		wide = make2vect (-wide.y, wide.x);
//		slth = make2vect (-slth.y, slth.x);
//	    }
	}	// End of switch


}

double	smdist (objq* a, objq* b)
/* Return distance between the two objects */
{
double	d;
//    d = vect2norm ( sub2vect (a->smpos, b->smpos) );
    d = vect2norm ( sub2vect (a->smp2, b->smp1) );
    return  d;
}

double	totdist (objq* ob)
/* Compute the total distance in the object queue */
{
double	d = 0.0;
objq	*p, *q;
    if (ob == NULL) return 0.0;
    C_loop (p, ob) {
	q = p->next;
	d += smdist (p, q);
    }
    return  d;
}


static	void	revert (objq* ob, int n)
/* Reverse order of the next n nodes after node ob */
{
int	k;
objq	*p, *q, *r;

    if (ob == NULL) return;
    if (n == 0) return;

/* We assume n is positive.  Later, allow negative with
a similar set of code. */
    if (n < 2) return;		// Null case, effectively

/* Pop off n nodes, and put them in a queue in reverse order */
    r = NULL;
    for (k=0; k<n; k++) {
	p = ob->next;
	q = pop_obj (p);	// pop from left of a list; forget q
	if (q != NULL) q = NULL;  // nonsense for lint
	r = push_obj (p, r);	// Push onto left of new list
	if (p != NULL) r = p;	// Remember new leftmost
    }

/* Push the new list back just after ob (left of ob->next) */
    p = ob->next;
    r = push_obj (r, p);
}

objq*	reorder (objq* ob, int kn, int dblev)
/* Re-order object list to minimize travel distance between objects
  on the mask.  Return pointer to object following the largest
  distance between objects.  Argument kn is size of maximum substring
  to search for possible inversion. */
/* dblev is debug level of this call.  Needs to be over 1 to report */
{
int	n = 0;
int	j;
int	m;
//int	bad;
int	changed;
int	change2;
double	ab, cd;
double	ac, bd;
double	rc;	// reversion cost
double	dr, di;
double	big;
double	befo, aftr;
objq	*p, *q;
objq	*a, *b, *c, *d;
//int  kk=0;	// debug only

// debug things
int	snc=0;
int	tnc=0;	// total number of changes

/* We assume that the object list is circular.  If it is not, we take
  the liberty of joining its two ends to each other. */

/* Take care of any null cases in the usual way... */
    if (ob == NULL) return ob;

/* Scan the list and count its members.  Correct any "ends" */
    G_circ (ob);
    n = listcount (ob);

/* If we have fewer than 4 nodes, we can't do the algorithm... */
    if (n < 4) return ob;

/* If kn is less than 2, the algorithm is effectively null */
    if (kn < 2) return ob;

/* Limit kn based on n... */
#ifndef  RLSFULL
    m = n/2;
    if (kn > m) kn = m;
#endif
// If inverse order has different length, we must allow this.

/* Set up the to and from slew locations in the list. */
    C_loop (p, ob) {
// Find start and end points smp1,smp2 much like maskcut does it...
//#ifdef  RLSFULL
//	endpoints (p);
//#else
//	p->smp1 = p->smp2 = p->smpos;
//#endif
	endpoints (p);
#ifndef  RLSFULL
	p->smp2 = p->smp1;
#endif
    }

/* The before picture; find total length of list... */
    befo = totdist (ob);
if (dblev+DEBUG_UTL > 0) {
    printf ("  Start reorder with %d nodes and total distance of %5.2f mm.\n",
	n, befo);
    fflush (stdout);		// debug only
}
// if (befo > 1.e6) return ob;	// debug exit

#ifdef  NNA
/* The nearest neighbor pre-qualify algorithm */

// Find nearest neighbor of each node
    C_loop (p, ob) {
double	x, s, z;
	z = s = 2.0e20;
	p->nn = p->nx = NULL;
	C_loop (q, p) {
	    if (q == p) continue;
	    x = smdist (p, q);
	    if (x < z) {	// new smallest or next smallest
/* s is smallest known distance (nn), z is next smallest (nx) */
		if (x < s) {	// Replace s by x and z by s
		    z = s;
		    p->nx = p->nn;
		    s = x;
		    p->nn = q;
		} else {	// Replace z by x
		    z = x;
		    p->nx = q;
		}
	    }
	}
    }
//    printf ("Start NNA algorithm.\n");
//    fflush (stdout);

/*  NOTE -- Suggestion to enhance this algorithm.
  The test case for very large masks had pairs of nearby slits, so
  each node would have a very near neighbor, and these were usually
  already located adjacent.  So this had little effect with very
  few neighbor nodes being moved.  An enhancement, while still general,
  would be to track also the second nearest neighbor, and if the
  nearest is adjacent to the node, check to see if the second nearest
  is also adjacent, and if not to move it to the side on which the
  nearest neighbor is not situated.  This should greatly change the
  node list order in a relatively few passes, and perhaps do the speed
  up of convergence which we want but did not achieve.  */

// Loop to process this some number of times
    m = 2 + (int)(log((double)n));
    for (tnc=0; tnc<m; tnc++) {
// Clear the flags on all nodes
	C_loop (p, ob) p->flag &= ~OBJECT_FLAG;

// Loop to insert nearest neighbor if not already done
	changed = 1;
	change2 = 0;
	j = 0;
	while (changed) {
	    changed = 0;
	    C_loop (p, ob) {
		if (p->flag & OBJECT_FLAG) continue;
/* Skip any flagged node, and flag this node */
		p->flag |= OBJECT_FLAG;
		if (p->next == p->nn) {
/* nearest is next, put next nearest in last location */
		    if (p->last == p->nx) continue;
		    q = p->nx;
		    d = pop_obj (q);
		    p = push_obj (q, p);
		    changed++;
		} else if (p->last == p->nn) {
/* nearest is last, put next nearest in next location */
		    if (p->next == p->nx) continue;
		    q = p->nx;
		    d = pop_obj (q);
		    q = push_obj (p, q);
		    changed++;
		} else {
/* nearest not adjacent, put it just before p */
		    q = p->nn;
		    d = pop_obj (q);
		    p = push_obj (q, p);
		    changed++;
		}
	    }
	    j++;
	    change2 += changed;
	}
if (dblev+DEBUG_UTL > 1) {
	printf (" NNA %d in %d loops changed %d; %5.2f\n",
		tnc+1, j, change2, totdist(ob) );
	fflush (stdout);
}
    }  // total loop on tnc
    tnc = 0;		// reset the count
#endif  // on NNA

/* Outer loop for various change strategies */
    changed = change2 = 1;
    while (changed) {		// We leave when we don't change anything
	changed = 0;	// Increment this on any change...
/* Loop around the list looking for things to switch.  End when we
  have gone around without finding anything. */

/* Within the loop, we shall remember these things:
	The last pivot node, or where we started
	The current nodes and their distances  .. */

	a = b = c = d = NULL;
/* It appears that doing the interchange loop at the start increases
  the rate of convergence, so we drop that debug stuff below. */
	for (j=0,a=ob; j <= n && changed < n; j++) {
/* In this loop, as a debug thing, we test on tnc, so that we only do
  this interchange after some changes by moving nodes... */
//	for (j=0,a=ob; j <= n && tnc > 0 && changed < n; j++) {
//	for (j=0,a=ob; j <= n && tnc > 0; j++) {
	    if (b == NULL) {	// (re-)initialize
		b = a->next;
		ab = smdist (a, b);
		q = a;		// the deja-vu node
	    }
	    c = b->next;
	    rc = smdist (c,b) - smdist (b,c);  // increase for inversion

//#ifdef  NNA
#ifndef  DOITTHEBADWAY
	    for (m=2; m<= kn; m++) {	// The inner loop
#else
	    for (m=2; m<= kn && m <= tnc+2; m++) {	// The inner loop
#endif  // on NNA

		d = c->next;
/* Here we test the next several strings for reversal */
		cd = smdist (c, d);
		ac = smdist (a, c);
		bd = smdist (b, d);
//		if (ac+bd < ab + cd) {	// Found reversion
		if ((ac+bd + rc) < (ab + cd)) {	// directed reversion
		    revert (a, m);
		    j = -1;		// Will increment at bottom
		    changed++;
		    q = a;		// Where we last reverted
		    b = a->next;	// should be c
		    ab = smdist (a, b);	// Should be = ac
		    rc = -rc;
		    m=0; break;	// Using m as a signal for reversion
		}
		rc += smdist (d,c) - smdist (c,d);
		c = d;	// Next longer substring
	    }
	    if (m == 0) continue;	// Did reversion, test again

    /* Move to the next item... */
	    a = a->next;
	    if (a == q) break;	// been there; location tracking
	    b = a->next;
	    ab = smdist (a, b);
	}	// End basic loop for pivot node (a)

/* See if the second loop has been run before with no changes
  having taken place... */
	if (changed + change2 == 0) break;
	change2 = 0;

/* Second change method.  For every node, see if it could be placed
  in a different place to reduce distance... */
	for (j=0,p=ob; j<=n && change2 < n; p = p->next,j++) {
//	for (j=0,p=ob; j<=n; p = p->next,j++) {   // look at all pivot nodes
// save distances which result if p is popped from a-p-b
	    a = p->last;
	    b = p->next;
	    ab = smdist (a, b);
	    dr = smdist (a, p) + smdist (p, b);	// removed
// Inner loop searches locations into which to put p, c-p-d
	    for (c=b; c != a; c = c->next) {
// Find distances resulting from p being pushed
		d = c->next;
		if (d == p) break;
		cd = smdist (c, d);
		di  = smdist (c, p) + smdist (p, d);
// See if improvement results.
// The old distance is dr and cd; the new would be ab and di
		if (dr+cd > ab+di) {   // Improvement exists
		    change2++;
		    j = -1;
		    q = pop_obj (p);
		    p = push_obj (p, d);
		    p = q;	// for outer loop location
#ifdef  NNA
		    a = p->last;
		    b = p->next;
		    if (b == ob) break;
		    ab = smdist (a, b);
		    dr = smdist (a, p) + smdist (p, b);	// removed
		    c = b;
#else
		    break;	// exit the inner loop
#endif
		}
	    }  // End inner loop on c coming around to a
// Outer loop stops when j counts out; p going around once with no change
	}  // End outer loop on p counted by j

	changed += change2;	// Now means ANY change
	if (changed) {		// debug output of outer loop
	    tnc++;
if (dblev+DEBUG_UTL > 1) {
	printf ("  Iteration %d made %d changes, %d in secondary, %5.2f\n",
	    tnc, changed, change2, totdist(ob) );
}
	}
	snc += changed;

    }	// End loop on change methods

/* At the end of the loop, we return the right node of the longest
  distance between adjacent nodes as the new start of list. */
    for (big=0.0,j=0,b=q=ob; j <= n; j++) {
	p = q->next;
	ab = smdist (q, p);
	if (ab > big) {
	    b = p;
	    big = ab;
	}
	q = p;
    }

/* The "after" picture; find new length, and how much we improved */
if (dblev+DEBUG_UTL > 0) {
    aftr = totdist (b);
    printf ("  End reorder with new distance %5.2f mm.\n", aftr );
//    printf ("  Used %d switches, total %d nodes.\n", kr, nn);
    printf ("  In %d iterations, made %d changes.\n", tnc, snc);
    ac = befo - aftr;
    bd = 60.0 * ac / (50.0 * Inch);	// Nominal slew rate
    printf ("  Saved %5.2f mm., or %5.1f seconds.", ac, bd);
    dr = 100.0 * ac / befo;
    printf ("  Improvement %5.2f%%.\n", dr);
    fflush (stdout);
// all debug only print above...
}

    return  b;
}


	/* - - - - - - - - - - *\
	|			|
	|    SMDF  reading	|
	|			|
	\* - - - - - - - - - - */

/*  ..--== NOTE:
  smdf reading program was previously an include item, it is put
 here to get it into the library.  It is not needed by mask generation,
 but is needed by many other applications which use this library.
 Notably, it is needed by maskcut.c and by smfix.c among others.  */
/* From smdfutil.c, now in mgutils.c */

/*  =====  Storage utilities  =====  */

/* kill_litlist was used only by read_smdf; now obsolete. */

static	element*  grab_element (element* gs, char* name)
/* Obtain element from gs queue, or if no gs queue, create an
  element structure to hold the name instead.  Used only in
  the reading of SMDF with a NULL gs value possible.  */
{
element	*e;

/* Some programs (such as maskcut) call read_smdf with a NULL value
  for the optical queue pointer.  They only need to use the name
  of the elements which are present in the file, so we use an otherwise
  NULL element queue item to store that name when gs is NULL, but
  provide a pointer to the real thign if gs is present.  */

    if (gs == NULL) e = newelement       (name);
    else            e = find_element (gs, name);
    return  e;
}

/*  =====  SMDF Reading subroutines  =====  */

#define  SMDFREC  127

/* SMDF reading subroutine */
Obs*	read_smdf (char* filename, element *gs, int special)
/* Argument is file name; return pointer to filled obs structure */
{
FILE*	in;
Obs	*obs;
objq	*ob = NULL;	/* Current or last object */
// char	*a;
char	buffer[SMDFREC+1];
char	tbuf[SMDFREC+1];
char	*key;
char	*rest;
char	*st;	/* scratch */
char	*fspec=NULL;
int	l;
int	sm;	/* Sub-Mask, Not yet used... */
int	k;	/* scratch */
//int	shape;	// using k to read this
namlist*	comq = NULL;
namlist*	preq = NULL;
namlist*	postq = NULL;
double	val;
int	yr, mo, dy;
char	*c;
char	buf[32];
//vect3	r, a, v;


/* Open the given file; error if none */
    findfile (&fspec, "r", filename, NULL);
    sprintf (buffer, "%s.SMF", filename);
    findfile (&fspec, "r", buffer, NULL);
    sprintf (buffer, "%s.smf", filename);
    findfile (&fspec, "r", buffer, NULL);

    if (fspec == NULL) {
	printf (" *Can't open file %s\n", filename);
	return  NULL;
    }
    in = iopen (fspec, "r");
    free (fspec);		// No longer need that
    if (in == NULL) {
	return  NULL;
    }

// Output message of processing the input file
//    printf ("Processing \"%s\".\n", filename);
//    if (special) printf (" -- Special input\n");	// debug only

/* Get a default obs structure */
    obs = defobs (gs);

/* Loop to read all records in the file */
    while ( fgets (buffer, SMDFREC, in) != NULL) {
/* Clean out terminating newline from buffer */
	trimend (buffer, NEWLINES);
	l = strlen (buffer);
	if (l < 1) continue;	/* Drop null lines */

/* Test for and deal with comment records */
	if (*buffer == '!' || *buffer == '#') {	/* comment */
/* Separate pre-post-global comments into separate queues */
	    if (buffer[1] == '.') {
		preq = pushtxt (buffer+2, preq);
	    } else if (buffer[1] == '-') {
		postq = pushtxt (buffer+2, postq);
	    } else if (buffer[1] == '=') {
		comq = pushtxt (buffer+2, comq);
	    } else {
/* Ignore me */
	    }
	    continue;		/* Otherwise ignore the comment */
	}

// Copy the buffer before destructive parsing occurs;
// this will be used if needed in special parsing of hole/slit
	strcpy (tbuf, buffer);

/* If special input, we can have a leading blank which indicates
  a special record indicating an inverse transform of an object. */
	if (special && (*buffer == ' ')) {	// special edit record
/* Fill in post comment queue for last object */
	    if (ob == NULL) obs->comments = push_name (postq, obs->comments);
	    else ob->dat.postcom = postq;
	    postq = NULL;

// Read the things we need for the object;
// flag the object as special
// obtain ra/dec using the inverse transform?

// this record is like an object record except:
// no hole/slit
// name exitst
// no ra/dec field
// other fields present.

/* Put the special object in the prequeue.. */
	    preq = pushtxt (buffer+1, preq);

/* Special object:  name, (nora/dec) x, y, width, [shape?] alen, blen, angle */
// non-destructive parsing of the components...
	    ob = xalloc (objq);
	    ob->next = ob->last = ob;
	    rest = parse (buffer, buf, WHITESPACE);
	    ob->dat.name  = dynamstr (buf);
// printf (" --Found special object %s\n", buf);	// debug only
// ra, dec will be obtained by inversion shortly
	    c = parse (rest, buf, WHITESPACE);
	    ob->smpos.x  = floatvalue (buf);
	    c = parse (c, buf, WHITESPACE);
	    ob->smpos.y  = floatvalue (buf);

	    unmaskvect (obs, ob->smpos, &(ob->dat.ra), &(ob->dat.dec) );

// Now, continue with width [shape?] alen, blen, angle
	    c = parse (c, buf, WHITESPACE);
	    ob->slit.width = floatvalue (buf);
	    c = parse (c, buf, WHITESPACE);

// Now, if buf is integer, it is the shape; otherwise it is
// the alen part of a slit, so we set shape appropriately.
// decode using sscanf on int, and if found is int...
//	    k = sscanf (buf, "%d", &shape);
// that above won't make it; 13.55 is decoded as 13 with no error!
// Search instead for a decimal point in buf...
	    if (strchr (buf, '.')) {	// it is alen...
		ob->slit.shape = RECTANGLE;
		ob->type = OBJ_OBJECT;
	    } else {
		ob->slit.shape = intvalue (buf);
		c = parse (c, buf, WHITESPACE);	// decode some more (alen)
		ob->type = OBJ_REFERENCE;
	    }

// Continue on with decoding alen, blen, angle
//	    c = parse (c, buf, WHITESPACE);
	    ob->slit.alen  = floatvalue (buf);
	    c = parse (c, buf, WHITESPACE);
	    ob->slit.blen  = floatvalue (buf);
// Those values are read as mm, no scaling is needed here...
	    c = parse (c, buf, WHITESPACE);
	    ob->slit.angle = floatvalue (buf) * Degree;

/* Put in the pre-queue and re-initialize it. */
	    ob->dat.precom = preq;
	    preq = NULL;
	    ob->val = NULL;

/* Put the object into the queue */
	    ob->flag = OBJECT_ACTIVE | OBJECT_SPECIAL;
	    obs->ob = push_obj (ob, obs->ob);

	    continue;	// Omit the decoding of the buffer again
	}	// end special edit record (leading blank)

/* Data record, parse the keyword and rest */
	key = strtok (buffer, WHITESPACE);	/* keyword */
	rest = strtok (NULL, "");
	if (key == NULL) continue;	/* ? may want to debug print? */
	if (rest == NULL) {
	    l = strlen (key);
	    rest = key + l;
	}

/* Trim whitespace from start and end of rest */
	trims (rest, WHITESPACE);

/* We uppercase the "key" string at this time,
  so the compares will be non-case-sensitive */
	for (st=key; *st != '\0'; st++) *st = toupper(*st);

/* Decision tree on keyword... */
	if        (strcmp (key, "OBSERVER") == 0) {
/* OBSERVER  put 2-16 characters into observer name  */
	    strncpy (obs->oname, strtok(rest, WHITESPACE), 16);
/* -NOTE- need to enforce 2 char. minimum! */
	} else if (strcmp (key, "NAME") == 0) {
	    strncpy (obs->fname, rest, 8);
	} else if (strcmp (key, "TITLE") == 0) {
	    if (obs->title != NULL) free (obs->title);
	    obs->title = dynamstr (rest);
	} else if (strcmp (key, "MADE") == 0) {
// The MADE entry is taken as commentary
	    ;
	} else if (strcmp (key, "TELESCOPE") == 0) {
/* focal-length, curvature, reflections, name */
// Get fl, curv, reflections, and name here...
	    strtok(rest, WHITESPACE);  // focal length
	    strtok(NULL, WHITESPACE);  // curvature
	    strtok(NULL, WHITESPACE);  // reflections
	    strcpy (buf, strtok (NULL, WHITESPACE) );  // name
	    if (!strcmp(buf, "Mag2NoADC")) strcpy (buf, "Magellan2");
/* The "maskcut" program reads a SMDF file without an optical element
  queue having been constructed; it is (so far) not needed in that
  application.  The name of the telescope is needed, so it is put into
  an otherwise blank element structure for the telescope.  */
	    obs->Tel_scope = grab_element (gs, buf);
	    if (obs->Tel_scope == NULL)
		printf (" **Telescope \"%s\" not found.\n", buf);

	} else if (strcmp (key, "INSTRUMENT") == 0) {
/* primary name, secondary name */
/* AS YET, NO PLACE TO STORE SECONDARY NAME */
// Treat instrument name same as telescope name; dummy element if the
// element queue is not present.
/* The "maskcut" program also uses instrument name, at least in
  some debug output... */
	    strcpy (buf, strtok (rest, WHITESPACE) );
	    obs->Instr_mt = grab_element (gs, buf);

	} else if (strcmp (key, "GRATING") == 0) {
/* name of grating - angle - order - lpmm */
/* NOTE that this is obsolete, and we should store disperser stuff */
// change to set disperser by the name, store order,
//	    printf (" **Obsolete GRATING item in SMDF input.\n");
	    strcpy (buf, strtok (rest, WHITESPACE) );
	    obs->Disperser = find_element (gs, buf);
	    obs->D_angle = floatvalue (strtok (NULL, WHITESPACE))*Degree;
	    obs->order   = intvalue   (strtok (NULL, WHITESPACE));
// We forget about the density in the input -- it is redundant.

	} else if (strcmp (key, "DISPERSER") == 0) {
/* disperser element -- name, order, grating angle  */
	    strcpy (buf, strtok (rest, WHITESPACE) );
	    obs->Disperser = find_element (gs, buf);
	    obs->order   = intvalue   (strtok (NULL, WHITESPACE));
	    obs->D_angle = floatvalue (strtok (NULL, WHITESPACE))*Degree;

	} else if (strcmp (key, "POSITION") == 0) {
/* sub-mask, ra (sexigis), dec (secigism), equinox, p.a. degrees,
  field center position on mask x mm, y mm. */
	    sm = intvalue (strtok (rest, WHITESPACE));
	    if (sm > 100) k = sm;	// Nonsense for lint...
/* --NOTE-- submask is not currently used. */
	    obs->ra  = posvalue (strtok (NULL, WHITESPACE)) * Hour;
	    obs->dec = posvalue (strtok (NULL, WHITESPACE)) * Degree;
	    obs->equinox = floatvalue (strtok (NULL, WHITESPACE));
	    obs->angle = floatvalue (strtok (NULL, WHITESPACE)) * Degree;
/* --NOTE-- we have also x,y field center and its not stored. */

/* construct transform matrix too */
	    get_mpos (obs);

	} else if (strcmp (key, "GUIDER") == 0) {
/* sub-mask, guider #, angular, radial, star X, star Y */
		;
	} else if (strcmp (key, "DATEOBS") == 0) {
/* year, month, day, hour angle sexigisimal */
	    yr = intvalue (strtok (rest, WHITESPACE));
	    mo = intvalue (strtok (NULL, WHITESPACE));
	    dy = intvalue (strtok (NULL, WHITESPACE));
/* --NOTE-- hour angle is not present in obs struvcture */
	    obs->epoch = Jyear (jdate (yr, mo, dy));
/* -NOTE- Add time of day or hour angle. this jyear is for 0 gmt on
  the given date; we want to add position and local sidereal time
  to obtain time of day in gmt for transit of position */
/* ==NOTE== for this, we also need a local sidereal time subroutine,
  given date and longitude.  Thus we need for telescope, by name,
  the location in longitude and latitude! */

	} else if (strcmp (key, "FILTER") == 0) {
/* name, blue limit, red limit */
/* NOTE -- should have standard filter names in optic data; so we
  should just read the name.  However, we also read the limits. */
		;
	} else if (strcmp (key, "WAVELENGTH") == 0) {
/* center wavelength in angstroms; defines grating angle */
/* -NOTE- should define grating angle, not center wavelength */
/* for a while, this is commentary */
/* -- need to put it into obs anyway... */
	    sscanf (rest, "%lf", &(obs->cw));
	} else if (strcmp (key, "WLIMIT") == 0) {
/* The passband, needed for inverse mapping stuff */
	    sscanf (rest, "%lf %lf", &(obs->pb.blue), &(obs->pb.red));
	} else if (strcmp (key, "DLIMIT") == 0) {
	    sscanf (rest, "%lf %lf", &(obs->pb2.blue), &(obs->pb2.red));
	    obs->pb2act = 1;
	} else if (strcmp (key, "TEMPERATURE") == 0) {
/* Expected temperature in degrees */
//	    sscanf (rest, "%lf", &(obs->temp));
	    sscanf (rest, "%lf", &val);
/* Limit temperature to a resonable range here... */
	    if (fabs(val) <= 40.0) obs->temp = val;
	} else if (strcmp (key, "EPOCH") == 0) {
/* Epoch years -- alternative to dateobs */
/* -NOTE- epoch and dateobs are same thing ?? */
	    sscanf (rest, "%lf", &(obs->epoch));
	} else if (strcmp (key, "HANGLE") == 0) {
	    if ( (sscanf (rest, "%lf", &(obs->ha)) != 1)) obs->ha = -25.0;
//	} else if (strcmp (key, "XXX") == 0) {


/* Handle the object keywords here -- these add entries to object queue */
	} else if (strcmp (key, "SLIT") == 0) {
/* Fill in post comment queue for last object */
	    if (ob == NULL) obs->comments = push_name (postq, obs->comments);
	    else ob->dat.postcom = postq;
	    postq = NULL;

/* What's supposed to happen:
	if special is 1, and if the cannonical ra/dec position does
	not give the x,y position to within 0.01 arc sec, we will
	need to set the special flag OBJECT_SPECIAL and also add
	the buffer to the pre-comment queue.
	Here, we print the .buffer stuff to tbuf.
	Below, before preq is put into precom, we test for
	the correct position.  But only if special is 1.
.. */

/* object, ra, dec, x, y, width, length -, length +, angle, magnitude */
/* Assign any pre/post comment queues */
	    ob = xalloc (objq);
	    ob->next = ob->last = ob;
	    ob->flag = 0;
	    ob->dat.name   = dynamstr (strtok (rest, WHITESPACE));
// NOTE -- HERE AND IN HOLE CASE BELOW, IF SPECIAL IS 1, WE NEED
// TO CHECK WHETHER THE FOLLOWING FIELD IS ACTUALLY RA AND DEC, OR IF
// THEY MIGHT BE MISSING.  IF MISSING, SUBSTITUTE FIELD CENTER FOR
// THE SAKE OF COMPUTATIONS, AND SET SPECIAL FLAG.
// IF NOT MISSING, DECODE THE RA AND DEC.
// OTHERWISE, (SPECIAL NOT 1) WE DECODE RA/DEC AS BEFORE.
	    if (special == 1) {
		st = strtok (NULL, WHITESPACE);
		if (strchr (st, ':')) {	// Have ra/dec
		    ob->dat.ra  = posvalue (st) * Hour;
		    ob->dat.dec = posvalue (strtok (NULL, WHITESPACE))*Degree;
		    ob->smpos.x = floatvalue (strtok (NULL, WHITESPACE));
		} else {		// ra/dec missing
		    ob->dat.ra  = obs->ra;
		    ob->dat.dec = obs->dec;
		    ob->flag |= OBJECT_SPECIAL;
		    ob->smpos.x    = floatvalue (st);
// and force precomment
		}
	    } else {
	    ob->dat.ra     = posvalue (strtok (NULL, WHITESPACE)) * Hour;
	    ob->dat.dec    = posvalue (strtok (NULL, WHITESPACE)) * Degree;
	    ob->smpos.x    = floatvalue (strtok (NULL, WHITESPACE));
	    }
	    ob->smpos.y    = floatvalue (strtok (NULL, WHITESPACE));
	    ob->slit.width = floatvalue (strtok (NULL, WHITESPACE));
	    ob->slit.alen  = floatvalue (strtok (NULL, WHITESPACE));
	    ob->slit.blen  = floatvalue (strtok (NULL, WHITESPACE));
// Those values are read as mm, no scaling is needed here...
	    ob->slit.angle = floatvalue (strtok (NULL, WHITESPACE)) * Degree;
	    ob->slit.shape = RECTANGLE;  /* from name slit */
// if special, test position and see if
// we need to flag the object and remember buffer...
	    if (special == 1) {
vect2	d;
/* Find difference between position from ra/dec and from record */
		d = sub2vect (ob->smpos,
				maskvect (obs, ob->dat.ra, ob->dat.dec));
		if (vect2norm(d) > 0.003) {	// wrong position
// The .003 mm is about .01 arc seconds for Magellan
// see if it is close to the record's position
// if not:
//   find true ra/dec -- no
//   put those in record -- no
//   put tbuf into precomment queue
		    preq = pushtxt (tbuf, preq);
//   mark it special
		    ob->flag |= OBJECT_SPECIAL;
		}
	    }
	    ob->dat.precom = preq;
	    preq = NULL;
	    ob->val = NULL;

/* -- NOTE -- we don't yet have these quantities:
	priority/magnitude is not properly used in structure
	wavelength limits in object.dat
	slit shape code is conflict
	type and flag are not done
.. */

/* Put the object into the queue */
	    ob->type = OBJ_OBJECT;
	    ob->flag |= OBJECT_ACTIVE;
	    obs->ob = push_obj (ob, obs->ob);

	} else if (strcmp (key, "HOLE") == 0) {
/* Fill in post comment queue for last object */
	    if (ob == NULL) obs->comments = push_name (postq, obs->comments);
	    else ob->dat.postcom = postq;
	    postq = NULL;

/* object, ra, dec, x, y, size, shape, (angle?,) magnitude */
/* Assign any pre/post comment queues */
	    ob = xalloc (objq);
	    ob->next = ob->last = ob;
	    ob->flag = 0;
	    ob->dat.name   = dynamstr (strtok (rest, WHITESPACE));
	    if (special == 1) {		// Test next token for colon
		st = strtok (NULL, WHITESPACE);
		if (strchr (st, ':')) {	// Have ra/dec
		    ob->dat.ra  = posvalue (st) * Hour;
		    ob->dat.dec = posvalue (strtok (NULL, WHITESPACE))*Degree;
		    ob->smpos.x = floatvalue (strtok (NULL, WHITESPACE));
		} else {		// ra/dec missing
		    ob->dat.ra  = obs->ra;
		    ob->dat.dec = obs->dec;
		    ob->flag |= OBJECT_SPECIAL;
		    ob->smpos.x    = floatvalue (st);
		}
	    } else {
	    ob->dat.ra     = posvalue (strtok (NULL, WHITESPACE)) * Hour;
	    ob->dat.dec    = posvalue (strtok (NULL, WHITESPACE)) * Degree;
	    ob->smpos.x    = floatvalue (strtok (NULL, WHITESPACE));
	    }
	    ob->smpos.y    = floatvalue (strtok (NULL, WHITESPACE));
	    ob->slit.width = floatvalue (strtok (NULL, WHITESPACE));
// SET K DIRECTLY INTO SLIT.SHAPE
// KDC FIX
	    k = intvalue (strtok (NULL, WHITESPACE));
//	    ob->slit.shape = (k == 1) ? SQUARE : CIRCLE;
	    ob->slit.shape = k;
// ADD DECODE OF ALEN, BLEN, ANGLE AS ABOVE
	    ob->slit.alen  = floatvalue (strtok (NULL, WHITESPACE));
	    ob->slit.blen  = floatvalue (strtok (NULL, WHITESPACE));
	    ob->slit.angle = floatvalue (strtok (NULL, WHITESPACE)) * Degree;
	    if (special == 1) {  // special, test for different
vect2	d;
/* Find difference between position from ra/dec and from record */
		d = sub2vect (ob->smpos,
				maskvect (obs, ob->dat.ra, ob->dat.dec));
		if (vect2norm(d) > 0.003) {	// wrong position
// The .003 mm is about .01 arc seconds for Magellan
		    preq = pushtxt (tbuf, preq);
		    ob->flag |= OBJECT_SPECIAL;
		}
	    }
	    ob->dat.precom = preq;
	    preq = NULL;
	    ob->val = NULL;

/* Put the object into the queue */
	    ob->type = OBJ_REFERENCE;
	    ob->flag |= OBJECT_ACTIVE;
	    obs->ob = push_obj (ob, obs->ob);

/* End with a default case, here } else { ... with error message */
	} else {		/* Default unknown keyword */
	    printf (" ** Unknown keyword \"%s\" ignored.\n", key);
	}	/* End of keyword parsing */

/* Copy from examples in reading obs and obf files in the mgutils.c
  source code.  */



    }	/* End of reading loop */

//printf (" End of reading SMDF file.\n"); fflush (stdout); //debug

/* Fix up the ending comment queues, and any processing depending
  on the data from the whole file */
    if (ob == NULL) obs->comments = push_name (postq, obs->comments);
    else ob->dat.postcom = postq;
    obs->comments = push_name (comq, obs->comments);

/* Close the file when all read. */
    fclose (in);
//printf (" Closed input file.\n"); fflush (stdout); //debug

//printf ("End read_smdf, telescope %s.\n", obs->Tel_scope->name); //debug

/* We do not do normalize_obs here, but do catch a few things which
  would normally be done there... */

/* Try to fill in the detector if possible */
    kill_ptr ((void*)&(obs->Detect));
    obs->Detect = get_detector (obs->Instr_mt, gs, obs->dlevel);

/* Return the structure when file is done */
//printf (" Returning from read_smdf.\n"); fflush (stdout); //debug
    return  obs;
}


	/* - - - - - - - - - - *\
	|			|
	|   Wing Chip Support	|
	|			|
	\* - - - - - - - - - - */

/*  ====  Wing Chip Support programs  ====  */

objq*	wing_holes (double x, double y, double d, int f)  {
/* Create a linked list of new object holes with diameter d, and
  coordinates in x, y and +/- variations of x,y as provided by the
  flag f:  1 for x-reflected, and 2 for y-reflected, 3 for both. */
objq	*b;
objq	*h;
char	buf[8];
char	sc='x';
char	oc='l';
int	hr;
static	int	kh=0;

/* Set up some local variables */
    hr = (x != 0.0) && (f & 1);		// Horizontal reflection
    sc = hr ? 'r' : 'x';		// Literal for right side
    oc = 'l';				// Literal for other side

/* Get an object data structure and clean it */
    b = xalloc (objq);
    memset (b, 0, sizeof(objq));	/* CLEAR TO ZERO */
    b->next = b->last = b;

/* Set up the data for the base hole */
    b->smpos.x = x;
    b->smpos.y = y;
    b->type = OBJ_REFERENCE;
    b->flag |= (OBJECT_ACTIVE | OBJECT_SERVICE);
    b->slit.shape = CIRCLE;
    b->slit.width = d;
    b->slit.alen = b->slit.blen = b->slit.width / 2.0;
    sprintf (buf, "wc%c%02d", sc, kh);
    b->dat.name = dynamstr (buf);

/* Get the x-reflected version too, if flag 1 is on */
    if (hr) {
	h = xalloc (objq);
	memcpy (h, b, sizeof(objq));
	h->next = h->last = h;
	h->smpos.x = -x;
	sprintf (buf, "wc%c%02d", oc, kh);
	h->dat.name = dynamstr (buf);
	h = push_obj (b, h);
    }
    kh++;

/* If y is not zero, we need reflected up/down versions, too. */
/* Consider using (fabs(y) > 5.0) here instead... */
    if (y != 0.0 && (f & 2) ) {
	h = xalloc (objq);
	memcpy (h, b, sizeof(objq));
	h->next = h->last = h;
	h->smpos.y = -y;
// Note that hr and sc were set above
	sprintf (buf, "wc%c%02d", sc, kh);
	h->dat.name = dynamstr (buf);
	h = push_obj (h, b);

	if (hr) {
	    h = xalloc (objq);
	    memcpy (h, b, sizeof(objq));
	    h->next = h->last = h;
	    h->smpos.x = -x;
	    h->smpos.y = -y;
	    sprintf (buf, "wc%c%02d", oc, kh);
	    h->dat.name = dynamstr (buf);
	    h = push_obj (h, b);
	}

	kh++;
    }

    return  b;
}


objq*	wing_set (Obs* obs)  {
/* Derive a wing chip set appropiate to the observation */
objq	*q;
objq	*p;
//vect3	a, r, v;
//double	ra, dec;
double	wcd;	// Wing chip hole size
int	k=0;	// debug

    q = NULL;
// check for nod/shuffle, null if so (dca or dcb non-zero)
    if (obs->dca != 0.0) return NULL;
    if (obs->dcb != 0.0) return NULL;
// if ifu mode, null (we're not here then anyway...)
    if (obs->IFUmode) return NULL;
// check for inhibit, null if so
// Also return NULL if using MOE, the echelle...
//    if (obs->Flag & OFmoe) return NULL;
// Possibly for MOE, we cut one hole, or at most 2 holes
// left and right in special locations.

/* If not the IMACS instrument, we skip the wing chip holes */
    if (!strncmp(obs->Instr_mt->name, "LDSS", 4)) return NULL;

// check Det_mode; if Direct or LongCam to long; if ShortCam do short
// change to look at instrument name...
// Actually, direct and longcam are 0, ShortCam is 2, Grating is 1
// Decide which instrument, to select the proper set; either
// long or grating is one case, shortcam is the other...
    if (obs->Flag & OFmoe) {			// Holes for MOE
	wcd = 0.345;	// About 1.0 arc seconds
	q =           wing_holes ( 176.0, 8.0, wcd, 1);
    } else if (!strncasecmp(obs->Instr_mt->name, "IMACS_s", 7)) {  // short
// IMACS short camera wing holes
// create queue of matching holes left/right up/down
	wcd = 0.5175;	// About 1.5 arc seconds
	q =           wing_holes (310.0,  0.0, wcd, 1);
	q = push_obj (wing_holes (308.2, 15.0, wcd, 3), q);
	q = push_obj (wing_holes (306.4, 30.0, wcd, 3), q);
	q = push_obj (wing_holes (304.6, 45.0, wcd, 3), q);
	q = push_obj (wing_holes (302.8, 60.0, wcd, 3), q);
	q = push_obj (wing_holes (301.0, 75.0, wcd, 3), q);

    } else if (!strncasecmp(obs->Instr_mt->name, "IMACS_l", 7)) {  // long
// IMACS long camera wing holes
// create queue of matching holes left/right up/down
	wcd = 0.345;	// About 1.0 arc seconds
	q =           wing_holes (176.0,  0.0, wcd, 1);
	q = push_obj (wing_holes (175.0,  8.5, wcd, 3), q);
	q = push_obj (wing_holes (174.0, 17.0, wcd, 3), q);
	q = push_obj (wing_holes (173.0, 25.5, wcd, 3), q);
	q = push_obj (wing_holes (172.0, 34.0, wcd, 3), q);
	q = push_obj (wing_holes (171.0, 42.5, wcd, 3), q);

    } else {	// Generate the standard hole set for other cam.
// No hole set for this yet...
	;
    }

// Put in source positions -- a loop to create them...
    C_loop (p, q) {
// Find the ra/dec for the hole position...  (inv Op_transform)
	unmaskvect (obs, p->smpos, &(p->dat.ra), &(p->dat.dec) );
	k++;
    }
    if (obs->dlevel > 0) printf (" Added %d wing chip holes.\n", k);

    return  q;
}


void	add_wings (Obs* obs)  {
/* Add wing chip support into the obs structure */
objq	*wco;	// Wing chip objects

/* Call to get the wing chip hole set queue */
	wco = wing_set (obs);

/* Add the wing chip queue to end of object queue. */
	obs->ob = push_obj (wco, obs->ob);

}


	/* - - - - - - - - - - *\
	|			|
	|    Slit  Extension	|
	|			|
	\* - - - - - - - - - - */

/*  ====  Slit Extension subroutines  ====  */

/* node push and pop standards for the sleg type. */
/* Note that the sleg queue is NOT circular */

sleg	*head_sleg (sleg* a) {
    return  G_head (a);
}

sleg	*head_sleg2 (sleg* a) {
/* Return the head of a; NULL if a is NULL, the preceeding item pointing
  to NULL if any, otherwise it is a itself */
sleg	*p;
    if (a == NULL) return a;
    for (p=a; p->last != NULL; p = p->last) {
	if (p->last == a) return a;
    }
    return  p;
}

sleg	*tail_sleg (sleg* a) {
    return  G_tail (a);
}

sleg	*tail_sleg2 (sleg* a) {
/* Return the tail of a; NULL if a is NULL, the following item pointing
  to NULL if any, otherwise it is the item pointing to a */
sleg	*p;
    if (a == NULL) return a;
    for (p=a; p->next != NULL; p = p->next) {
	if (p->next == a) break;
    }
    return  p;
}

static	sleg	*push_sleg (sleg* a, sleg* b) {
    return  G_push_l ( a,  b );
}

static	sleg	*pop_sleg (sleg* a) {
    return G_pop_l (a);
}

/* The above push and pop entries are declared static, and used only
  within this module.  Restore to the header file and remove the
  static declarations if they are to be used outside this module.  */

sleg	*kill_sleg (sleg* a) {
/* Remove the object, and return the next one */
sleg	*r;
    if (a == NULL) return a;
    r = pop_sleg (a);	// The one to return
/* We have no dynamic allocated storage pointers in sleg yet */
    free (a);
    return  r;
}


// Add sort decision routine for free length here;
// copy from the decision routines for object sorting.

static	int	sleg_decide (sleg* a, sleg* b) {
/* Return -1, 0, or 1 for a>b, a=b, a<b */
double	ac, bc;
int	as, bs;

/* If either is NULL, return 0 */
    if (a == NULL) return 0;
    if (b == NULL) return 0;

/* Make first found fixed node come first */
/* Also, put free nodes ahead of provisional */
    as = a->objec->val->stat[a->side] & (Stat_Free | Stat_Prov);
    bs = b->objec->val->stat[b->side] & (Stat_Free | Stat_Prov);
    if (as < bs) return  1;
    if (as > bs) return -1;

/* Basically, we decide on basis of free length.  However,
  it is good for a fixed node to come first. */
    if (a->free < b->free) return  1;
    if (a->free > b->free) return -1;

/* Also choose by Bx position of edge, and then of center */
    ac = (a->side)?a->objec->val->sod->e2.p.x:a->objec->val->sod->e1.p.x;
    bc = (b->side)?b->objec->val->sod->e2.p.x:b->objec->val->sod->e1.p.x;
    if (ac < bc) return -1;
    if (ac > bc) return  1;
    ac = a->objec->val->sod->center.x;
    bc = b->objec->val->sod->center.x;
    if (ac < bc) return -1;
    if (ac > bc) return  1;

    return  0;
}


static	sleg*	sort_merge_eq (sleg *a, sleg *b, int(decision)(sleg*,sleg*))
/* Return pointer to head of merged queue with all nodes
in sort order. */
/* ==NOTE== the decision function is implicitly global here */
{
sleg	*p;
sleg	*q;

/* Binary merge operation:
  Start a new queue at a null head.
  Loop while both a, b are not null:
	Test element at head of a and b queues
	If a comes first, pop that element from a and
	push it onto the result queue.
	Otherwise, b comes first, pop that element and
	push it onto the result.
  The loop exited when a or b became null.  Push each of these
  onto the result queue, a null being pushed is a null operation.

.. */
    p = NULL;
    while ( (a != NULL) && (b != NULL) ) {
/* Binary merge, decide which one is pulled */
	if (decision (a, b) < 0)	b = pop_sleg (q=b);
	else				a = pop_sleg (q=a);
// Can also use the general routine directly:
//	else				a = G_pop_l (q=a);
	p = push_sleg (p, q);
    }
    p = push_sleg (p, a);
    p = push_sleg (p, b);
    return  p;
}

static	sleg*	sort_sleg (sleg *head, int(decision)(sleg*,sleg*))
/* Return pointer to head of sorted queue */
{
sortq	*sqh, *sqa, *sqb;
sleg	*cur;

/* Allocate meta-queue heads in a list, each pointing to a
sort-order queue.  Split the queues between these, probably
with NULL split ends. */

/* How to start:
  If the object queue is null, return null.

  Create outer loop - initialize by obtaining a current object by
  popping it from the object queue, and set sort queue head null.
  Outer loop (while current object not null):
    Obtain a new sort queue element
    Do the inner loop as follows--
    Inner loop (forever):
	Put the current object into the sort object queue.
	Obtain next object from object queue; if null, break.
	Test current object w.r.t. last item in sort object q;
	if correct order, push into queue, if not, break.
    Push the sort queue element into the sort queue.
  -- end of outer loop.
  Result is the sort queue.

.. */

    sqh = NULL;
    head = pop_sleg (cur=head);
    while (cur != NULL) {
/*	sqa = (sortq*) malloc (sizeof(sortq)); */
	sqa = xalloc (sortq);
	sqa->q = NULL;
	sqa->next = sqa->last = NULL;
//	while (1) {
	for (ever) {
	    sqa->q = push_sleg ( (sleg*)(sqa->q), cur);
// NOTE (sleg*)sqa->q also works, at least in solaris.
	    head = pop_sleg(cur=head);
	    if (cur == NULL) break;
	    if (decision (tail_sleg((sleg*)(sqa->q)), cur) < 0) break;
// PERFORMANCE ENHANCEMENT
// The tail_sleg slows us down; remember the tail (former cur)
// and use that; set it null at the right time, when new sqa is made.
	}
	sqh = push_sq (sqa, sqh);
    }

/* The meta-queue list is circular.  Traverse it, popping off
two items at a time, merging them, and pushing the result back onto
the meta-queue.  When second pop gives NULL, the first one
is the result, which is returned. */
    if (sqh == NULL) return NULL;
    while (sqh != NULL) {
	sqh = pop_sq (sqa=sqh);
	if (sqh == NULL) break;
	sqh = pop_sq (sqb=sqh);
	sqa->q = sort_merge_eq ((sleg*)(sqa->q), (sleg*)(sqb->q), decision);
	sqb = kill_sq (sqb);
	sqh = push_sq (sqa, sqh);
    }
    cur = (sleg*)(sqa->q);
// NOTE  (objq*)sqa->q; also works, at least in solaris
    kill_sq (sqa);

    return  cur;
//    return  sqa->q;
}


static	void	augment_side (Obs* obs, objq* ob, int e, double a) {
/* Augment side e of object ob by a, creating a new extended queue
  if necessary.  Change parameters in ob to reflect all the new
  dimensions.  Similar computations to getspect, but operating on
  an existing structure.  */
double	olap;
spect	*sod;
int	x;

/* Test for trivialities */
    if (ob == NULL) return;
    if (a == 0.0) return;	// Nothing to do

//	TX_beg (TS_autst);
    sod = ob->val->sod;		// To reduce expression size

/* Change the indicated side of the spectrum */
    olap = 0.5 * obs->Detect->psize.x * obs->minsep;
    if (e == 0) {		// left edge is 0, right is 1
	ob->slit.alen += a;
	sod->e1 = getedge (obs, ob, lslit(ob->slit),
		olap, sod->center.x, obs->order,
		obs->pb );
    } else {
	ob->slit.blen += a;
	sod->e2 = getedge (obs, ob, rslit(ob->slit),
		olap, sod->center.x, obs->order,
		obs->pb );
    }

//	TX_beg(TS_augsb);
/* Update the bounds and on detector value; and secondary bound */
//    sod->bb = boundor (bbedge(sod->e1), bbedge(sod->e2));
    sod->bb = bbspec (sod);
    sod->on_det = (bbover (sod->bb, obs->Detect->bb) & 7);
// And, the second bound box logic...
    if (obs->pb2act) {	// Secondary
sp_edge	e1, e2;
	e1 = getedge (obs, ob, lslit(ob->slit),
	    olap, sod->center.x, obs->order, obs->pb2);
	e2 = getedge (obs, ob, rslit(ob->slit),
	    olap, sod->center.x, obs->order, obs->pb2);
// Use similar code to find bounds...
	sod->bd = boundor (bbedge(e1), bbedge(e2));
    } else {
	sod->bd = sod->bb;
    }
//	TX_enm(TS_augsb, "augment_side bounds");

/* If extra orders are enabled, replace the image queue with a new one */
#ifdef  EXORDX
    if (ob->type == OBJ_OBJECT) x = 2;
    else if (ob->type == OBJ_REFERENCE) x = 1;
    else  x = 0;
#else
    x = 3;	// compatable with old methods
#endif
/* Above sets up the extra order selection bit. */
    if (obs->ex_order & x) {
/* Enhance this use just as the special ex_order flag is done in the
  regular implementation; decide reference vs. object and use the
  correct bit in the ex_order flag for compares. */
int	ord;
spect	*sp;
//	TX_beg (TS_augxo);
// Kill any existing queue
	while (ob->val->img != NULL) ob->val->img = kill_spect (ob->val->img);
// go up in order
	for (ord=obs->order + 1; ord < XordLim; ord++) {
	    sp = getspect (obs, ob, ord, obs->pb);
	    if (sp->on_det) {
		sp->next = ob->val->img;
		ob->val->img = sp;
	    } else {
		sp = kill_spect (sp);
		break;
	    }
	}
// and down
// NOTE - should start from an int baseorder, and test
// to baseorder-XordLim; then set the limit value at top to ~10 or ~15
	for (ord=obs->order - 1; ord > -XordLim; ord--) {
	    sp = getspect (obs, ob, ord, obs->pb);
	    if (sp->on_det) {
		sp->next = ob->val->img;
		ob->val->img = sp;
	    } else {
		sp = kill_spect (sp);
		break;
	    }
	}
// Extra orders vastly increase the run time here.
//	TX_enm (TS_augxo, "augment_side extra-ord");
    }
/* That's it.  Void routine needs no explicit return */
//	TX_enm (TS_autst, "augment_side total");
}


static	objq	*conflict_x (Obs* ob, objq* b, int e, double a, objq* f) {
/* Find and return conflict to edge e of object b, when augmented
  by length a.  (Mark if m nonzero.)  Use fast algorithm if f non-null,
  i.e. check for conflict with f first, and return it if found, if
  not, then search and return first conflict found.  When f is null,
  do a complete search.  Use self (b) to indicate conflict with the
  detector edge.  The Obs strcuture is also needed here.
  This is the heart of slit extension.    */
objq	*r=NULL;	// Return value
objq	*q;	// scratch pointer
objq	tob;	// Local test object
objval	xval;	// Local test values
spect	sod;	// Local spectrum
spect	*img=NULL;	// Local image queue
int	cx;	// conflict status
int	cs = 0;	// complete search flag

/* Take care of trivial cases */
    if (b == NULL) return b;
    if (a == 0.0) return NULL;	// May need to really test this case.
    TX_beg (TS_conxs);
    Tobj = &tob;
/* The above suppresses the creating of conflict queue entries when
  find_conflict is called using &tob as the primary node. */

/* Construct a local copy of the node b in tob, xval, sod */
// copy b to tob
//	TX_beg(TS_initc);
//	TX_beg(TS_copsd);
    tob = *b;
// change all dynamic things in tob (links, objval)
    tob.last = tob.next = NULL;
    xval = *(b->val);
    tob.val = &xval;
// change all dynamic things in objval (cfq, sod, img)
// use a static sod to save time
    xval.cfq = NULL;	// maybe unnecessary
    sod = *(xval.sod);
    xval.sod = &sod;
    sod.next = NULL;
    xval.img = NULL;

/* Augment the local copy by the requested amount a, making a test object */
//	TX_beg(TS_augms);
    augment_side (ob, &tob, e, a);
//	TX_enm(TS_augms, "augment side in conf-x");
    img = xval.img;		// Save it for later

/* The local copy in tob is now created, and resembles the actual
  object for which conflict status is desired.  */
//	TX_enm(TS_copsd, "copy sod");
/*  Set the scan delta for find_conflicts;
  This is needed whether or not we do a full scan here */
    DeltaX = 1.1 * (bigest_bounding_box (ob->ob) + a);
// PERFORMANCE ENHANCEMENT
// NOTE that this computes a full scan of all objects
// for every conflict test; we should globally store that bbb comp!
/* Quick conflict finding feature -- if f present, test first for
  a conflict there.  If no f, or no conflict, a full search is needed. */
    cs = 0;		// No search unless needed
    if (f != NULL) {

//		TX_beg(TS_inifc);
// Select the secondary object; either ourselves or f...
	q = (f == b) ? &tob : f;
	cx = ( find_conflict (&tob, q, ob) > 0 );

	if (cx) r = f;	// found conflict and return it
	else    cs = 1;	// no conflict, do search
//		TX_enm(TS_inifc, "Initial find conf.");
    } else  cs = 1;	// The null f case comes here...

/* The find_conflict on self now checks for slit within the IMACS radius */
//	TX_enm (TS_initc, "Initial conf-x (2-4)");
    if (cs) {	// Need a complete search
objq	*p;
int	i;

	TX_beg (TS_confx);

/* Test for self conflict first... */
	if (find_conflict (&tob, &tob, ob)) r = b;

/* Determine the scan direction for this edge.  If edge center less
  object (in X) is positive, we scan to next links, otherwise last */
	if ( (e?sod.e2.p.x:sod.e1.p.x) > sod.center.x) {	// next
	    for (p=b->next; p != NULL && r == NULL;
		p = (p->next == b) ? NULL : p->next) {
		if ((i = find_conflict(&tob, p, ob)) < 0) break;
		if (i > 0) r = p;	// Found the conflict
//		if (p->next == b) break;
// BAD CONSTRUCT BREAK
	    }	// The search loop.
	} else {						// last
	    for (p=b->last; p != NULL && r == NULL;
		p = (p->last == b) ? NULL : p->last) {
		if ((i = find_conflict(&tob, p, ob)) < 0) break;
		if (i > 0) r = p;	// Found the conflict
//		if (p->last == b) break;
// BAD CONSTRUCT BREAK
	    }	// The search loop.
	}
	TX_enm (TS_confx, "Conf-x Fullsearch");
    }	// End of complete search
/* When done, return r as our result */

/* Mark the conflicted node... */
//     if ( (r != NULL) && (r != b) && m) {
// // Mark the oposite side of the conflicted node
// 	r->val->stat[1-e] |= Stat_Mark;
//     }

/* Kill any local image queue which happens to be left over */
// This can cause elevated run times...
    while (img != NULL) img = kill_spect (img);
    TX_enm (TS_conxs, "Conf-x Search");
    return  r;
}


static	double	free_space (Obs* obs, sleg* eg, int cc) {
/* Return the free space available to slit edge eg */
/* cc is convergence criteria; 0 = low, 2 = high accuracy
 -- this is a debug try, to minimize time expended -- */
double	s = 0.0;
double	lv, hv;
double	df, tv;
objq	*ob;
objq	*cf;
objq	*pc;
int	e;
double	lvm[5] = {0.3, 0.4, 0.5, 0.5, 0.5};
double	dfm[5] = {3.0, 2.5, 2.0, 1.8, 1.8};
double	cvg[5] = {1.5e-2, 6.0e-3, 4.0e-3, 2.0e-3, 1.0e-3};

/* Usual trivial case tests... */
    if (obs == NULL) return s;
    if (eg == NULL) return s;
    ob = eg->objec;
    if (ob == NULL) return s;

/* Keep cc within bounds here (by calling conventions, cc is
  here a local variable) by limiting it to array dimensions */
    if (cc < 0) cc = 0;
    if (cc > 4) cc = 4;
// This should reduce some panics, perhaps.

/* If the edge is not free or provisional, it has defined no
  free space, so return zero... */
    e = eg->side;
    if ( !(ob->val->stat[e] & (Stat_Free | Stat_Prov)) ) return s;

    TX_beg (TS_freet);

/* Preliminary -- Find the interval; a low and high value for
  the free length, with the low value resulting in no conflict,
  and the high value resulting in a known conflict. */

/* First part, find an acceptable low value.  Try the current
  value, it should work.  If not, try half that.  If the current
  value is zero, and it conflicts, we have a serious error.  Any
  conflicting value found is saved as a high value.  */
    pc = NULL;

    lv = eg->free - 0.001;	// Reduce it a little to avoid panic?
// The reduction is needed; without it, execution time increases a lot.
// IDEA -- reduce it if it is greater than zero only??
//    lv = eg->free;
//    if (lv > 0.0) lv -= 0.001;
// The above method is observed to result in some panics at lv = 0.0

    hv = -1.0;		// An impossible value as a flag
    cf = conflict_x (obs, ob, e, lv, pc);

    while (cf != NULL) {
	if (lv < 1.e-4) {
/* PANIC - conflict exists at zero  extension! */
// Sometimes can find this if we jam two slits together; to suppress
// extra messages, we will allow a low of -1 microns for this test...
// Funny thing - the report does not show the -0.001 of lv above!
/* Message only if lv is not 0 and not the -0.001 fudge value.
  We kill all low lv values anyway; warn only if we can't find lv
  from non-zero free value.   */
/* Report the node name, and return zero */
//	    if (lv < -0.001) return s;	// skip the panic messages
	    if (lv > -0.001 && lv != 0.0) {  // Skip if lv is 0 also...
// Otherwise we get messages for valid non-extend condition...
//	    if (lv <= 0.0) {
		printf (" **Panic: Conf-zext: %s (%d) Free %.3f stat %d",
		    ob->dat.name, e, lv, ob->val->stat[eg->side]);
		printf (" with %s st %d.\n",
		    cf->dat.name, cf->val->stat[1-eg->side]);
	    }
	    TX_end (TS_freet);
// may want to set confl to cf here?
	    if (eg->confl == NULL) eg->confl = cf;	// to not return null
	    return  s;
	}
// Remember the current conflict in the structure, and set high
	eg->confl = cf;		// Put above before the conditional?
	hv = lv;
// Try half the low value
//	lv *= 0.5;
// Be more aggressive at low cc
	lv *= lvm [cc];
// double	lvm[5] = {0.3, 0.4, 0.5, 0.5, 0.5};
	pc = cf;
	cf = conflict_x (obs, ob, e, lv, pc);
    }
// Upon exit from loop, lv is an acceptable non-conflict value.

/* Second part, find an acceptable high value.  If a conflict was
  found in the first part, it is used.  Otherwise, we need to try
  a value slightly higher than low.  If low is zero, use 1, and if
  not zero, use low/4.  Then double the addition until a conflict
  is found.  If any non-conflict is found here, keep as a new low. */
    df = (lv > 0.0) ? lv / 4.0 : 1.0;
    while (hv < 0.0) {
	tv = lv + df;
	cf = conflict_x (obs, ob, e, tv, pc);
	if (cf == NULL) {	// Try again, and update low
	    if (fabs(df) > 2.0*D_maskrad) {	// whatever; a large value
/* PANIC -- too great a test extension */
/* Report node name, and return zero */
		printf (" **Panic: No conflict extending %s (%d)",
			ob->dat.name, e);
		printf (" by %11.3e.\n", tv);
		TX_end (TS_freet);
		return  s;
	    }
	    lv = tv;
//	    df *= 2.0;	// Maybe not needed, but faster
// Change increment with cc
//double	dfm[5] = {3.5, 3.0, 2.0, 1.8, 1.8};
	    df *= dfm[cc];
	} else {
	    hv = tv;	// And, we terminate loop
	    eg->confl = pc = cf;
	}
    }
// Upon exit from loop, hv is a known conflict value

/* Refinement of interval; a binary search between the low and
  high lengths, until the distance between them is reduced to
  a determined tolerance.  The low side is the free length found. */
//    while ( fabs (df = hv - lv) > 0.75e-3 ) {
//    while ( fabs (df = hv - lv) > 4.0e-3 ) {
// Change convergence with cc
//double	cvg[5] = {1.5e-2, 6.0e-3, 4.0e-3, 2.0e-3, 1.0e-3};
    TX_beg (TS_freeb);
    while ( fabs (df = hv - lv) > cvg[cc] ) {
	tv = 0.5 * (lv + hv);	// halfway between
//	cf = conflict_x (obs, ob, e, tv, 1, pc);
TX_beg (19);
	cf = conflict_x (obs, ob, e, tv, pc);
TX_enm (19, "bsearch confl-x");
	if (cf == NULL) {	// a new low
	    lv = tv;
	} else {		// a new high
	    hv = tv;
	    eg->confl = pc = cf;
	}
    }
    TX_enm (TS_freeb, "Free space bsearch");
// At end of loop, lv is a safe value
    TX_end (TS_freet);
    s = lv;
/* if s is not greater than 0, we should fix the edge here? */
    return  s;
}


sleg	*scan_free (Obs* obs) {
/* Scan object queue, return slit edge queue */
sleg	*r=NULL;	// Returned
sleg	*s;
objq	*ob;
objq	*p;	// had.. , *q;
objval	*val;
int	good;	// test flag
int	e;	// edge

/* Check for real data */
    if (obs == NULL) return NULL;
    ob = obs->ob;
    if (ob == NULL) return NULL;

/* Main loop over the object queue... */
    C_loop (p, ob) {

/* The candidate object must meet several criteria; must be
  active, an object rather than alignment, have a slit shape, and
  perhaps others.  Collect these here into a single decision variable */
	good = p->flag & OBJECT_ACTIVE;
	good = good && (p->type & OBJ_OBJECT);	// No gaps or reference
	good = good && (p->slit.shape == RECTANGLE);
	val = p->val;
	good = good && (val != NULL);

/* The slit extension will apply to both left and right sides
  of all good slits... */
	for (e=0; good && (e<2); e++) {
	    s = xalloc (sleg);
	    s->next = s->last = NULL;
	    s->objec = p;
	    s->confl = NULL;
	    s->side = e;
//	    s->xtnd = 0.0;	// may not be needed
	    s->free = 0.0;	// Eliminate that nan junque
	    val->stat[e] |= Stat_Free;		// Set edge free
//TX_beg (15);
	    s->free = free_space (obs, s, 0);
//TX_enm (15, "Free space scan");
	    r = push_sleg (r, s);
	}	// Node tested good; do both sides
    }	// End loop on p over object queue
    return  r;
}

sleg	*scan_sle (sleg* a, Obs* obs) {
/* Scan the list, update any marked nodes for extension length,
  and update any status needed.  */
//sleg	*r=NULL;	// Returned
sleg	*p;
objq	*ob;
double	fl;
int	e;
int	sx;

    if (a == NULL) return a;

    for (p=head_sleg(a); p != NULL; p = p->next) {
	e = p->side;
	ob = p->objec;
	sx = ob->val->stat[e];
	if (sx & Stat_Mark) {
//TX_beg (16);
	    fl = free_space (obs, p, 1);
//TX_enm (16, "Free scan marked");
	    p->free = fl;
//	    ob->val->stat[e] &= (!Stat_Mark);
// The above seems to set status wrongly!
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  The problem is that the ! operator is a logical negation,	 */
/*  while what is desired is a bitwise inversion.  That is	 */
/*  done by the ~ operator instead!!!  Fix all instances.	 */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
// if fl is zero, we might set prov status??
//	    ob->val->stat[e] = sx ^ Stat_Mark;
	    ob->val->stat[e] &= ~Stat_Mark;
	}
    }

    return  head_sleg(a);
}

sleg	*sort_sle (sleg* a) {
/* Sort the list in extension length order, all fixed nodes
  come first.  If a tie, choose by object Bx position, and if
  still tied, choose edge center Bx position */
sleg	*r;	// Returned

    r = sort_sleg (a, sleg_decide);
    return  head_sleg (r);
//    return  a;
}

sleg	*fix_sle (sleg* a, Obs* obs) {
/* Scan the list, removing all fixed nodes, setting any nodes
  to fixed which cannot move.  First node found with free status
  will be set to provisionally fixed.  Any further fixed nodes
  are also removed.  Return the resulting list head. */
sleg	*r=NULL;	// Returned
sleg	*p;
objq	*ob;
objq	*obc;	// conflicting object
int	e;
int	stx;
int	stc;
double	fl;
#ifdef  debug_slex
static int dbk=0;	// debug counter
#endif

/* Trivial cases */
    if (a == NULL) return a;

/* Scan in a simple loop.  List is not circular */
    for (p=head_sleg(a); p != NULL; p = p->next) {

/* Remove all fixed nodes we find... */
	while (p != NULL) {	// note the break condidions
	    if (p->objec == NULL) { p = kill_sleg(p); continue; }
	    e = p->side;
	    stx = p->objec->val->stat[e];
	    if ( stx & (Stat_Free | Stat_Prov) ) break;
	    p = kill_sleg (p);	// Advances to next node here
	}
	if (p == NULL) break;	// Could end early here...

/* Find current free length; then if the conflicting node is fixed,
  we fix this node and go on. */
//TX_beg (17);
	fl = free_space (obs, p, 2);
//TX_enm (17, "Free fix");
	e = p->side;
	ob = p->objec;
	stx = ob->val->stat[e];
	obc = p->confl;
// TRY THIS -- if obs is null, it has no conflict;
// then stc is undefined;  If so, we set stc to 0
	if (obc != NULL) {
	    if (obc->val != NULL) {
		stc = obc->val->stat[1-e];
	    } else {
		stc = 0;	// Conflict has no val field.
	    }
	} else {
	    stc = 0;	// No conflict known
	}

/* NOTE -- the actual conflict is NOT ALWAYS with the other edge... */
	if (obc == ob) stc = 0;		// self conflict is fixed

	if ( !(stc & (Stat_Free | Stat_Prov)) ) {
#ifdef  debug_slex
if (dbk++ < 30) {
	printf (" 1 Extend %s (%d) by %.3f; free = %.3f, stat=%d",
		ob->dat.name, e, fl, p->free, stx);
	printf (" Conf: %s stat=%d.\n",
		obc->dat.name, stc);
}
#endif
	    augment_side (obs, ob, e, fl);
//	    ob->val->stat[e] &= ~(Stat_Free | Stat_Prov);
//	    ob->val->stat[e] |= Stat_Mark;
	    ob->val->stat[e] = 0;
	    p->free = 0.0;
//	    continue;		// try break?
	    break;		// no, that's not it...
	}

// other free, this free -> move 1/2 way and set prov; break
// other free, this prov -> stay put.
// other prov, this free -> move all the way and set prov
// other prov, this prov -> move all the way and fix; break


/* Node is not fixed.  Treatment depends on status of conflict and p */
	if (stc & Stat_Prov) {
/* Against a provisional edge, move our total free length. */
#ifdef  debug_slex
if (dbk++ < 30) {
	printf (" 2 Extend %s (%d) by %.3f; free = %.3f, stat=%d",
		ob->dat.name, e, fl, p->free, stx);
	printf (" Conf: %s stat=%d.\n",
		p->confl->dat.name, stc);
}
#endif
	    augment_side (obs, ob, e, fl);
	    p->free = 0.0;
	    if (stx & Stat_Free) {
		ob->val->stat[e] &= ~Stat_Free;
//		ob->val->stat[e] |= (Stat_Prov | Stat_Mark);
		ob->val->stat[e] |= Stat_Prov;
//		break;
	    }
	    if (stx & Stat_Prov) {
		ob->val->stat[e] &= ~(Stat_Free | Stat_Prov);
//		ob->val->stat[e] |=  Stat_Mark;
//		break;
	    }
	    break;		// See if this helps (no, but doesn't hurt)
	}

	if (stc & Stat_Free) {
/* Against a free edge.  Move if we are free... */
	    if (stx & Stat_Free) {
		fl *= 0.5;
#ifdef  debug_slex
if (dbk++ < 30) {
	printf (" 3 Extend %s (%d) by %.3f; free = %.3f, stat=%d",
		ob->dat.name, e, fl, p->free, stx);
	printf (" Conf: %s stat=%d.\n",
		p->confl->dat.name, stc);
}
#endif
		augment_side (obs, ob, e, fl);
		p->free -= fl;
		ob->val->stat[e] &= ~Stat_Free;
//		ob->val->stat[e] |= (Stat_Prov | Stat_Mark);
		ob->val->stat[e] |= Stat_Prov;
		p->confl->val->stat[1-e] |= Stat_Mark;
		break;
	    }
//	    break;		// See if this helps (no, but may speed it)
	}

    }	// End of scan of list
    r = head_sleg (p);
    return  r;
}




void	extend_slits (Obs* ob) {
/* Do slit extension for active objects */
sleg	*r;	// Slit edge queue
clock_t	bt, et;
double	tu;
#ifdef  debug_slex
int	dbk=0;	// debug counter
#endif

// DEBUG output of MOE horizontal bounding box values
// #undef DMOE
#if (DEBUG_GEN > 0)
#ifdef DMOE
    if ( ob->Flag & OFmoe ) {
objq	*p;
objval	*v;
spect	*s;
//	for (p=ob->ob; p != NULL; p = (p->next == ob->ob) ? NULL : p->next) {
	C_loop (p, ob->ob) {
	    if (!(p->flag & OBJECT_ACTIVE)) continue;
	    v = p->val;
	    if (v == NULL) continue;
	    s = v->sod;
	    if (s == NULL) continue;
	    printf (" Obj: %04s, Mask: %8.3f %8.3f; Det: %8.3f %8.3f\n",
		p->dat.name, p->smpos.x, p->smpos.y,
		s->bb.x.lo, s->bb.x.hi);
	}
    }
#endif  // DMOE
#endif  // DEBUG_GEN

/* --ENHANCEMENT-- Should be able to return for report either or both
  the number of slit edges extended and the total length added ?? */

/* See if enabled */
    if (! ob->slext) return;

/* Cancel and give message if incompatable mode detected... */

/* Check for N&S mode, no extensions given if enabled */
    if ( (ob->dca != 0.0) || (ob->dcb != 0.0) ) {
	printf (" -- Slit extension canceled due to Nod & Shuffle mode.\n");
	return;
    }

/* Will also return if we have echellette (MOE) */
    if ( ob->Flag & OFmoe ) {
	printf (" -- Slit extension canceled for MOE disperser.\n");
	return;
    }

/* And, IFU mode is a show-stopper too */
    if (ob->IFUmode) {	// Change to ob->Flag & OFifu
	printf (" -- Slit extension canceled due to IFU mode.\n");
	return;
    }

/*  If we get this far, we must be able to actually extend
  some slits.  So we start the process */

    if (ob->dlevel > 0) printf ("\n  -- (Extending slits... ) --\n");
    bt = clock();
    fflush (stdout);
    ob->slext |= 2;	// Mark as being done

#ifdef  TS_dotime
    TX_rep (-1);	// Any which had been started...
    TX_con (TX_size);
#else
    TX_con (0);
#endif

    TX_beg (TS_whole);
    TX_nam (TS_freet, "Free Space - Total.");

/*  (1)  Sort the object list by x position (should be already so) */
    TX_beg (TS_sorts);
    ob->ob = sort_obj (ob->ob, sort_decide);
    TX_enm (TS_sorts, "Slit Ext. Sorts.");

/*  (2)  Scan the list and form a queue of free edges */
    r = scan_free (ob);


/*  (3)  Loop while the free edge list is non-null */
    while ( r != NULL ) {

/*  (4)  Call to scan the list to update free lengths */
	r = scan_sle (r, ob);
// The scan is not really needed, and few nodes are marked anyway

/*  (5)  Call to sort the list by free length */
	TX_beg (TS_sorts);
	r = sort_sle (r);
	TX_end (TS_sorts);

/* Debug section - print the queue  */
#ifdef  debug_slex
	if (dbk++ < 5) {
int	k=0;
int	n=0;
sleg	*p;
int	e;
int	stx;
objq	*ob;
// #define	ever  (;;)	// already done in this module
	    p = head_sleg(r);
	    for (ever) {
		if (p == NULL) break;
		n++;
		e = p->side;
		ob = p->objec;
		stx = ob->val->stat[e];
		if (
//			(stx & (Stat_Free | Stat_Prov)) &&
//			(p->free != 0.0) &&
			(k < 4)
		) {
		    printf (" %3d Slit %s (%d) Stat %d free = %.3f\n",
			n, p->objec->dat.name, p->side, stx, p->free);
		    k++;
		}
		p = p->next;
	    }	// loop forever over queue r
	    printf ("  (%d total nodes.)\n", n);
	}
/* End debug section */
#endif

/*  (6)  Call to drop node(s) from top of list */
	r = fix_sle (r, ob);

    }	// End of algorithm on free edge list
/*  (7)  End of loop, list complete */
/* (No list remaining, so our work here is done... */

    TX_enm (TS_whole, "Whole Slit Extension");

#ifdef  TS_dotime
    TX_rep (-1);	// Report and clear all
#else
    TX_con (0);
#endif

    et = clock();
    tu = (double)(et - bt) / (double)CLOCKS_PER_SEC;
    if (ob->dlevel > 0) printf ("  -- Took %.2f Seconds --\n", tu);

}


/*  ===  End Subroutines  ===  */

/*  ===  Revsion History  ===

2006/02-22 -- Version 1.46
	-- Include code for OBJCHECK feature.

2006/02-13 -- Version 1.45
	-- Drop a debug exit to reorder which was inhibiting the
	ordering of a very large set.
	-- Trying to speed up convergence in reorder.

2006/01-19 -- Version 1.44
	-- Drop read of integer from guide star data in read obs file.
	-- Also supports guide star data.

2005/12-23 -- Version 1.43
	-- Add check_temp to set ldss temperature default

2005/12-20 -- Version 1.42
	-- Add support for "date" field in .obs file to reading.

2005/12-08 -- Version 1.41
	-- Add a type name (object/alignment/unknown) on several messages
	about objects being removed.
	-- Got slit extend errors reduced by defining all elements of
	the newly created sleg structure.

2005/12-06 -- Version 1.40
	-- Change dodebug to obs->dlevel and compare to a value rather
	than use as a boolean.
	-- Add debug level to get_optic_data and get_cut_data calls.
	-- get_detector is now a static (local) program, with debug level

2005/11-28 -- Version 1.38
	-- Change reading of maskcut data to use get_cut_data, and that
	and get_optic_data to use common file finding routine.

2005/11-21 -- Version 1.37
	-- Refine storage of maskcut data

2005/11-18 -- Version 1.36
	-- Add new storage method for the maskcut.cfg data.
	-- Add comment giving observing catalog data
	-- Add support of reflimit parameter in obs (read/write/init)

2005/11-09 -- Version 1.35
	-- Adding better computation of object being in or out of a
	circle.  To use for on mask, and out of avoids.
	-- Move reading of maskcut.cfg data here from various.

2005/10-20 -- Version 1.34
	-- Add computation of direct image of object and storing it
	into objq.val->dimg, a bounding box on detector plane.
	This adds the "direct_loc" subroutine, computing a bbox.
	-- Add warning message for reference object found in conflict
	scan to have a direct image off the detector.
	-- Add code to remove (optionally) conflicting alignment objs.
	-- Add report of remaining conflicts (repconf program)

2005/10-07 -- Version 1.33
	-- The DGrism data type is no longer supported, the EGrism
	data type has replaced it.

2005/08-25 -- Version 1.32
	-- The optutils.c is supporting extended grism data type.
	-- MOE local routines now declared static.
	-- Some dodebug uses have been assigned to DEBUG_UTL
	-- Support for echelle center wavelength computation being
	dynamic is included now.

2005/08-24 -- Version 1.31
	-- Write DLIMIT if passband data differ as well as when
	a DLIMIT input has been written.  This fixes intgui error
	of not recognizing the dlimit input.
	-- Wrong - that was for write smdf, not write obs...
	-- We DO make default telescope depend on instrument in normalize.

2005/06-30 -- Version 1.30
	-- Includes support for virtual gap feature
	-- Removes some debug and unused code
	-- More subroutines are made static to this module

2005/06-16 -- Version 1.29
	-- Support for detector information in opticdef.dat

2005/06-14 -- Version 1.28
	-- Removal of grating structure from Obs structure
	-- Removal of instrument pointer from Obs structure
	-- Removal of telescope pointer from Obs structure
	-- Putting telescope name and instrument name into an
	otherwise null element structure if no optical queue during
	read of SMDF file input (for maskcut).
	-- Removal of some initialization in normalize_obs
	-- Removal of much commented code
	-- Should be functionally equivalent to previous versions

2005/06-10 -- Version 1.27
	-- Removal of Optic_element structures from Obs structure

2005/05-26 -- Version 1.26
	-- Implements bit selected extra order feature.
	-- Code for full reordering using proper endpoints included

2005/05-24 -- Version 1.25
	-- Add "genlist.c" as a literal inclusion section
	-- Removing various old, unused code sections (also from mgutils.h)

2005/05-20 -- Version 1.24
	-- Put reorder things on a debug level

2005/05-11 -- Version 1.23
	-- Obsolete Gwavl program deleted.
	-- Addition of pointers in grism types for glass element improves
	speed of slit extension by factor of 2 to 3.

2005/05-10 -- Version 1.22
	-- Support for 2-grism type
	-- Includes the count conflicts debug routine

2005/04-20 -- Version 1.21
	-- Introduce pbmean function, local.
	-- Use mean wavelength for spectrum center wavelength computation.
	-- Use obs.sw as either cw or mean of passband for maskvect.
	-- Drop the wlimit member of object data structure, use obs->pb
	-- Write, read WLIMIT and DLIMIT in .SMF read/write

2005/04-12 -- Version 1.20
	-- New smdf read/write routines including global optical
	structure pointer, special activity flags.
	-- Add pre/post comments to smdf i/o
	-- This is support of smfix

2005/03-23 -- Version 1.15
	-- Move features to mgfeats.h to allow globalization.

2005/03-22 -- Version 1.14
	-- Fixing a bug in slit extension which can extend a slit past
	the maximum cut radius of the mask.

2005/03-15 -- Version 1.13
	-- Change many instances of " " as search string in parse and strtok
	to WHITESPACE so as to recognize tabs too.  This should bring the
	behavior to agree with advertised delimiters.

2005/03-02 -- Version 1.12
	-- Support differential refraction

2005/02-01 -- Version 1.11
	-- More support for new LDSS

2005/01-13 -- Version 1.10
	-- Add LDSS_NEW symbol and support

2004/11-19 -- Version 1.09
	-- Some debug print output for MOE support.

2004/11-12 -- Version 1.08
	-- New subroutine for adding wing chip holes, and add wing
	chip holes to echellette support.

2004/11-09 -- Version 1.07
	-- Adding support programs for MOE - echellette

2004/11-02 -- Version 1.06
	-- Changing some loop end logic to new method.
	-- Read both #! and !# as pass-through comments in .obs file
	-- Add comment in SMDF that slits have been extended.

2004/11-01 -- Version 1.05
	-- Change order of node counting routines.
	-- Drop unused "nxbyn" subroutine - debug only.

2004/10-29 -- Version 1.04
	-- Add symbol PRIORDER to set order priority; initially set
	on, to give earlier objects the priority over later ones.
	-- Using fill_objects in maskgen instead of inline code.
	-- Change fillobject to return activity status of the object,
	and object_setup to return count of active objects found rather
	than total objects.
	-- Fix use of !OBJECT_ACTIVE to ~OBJECT_ACTIVE and tests
	of equality to activity flag to simply & operation.  This will
	correctly implement the flag bit concept.
	-- Global search for ! and & shows proper use in all cases.
	-- Add global MESCOUNT and message limits

2004/10-28 -- Version 1.03
	-- Improve again the duplicate removal in single_confs;
	remove earlier order if equal in priority
	-- Add removal of conflicts with musthave before any other
	conflict resolution...
	-- In priority sort, when equal sort by (inverse) order
	-- Remove conflictedest as being obsolete.
	-- In conflict count, for otherwise equal nodes, break a
	tie with order, earlier node is dropped.

2004/10-27 -- Version 1.02
	-- Removing some redundant code regarding conflict removal.
	-- Change to consider priority when removing duplicate objects.
	-- Improve and extend the single conflict finding program.
	-- Remove predeconflict as single conflict does it all now.
	-- Put single conflict tests of musthave conflicts under a
	symbol and let them be resolved in full resolution case.

2004/10-26 -- Version 1.01
	-- Add comment to .SMF file about grating angle when we have
	a grating (long camera).  (GratingAngle)
	-- Requires pcom field in obs data, needs to be initialized
	and killed properly.  Added to read_obs, clean_obs, and others,
	similarly to the "comments" field of obs structure.

2004/10-14 -- Version 1.00
	-- Introduced the MGUVERS string, basically global to be used by
	routines here, and by  including programs.  Also MGUNAME is defined
	just in case we want the include module name to be reported.
	-- Can track the utility changes this way, too.  The MGUVERS will
	be appended to a new form version in main maskgen, and others;
	eg ( ... %s%s ... VERSION, MGUVERS )
	-- For version numbering of just the utilities, a high order
	digit of utility version is in MGHIVD.  (This is NOT used yet)
	-- Symbol TS_dotime is used to control reporting of detail timing
	for the slit extension.  Turned on during debug work, and usually
	off for distributed versions.
	-- Some slit extension subroutines declared static, private to
	this module. (augment_side, conflict_x, free_space) -- We recognize
	that other subroutines should be static also...


... */

/* ----- Things "removed" are here... ----- */

/*  ===  Object queue manipulation  ===  */

/* How about a general push, pop for any queue whose starting
  pointers are last, next ?? */
