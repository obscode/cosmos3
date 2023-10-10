/*  optutils.c  ==  Generalized Optical Transform utilities  */

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

/* Define for proper external use */
// #define  IS_OPT_C

/* Standard includes */
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

//#include	"mgfeats.h"
// features not included by maskdat.h

/* Required Includes */
// #include	"optutils.h"
// #include	"../../utils/ioutils.h"
// LATER, MAKE INTO A SYMBOLIC INCLUDEs LIKE OTHER LIBRARIES
#include	"mgutils.h"
// mgutils.h is needed for queue utilities (listcount, etc.)
//#include	OPTUTILH
// optutils.h is included by mgutils.h
//#include	IOUTILH
// ioutils.h included by maskdat.h, included by optutils.h


/* Include the general reading code */
// #include	"../../gmap/datafile.c"
// the above will later be placed in ioutils.c and ioutils.h

/*  -- Globals --  */

int	mdebug = 1;	// global debug stuff for optutils.c

/* The following are present here, and are defined as extern
  in optutils.h which is then included in many places. */

int	dodebug = 0;	// global for compatability.
double*	GRangle = NULL;		// Grating angle[s]
int*	GRorder = NULL;		// Grating or Grism order[s]

element	*Cur_grating = NULL;
element *Cur_grism = NULL;

#ifdef  OTDEBUG
double	Cam_angle = M_PI / 4.0;	// special debug value...
int	llong = 1;		// Another debug thing
int	direct = 0;		// Another debug thing
double	cur_den = 300.0;	// Also debug for ntest
double	grindex = 1.0;		// Also debug for ntest
#endif

/* INITIALIZED GLOBALS MUST BE PUT HERE RATHER THAN IN THE HEADER FILE */

/*  ===  Keyword decoding  ===  */

Keyword Decode_Keyword (char* key)
{
Keyword kt;
 
    kt = VOID;
    if      (!strcasecmp (key, "END"))		kt = END;
    else if (!strcasecmp (key, "telescope"))	kt = Telescope;
    else if (!strcasecmp (key, "focuser"))	kt = FOCuser;
    else if (!strcasecmp (key, "camera"))	kt = FOCuser;
    else if (!strcasecmp (key, "colimator"))	kt = Colimator;
    else if (!strcasecmp (key, "collimator"))	kt = Colimator;
    else if (!strcasecmp (key, "reimager"))	kt = Reimager;
    else if (!strcasecmp (key, "grating"))	kt = GRAting;
    else if (!strcasecmp (key, "grism"))	kt = GRIsm;
//    else if (!strcasecmp (key, "dgrism"))	kt = DGRIsm;
    else if (!strcasecmp (key, "egrism"))	kt = EGRIsm;
    else if (!strcasecmp (key, "echelle"))	kt = ECHelle;
    else if (!strcasecmp (key, "reflection"))	kt = Reflection;
    else if (!strcasecmp (key, "angle"))	kt = ANgle;
    else if (!strcasecmp (key, "filter"))	kt = Filter;
    else if (!strcasecmp (key, "detector"))	kt = DETector;
    else if (!strcasecmp (key, "instrument"))	kt = INSTrument;
    else if (!strcasecmp (key, "foclen"))	kt = FOCLEN;
    else if (!strcasecmp (key, "scale"))	kt = SCale;
    else if (!strcasecmp (key, "foccurv"))	kt = FOCCURV;
    else if (!strcasecmp (key, "field"))	kt = Field;
    else if (!strcasecmp (key, "pixels"))	kt = Pixels;
    else if (!strcasecmp (key, "axis"))		kt = AXIS;
//    else if (!strcasecmp (key, "rotation"))	kt = ROTATION;
    else if (!strcasecmp (key, "lines"))	kt = Lines;
    else if (!strcasecmp (key, "glass"))	kt = GLASS;
    else if (!strcasecmp (key, "RMOD"))		kt = RMOD;
    else if (!strcasecmp (key, "TMOD"))		kt = TMOD;
    else if (!strcasecmp (key, "WMOD"))		kt = WMOD;
    else if (!strcasecmp (key, "ALIGNROT"))	kt = ALIGNROT;
    else if (!strcasecmp (key, "Prism1"))	kt = PRIS1;
    else if (!strcasecmp (key, "Prism2"))	kt = PRIS2;
    else if (!strcasecmp (key, "Gnormal"))	kt = GNORM;
    else if (!strcasecmp (key, "Gdispv"))	kt = GDISP;
    else if (!strcasecmp (key, "Gtilt"))	kt = GTILT;
    else if (!strcasecmp (key, "Gap"))		kt = GAP;
//    else  fprintf (stdbugf,
    else  fprintf (stderr,
        " ** Unrecognized keyword \"%s\" found; set to Void.\n", key );

    return  kt;
}    


char*	Keyname (Keyword k)
/* Return pointer to literal name of the keyword */
{
char	*c;
    switch (k) {
        case VOID:		c = "Void";		break;
        case END:		c = "END";		break;
        case Telescope:		c = "Telescope";	break;
        case FOCuser:		c = "Camera";		break;
        case Colimator:		c = "Colimator";	break;
        case Reimager:		c = "Reimager";		break;
        case GRAting:		c = "Grating";		break;
        case GRIsm:		c = "Grism";		break;
//	case DGRIsm:		c = "DGrism";		break;
        case EGRIsm:		c = "EGrism";		break;
        case ECHelle:		c = "Echelle";		break;
        case Reflection:	c = "Reflection";	break;
        case ANgle:		c = "Angle";		break;
        case Filter:		c = "Filter";		break;
        case DETector:		c = "Detector";		break;
        case INSTrument:	c = "Instrument";	break;
        case FOCLEN:		c = "FocLen";		break;
        case SCale:		c = "Scale";		break;
        case FOCCURV:		c = "FocCurv";		break;
        case Field:		c = "Field";		break;
        case Pixels:		c = "Pixels";		break;
        case AXIS:		c = "Axis";		break;
//        case ROTATION:		c = "Rotation";		break;
        case Lines:		c = "Lines";		break;
        case GLASS:		c = "Glass";		break;
        case RMOD:		c = "RMOD";		break;
        case TMOD:		c = "TMOD";		break;
        case WMOD:		c = "WMOD";		break;
	case ALIGNROT:		c = "ALIGNROT";		break;
	case PRIS1:		c = "Prism1";		break;
	case PRIS2:		c = "Prism2";		break;
	case GNORM:		c = "Gnormal";		break;
	case GDISP:		c = "Gdispv";		break;
	case GTILT:		c = "Gtilt";		break;
	case GAP:		c = "Gap";		break;
        default:		c = "-\?\?-";		break;
// The \? is used lest Linux interpret it as a trigraph (x?y:z)
    }
    return  c;
}


char*	typename (Etype k)
/* Return pointer to literal name of type */
{
char	*c;
    switch (k) {
	case Focuser:		c = "Focuser";		break;
	case Grism:		c = "Grism";		break;
//	case DGrism:		c = "2-Grism";		break;
	case EGrism:		c = "E-Grism";		break;
	case Glass:		c = "Glass";		break;
	case Echelle:		c = "Echelle";		break;
	case Instrument:	c = "Instrument";	break;
	case FloatValue:	c = "Float_Value";	break;
	case Secondary:		c = "Secondary";	break;
	case Detector:		c = "Detector";		break;
	case Gap:		c = "Gap";		break;
	case Other:		c = "Other";		break;
	default:
	case Void:		c = "Void";		break;
    }
    return c;
}



// void    showkey (FILE* f, Keyword k)
// /* Show keyword as a name; inverse of decode... */
// {
// char	*c;
//     c = keyname (k);
//     fprintf (f, "%s", c);
// }
// Above written to replace a similar function in an earlier program;
// is not actually used here, so it is commented.  It is retained as
// a comment in case it should become useful in the future.


Etype	type_of (Keyword k)
/* Find element (data) type for given keyword... */
{
Etype	v;

    switch (k) {
        case Telescope:
        case FOCuser:
        case Colimator:
        case Reimager:
	    v = Focuser;
	    break;
        case GRAting:
        case ANgle:
	    v = FloatValue;
	    break;
        case GRIsm:
	    v = Grism;
	    break;
        case EGRIsm:
	    v = EGrism;
	    break;
	case ECHelle:
	    v = Echelle;
	    break;
        case GLASS:
	    v = Glass;
	    break;
        case INSTrument:
	    v = Instrument;
	    break;
        case DETector:
	    v = Detector;
	    break;
	case GAP:
	    v = Gap;
	    break;
	case Reflection:
        case Filter:
	    v = Other;
	    break;
/* Following are secondary keywords having no valid e-type */
        case FOCLEN:
        case SCale:
        case FOCCURV:
        case Field:
	case Pixels:
        case AXIS:
//        case ROTATION:
        case Lines:
        case RMOD:
        case TMOD:
        case WMOD:
	case ALIGNROT:
	case PRIS1:
	case PRIS2:
	case GNORM:
	case GDISP:
	case GTILT:
        default:	// Any others must be secondary
	    v = Secondary;
	    break;
/* Special cases... */
        case VOID:
        case END:
	    v = Void;
	    break;
    }
    return  v;
}


/* Data storing subroutine(s) */

	/* - - - - - - - - - - *\
	|			|
	|    Data Storage	|
	|			|
	\* - - - - - - - - - - */

element*  newelement (char* name)
/* Obtain, initialize and return a element structure */
/* Put (copy of) name given into the new structure */
{
element	*p;	// return value

    p = xalloc (element);
    p->last = p->next = p;	// rather than NULL, make sure it's circular
    p->name = dynamstr (name);
    p->flag = 0;
    p->kw = 0;
    p->data = NULL;	// The data pointer, currently void
    p->head = NULL;

    return  p;
}


focdat*   newfdat (void)
{
focdat	*s;
    s = xalloc (focdat);
    s->scale = 0.0;
    s->rmod = NULL;
    s->tmod = NULL;
    s->wmod = NULL;
    s->curv = NULL;
    return  s;
}



datlist*	newdat (char* line)
/* Return a new datlist structure containing COPY OF line */
{
datlist	*d;
char	*cl;
int	j;

    j = strlen (line) + 1;
    cl = malloc (j);
    strcpy (cl, line);	// Includes terminating null
    d = xalloc (datlist);
    d->next = NULL;
    d->data = cl;
    d->eg = NULL;
    return  d;
}


/*  ====  Data Storage cleanup (destructor) routines  ==== */

datlist*	freedat (datlist* p)
/* Free the pointed datlist, and return its successor */
{
datlist *d;
    if (p == NULL) return NULL;
    d = p->next;
if (mdebug > 0) {
    printf (" d");
    fflush (stdout);
}
    free (p->data);
    free (p);
if (mdebug > 0) {
    printf (".");
    fflush (stdout);
}
    return  d;
}

void	freedatlist (datlist* p)
/* Free the total data list at p */
{
    if (p == NULL) return;
if (mdebug > 0) {
    printf (" dl ");
    fflush (stdout);
}
    while (p != NULL) p = freedat (p);
}

void	freepoly (poly* p)
/* Free the storage for polynomial */
{
    if (p == NULL) return;
if (mdebug > 0) {
    printf (" p");
    fflush (stdout);
}
    if (p->coef != NULL) free (p->coef);
if (mdebug > 0) {
    printf (".");
    fflush (stdout);
}
    free (p);
if (mdebug > 0) {
    printf (".");
    fflush (stdout);
}
}

void	freestgrism (stgrism* p)
{
    if (p == NULL) return;
if (mdebug > 0) {
    printf (" G ");
    fflush (stdout);
}
    free (p->glass);
    free (p);
}


void	freexgrism (egrism* p)
{
    if (p == NULL) return;
if (mdebug > 0) {
    printf (" EG ");
    fflush (stdout);
}
    free (p->inglass);
    free (p->eglass);
    free (p);
}


void	freestechel (stechel* p)
{
    if (p == NULL) return;
if (mdebug > 0) {
    printf (" E ");
    fflush (stdout);
}
    free (p->glass);
    free (p);
}

void	freefdat (focdat* p)
{
    if (p == NULL) return;
if (mdebug > 0) {
    printf (" F");
    fflush (stdout);
}
    freepoly (p->rmod);
    freepoly (p->tmod);
    freepoly (p->wmod);
    freepoly (p->curv);
if (mdebug > 0) {
    printf (":");
    fflush (stdout);
}
    free (p);
}


gapdat*	newgap (void) {
gapdat	*p;
    p = xalloc (gapdat);
    p->last = p->next = p;	// rather than NULL, make sure it's circular
    p->name = NULL;
    p->instr = NULL;
    p->ns = 0;
    p->x1u = p->x1b = 0.0;
    p->x2u = p->x2b = 0.0;
    return  p;
}

gapdat*	pop_gap (gapdat* p) {
    return (gapdat*)G_pop_c (p);
}

gapdat*	push_gap (gapdat* a, gapdat* b) {
    return (gapdat*)G_push_cl (a, b);
}

void	freegap (gapdat* p)
{
gapdat	*q;
    while (p != NULL) {
	q = pop_gap (p);
	free (p->name);
	free (p->instr);
	free (p);
	p = q;
    }
}

element*  pop_element (element* p)
/* Pop the element and return the "next" one */
{
element	*l, *r;
    if (p == NULL) return NULL;
    l = p->last;
    if (l == p) l = NULL;
    r = p->next;
    if (r == p) r = NULL;
    if (l != NULL) l->next = r;
    if (r != NULL) r->last = l;
    p->last = p->next = p;	// For safety
    if (r == NULL) r = l;	// That's what we return
    return  r;
}

element*	push_element (element* b, element* a)
/* Put the queue or element b into queue a, so that it
  preceeds a.  Return the new head (a) pointer. */
{
element	*p;
element	*r;

/* Take care of null cases */
    if (b == NULL) return a;
/* Find queue ends such that a...p and b...r are merged, into
  something like  ...p - b...r - a... ; this replaces the links
  ...p - a... and ...r - b... with the p - b and r - a links. */
/* If a has null predicessor, we link ...0 - b...r - a...; if
  b has null predicessor, find the next end and use that. */
/* Link it up */
    if (a == NULL) return b;
    p = a->last;		/* ...p - a... */
    if (p == NULL) {	/* Scan to end */
	for (p=a; p->next != NULL; p = p->next) {
	    if (p->next == a) break;
// BAD CONSTRUCT BREAK
	}
    }
    if (p == NULL) p = a;	/* maybe a pointers got nulled ? */

    r = b->last;		/* ...r - b... */
    if (r == NULL) {	/* Scan to end */
	for (r = b; r->next != NULL; r = r->next) {
	    if (r->next == b) break;
// BAD CONSTRUCT BREAK
	}
    }
    if (r == NULL) r = b;	/* Maybe b pointers got nulled ? */
//    if (p != NULL) p->next = b;
    p->next = b;
    b->last = p;
    r->next = a;
    a->last = r;		/* ...p - b...r - a... */
    return a;
}



element*  free_element (element* p)
/* Remove and free storage of element at p; return remainder */
{
element	*r;
    if (p == NULL) return NULL;
if (mdebug > 0) {		// debug
    printf ("  - %s (%s) ->", p->name, Keyname (p->kw) );
    fflush (stdout);
}
    r = pop_element (p);
    free (p->name);
    freedatlist (p->head);
/* The hard part -- free the type-dependent data under the
  void (many possible types) pointer at "data"...  */
    switch (type_of(p->kw)) {
	case Focuser:
	    freefdat ( (focdat*) p->data );
	    break;
	case Grism:
	    freestgrism ( (stgrism*) p->data );
	    break;
//	case DGrism:
//	    freesdgrism ( (sdgrism*) p->data );
//	    break;
	case EGrism:
	    freexgrism ( (egrism*) p->data );
	    break;
	case Echelle:
	    freestechel ( (stechel*) p->data );
	    break;
	case Instrument:
/* The data points to head, which is already freed... */
	    p->data = NULL;
	    break;
	case Glass:
	    freepoly ( (poly*) p->data);
	    break;
	case FloatValue:
	    free ( (double*) p->data);
	    break;
	case DETector:
	    free ( (detector*) p->data);
	    break;
	case Gap:
	    freegap ( (gapdat*) p->data);
	    break;
	default:
	    free (p->data);
	    break;
    }
if (mdebug > 0) {
    printf (" - ");
    fflush (stdout);
}
    if (p != NULL) free (p);
if (mdebug > 0) {
    printf (" Done.\n");
    fflush (stdout);
}
    return  r;
}


/* REMOVING FLOATVALUE, INTVALUE; THEY ARE NOW IN IOUTILS */


double*	newfv (char* c)
{
double	*p;
/* Get a double storage element, and read c data into it */
    p = xalloc (double);
    *p = floatvalue (c);
    return  p;
}



/*  ===  Parsing Subroutines  ===  */

// now in ioutil -- include ioutil.h and link with imacs library...


vect3	readv3 (char* text) {
/* Read the text as a 3-vector and return... */
vect3	result;
int	k;

/* Read up to 3 float values... */
    k = sscanf (text, "%lf %lf %lf",
	&result.x, &result.y, &result.z);

/* Clear any values not found on reading */
    if (k < 3) result.z = 0.0;
    if (k < 2) result.y = 0.0;
    if (k < 1) result.x = 0.0;
    return  result;
}


	/* - - - - - - - - - - *\
	|			|
	|  Polynomial Support	|
	|			|
	\* - - - - - - - - - - */

poly*	read_poly (char* text)
/* Read the polynomial text, and generate an approparite
  polynomial structure, returning a pointer to the new structure. */
/* Use non-destructive parsing, from above */
{
char	token[64];
char	*c;
int	k;
int	ord=0;
poly	*p; // was.. = NULL;
double	*d; // was.. = NULL;
int	i, j;

/* Separate the first token, which is the order */
    c = parse (text, token, WHITESPACE);

/* Read the order */
    k = sscanf (token, "%d", &ord);
    if (k < 1) return NULL;

/* Get the basic structure, and also an array to hold
  the coeficients. */
    p = xalloc (poly);
    p->order = ord;
    k = ord + 1;
    d = (double *) malloc (k * sizeof (double));
    p->coef = d;

//    for (i=0; i < k; i++) d[i] = 0.0;

/* Loop to get tokens, read them, and fill the coeficients. */
    for (i=0; i < k; i++) {
	c = parse (c, token, WHITESPACE);
	j = sscanf (token, "%lf", &(d[i]) );
	if (j == 0) d[i] = 0.0;
    }

/* Drop any null or zero polynomial... */
    if (ord <= 0) {
	if (*d == 0.0) {
	    free (d);
	    free (p);
	    p = NULL;
	}
    }

    return  p;
}

/* Display of polynomial - used in debug routines */

void	showpoly (poly* p, FILE* f)
/* Display debug version of polynomial on file f */
{
int	i;
double	*c;
/* Write the order, and ready for coeficients */
    fprintf (f, "(%d){", p->order);
    for (i=0, c = p->coef; i <= p->order; i++) {
	fprintf (f, "%14.7e", *(c++));
	if (i < p->order) fprintf (f, ", ");
	if (i > 32) break;	// debug - emergency halt
    }
    putc ('}', f);
}

/*  ====  Polynomial Evaluation  ====  */

double  eval (poly*  p, double arg)
/* Evaluate the polynomial at argument arg */
{
double  v;
int     i;

    if (p == NULL) return 0.0;
    for (i=p->order,v=p->coef[i],i--; i>=0; i--) {
        v *= arg;
        v += p->coef[i];
    }
    return  v;
}    



	/* - - - - - - - - - - - - - - *\
	|				|
	|   Optical Data File Reading	|
	|				|
	\* - - - - - - - - - - - - - - */


/*  NOTE -- Following is a standard "user program" for the generalized
  configuration file reading code in ioutils library.
  It is passed by the call to do the reading, which originates outside
  this module...   */


void*	optstore (char* line, void* head)
/* Imitation data storage for optics - for debug of reading */
{
#define  KEYSIZE  64
static  element*  curelm=NULL;
static  datlist*  curend=NULL;
char	keyword[KEYSIZE];
char	name[KEYSIZE];
char	*c;
// element	*p;
datlist	*d;
stgrism	*sg;
//sdgrism	*sdg;
egrism	*sxg;
stechel	*se;
detector *dp;
static	int	dc = 0;		// debug count -- debug only

Keyword	kw;
Etype	ke;
Etype	ce;
double	v;
int	j;
focdat	*fd;
poly	*py;
int	ks;  // debug flag for secondary data storing
/* "ks" is used to indicate special secondary data lines; it is incremented
  when secondary data is stored, so that the secondary line is then not
  placed in the list of secondary data lines.  The goal is to store all
  valid information this way, so secondary lines with ks==0 should 
  eventually be flagged as errors.  */

/* debug reading version.
  if a ".", add line to element; otherwise, keyword names
  the next element.  */

/* Report the line without a dot. */
/* First token is the type of thing; second token is the name, use
  for name of element. */

/* Later, we extract polynomials from those things having them
  as rest of line, and report them. */

/* Also later, we extract such things as focuser, colimator, etc., and
  save their polynomials.  Also save instrument, angle, and such things
  in separate stacks. */

/* Devise a storage scheme, and a debug program to report it; and a
  lookup program to compute things from it. */

/* If line is NULL, return head unconditionally */
    if (line == NULL) {
	if (mdebug > 0) printf (" (Processed %d lines.)\n", dc);  //debug
	return  head;
    }

/* If head is null, that means we just get a new one and have it
  linked to itself rather than to the head pointer.  Once a head
  pointer is made, we link to the predecessor of that which is pointed
  to by the head pointer. */

/* The void* is actually a gdata* pointer here.  It is "void" in the
  calling program, the general reading program.  */

/* For DEBUG, we just echo the input line, with some stuff around
  it to show what we have to work with... */
	dc++;	// debug
//	if (dc < 15) printf (" >\"%s\"\n", line);	//debug

/* Remember where we are - current element is static storage;
  start with NULL, and if NULL, ignore as being not in a element. */
/* Also, remember the current last entry in the element datlist */


/* Parse (non-destructive) the initial keyword in the line and find the
  location of the follow-on part of the line.  The "parse" subroutine
  will set *c to the "rest" of line, and copy the keyword into the
  provided "keyword" buffer; the delimiter set here is WHITESPACE */

    c = parse (line, keyword, WHITESPACE);


// Put the result in last place in the appropriate head structure...
// use the push-element for that.

// Don't yet have push/pop element; also need a kill element code.


/* If an optional dot is present at start of line, and has been parsed
  as the keyword, it is now removed, and the rest of the line re-parsed,
  with the "line" pointer reset to the new start of line.  If later the
  line pointer is found to be needed, a copy should be made. */
    if (!strcmp (keyword, ".") ) {
	line = c;			// Drop the dot, effectively
	c = parse (line, keyword, WHITESPACE);	// re-parse
    }

/* Now, decode the keyword and find the keyword level from special purpose
  subroutines above. */
    kw = Decode_Keyword (keyword);
    ke = type_of (kw);	// Looking for Secondary keywords
    ks = 0;  // not special (yet)


/* Store the data in the line, depending on the keyword; we do this
  differently for primary and secondary keywords...
  This section is customized for the storage design...
..*/

	/* - - - - - - - - - - *\
	|			|
	|    Void  Keywords	|
	|			|
	\* - - - - - - - - - - */

    if (ke == Void) {		// Other or void keyword
/*	if (kw != END)  */
	printf ("  Void Keyword \"%s\", rest \"%s\".\n", keyword, c);

	/* - - - - - - - - - - *\
	|			|
	|  Secondary  Keywords	|
	|			|
	\* - - - - - - - - - - */

    } else if (ke == Secondary) {		// sub-keyword
/* Data from secondary keyword is stored in the structure of the
  current primary keyword, if the type is correct.  Otherwise, ignore
  any out-of-place secondary keywords. */

	if (dc < 0) {	// debug section
	    printf (" ..Sub-keyword \"%s\", rest \"%s\"\n", keyword, c);
	    fflush (stdout);	// debug
	}

	switch (kw) {
	    case FOCLEN:	// recip. if colimator
	    case SCale:	// directly into focuser
		v = floatvalue (c);
//		if (kw == FOCLEN && type_of(curelm->kw) == Colimator) v = 1.0/v;
// NEED TO KEEP MASTER ELEMENT KEYWORD RATHER THAN TYPE TO DECIDE
// IF IT WAS COLIMATOR!!
		if (type_of(curelm->kw) == Focuser) {
		    fd = curelm->data;
		    fd->scale = v;
		}
		ks++;  // special, debug
		break;

/* Modifier "AXIS" applies to grating, grism and angle data */
	    case AXIS:	// integer value
		if (curelm->kw == GRAting ||
		    curelm->kw == GRIsm ||
//		    curelm->kw == DGRIsm ||
		    curelm->kw == EGRIsm ||
		    curelm->kw == Reflection ||
		    curelm->kw == ANgle ) {
		    j = intvalue (c);
		    if (j & 1) curelm->flag |= 2;
		}
		break;

/* Modifier "LINES" can apply to grating or grism data */
	    case Lines:	// float value
		if (curelm->kw == GRAting) 
		    *((double*)(curelm->data)) = floatvalue (c);
		else if (curelm->kw == GRIsm)
		    ((stgrism*)curelm->data)->lpmm = floatvalue (c);
		else if (curelm->kw == EGRIsm)
		    ((egrism*)curelm->data)->lpmm = floatvalue (c);
//		else if (curelm->kw == DGRIsm)
//		    ((sdgrism*)curelm->data)->grism.lpmm = floatvalue (c);
		else if (curelm->kw == ECHelle)
		    ((stechel*)curelm->data)->lpmm = floatvalue (c);
		break;

/* Echelle vector secondary keywords... */
	    case PRIS1:
		if (curelm->kw == ECHelle)
		    ((stechel*)curelm->data)->pris1 = readv3 (c);
		break;

	    case PRIS2:
		if (curelm->kw == ECHelle)
		    ((stechel*)curelm->data)->pris2 = readv3 (c);
		break;

	    case GNORM:
		if (curelm->kw == ECHelle)
		    ((stechel*)curelm->data)->grnor = readv3 (c);
		break;

	    case GDISP:
		if (curelm->kw == ECHelle)
		    ((stechel*)curelm->data)->grdis = readv3 (c);
		break;

#ifdef  GTILTMOE
	    case GTILT:
		if (curelm->kw == ECHelle)
		    ((stechel*)curelm->data)->gtilt = floatvalue (c);
		break;
#endif


// Remove all reference to ROTATION

	    case RMOD:	// polynomials for focuser...
	    case TMOD:
	    case WMOD:
	    case FOCCURV:
		if (type_of(curelm->kw) == Focuser) {
		    py = read_poly (c);
		    fd = curelm->data;
		    if      (kw == RMOD) fd->rmod = py;
		    else if (kw == TMOD) fd->tmod = py;
		    else if (kw == WMOD) fd->wmod = py;
		    else if (kw == FOCCURV) fd->curv = py;
		}
		ks++;  // special, debug
		break;

/* ALIGNROT = alignment rotation value, in degrees; store as a
  special line in curelm->head as a datlist entry */
	    case ALIGNROT:	// Special entry
		ce = type_of (curelm->kw);
		if (ce == Instrument) break;
/* We don't do this for an entire instrument - data conflict */
		if (ce == Void || ce == Focuser ||
			ce == Glass || ce == Other) break;
/* We don't do this for those types either; intended to do this
  only for Grism, FloatValue (i.e. Grating) only */
		v = floatvalue (c);	// The rotation quantity
		v *= Degree;		// Convert to radians
		d = xalloc (datlist);	// A datlist block
		d->next = NULL;
		d->data = malloc (16);	// 16 characters in a block
		d->eg = NULL;		// Should never be used!
		memset (d->data, 0, 16);
		strcpy (d->data, "ARot");	// Special identifier
		memcpy (&(d->data[4]), &v, 8);	// The angle in radians
		if (curend == NULL) curelm->head = d;
		else		curend->next = d;
		curend = d;
		break;

/* To use this (below in transform) see if the dispersion element has
  a "head" pointer, and if that pointer has a non-NULL data field, and
  if that data field has the special key "ARot" in its first 4 characters.
  If so, the rotation error value is the next 8 characters, use memcpy
  to extract the value.  Scan the "head" queue to find the first such
  value, and if none found the rotational error value is zero. */

/* Field values for detector are stored in the field part of
  the detector structure pointed to by the main element */
	    case Field:	// unknown yet
		if (curelm->kw == DETector) {
double	v1, v2;	// scratch
//intval	w;	// scratch
bbox	*bb;	// pointer to bounds to set
/* The detector field entries are of this form:
  FIELD  n  x1 x2  y1 y2
  where n is either N or S for normal or N&S orientation, or possibly
  D for default orientation.  The x1 x2 is an intval, as is y1 y2.
  Those interval values are stored in the appropriate bounding box
  area for the field description.  */

/* Decode the n x1 x2 y1 y2 bounding box descriptor */
		    c = parse (c, name, WHITESPACE);
		    if (strncasecmp (name, "N", 1) == 0)
			bb = &( ((detector*)(curelm->data))->bbn );
		    else if (strncasecmp (name, "S", 1) == 0)
			bb = &( ((detector*)(curelm->data))->bbs );
		    else if (strncasecmp (name, "D", 1) == 0)
			bb = &( ((detector*)(curelm->data))->bbd );
		    else  bb = NULL;
		    if (bb != NULL) {
			c = parse (c, name, WHITESPACE);
			v1 = floatvalue (name);
			c = parse (c, name, WHITESPACE);
			v2 = floatvalue (name);
			bb->x = interval (v1, v2);
			c = parse (c, name, WHITESPACE);
			v1 = floatvalue (name);
			c = parse (c, name, WHITESPACE);
			v2 = floatvalue (name);
			bb->y = interval (v1, v2);
		    }
		}	// end primary keyword selections

		break;

/* Pixels values for detector are stored in the pcount part
  of the detector structure pointed to by the main element */
	    case Pixels:	// unknown yet
		if (curelm->kw == DETector) {
// Modifiers X (value) and Y (value)
//printf (" >>field(detector) parsing: \"%s\".\n", c);
//printf (" .. detector \"%s\".\n", curelm->name);
//fflush (stdout);	// debug
		    while (strlen(c) > 0) {
			c = parse (c, name, WHITESPACE);
			if (!strcasecmp(name, "X")) {
			    c = parse (c, name, WHITESPACE);
			    v = floatvalue (name);
			    ((detector*)(curelm->data))->pcount.x =  v;
			    ks++;	// Flag as translated
			} else if (!strcasecmp(name, "Y")) {
			    c = parse (c, name, WHITESPACE);
			    v = floatvalue (name);
			    ((detector*)(curelm->data))->pcount.y =  v;
			    ks++;	// Flag as translated
			}	// end x/y detection
		    }	// end loop on parameters
		}	// end primary keyword selections

		break;

//	    case FOCCURV:	// forget it
	    default:
		break;
	}

// for debug, we still add the line as a data list value;
// in the real thing, this is not used.
	if (curelm != NULL && ks == 0) {  // Only add unused lines
// can't add to unknown element
	    d = newdat (line);
	    if (curend == NULL) curelm->head = d;
	    else		curend->next = d;
	    curend = d;
	}

	/* - - - - - - - - - - *\
	|			|
	|       Gap  Data	|
	|			|
	\* - - - - - - - - - - */

    } else if (ke == Gap) {		// special data
// At this point, we have a gap line; keyword is "gap" and
// c points to the rest of the line.
// line contains:  inst ns name x1u x1b x2u x2b
// Do this:
//  decode the rest of line into a new gapdat structure element
//  use newgap function to create it.
//  use dynamstr to get dynamic copies of string elements
//  find the special element "GapHead", if none, make one
//  data field of special element is head to gap structure queue
//  push new data block into the element's data block.
gapdat	*pg;
element	*gx;
/* Read and store the GAP line data here */
	pg = newgap ();
	c = parse (c, name, WHITESPACE);	// Instrument name
	pg->instr = dynamstr (name);
	c = parse (c, name, WHITESPACE);	// N&S state
	pg->ns = intvalue (name);
	c = parse (c, name, WHITESPACE);	// Gap name
	pg->name = dynamstr (name);
	c = parse (c, name, WHITESPACE);	// X-value
	pg->x1u = floatvalue (name);
	c = parse (c, name, WHITESPACE);	// X-value
	pg->x1b = floatvalue (name);
	c = parse (c, name, WHITESPACE);	// X-value
	pg->x2u = floatvalue (name);
	c = parse (c, name, WHITESPACE);	// X-value
	pg->x2b = floatvalue (name);
	ks++;		// flag as translated

/* Store gap structure into the special gap queue header */
	gx = find_element (head, GAPHEAD);
	if (gx == NULL) {	// Make special element
	    gx = newelement (GAPHEAD);
	    gx->kw = kw;
	    head = push_element (gx, head);
	}
	gx->data = push_gap (pg, (gapdat*)gx->data);
// We do not store curelm for gap data, that's why it is special

	/* - - - - - - - - - - *\
	|			|
	|   Primary Keywords	|
	|			|
	\* - - - - - - - - - - */

    } else {		// Other type, must be actual keyword...
	c = parse (c, name, WHITESPACE);
	curelm = newelement (name);
	curelm->kw = kw;
	curend = NULL;		// start a new element queue
    if (dc < 0) {	// debug
	printf (" ..Actual keyword \"%s\", name \"%s\".\n", keyword, name);
	fflush (stdout);	// debug
    }

	head = push_element (curelm, head);

/* Decode some data following the name... */

/* Obtain storage for the element data; depends on element type */
	switch (ke) {		// On type of storage
	    case Instrument:
		while (*c != AN) {
		    c = parse (c, name, WHITESPACE);
		    d = newdat (name);
		    if (curend == NULL) curelm->head = d;
		    else		curend->next = d;
		    curend = d;
		}
		curelm->data = curelm->head;	// The list...
		break;

	    case Focuser:
/* Obtain a focuser data structure for these */
//		curelm->data = newfdat();
		fd = newfdat ();
		curelm->data = fd;
		v = floatvalue (c);
// Scale is v, unless we have collimator, when it is 1/v
		fd->scale = (kw == Colimator) ? 1.0/v : v;
		break;

	    case FloatValue:
/* These use just a floating value */
// Currently float values are GRATING, ANGLE
		curelm->data = newfv (c);
		break;

	    case Grism:
// was  3 things -- lpmm, name of glass, angle (deg) in that order...
// Read 3 things -- name of glass, angle (deg), lpmm in that order...
		sg = xalloc (stgrism);
		curelm->data = sg;

		c = parse (c, name, WHITESPACE);
		sg->glass = malloc (strlen(name) + 1);
		strcpy (sg->glass, name);
		sg->eg = find_element (head, name);

		c = parse (c, name, WHITESPACE);
		sg->angle = floatvalue (name);

		sg->lpmm = floatvalue (c);

		break;
// Order above changed due to "lines" modifier being allowed; primary
// data of glass and grism angle is needed first.

// We may want to define secondary keyword modifiers for the Grism type,
// To give the grism angle, and the glass type used in the grism.
// These can't be same as primary keywords Glass or Angle which have
// different meanings in this file.


	    case EGrism:
// Add reading of extended grism information...
// Left to right:  angle, glass, angle, lpmm, glass, angle
// Things left off right end are null
// glass name can be "null" -- case insensitive
		sxg = xalloc (egrism);
		curelm->data = sxg;

		c = parse (c, name, WHITESPACE);
		sxg->inangle = floatvalue (name);

		c = parse (c, name, WHITESPACE);
		if (!strcasecmp(name, "null") || strlen(name) == 0) {
		    sxg->inglass = NULL;
		    sxg->ig = NULL;  // or sxg->ig = &sxg;
		} else {
		    sxg->inglass = malloc (strlen(name) + 1);
		    strcpy (sxg->inglass, name);
		    sxg->ig = find_element (head, name);
		}

		c = parse (c, name, WHITESPACE);
		sxg->gangle = floatvalue (name);

		c = parse (c, name, WHITESPACE);
		sxg->lpmm = floatvalue (name);

		c = parse (c, name, WHITESPACE);
		if (!strcasecmp(name, "null") || strlen(name) == 0) {
		    sxg->eglass = NULL;
		    sxg->eg = NULL;  // or sxg->eg = &sxg;
		} else {
		    sxg->eglass = malloc (strlen(name) + 1);
		    strcpy (sxg->eglass, name);
		    sxg->eg = find_element (head, name);
		}

		c = parse (c, name, WHITESPACE);
		sxg->eangle = floatvalue (name);

		break;


	    case Echelle:
// Just the name of glass, as done for grism above.
// the vector values are done by secondary keywords.
		se = xalloc (stechel);
		curelm->data = se;

		c = parse (c, name, WHITESPACE);
		se->glass = malloc (strlen(name) + 1);
		strcpy (se->glass, name);
		se->eg = find_element (head, name);
#ifdef  GTILTMOE
		se->gtilt = 0.0;
		se->gtnor.x = -999.0;	// special flag value
#endif
		break;

	    case Glass:
// Read the polynomial in "c" into a poly, and link to it...
		py = read_poly (c);
		curelm->data = py;
		break;

	    case Detector:
// Read detector: array size in mm, pixel count on side, vignetting
		dp = xalloc (detector);
		curelm->data = dp;
		memset (dp, 0, sizeof(detector));

		c = parse (c, name, WHITESPACE);
		v = floatvalue (name) / 2.0;	// bounding box
		dp->bb.x.lo = dp->bb.y.lo = -v;
		dp->bb.x.hi = dp->bb.y.hi =  v;
// Set bbd, bbn and bbs default bounding boxes too.
// Change that to setting bbd as default bound box.
		dp->bbs = dp->bbn = dp->bbd = dp->bb;
// Add "field" modifiers to store bbn and bbs bounds too.
// Initialize all bb values when detector structure is obtained.

		c = parse (c, name, WHITESPACE);
		v = floatvalue (name);		// pixels
		dp->pcount = make2vect (v, v);

		c = parse (c, name, WHITESPACE);
		dp->vrad = floatvalue (name);

		break;

// For "REFLECTION" - leave data pointer null; we apply axis modifier
// to element structure.  Otherwise use grating with 0.0 lines/mm.

// ITEMS BELOW NEED TO BE PROPERLY STORED. (keyword values commented)
//	    case Reflection:
//	    case Filter:
	    default:
		break;
	}


    }	// End of keyword level if nest

// See ../../gmap/dtest.c for example of keyword decoding
// from which some of this code is derived...

/* Return the head pointer */
    return  head;
}



/*  ====  Optics Structure utility code  ====  */

	/* - - - - - - - - - - *\
	|			|
	|    Element Queue	|
	|			|
	\* - - - - - - - - - - */


element*	find_element (element* head, char* name)
/* Find the named element in the list */
{
element*  p = NULL;
    if (head == NULL) return p;
    C_loop (p, head) { if (! strcasecmp (p->name, name)) break; }
    return  p;
}

int	count_elements (element* head)
{
    return  listcount (head);
}

int	count_types (element* head, Keyword kw)
{
element*  p;
int	k = 0;
    if (head == NULL) return k;
    C_loop (p, head) { if (p->kw == kw) k++; }
    return  k;
}



/*  ====  Optic Structure error handling  ====  */

	/* - - - - - - - - - - - - - - *\
	|				|
	|     Element Queue Errors	|
	|				|
	\* - - - - - - - - - - - - - - */


void	error_title (int te)
{
    if (te > 0) return;
    printf ("\n ** Errors in Optics Data found -- \n\n");
}


int	find_errors (element* head)
/* Find (and report) errors in element queue */
{
element	*p;
stgrism	*sg;
//sdgrism	*sdg;
egrism	*sxg;
stechel	*se;
datlist	*d;
element *x;
int	j;
int	k;
Keyword	t;
int	te = 0;	// error count

    if (head == NULL) {
	error_title (te);
	printf (" ** No element queue present.\n");
	return 1;
    }

    C_loop (p, head) {
//	te += display_element (p, 1);
// test the element, just as in display_element below...

	if (p->head == NULL && p->data == NULL) {
	    error_title (te);
	    printf ("\n No data in element \"%s\" (%s).\n", 
		    p->name, typename(type_of(p->kw) ) );
	    te++;
	}

	switch (type_of(p->kw)) {

	    case Grism:
	    sg = (stgrism*) p->data;
// VERIFY THAT THE GLASS EXISTS
	    x = find_element (p, sg->glass);
	    if (x == NULL) {
		error_title (te);
		printf ("  ** For Grism \"%s\", glass \"%s\" does not EXIST!\n",
			p->name, sg->glass);
		te++;
	    } else {
		if (sg->eg == NULL) sg->eg = x;
	    }
	    break;

	    case EGrism:
	    sxg = (egrism*) p->data;
// Verify that glasses exist as above...
	if (sxg->inglass != NULL) {
	    x = find_element (p, sxg->inglass);
	    if (x == NULL) {
		error_title (te);
		printf (
		    "  ** For EGrism \"%s\", glass \"%s\" does not EXIST!\n",
			p->name, sxg->inglass);
		te++;
	    } else {
		if (sxg->ig == NULL) sxg->ig = x;
	    }
	}
	if (sxg->eglass != NULL) {
	    x = find_element (p, sxg->eglass);
	    if (x == NULL && strcmp(sxg->inglass, sxg->eglass)) {
		error_title (te);
		printf (
		    "  ** For EGrism \"%s\", glass \"%s\" does not EXIST!\n",
			p->name, sxg->eglass);
		te++;
	    } else {
		if (sxg->eg == NULL) sxg->eg = x;
	    }
	}
	    break;

	    case Echelle:
	    se = (stechel*) p->data;
// VERIFY THAT THE GLASS EXISTS
	    x = find_element (p, se->glass);
	    if (x == NULL) {
		error_title (te);
		printf ("  ** For Echelle \"%s\", ", p->name);
		printf ("glass \"%s\" does not EXIST!\n", se->glass);
		te++;
	    } else {
		if (se->eg == NULL) se->eg = x;
	    }
	    break;

	    case Instrument:
/* Verify existance of each element in the datlist */
	    for (j=1, d = p->head; d != NULL; d = d->next, j++) {
//		printf ("  %d) %s", j, d->data);
/* If we have GRISM or GRATING, count number of such in element list */
		if      (! strcasecmp (d->data,  "GRISM" )) t = GRIsm;
//		else if (! strcasecmp (d->data,  "DGRISM")) t = DGRIsm;
		else if (! strcasecmp (d->data, "GRATING")) t = GRAting;
		else  t = VOID;

		if (t != VOID) {
		    k = count_types (head, t);
		    if (k < 1) {
			error_title (te);
			printf (
	"  ** For Instrument \"%s\", no elements of type \"%s\" EXIST!\n",
	p->name, Keyname (t) );
			te++;
		    }
		} else {
/* Should be actual element name, verify its existance */
// Special case:  skip alignment angle  ARot
		    if (strncmp (d->data, "ARot",4)) {
			x = find_element (p, d->data);
			if (x == NULL) {
			    error_title (te);
			    printf (
	"  ** For Instrument \"%s\", no element named \"%s\" EXISTS!\n",
		p->name, d->data );
			    te++;
			}
		    }	// End special case
		}	// End data name test
	    }		// End loop on datlist for Instrument case
	    break;
/* For lint and its clones... */
	    default:
	    break;
	}	// end switch
    }		// end loop on elements in queue
    return  te;
}



	/* - - - - - - *\
	|		|
	|   Refraction	|
	|		|
	\* - - - - - - */


vect3	snell (vect3 a, vect3 s, double ndx) {
/* Refract ray a at surface normal s into medium of relative index ndx */
vect3	A;	// Input normalized
vect3	S;	// Input normalized;
vect3	P;	// Input perpendicular vector
vect3	U;	// Input parallel vector
vect3	V;	// Output parallel vector
vect3	Q;	// Output perpendicular vector

double	d;	// Perpendicular vector, fraction of S
double	x;	// Normalization


/* Normalize input and surface normal vectors to unit length */
    A = norm3vect (a);
    if (ndx == 1.0) return A;
    if (ndx == 0.0) return A;
    S = norm3vect (s);

/* Compute components of input vector normal to and parallel to surface */
    d = dot3vect (A, S);
    P = mul3vect (S, d);
    U = sub3vect (A, P);

/* Find components of result vector, normal and parallel */
    V = mul3vect (U, 1.0/ndx);
    x = 1.0 - dot3vect (V, V);
    if (x < 0.0) {	// Total internal reflection case
	Q = mul3vect (P, -1.0/ndx);
    } else {
	x = sqrt (x);
	if (d < 0.0) x = -x;
	Q = mul3vect (S, x);
    }

/* Return the combined result vector */
    return norm3vect (sum3vect (V, Q));
// The normalization is needed in the internal reflection case to
// obtain a unit vector.  Could use Q = -P, V = U and get it right.
// Normalization of Q + V is then not needed.
}

vect3	diffract (vect3 a, vect3 s, vect3 g, double v) {
/* Compute diffracted output vector; normal to grating is s, dispersion
  vector in grating plane is g; input vector is a; value v is the
  product of wavelength x ruling_density x order.  Output is in the
  direction of s away from grating, and (g, v) define the direction
  of dispersion.  */
vect3	A;	// Input normalized
vect3	S;	// Grating normal unit vector;
vect3	G;	// Dispersion unit vector;
vect3	P;	// Input perpendicular vector
vect3	U;	// Input parallel vector, dispersion direction
vect3	V;	// Input parallel vector, grove direction

double	d;	// Perpendicular vector, fraction of S
double	x;	// Normalization
double	z;	// normalization of output

/* Normalize the input vectors to unit length */
    A = norm3vect (a);
    S = norm3vect (s);
    G = norm3vect (g);

/* Resolve input vector along grating coordinate directions */
    d = dot3vect (A, S);
    P = mul3vect (S, d);
    V = sub3vect (A, P);
    U = mul3vect (G, dot3vect (V, G));
    V = sub3vect (V, U);

/* Change the resolved inputs to output normal vectors */
    U = sum3vect (U, mul3vect (G, v));
    x = 1.0 - dot3vect (U,U) - dot3vect (V,V);
    z = (x < 0.0) ? fabs(d) : sqrt (x);
    P = mul3vect (S, z);
    return norm3vect (sum3vect (sum3vect(U, V), P));
}


	/* - - - - - - - - - - *\
	|			|
	|    Grating  Angle	|
	|			|
	\* - - - - - - - - - - */

double  Gr_Angle ( double lambda, double density, double Camangle, int order)
/* Compute the grating angle, radians */
/* Return 0.0, an otherwise out of range value, if no
  solution can be found. */
{
double	warg;
double	cx;

double	fa, fb;
double  ws;
double	cv, dv;
double	qc;

double	gx;
double	px;

/* Consolidate the arguments for wavelength, grating density, order */
    warg = lambda * Angstrom * order * density;

/* The IMACS camera angle is 180 + the angle through which the
  reflection is done... */
    cx = Camangle - M_PI;

/* Now, the equation which needs to be solved is:
  warg = sin (gx) - sin (cx - gx)
  and we note:  sin (cx - gx) = sin(cx)*cos(gx) - sin(gx)*cos(cx)
  and:  cos^2(gx) + sin^2(gx) = 1.0  */
/* Resulting quadratic equation is:
  2*sin^2(gx) - 2*warg*sin(gx) + warg^2/(1+cos(cx)) - (1-cos(cx)) = 0
  which is then solved for sin(gx).  */

    fa = 0.5 * warg;
    qc = cos (cx);
    ws = warg * warg;
    cv = qc - 1.0 + ws / (1.0 + qc) ;
    dv = ws - 2.0 * cv;
    if (dv < 0.0) return 0.0;	// Error condition
    fb = 0.5 * sqrt (dv) ;
// ALSO RETURN 0.0 IF ARG. TO ASIN IS OUTSIDE -1 TO +1
    gx = asin (fa+fb);
    px = asin (fa-fb);

// gx will usually be the angle; however px is also a solution;
// the one in range 22.5 to 67.5 degrees should be selected.
// px will be (gx - cx), and in range -22.5 to +22.5, usually.
    if (px > 22.5*Degree && px < 67.5*Degree && gx > 67.5*Degree)
	gx = px;

    return  gx;
}


	/* - - - - - - - - - - *\
	|			|
	|  Index of Refraction	|
	|			|
	\* - - - - - - - - - - */


static	double	indxref (char* gn, element** ge, element* gs, double rls)
/* Function declared static as it is only used in general transform */
/* Find index of refraction of glass.
  Argument gn is a glass name, and *ge is a pointer to the optical
  element of that glass.  The rls value is 1.0/(wavelength squared)
  with the wavelength in microns, and gs is an element queue pointer.
  The glass element will be found if *ge is NULL, and stored in *ge.
  The index will return as 1.0 if the glass name is null, or if the
  rls value is 0.0, or if the glass is not in the element queue.
  Otherwise, compute the index from the polynomial data present in
  the valid element queue member using rls.   */
{
double	nx = 1.0;
element	*gx;

/* Take out some null cases first */
    if (gn == NULL) return  nx;
    if (rls == 0.0) return  nx;
/* Allow ge to be NULL if we don't want to keep glass pointer */
    if (ge == NULL) {
	gx = find_element (gs, gn);
    } else {
	gx = *ge;
	if (gx == NULL) gx = *ge = find_element (gs, gn);
    }
    if (gx != NULL) nx = sqrt(eval((poly*)(gx->data), rls)/rls);
// YA KNOW, WE REALLY SHOULD CHECK THAT GX ELEMENT FOR BEING
// OF THE TYPE GLASS BEFORE DOING THAT...
    return  nx;
}


	/* - - - - - - - - - - - - - - *\
	|				|
	|    Generalized  Transform	|
	|				|
	\* - - - - - - - - - - - - - - */

/*  ----  Generalized Transform Code  ----  */


vect3	Op_transform (vect3 A, element* p, double wavl, double temp)
/* Using the element p, transform vector A into the output vector
  returned to the user... */
{
vect3	B;	// The output
vect3	C;	// Intermediate
pvec3	P;	// Polar vector intemediate
Etype	type;
Keyword	kw;
datlist	*d;
element	*ex;
stgrism	*gp;
//sdgrism	*dg;
egrism	*xg;
stechel	*ep;
focdat	*fd;
double	v;
double	ga;
double	nx;
double	rls;
double	wm;
double	ordr;
double	rot;

/*  Note that this routine is basically recursive for the Instrument
 element, wherein it will call itself.  For defined elements, we do
 the defined transform based on the element data.  We need to use the
 header for the data types used in the data reading code.  */

/* Set up the output, and return any null case */
    B = A;
    if (p == NULL) return B;

/* Debug print of the input, if desired */
    if (dodebug > 1) {
	printf ("  Xform %s (%s) %6.3e %6.3e\n",
		p->name, Keyname(p->kw), A.x, A.y);
	printf ("  (= %4.3f %4.3f Inches;  %4.3f %4.3f Degrees)\n",
		A.x/Inch, A.y/Inch,  A.x/Degree, A.y/Degree);
    }

/* Establish keyword and type, etc. for a switch based on data type */
    kw = p->kw;
    type = type_of (kw);
    wm = wavl * Angstrom / Micron;		// wavelength microns
    rls = (wm == 0.0) ? 0.0 : 1.0 / (wm*wm);	// recip. lambda squared

/* 
  Find the rotational error value (if any) for any dispersor element
  which has such a value.  The type will be Grism or FloatValue for
  these things.  Scan any p->head datalist queue for an entry with
  a non-NULL data field having the special key "ARot" in the first 4
  characters.  If found, the rotational value is the following 8 bytes.
  Otherwise, the rotational error value is set to zero.
*/
    rot = 0.0;		// The usual value here.
    if (type == Grism || type == Echelle ||
//	type == DGrism ||
	type == EGrism ||
	type == FloatValue) {	// disperser alignment
	C_loop (d, p->head) {
	    if (d->data != NULL) {
		if (!strncmp(d->data, "ARot", 4)) {
		    memcpy (&rot, &(d->data[4]), 8);
		    break;
		}
	    }
	}	// Loop over data values
    }				// disperser alignment

    switch (type) {

// focuser - apply the scale to the vector, then apply the
// distortions if needed.
	case Focuser:
	    fd = (focdat*) p->data;
	    P = polar3 (A);
// If input is direction, convert to angle
	    if (P.z != 0.0) P.r = atan (P.r/P.z);
	    v = P.r * P.r;
	    P.r *= fd->scale;
// oldway:	    if (P.z != 0.0) P.r /= P.z;
// The scale/z is used in camera code old style;
// HOWEVER, we must be careful, as sometimes the z value is not set!

	    P.r *= 1.0 + eval (fd->rmod, v);
	    P.r *= 1.0 + eval (fd->wmod, rls);
	    P.r *= 1.0 + eval (fd->tmod, temp);
// If output is an angle, convert to vector; if not, clear z
// Add sign change in vector conversion as image is before lens
	    if (P.z == 0.0) sincos (-P.r, &(P.r), &(P.z));
	    else  P.z = 0.0;

	    B = rect3 (P);
//	    B.z = znorm (B);
// SHOULD PROBABLY DO SCALE MULTIPLY IF CAMERA AND ZNORM IF COLLIMATOR!!
// oldway:	    B.z = (P.z == 0.0) ? znorm (B) : 0.0;
/* We use z = 0 for all focal (distance) locations, and non-zero for
  angular (normalized) vectors.  NOTE that this will not work correctly
  if the transform is from focus to focus or angle to angle... */
	    break;

// The temperature will probably be stored as a global...

// SEPARATE ALL THE FloatValue TYPE THINGS WITH A SECOND SWITCH.
	case FloatValue:
		v = *( (double*)p->data );
	    switch (kw) {
		case GRAting:
//		if (dodebug>0) printf ("  Lines/mm = %7.3f\n", v);
// grating - apply global angle, order and element's lpmm
//		ga = (v == 0.0) ? 0.0 : *GRangle;
/* Only use of GRangle is to apply a grating positioning anble here */
		if (GRangle != NULL) ga = *GRangle;
		else ga = 0.0;
#ifdef  OTDEBUG
		if (v == 0.0) ga = (Cam_angle - M_PI)/2.0;  // try this...
#endif
		C = rotate (A, ga, 1);
		if (GRorder != NULL) ordr = *GRorder;
		else  ordr = 1.0;
// find, rotate about z axis by rotational error angle...
		C = rotate (C, rot, 3);
//		C.y += wavl * Angstrom * (*GRorder) * v;
		C.y += wavl * Angstrom * ordr * v;
// rotate back about z axis by negative rotational error
		C = rotate (C, -rot, 3);
		C.z = -znorm (C);
		B = rotate (C, -ga, 1);

#ifdef  OTDEBUG
		if (dodebug > 0 && v == 0.0) {
		    B = rotate (B, Cam_angle, 1);
		    B.y *= -1.0;
		    B = rotate (B, -Cam_angle, 1);
		}
// The above compensates for an error in the older test program!
// The error was that instead of doing the reflection and rotation by
// camera angle, the direct code would simply copy the input.  This
// resulted in a non-reflection about y in that code.  So, after
// debug, we take out the correction section above...
#endif

// Will need a separate section if flag & 2 indicates we use y axis
// instead of x axis as the invariant...
		break;

		case ANgle:
// rotation/angle -- transform vector by rotation
// SET AXIS BASED ON FLAG MASK 2 VALUE
// NOTE Angle value in the optics data structure is Degrees!
		if (dodebug>1) printf ("  Angle = %7.3f Deg.\n", v);
		B = rotate (A, v * Degree, 1);
		break;
/* For lint and its clones.. */
		default:
		break;
	    }		// end of subsidary keyword switch
	    break;


	case Grism:
// grism - apply the global order and do usual grism xform

// should also set the axis based on flag & 2 value
	    gp = (stgrism*) p->data;
//	    ex = gp->eg;
//	    if (ex == NULL) {
//		gp->eg = ex = find_element (p->next, gp->glass);  // The glass
//	    }
// using p->next avoids problem with faked grism from dgrism case.

	if (gp->lpmm != 0.0) {
//	    nx = (rls == 0.0) ? 1.0 : sqrt(eval((poly*)ex->data, rls)/rls);	// It's index

	    nx = indxref (gp->glass, &(gp->eg), p->next, rls);
// using p->next avoids problem with faked grism from dgrism case.

//        The argument is a = wavl^-2, and the index will be
//        sqrt (eval (poly, a) / a);

//	    if (dodebug>0) printf ("  Index = %6.4f\n", nx);
	    ga = gp->angle * Degree;
#ifdef  OTDEBUG
	if (GRangle != NULL) ga = *GRangle;	// To mimic the old way error!
#endif
// Entrance at face normal to coordinates here...
	    C.x = A.x / nx;
	    C.y = A.y / nx;
	    C.z = znorm (C);
// rotational alignment compensation; could also rotate A above
	    C = rotate (C, rot, 3);
// rotation by grism angle for y-z axis here
	    C = rotate (C, ga, 1);
// refraction coming out at angled surface,
// and angled surface has grism grating on it...
	    B.x = C.x * nx;
	    if (GRorder != NULL) ordr = *GRorder;
/*	    else  ordr = 1.0;
  ... Above is wrong; order must agree in sign with angle ga here */
	    else  ordr = (ga < 0.0) ? -1.0 : 1.0;
//	    B.y = C.y * nx + wavl * Angstrom * (*GRorder) * gp->lpmm;
	    B.y = C.y * nx + wavl * Angstrom * ordr * gp->lpmm;
	    B.z = znorm (B);
// un-rotate by the grism angle following the grating
	    B = rotate (B, -ga, 1);
// rotational alignment un-compensation
	    B = rotate (B, -rot, 3);
	}	// Only do that stuff if we have non-direct case!

	    break;

	case EGrism:
/* Extended and generalized grism has a prism with an entrance angle
  and glass, followed by a diffraction grating at another angle,
  followed by another prism with glass and exit angle.  May be used
  to model a simple grism (entrance angle zero, final prism null) or
  a grism with two prisms, or with a leading grating, or simply a
  grating itself, or even a prism without a grating.  */
	{	// Local section for some variables...
double	ni, ne;	// Index values for in, exit prisms
double	a1, a2, a3;	// Orientation angles

/* First, find all the index values and angles in radians */
	    xg = (egrism*) p->data;
	    a1 = xg->inangle * Degree;
	    ni = indxref (xg->inglass, &(xg->ig), p, rls);
	    a2 = xg->gangle * Degree;
	    ne = indxref (xg->eglass, &(xg->eg), p, rls);
	    a3 = xg->eangle * Degree;

/* Do an initial alignment rotation about optic axis, may be zero */
	    C = rotate (A, rot, 3);

/* Go to the input prism normal and refract into the prism. */
	    C = rotate (C, a1, 1);
	    B.x = C.x / ni;
	    B.y = C.y / ni;
	    B.z = znorm (B);

/* Go to grating normal coordinates, and do the refraction between
  prisms and the diffraction at the grating */
	    B = rotate (B, a2-a1, 1);
	    C.x = B.x * ni/ne;
	    if (GRorder != NULL) ordr = *GRorder;
	    else  ordr = 1.0;	// No special angle stuff allowed
	    C.y = (B.y * ni + wavl * Angstrom * ordr * xg->lpmm)/ne;
	    C.z = znorm (C);

/* Go to exit prism normal coordinates and refract out of prism */
	    C = rotate (C, a3-a2, 1);
	    B.x = C.x * ne;
	    B.y = C.y * ne;
	    B.z = znorm (B);

/* Finally, remove coordinate rotations for final face and alignment */
	    B = rotate (B, -a3, 1);
	    B = rotate (B, -rot, 3);

	}	// End local section
	    break;

// Add echelle case here:	case Echelle:
	case Echelle:
/* Apply the echelle transforms to obtain reflected output vector */
	    ep = (stechel*) p->data;
/* Find our index and grating constants */
//	    ex = ep->eg;
//	    if (ex == NULL) {
//		ep->eg = ex = find_element (p, ep->glass);  // The glass
//	    }
//	    nx = (rls == 0.0) ? 1.0 : sqrt(eval((poly*)ex->data, rls)/rls);	// It's index
	    nx = indxref (ep->glass, &(ep->eg), p, rls);

	    ordr = (GRorder == NULL) ? 1.0 : *GRorder;
	    v = wavl * Angstrom * ordr * ep->lpmm;

/* Start with adjustment rotation */
	    C = rotate (A, rot, 3);

/* Two prism surfaces going in... */
	    C = snell (C, ep->pris1, nx);
	    C = snell (C, ep->pris2, 1.0/nx);

/* The grating itself */
#ifdef  GTILTMOE
/* ---- Special debug section start ---- */
/* Add a tilt here, dynamically */
//		se->gtnor.x = -999.0;	// special flag value
	if (ep->gtnor.x == -999.0)  {	// local variable section
vect3	tnor, tdis;
double	st, ct;
    /* Find amount of tilt */
	    sincos (ep->gtilt * Degree, &st, &ct);
	    tnor = sum3vect (mul3vect(ep->grnor, ct), mul3vect(ep->grdis, st) );
	    ep->gtnor = tnor;
	    tdis = sub3vect (mul3vect(ep->grdis, ct), mul3vect(ep->grnor, st) );
	    ep->gtdis = tdis;
//	B = diffract (C, tnor, tdis, v);
	}	// end local section
	    B = diffract (C, ep->gtnor, ep->gtdis, v);
/* ---- Special debug section end ---- */
#else
	    B = diffract (C, ep->grnor, ep->grdis, v);
#endif


/* Two prism surfaces coming back out... */
	    B = snell (B, ep->pris2, nx);
	    B = snell (B, ep->pris1, 1.0/nx);

/* End with adjustment un-rotation */
	    B = rotate (B, -rot, 3);
	    break;


// detector -- no transform needed, return A

// instrument  --  start a loop over the element names
//    for each name, find the element, and call to transform
//    the input to a new output, recursively
//   The last output gets returned as the result.
//   this allows instruments to nest naturally
	case Instrument:
/* Debug print of mode */
#ifdef  OTDEBUG
	    if (dodebug > 0) {
		printf (" Mode: ");
		if (llong) printf ("Long, ");
		else printf ("Short, ");
		if (direct) printf ("Direct.");
		else {
		    if (llong) printf ("Grating");
		    else printf ("Grism");
		}
		if (direct) printf ("\n");
		else {
		    printf ("\n Lines/mm: %8.6f  Angle: %8.5f  order: %d",
			cur_den, *GRangle/Degree, *GRorder);
		    if (llong) printf ("\n");
		    else printf (" index: %8.6f\n", grindex);
		}
	    }
#endif
	    C_loop (d, p->head) {	// All instrument members
		if (!strcasecmp(d->data, "GRATING"))
		    ex = Cur_grating;
		else if (!strcasecmp(d->data, "GRISM"))
		    ex = Cur_grism;
		else
		    ex = d->eg;
		if (ex == NULL) ex = d->eg = find_element (p, d->data);
//		    ex = find_element (p, d->data);

#ifdef  OTDEBUG
/* Special debug section.  Report inputs to colimator, grating and
  camera within an instrument.  To match output of older test program. */
		if (dodebug > 0) {
		    switch (ex->kw) {
		    case Colimator:
			printf ("S(in): %8.6f %8.6f\n", B.x/Inch, B.y/Inch);
			break;
		    case GRAting:
		    case GRIsm:
			printf ("A(dg): %8.6f %8.6f\n", B.x/Degree, B.y/Degree);
			break;
		    case FOCuser:
			printf ("B(dg): %8.6f %8.6f\n", B.x/Degree, B.y/Degree);
			break;
		    }
		}
#endif

		B = Op_transform (B, ex, wavl, temp);  // Recursive.
	    }	// Loop over instrument components
	    break;

// other -- maybe an error?  we need to transform all possible
// element types.
	default:
	    break;
/* It is known that if above 2 lines are not present, the enum
  values not handled are: Void Glass Secondary Other   */
    }  // end of keyword switch
// Return the created transform
    return  B;
}


/*  ====  Inverse Transform Code  ====  */


	/* - - - - - - - - - - *\
	|			|
	|    Inverse Support	|
	|			|
	\* - - - - - - - - - - */

static	double	zval (element* p) {
double	z = 1.0;
Etype	type;
Keyword	kw;
datlist	*d;
focdat	*fd;
element	*ex;

/* Check for any NULL case */
    if (p == NULL) return z;

/* Obtain proper kind of input vector, depending on element... */
    kw = p->kw;
    type = type_of (kw);

/* We have z = 0 for all focal (distance) locations and non-zero for
  angular (normalized) vectors.  The transform (forward) code relies
  on its input vector for this information.  It really assumes that
  this is correct, and can fail badly if it is not.  We do not have
  that luxury (of a defined input).  So, we must set our initial input
  z value from an analysis of the element transform.
  If we have a grism, glass, float value, or other we probably have
  an angular input.  If a focuser, we check the scale.  If it is greater
  than 1, it is an image forming element with angular input and a
  focal output.  If less than 1 it is a collimator with focal input
  and angular output.
  If it is an instrument, we check the type of the first element
  recursively.  This is similar to the way the forward transform
  operates.  */

    switch (type) {
	case Focuser:
	    fd = (focdat*) p->data;
/* If scale > 1 we have focuser with angular input (z=1) otherwise
  we have colimator with focal input (z=0)  */
	    z = (fd->scale > 1.0) ? 1.0 : 0.0;
	    break;
	case Instrument:
	    d = p->head;	// Looking at first element...
	    if (!strcasecmp(d->data, "GRATING")) z = 1.0; // ex = Cur_grating;
	    else if (!strcasecmp(d->data, "GRISM")) z = 1.0; // ex = Cur_grism;
	    else {
		ex = d->eg;
		if (ex == NULL) ex = d->eg = find_element (p, d->data);
//		ex = find_element (p, d->data);
		z = zval (ex);	// Recursive
	    }
	    break;
	case FloatValue:	// Angular input
	case Grism:		// Angular input
//	case DGrism:		// Angular input
	case EGrism:		// Angular input
	case Echelle:		// Angular input
	default:		// Angular input
	    z = 1.0;
	    break;
    }  // end of keyword switch
    return  z;
}


static	vect3	zfix (vect2 v, double z) {
/* Special input vector generation.  Obtain a 3 vector with the z
  value properly set.  If z is 0.0, use that as the z of the result;
  otherwise, set it to a normalized value. */
vect3	result;
double	N;
    if (z == 0.0) {
	result = get3de2v (v, z);
    } else {
	result.x = v.x;
	result.y = v.y;
	N = 1.0 - v.x*v.x - v.y*v.y;
	result.z = (N < 0.0) ? 0.0 : sqrt(N);
    }
    return  result;
}


	/* - - - - - - - - - - - - - - - - - - *\
	|					|
	|     The General Inverse Transform	|
	|					|
	\* - - - - - - - - - - - - - - - - - - */




vect3   Inv_transform (vect3 B, element* p, double wavl, double temp) {
/* Using the element p, transform output vector B into the approximate
  input vector A which will produce that output from Opt_transform */
vect3	A;	// Our result
matrix2	pdf;	// Forward partial derivatives
matrix2	trm;	// Inverse transform matrix
vect2	a0, a1;	// Candidate input vectors
vect2	b0, b1;	// output vectors found
vect2	b;	// Desired result
vect2	e;	// Current error vector
vect2	s;	// Best solution movement (a - a0)
double	z;	// Zvalue; 1 for angular input, 0 for dimensional
double	small;	// Convergence criterion
double	d;	// Input difference
int	nx=0;
int	ny=0;	// debug iteration counters...

/* Check for any NULL case */
    A = B;
    if (p == NULL) return A;
    b = get2de3v (B);

/* Obtain proper kind of input vector, depending on element... */
    z = zval (p);
/* If z is 0.0 we have a focal (distance) input desired, and if
  1.0 we have an angular (normalized) input.  */
/* For the z == 0.0 case, our initial dx and dy will be 1.0 (mm), and
  we will set our input vector's z to 0.0; the convergence criteria
  is then 10 microns.
  For the z == 1.0 case, our initial dx and dy will be 0.001 (radian),
  and we will normalize our input vector's z values; the convergence
  will be set to 0.1 arcsecond.  */
/* All that is done by zfix (A, z);  */
    if (z == 0.0) {		// Distance input
	d = 0.1;
	small = 7.07e-3;	// 0.01 mm * sqrt(0.5)
    } else {			// Angular input
	d = 1.0e-4;
	small = 3.428e-7;	// 0.1 Arcsecond * sqrt(0.5)
    }
// Note that d must exceed small at this point to iterate properly.


/* Initialize:  Obtain 2 partial derivative vectors */
    a0 = make2vect (0.0, 0.0);	// Initial (center) a
    b0 = get2de3v (Op_transform (zfix (a0, z), p, wavl, temp));

    a1 = sum2vect (make2vect(d,0.0), a0);	// d in x direction
    b1 = get2de3v (Op_transform (zfix (a1, z), p, wavl, temp));
    pdf.x = mul2vect (sub2vect (b1, b0), 1.0/d);

    a1 = sum2vect (make2vect(0.0,d), a0);	// d in y direction
    b1 = get2de3v (Op_transform (zfix (a1, z), p, wavl, temp));
    pdf.y = mul2vect (sub2vect (b1, b0), 1.0/d);
/* pdf is actually the transpose of the forward transform derivative */

/* Iterative improvement section */
//    while (1) {		// See break below to exit loop.
    for (;;) {		// See break below to exit loop.
	trm = minv2 (mtr2(pdf));	// Inverse transform...
	e = sub2vect (b, b0);	// Error vector
	s = mmul2v (trm, e);	// Solution vector
	if (fabs(d) < small) break;
/* Step in the coordinate with the largest error, and compute the
  partial derivative in that coordinate at the current location. */
	if (fabs(s.x) > fabs(s.y)) {	// Move in x
	    d = s.x;
	    a1 = sum2vect (make2vect(d,0.0), a0);
	    b1 = get2de3v (Op_transform (zfix (a1, z), p, wavl, temp));
	    pdf.x = mul2vect (sub2vect (b1, b0), 1.0/d);
	    nx++;
	} else {			// Move in y
	    d = s.y;
	    if (fabs(d) < small) break;	// For stability.
	    a1 = sum2vect (make2vect(0.0,d), a0);
	    b1 = get2de3v (Op_transform (zfix (a1, z), p, wavl, temp));
	    pdf.y = mul2vect (sub2vect (b1, b0), 1.0/d);
	    ny++;
	}
	b0 = b1;
	a0 = a1;
	if (nx+ny > 50) break;		// Just in case
    }

/* Add another iteration for full accuracy... */
    a1 = sum2vect (s, a0);
    b1 = get2de3v (Op_transform (zfix (a1, z), p, wavl, temp));
    e = sub2vect (b, b1);
    s = mmul2v (trm, e);
    a0 = a1;

/* Final computation */
//    a1 = sum2vect (s, a0);
//    A = zfix (a1, z);
    A = zfix (sum2vect (s, a0), z);
    return  A;
}


/*  ====  Disperser Utilities  ====  */

/* A generalized version which will also find the center wavelength
  of a given order for an echelle as well as any grism, dgrism
  or egrism currently exists in mgutils for the echelle.  */

	/* - - - - - - - - - - - - - - *\
	|				|
	|    Grism Center Wavelength	|
	|				|
	\* - - - - - - - - - - - - - - */


/* The following function finds center wavelength of a grism, given
  the grism element.  Uses current order (that's implicit to the
  general transform) and current temperature.  Consider putting this
  into the optutils subroutines.  */
/* Uses GRorder as a global pointer to double; defaults to 1 if none */
/* Consider adding order as argument, sub. it into GRorder and restore */

double	grwavc (element* te, double cur_temp)
/* Search wavelength space for center line input producing
  centerline output.  Return 0.0 if not found... */
/* Argument is a grism element */
{
vect3	inv={0.0, 0.0, 1.0};
double  ah, al;
double  yh, yl;
double  a, y;
double  b;
int     k;

/* Check for current element (should be grism) */
    if (te == NULL) {
	printf (" **No current grism.\n");
	return  0.0;
    }

/* Check for proper element type */
//    if (te->kw != GRIsm && te->kw != DGRIsm && te->kw != EGRIsm) {
    if (te->kw != GRIsm && te->kw != EGRIsm) {
	printf (" **Element %s is not a grism.\n", te->name);
	return 0.0;
    }

/* Transforming from on-axis vector inv through the grism te, yields
  a y-value:   Op_transform (inv, te, cur_wavl, cur_temp).y
  which should be zero when the center wavelendth is current.
  We vary cur_wavl here to find the zero y value.  */

 
/* Roll through solution space by 500 Angstrom increment, searching for
  the y value to change sign; then do binary search about the interval */
    al = 300.0;		/* Start truly low */
    yl = Op_transform (inv, te, al, cur_temp).y;
    if (yl == 0.0) return al;
    for (;;) {
        ah = al + 500.0;
	yh = Op_transform (inv, te, ah, cur_temp).y;

/* Debug report solution space */
/*      printf ("  ..Wavelength %f y = %f\n", ah, yh); */  /* debug */
        if (yh == 0.0) return  ah;
        if ( (yl < 0.0) && (yh > 0.0) ) break;
        if ( (yl > 0.0) && (yh < 0.0) ) break;
        if (ah > 15000.0) {
            printf ("  **No grism solution found for %s.\n", te->name);
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

	y = Op_transform (inv, te, a, cur_temp).y;
/* Debug trace of solution */
/*      printf ("  ...Wavelength %f y = %f\n", a, y); */ /* debug */
        if (y == 0.0) return  a;
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
    return  a;

}



