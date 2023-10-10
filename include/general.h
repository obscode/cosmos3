
/* General header information for data structure organization */

/* ------------------------------------------------------------ */
/*      General storage definitions for mask generation         */
/* ------------------------------------------------------------ */

#ifndef  INCLUDE_GENERAL_H
#define  INCLUDE_GENERAL_H



/* namlist -- structure; members are:  last, next, name (char*) */

/*  ====  */
#endif     /*  INCLUDE_GENERAL_H  */
/*  ====  */



/* ------------------------------------------------------------ */
/*   Examples of structures and typedefs from elsewhere         */
/* ------------------------------------------------------------ */
/*                   All below is commentary                    */
/* ------------------------------------------------------------ */
/* 
struct  vector3  { double  x, y, z; };
typedef struct  vector3 vect3;
*/

/* from mgutils.h: */


/* Object queue -- things per object from object data list */
/*
typedef  struct  objq {
        struct objq  *last, *next;
*/
/* Need to add - image coordinates, dispersed coordinates,
maybe associated slit structure... */
/* Also - slit width, length in 2 directions, slit angle;
  a list of conflicting slits.  If reference slit parameters
  are unioned with reference hole size and type. */
//        char*   name;           /* Name, dynamic string */
//        double  ra, dec;        /* Celestial Coordinates */
//        int     type;           /* Object vs. reference */
//        int     priority;       /* For optimization use */
//        int     flag;           /* In-use, etc. */
//        namlist *precom;        /* Comments before record */
//        namlist *postcom;       /* Comments after record */
//        } objq;


/* Object list -- Master structure containing all data from the
  object list file.  */
//typedef  struct  objectlist {
//        char*   title;          /* Observation title */
/* May also need a mask identification field, astronomer name field
and some explanation (comments) */
//        double  ra, dec;        /* Center field location */
/* Need:  dispersion, spectroscope constants such as grating angle.  */
/* Field rotation angle and preference */
/* Also - default slit width and minimum slit length; slit angle;
 scale, orientation, distortion.  Detector extent in the
  configured detector space.  Detector seams or each separate
  detector area.  */
/* Name of a FITS file for image input.  Option flag to
  plot mask output.  Name of SMDF output file. */
//        namlist*  comments;     /* Comments from input file */
//        objq*   head;           /* Per-object data listed */
//        } objectlist;


/* Hole queue -- details of slit mask hole, slit or alignment */
//typedef  struct  holeq {
//        struct holeq  *last, *next;
//        double  x, y;           /* Location in mm. */
//        double  a;              /* Slit angle */
//        double  sl, sr;         /* Slit lengths from x,y locaiton */
//        double  w;              /* Slit width */
/* Use "slit" structure in place of the above */
/* Also, need celestial coordinates of location, sizes in arcseconds,
  object name (names?) obj/ref type, ref hole kind */
//        namlist *precom;        /* Comments before record */
//        namlist *postcom;       /* Comments after record */
//        } holeq;
 
/* SMDF master data structure */
//typedef  struct  slitmask {
//        char*   title;          /* Mask title */
/* Also needs celestial coordinates of center, guide star data,
  observer name, field rotation angle, others? */
//        namlist  *comments;     /* Comments from all sources */
//        holeq*  head;           /* List of holes in mask */
//        } slitmask;


/* Some from ccbxlib.h: */

//typedef struct button_tag {
//  Window   win;
//  GC       gc;
//  int      w,h;
//  char     text[32];
//  int      flag;
//  int      status;
//  int      symbol;
//  u_long   color;      /* bg */
//} Button;

//typedef struct menuentry_tag {
//  char text[64];
//  int  value;
//  int  flag;
//} MenuEntry;

//typedef struct menu_tag {
//  Window     win;
//  GC         gc; 
//  char       title[128];
//  int        x,y,w,h;
//  MenuEntry  entry[32];
//  int        n;  
//  int        max;
//  int        sel;
//  int        flag;
//  u_long     bg; 
//} Menu;

