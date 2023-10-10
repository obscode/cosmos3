/*  mgutils.h == Mask Generation Utility Definitions */

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

/*  Utility (library) subroutines for mask generation and related
  application programs.  */

#ifndef  INCLUDE_MGUTILS_H
#define  INCLUDE_MGUTILS_H

#include	"optutils.h"	/* general optics structures */
/* optutils.h includes maskdat.h */
/* maskdat.h  includes ioutil.h, kdcutil.h */
/* We need ioutil.h, kdcutil.h, maskdat.h also */

/*  ===  Globals  ===  */
extern int  dodebug;	/* Defined in main program */
extern	FILE	*dbfile;

extern	int	Nfaroff;	/* Far off mask messages */
extern	int	Noutdet;	/* Outside detector messages */
extern	int	Nduppos;	/* Duplicate mask positions */
extern	int	Nconflx;	/* Conflicts reported */


/*  ====  Definitions  ====  */
#define  OBSLIST_ENTRY  257	/* Length of record */
#define  OBJECTLIST_ENTRY  163	/* Length of record */
#define  SMDF_ENTRY  255	/* Length of record */
#define  NEWLINES  "\n\r"

#define  MESCOUNT  6		/* Message count limit */


/*  ====  Some definitions to remove error messages...  ====  */
typedef  int  bool;
typedef  int  boolean;

/*  --NOTE-- We have to fix the above to real things!!! */


/*  ====  Typedefs  ====  */
/* namlist -- structure; members are:  next, last, name (char*) */
/* The above is an older version - see maskdat.h for new... */

/* Object list -- Master structure containing all data from the
  object list file.  */
typedef  struct  objectlist {
	char*	title;		/* Observation title */
// Only title (not needed) and comment and objq headers are active!
	namlist*  comments;	/* Comments from input file */
	objq*	head;		/* Per-object data listed */
	} objectlist;

typedef  struct  sort_tag {
	struct  sort_tag  *next, *last;
	void	*q;
//	objq	*q;
	}  sortq;

/*  -- Re-globalization of instrument index --  */

typedef	enum	Indx {
// Names for Instrument Index values (found by InstIndx)
	I_imacs,	// Any IMACS mask
// Experimentally switching these next two
	I_lds2,		// LDSS-2 masks
	I_cfg,		// Special centerfield mask
	I_lds3,		// LDSS-3 masks
	I_other
	} Indx;
#define  Instruments  4
// The Instruments should be same as I_other count.

//typedef	enum	IXndx {
//// Names for Instrument extended index (found by ExtIndx)
//	IX_ldss,
//	IX_imx2,
//	IX_imx4
//	} IXndx;

/* General storage structure; indx is application dependent index,
 type is a GS_types value, and data is a pointer to the data. */
typedef  struct  gsto {
	struct	gsto	*next, *last;
	int	indx;
	int	type;
	void	*data;
	} gsto;

/* General storage types, indicating the data; gsto is pointer to
  another gsto structure, int to a single int, double to a double,
  and strx to a general linked structure with no dynamic storage.
  -- This enum is extensible, when the supporting code is also
  extended to support all types. -- */
typedef	enum	GS_types {
	GS_gsto,	// Secondary queue
	GS_int,		// Integer value
	GS_double,	// Double value
	GS_strx		// Queue pointer
	}GS_types;

/* General storage for cut data; the enum subscripts the data
  which is stored in the base level gsto structure.  A function
  returns the design type given the index here. */
typedef	enum	CD_indx {
	CD_refr,	// any reference values
	CD_curve,	// curvature in diopters
	CD_beam,	// beam diameter in mm.
	CD_slew,	// slew rate mm/min
	CD_cut,		// cut rate mm/min
	CD_fine,	// fine motion mm/min
	CD_afoc,	// autofocus range, mm
	CD_dx,		// offset in x, mm
	CD_dy,		// offset in y, mm
	CD_dz,		// offset in z, mm
	CD_flange,	// flange radius, mm (IMACS masks)
	CD_aangle,	// adjustment angle, degrees (IMACS)
	CD_aclear,	// alignment clearance, mm (IMACS)
	CD_LongSlit,	// longest unsupported slit, mm
	CD_SlitGap,	// gap inserted, mm
	CD_tool,	// tool clearance, mm
	CD_cte,		// thermal expansion
	CD_tcut,	// Cutting temperature, Celsius
	CD_dmx,		// Maximum cut distance, mm
	CD_av,		// avoidance queue, pointer
	CD_mma		// millimeter mode allowed flag, int
	}CD_indx;


/*  ====  Default Value coordination  ====  */
#define  D_refsiz    5.803
#define  D_minsep    2.50
#define  D_maxradi 309.44
// maxradi is 15 arc minutes with nominal focus
#define  D_maxradl  88.0

/* Rotator offset angles here */
#define CHUECO (-46.150)
#define LDSSROA (-60.3)


/*  ====  Function Prototypes  ====  */

/*  ===  Special %#%#% genlist prototypes  ===  */

void*	G_pop_c (void* x);
void*	G_pop_l (void* x);
void*	G_pop (void* x, int c);

void*	G_head (void* a);
void*	G_tail (void* a);
void	G_circ (void* a);

void*	G_push_l (void* a, void* b);
void*	G_push_cl (void* y, void* x);
void*	G_push_c (void* y, void* x);

int	listcount (void* x);

/*  ===  Interval and Bounding Box support  ===  */

intval	interval (double a, double b);
int	inintv (double a, intval v);
int	ovr_lap (intval a, intval b);
int	bbover (bbox a, bbox b);
bbox	boundbox (intval x, intval y);
bbox	boundbv2 (vect2 a, vect2 b);
intval	intex (intval v, double p);
bbox	boundex (bbox b, vect2 p);
bbox	boundor (bbox a, bbox b);
bbox	bbedge (sp_edge e);
//bbox	bbspec (spect* s);
double	intlen (intval x);

/*  ===  Maskgen specific I/O  ===  */
void	data_file (char** fspec, char* filename, char* ev, int mf);

/*  ===  Subroutine to read optical data  ===  */
element*  get_optic_data (char* filename, int dl);

/*  ===  Globalization  ===  */
Indx	InstIndx (char* name);
char*	InstName (Indx i);


/*  ===  Object Queue manipulation utilities  ===  */
objq    *pop_obj (objq*);
spect	*kill_spect (spect*);
objval	*kill_objval (objval*);
objq	*kill_obj (objq*);
objq	*kill_objq (objq*);
objq	*push_obj (objq*, objq*);

/*  ===  Other dynamic storage  ===  */
void	kill_ptr (void**);
void	clean_obs (Obs*);
Obs*	kill_obs (Obs*);

/*  ===  Decoding data  ===  */
/*  ===  Encoding data  ===  */
// Those routines have been moved to ioutils.c

/*  --  Debug temporary  --  */
void	dbroq (char* title, Obs* ob);

/*  ====  Generalized Storage Support  ====  */
gsto	*new_gsto (void);
gsto	*kill_gsto (gsto* a);
int	get_CDtype (int item);
void*	GS_lookup (gsto* *head, int index, int item, int flag);


/* === Cutting parameters and avoidance support === */
gsto*	get_cut_data (char* filename, int dl);

/* -- Geometry of slit:  Find largest radius to slit edge */
int	Xcircle (objq* q, vect2 p, double r);
int	avoidanc (objq* ob, avoids* av, double tool);
// double	Rmax (objq*  p);

/* -- Grating angle computation */
// all obsolete or local (static) in mgutils

/*  ---  Optical utility (local to mgutils)  --  */
// double	foclen (Obs* ob);

/*  ====  Checking the object list warning flags  ====  */
#ifdef  OBJCHECK
int	object_check (objq* ob);
#endif

/*  =====  Filling in some structures  =====  */
objq*  get_object (char* data, int type, namlist* preq, Obs* ob);
//grating*	get_grating (char* name);
// detector*	get_detector (element* inst, element* gs, int dl);
// Above 3 routines may be local to mgutils and should not appear here.

/*  ===  Object queue manipulation  ===  */
// objq*	fill_object (char*, int, namlist*);
// Above is well and truly obsolete and is removed.

/*  ===  Reading files  ===  */
Obs*	defobs (element*);
void	get_mpos (Obs*);
void	normalize_obs (Obs* a);
void	check_temp (Obs* obs, gsto* Chead);

Obs*	read_obsfile (FILE*, element*);
objq*	read_objectlist (FILE*, Obs*);
objq*	read_objectqueue (Obs* ob);
void	fill_objects (Obs* ob);

/*  ===  Fill in object data from optics data  ===  */

/*  ===  Conflict geometry  ===  */
int	within (double, double, double);
bool	overlap (double, double, double, double);
#ifdef  DIFREF
static	vect3	drefractn (vect3 pos, vect2 zp, double rc);
static	vect3	urefractn (vect3 pos, vect2 zp, double rc);
#endif  // on DIFREF
vect2	maskvect (Obs*, double, double);
void    unmaskvect (Obs*, vect2, double*, double*);
void	pbnormal (passband*);
vect2	det_pos (Obs*, double, vect3, int);
sp_edge	getedge (Obs*, objq*, vect2, double, double, int, passband);
vect2	lslit (slit);
vect2	rslit (slit);
vect2	wslit (slit);
spect*	getspect (Obs*, objq*, int, passband);

/* === Multiple Object Echellette support === */
// static	double	moecw (int n);
// static	int	moeorder (double  wavl);
// static	passband	moepb (int n);
// above 2 may be declared static later?
spect*	moespect (Obs *obs, objq *ob);

// spect*	killspect (spect*);	// see kill_spect

/*  ===  Conflict geometry  ===  */
// int	spect_edge_intersect (vect2, vect2, double, sp_edge, double, int);
/* Add debug int to end of all these intersection programs... */
/* boolean	spect_intersect (vect2, vect2, vect2, spect*); */
// boolean	spect_intersect (vect2, vect2, vect2, spect*, int); /* debug only */
/* boolean	spect_overlap (spect*, spect*); */
// boolean	spect_overlap (spect*, spect*, int);  /* debug only */

/*  ===  Fill in object data from optics data  ===  */
// int	fillobject (objq*, Obs*);
// used only locally, fillobject is now static...
int	object_setup (Obs*);

/*  ====  Generating virtual objects for gaps  ====  */
#ifdef  GAPVO
// void	Gap_count (int, objq*);
// static	objq*	gap_queue (Obs* ob);
void	add_gaps (Obs* ob);
#endif

/*  ====  Sorting the object queue  ====  */

int	sort_decide (objq*, objq*);
int	sort_check (objq*, int(decision)(objq*,objq*)); // debug routine
objq*	sort_merge (objq*, objq*, int(decision)(objq*,objq*));
/* Need to typedef sortq too */
sortq*	push_sq (sortq*, sortq*);
sortq*	pop_sq (sortq*);
sortq*	kill_sq (sortq*);
objq*	sort_obj (objq*, int(decision)(objq*,objq*));


/*  ===  Conflict queue management  === */
cfl*	push_cfl (cfl*, cfl*);
cfl*	pop_cfl (cfl*);
cfl*	spop_cfl (cfl*);
cfl*	find_cfl (objq*);

/*  ===  Add and drop conflicts  ===  */
cfl*	add_conflict (objq*, objq*, int, bool);
cfl*	drop_conflict (cfl*, bool);

/*  ===  Find conflicts  ===  */
int	find_conflict (objq*, objq*, Obs*);
int	seek_conflict (objq*, bool, Obs*);
int	scan_conflicts (Obs*);

/* -- Debug only routine for maskgen -- */
void	check_cfcount (objq* ob, char* title);

/*  ===  Manage conflicts - drop, count, resolve  ===  */
double	hcpriority (objq*, double);  // use only by conflictedest
int	drop_single_confs (Obs*);
//int	count_conobj (objq*);	// no outside use
int	repconf (Obs*);

/*  ===  Object queue status -- set active/inactive  === */
int	set_active (objq*, Obs*);
int	set_inactive (objq*);
objq*	max_conflict (objq** loq, objq** hiq);
int	order_decide (objq* a, objq* b);
int	priority_decide (objq *a, objq *b);
int	deconflict (Obs *obs);

objq*	inact_obj (Obs *obs);

/*  ===  Writing SMDF file from object queue  ===  */
void	write_smdf (Obs* obs);
void	write_obf (Obs* obs);

/*  ===  Object queue utilities  ===  */
objq*	find_object (objq* a, char* name);
int	active_objs (objq*);
void	repobj (objq* ob);

/*  ===  Object queue ordering  ===  */
double	smdist (objq*, objq*);
double	totdist (objq*);
/* static void revert (objq* ob, int n);  */
objq*	reorder (objq*, int, int);

/*  =====  SMDF Reading subroutines  =====  */
/* static namlist*  kill_litlist (namlist*  n); */
Obs*	read_smdf (char* filename, element *gs, int special);

/*  ====  Wing Chip Support programs  ====  */
objq*	wing_holes (double x, double y, double d, int f);
objq*	wing_set (Obs* obs);
void	add_wings (Obs* obs);

/*  ====  Slit Extension subroutines  ====  */
sleg	*head_sleg (sleg* a);
sleg	*tail_sleg (sleg* a);
//sleg	*push_sleg (sleg*, sleg*);
//sleg	*pop_sleg (sleg*);
sleg	*kill_sleg (sleg*);
//void	augment_side (Obs* obs, objq* ob, int e, double a);
//objq	*conflict_x (Obs* ob, objq* b, int e, double a, objq* f);
// Above routines declared static
sleg	*scan_free (Obs* ob);
sleg	*scan_sle (sleg* a, Obs* obs);
sleg	*sort_sle (sleg* a);
sleg	*fix_sle (sleg* a, Obs* obs);
void	extend_slits (Obs* ob);

/*  ====  */
#endif     /*  INCLUDE_MGUTILS_H  */
/*  ====  */

