/*  mgfeats.h == Mask Generation Features Definitions */

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

/*  Defines the program features being supported.  Coordinates these
  features between all the source files which have conditional
  compilation based on feature set.  */

#ifndef  INCLUDE_MGFEATS_H
#define  INCLUDE_MGFEATS_H

/* This is now included by maskdat.h */

	/* - - - - - - - - - - - - - - *\
	|				|
	|     Application  Features	|
	|				|
	\* - - - - - - - - - - - - - - */


#define  OBJCHECK
//#undef  OBJCHECK
/* Symbol  OBJCHECK  is defined to implement a consistency check on
  all objects in maskgen.  When present the check is done, and a
  command line flag allows processing if it fails.  */

//#define  GSREQUIRED
#undef  GSREQUIRED
/* Symbol GSREQUIRED is defined to require the use of guide star data
  for mask generation step.  When present, conditionals requiring the
  presence of guide star data in the .obs file are enforced for all
  IMACS instruments, but not for LDSS instruments.  */

#define  GSDATA
//#undef  GSDATA
/* Symbol GSDATA is defined to implement the feature of designating and
  passing on the locations of up to 2 guide stars for the observing
  catalog comments written to the .SMF file.  */

#define  GUIDATE
//#undef  GUIDATE
/* Symbol GUIDATE is defined to enable a date field in place of the
  temperature window in the interface gui.  This date enables the
  display of a sun-twilight-night color bar below the graphic for
  rotator reversals.  The date information is stored as MJD in the
  .obs file, with a default value of 0 meaning no date and no
  sun display.  The windows are pull-down menus for date selection.
  When GUIDATE is defined, the maskcut program will use the temperature,
  probably the default value, in the cutting process as a default.  */

#ifdef  GUIDATE
#define  MOONLINE
//#undef  MOONLINE
#else  // on GUIDATE
#undef  MOONLINE
#endif  // on GUIDATE
/* Symbol  MOONLINE  requires GUIDATE be defined, and if defined it
  allows plotting of a moon rise/set graphic and moon % illumination
  value below the sun rise/set graphic line.  */

#define  REFLIM
//#undef  REFLIM
/* Symbol REFLIM is defined to establish the new method of limiting
  the number of reference stars.  It replaces the old Ref. Select
  window in intgui with Ref. Limit window, and requires the default
  value of -a in maskgen.  Old value is still supported on reading,
  but is not in the gui.  */

//#define  GAPVO
#undef  GAPVO
/* Symbol GAPVO is defined to support the virtual object creation and
  use to model detector gaps.  If not defined, no intgui controls are
  created, and maskgen does not create virtual objects.  When defined,
  full support is compiled.  */

#undef  GAPLDSS
//#define  GAPLDSS
/* Symbol GAPLDSS is defined with the LDSS instrument has a known
  gap value in the optical file; without it the menu is disabled. */

//#define  DIFREF
#undef  DIFREF
/* Symbol DIFREF is defined to support the differential refraction
  computation in mask positions.  If not defined, no change is made
  to the mask position.  Supports Hour Angle box in intgui. */

//#define  SLXB
#undef   SLXB
/* Symbol  SLXB is defined to support the slit extension feature,
  and mostly the button in intgui implementing it.  */

//#define  DMOE
#undef  DMOE
// DMOE is enabling flag for MOE support of extra disperser

#define  DWAV
//#undef  DWAV
/* Symbol DWAV is defined to allow entry of detector wavelength
  bounds separately from conflict wavelength bounds.  Required
  to be defined if DMOE is defined. */
#ifdef  DMOE
 #define  DWAV
#endif

#define  GTILTMOE
//#undef  GTILTMOE
/* Symbol GTILTMOE allows a grating tilt to be specified in opticdef
for the echellette specification.  */

//#define  RLSFULL
#undef  RLSFULL
/* Symbol RLSFULL allows the full computation of reordering using
  actual end points of each slit.  If not on, actual start point is
  used for both start and end of the slit. */
/* Actually used only in mgutils.c, should be coordinated between
  that use and cututils.c for reordering.  */

#define  GCENT
//#undef  GCENT
// GCENT is the global centering feature 
// Defined for skywin, and optionally used in intgui

#define  NNA
//#undef  NNA
/* Symbol NNA activates the Nearest Neighbor Algorithm, which pre-sets
  each node's nearest neighbor as next to that node to jump start the
  shortest total distance part of the algorithm, hoping to speed it
  up siginficantly in the case of many nodes.  Used in the reorder
  program in mgutils; used here to define needed items in the objq
  structure in maskdat.h . */


/* --- Following items are normally defined --- */

// #define  LDSS+NEW
//#undef  LDSS+NEW
/* Symbol  LDSS+NEW  is defined to support new LDSS spectrograph, 
  using its specifics instead of the older LDSS ones.   */
/* The LDSS+NEW symbol compiles special code for cutting multi-mask
 sets when the instrument type is the new LDSS.
 -- It requires the symbol LDSSOFFSET be defined in cututil.c */
/* The + is inserted here to remove this symbol without catching
  any global searches for the name without +.  */

#define  EXORDX
//#undef  EXORDX
/* Symbol EXORDX implements the feature to allow an extra order check
  on comparison star spectra, but not on object spectra.  This is a
  variable in the obs structure, and supported by a pull-down menu
  rather than a check box in the interface GUI.  */

#define  AXBUT
//#undef  AXBUT
/* Symbol AXBUT implements a button which will save the current
  data, run mask generation, and return to the window.  */

#define  PRIORDER
/* Symbol PRIORDER is used to drive the selection of otherwise similar
  nodes in conflict resolution activities.  When set, lower ordered
  objects will be selected over higher ordered ones; when unset the
  opposite happens, later objects override the earlier.  */



	/* - - - - - - - - - - - - - - *\
	|				|
	|        Debug  Features	|
	|				|
	\* - - - - - - - - - - - - - - */

/* Use these in a statement like:
  #if (DEBUG_XX > 1)
	printf ("bla bla bla\n");
  #endif
.. Set the level to 0 to drop the messages, and higher integers
get progressively more of the messages.  Then its easy to turn
the debug messages on/off.   */

#define  DEBUG_GUI  0
// For controlled messages in intgui

#define  DEBUG_GEN  0
// For controlled messages in maskgen

#define  DEBUG_CUT  0
// For controlled messages in maskcut

#define  DEBUG_UTL  0
// For controlled messages in general utilities


	/* - - - - - - - - - - - - - - *\
	|				|
	|     Obsolete  Features	|
	|				|
	\* - - - - - - - - - - - - - - */

/* Below are mostly unused or obsolete */

// #define  LDSS
#undef  LDSS
// LDSS is enabling flag for LDSS-2 instrument support


/*  ====  */
#endif     /*  INCLUDE_MGFEATS_H  */
/*  ====  */

/*  ===  Revsion History  ===  

2005/03-23	Started this feature.

... */
