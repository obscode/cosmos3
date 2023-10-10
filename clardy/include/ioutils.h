/*  ioutils.h -- Input/Output Utility Subroutines */

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

#ifndef  INCLUDE_IOUTILS_H
#define  INCLUDE_IOUTILS_H

#include	<stdio.h>

/*  ====  Definitions  ====  */
#define  AN  '\0'
#define  WHITESPACE " \t"

#ifdef IS_IOU_C
#define EXT_IOU
#else
#define EXT_IOU extern
#endif

/*  ====  Typedefs  ====  */
typedef  struct  namlist {
        struct namlist  *next, *last;
        char    *name;
        } namlist;

/*  ====  Globals  ====  */
EXT_IOU	char	LastXopen[256];

/*  ====  Function Prototypes  ====  */

/*  ====  String Manipulation  ====  */
void  trimend (char* a, char* space);
void  trimbeg (char* a, char* space);
void  trims (char* a, char* space);
char  *isotime (char *tbuf);
char*	upcase (char* s);
char*	lowcase (char* s);
char*	upcasen (char* s, int n);
char*	lowcasen (char* s, int n);

/*  ====  Queue Manipulation  ==== */
char	*dynamstr (char *a);
char    *pop_name (namlist **a);
namlist	*get_nam (char *a);
namlist *push_name (namlist *b, namlist *a);
namlist *pushtxt (char *text, namlist *a);
namlist	*kill_name (namlist *a);

/*  ====  Prompting  ====  */
void  askme (char* buf, int blen, char* prompt);
void  askmestr (char *buf, int blen, char *message);

/*  ====  I/O Routines  ====  */
FILE*	iopen (char* fs, char* tp);
FILE*	iopenx (char *fs, char *tp, namlist *exts);
// char*	findrfile (char* fname, namlist* nl);
void	findfile (char** spec, char* mode, char* name, char* dir);

/*  ====  Parsing Routines  ====  */

char*   skipover (char* line, char* delim);
char*   parse (char* line, char* token, char* delim);

char*	cfilename (char* fspec);


/*  ====  Parameter File Reading  ====  */

void*  get_data (char* filespec, void* (*userprog)(char*, void*),
	char* delim,	char* leadcom,
	char* anycom,	char* dxcom,
	char* xdcom,	char* dxdcom );

/*  ====  Sexigisimal decode/encode  ====  */

double	posvalue (char*  field);
int	intvalue (char*  field);
double	floatvalue (char*  field);
char*	decomp (int i, char s, int n, int r);
char*	sexig (double val, int np, int nd);
void	sexigw (char* buf, double val, int np, int nd);

/*  ====  Date decode and encode  ====  */

int     mjdate (int y, int m, int d);
void    date4mjd (int mjd, int* year, int* month, int* day);
int	weekofyear (int mj);

/*  ====  */
#endif     /*  INCLUDE_IOUTILS_H  */
/*  ====  */
