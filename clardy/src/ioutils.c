/* ioutils.c -- Input/Output Utility subroutines */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
**                                                                   **
**                   Copyright (C) 2002 by                           **
**     Ken Clardy, Carnegie Observatories, Pasadena California.      **
**                                                                   **
**      This software is proprietary, and may be used or copied      **
**          only with the written consent of the author.             **
**              All copies must retain this notice.                  **
**                                                                   **
** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*  ====  Needed Includes  ====  */
#include	<stdlib.h>
#include	<string.h>
#include	<time.h>
#include	<unistd.h>	// for bits in access
#include	<math.h>	// for sexigisimal formatting; exp, log
#include	<ctype.h>	// for toupper/tolower

#define   IS_IOU_C
#include	"ioutils.h"

/*  ====  Standard Definitions  ====  */

/*  -- Code, in order specified in utils.h file -- */


/*  ====  String Manipulation  ====  */

void  trimend (char* a, char* space)
/* Remove trailing whitespace from a string */
{
char*	p;
    if (a == NULL) return;
    if (*a == '\0') return;
    for (p=a; *p != '\0'; p++) ;
/* Points p at null termination of string a */
    for (p--; p >= a; p--) {
	if (strchr (space, *p) == NULL) break;
	*p = '\0';
    }
}

void  trimbeg (char* a, char* space)
/* Remove leading whitespace from a string */
{
char*	p;
    if (a == NULL) return;
    for (p=a; *p != '\0'; p++) {
	if (strchr (space, *p) == NULL) break;
    }
/* Points p at first non-whitespace character in string a */
    if (p == a) return;	/* null case */
    while (*p != '\0') *a++ = *p++;
    *a = '\0';
}

void  trims (char* a, char* space)
/* Remove leading and trailing whitespace from a string */
{
char	*b, *p;
    if (a == NULL) return;
    for (p=b=a; *p != '\0'; p++) {
	if (strchr (space, *p) == NULL) break;
    }
/* Points p at first non-whitespace character in string a */
    if (p != a) {	/* null case */
	while (*p != '\0') *b++ = *p++;
	*b = '\0';
    } else {
	for (b=a; *b != '\0'; b++) ;
    }
    if (*a == '\0') return;	/* If became empty */
/* From above, b points to current null end. */
    for (b--; b >= a; b--) {
	if (strchr (space, *b) == NULL) break;
	*b = '\0';
    }
}


char  *isotime (char *tbuf)
/* Provide the current time in ISO-8601 format.  Use the given
string if tbuf is not null, otherwise obtain a string dynamically.
Returns the string used.  Up to 20 characters are put into string.
Time format is:  yyyy-mm-ddThh:mm:ss
Requires 19 characters, and a 20 character buffer with null end. */
{

/* utime.h or time.h or sys/timeb.h
  see ~/source/lgsub.c or clock.c 
  utime is file access and modification times
  sys/timeb.h gives milliseconds of day; see man ftime
  time.h is usual time of day and date, see man ctime
  sys/time.h defines timeval structure, seconds + microseconds -
    see man gettimeofday
  Actual current time is from "man -s2 time" -- the time function.
  See man strftime - can convert data and time to string with formatting */
/* For the format see /home/draco/clardy/ISO8601.html, and 
  /home/draco/clardy/text/8601v3.pdf and ~/text/8601.pdf 
 also see  http://www.itl.nist.gov/fipspubs/fip4-1.htm
 or do a search on "iso-8601" for further. */

time_t	tloc;
struct tm *tm;
char	*buf;
/* Get the appropriate string pointer into local buf pointer */
    if (tbuf != NULL)	buf = tbuf;			/* Given */
    else		buf = (char*) malloc (20);	/* Gotten */

/* Obtain current time, put into structure, and format it to string */
    time (&tloc);		/* Gets the current time into tloc */
    tm = localtime (&tloc);		/* Find y-m-d-h-m-s+ in *tm */
    strftime (buf, 20, "%Y-%m-%dT%T", tm);	/* Format time to string */
    return  buf;
}


char*	upcase (char* s) {
/* Change string to uppercase; return pointer to argument */
char*	p;
    if (s == NULL) return s;
    for (p=s; *p != '\0'; p++)
	*p = toupper (*p);
    return  s;
}

char*	lowcase (char* s) {
/* Change string to lowercase; return pointer to argument */
char*	p;
    if (s == NULL) return s;
    for (p=s; *p != '\0'; p++)
	*p = tolower (*p);
    return  s;
}

char*	upcasen (char* s, int n) {
/* Change up to n characters of s to upper case, return pointer to arg. */
char*	p;
    if (s == NULL) return s;
    for (p=s; (*p != '\0' && n > 0); p++, n--)
	*p = toupper (*p);
    return  s;
}

char*	lowcasen (char* s, int n) {
/* Change up to n characters of s to lower case, return pointer to arg. */
char*	p;
    if (s == NULL) return s;
    for (p=s; (*p != '\0' && n > 0); p++, n--)
	*p = tolower (*p);
    return  s;
}




/*  ====  Queue Manipulation  ==== */


char	*dynamstr (char *a)
/* Return a dynamic string copy of given string a.  If a is NULL
  or empty, returns NULL */
{
char	*r;
int	l;
    if (a == NULL) return NULL;
    l = strlen (a);
    if (l < 1) return NULL;
    r = (char *) malloc (l + 1);
    strcpy (r, a);
    return  r;
}

char    *pop_name (namlist **a)
/* Remove **a from the queue, and return its name (value) pointer.
  Adjust *a to point to next node in queue, or NULL if emptied. */
{
char    *r;
namlist *h;
namlist *p, *f;
    if (a == NULL) return NULL;
    h = *a;
    if (h == NULL) return NULL;
    r = h->name;
    p = h->last;
    f = h->next;
    if (p != NULL) p->next = f;
    if (f != NULL) f->last = p;
    *a = (f == h) ? NULL : f;
    free (h);
    return  r;
}

namlist	*get_nam (char *a)
/* Create a namlist entry for character value a.  Returns NULL
  if a is NULL or empty.  Puts pointer a directly into structure,
  be sure to use dynamstr if you have a static buffer input. */
{
namlist	*r;
int	l;
    if (a == NULL) return NULL;
    l = strlen (a);
    if (l < 1) return NULL;
    r = (namlist *) malloc (sizeof(namlist));
    r->last = r->next = r;	/* Pointers never null */
    r->name = a;
    return  r;
}

namlist *push_name (namlist *b, namlist *a)
/* Put the queue or element b into queue a, so that it 
  preceeds a.  Return the new head (a) pointer.  */
{
namlist *p;
namlist *r;

/* Take care of null cases */
    if (b == NULL) return a;
/* Find queue ends such that a...p and b...r are merged, into
  something like  ...p - b...r - a... ; this replaces the links
  ...p - a... and ...r - b... with the p - b and r - a links. */
/* If a has null predicessor, we link ...0 - b...r - a...; if
  b has null predicessor, find the next end and use that. */
/* Link it up */
    if (a == NULL) return b;
    p = a->last;
    r = b->last;
    if (r == NULL) {	/* Scan to end */
	for (r = b; r->next != NULL; r = r->next) {
		if (r->next == b) break;
	}
    }
    if (p != NULL) p->next = b;
    b->last = p;
    if (r != NULL) r->next = a;
    a->last = r;
    return a;
}

namlist*  pushtxt (char* text, namlist* a)
/* Get a dynamic string and push it onto the a queue */
{
namlist *b;
char	*c;
    c = dynamstr (text);
    b = get_nam (c);
    return  push_name (b, a);
}

namlist	*kill_name (namlist *a)
/* Kill the entire queue from a... */
{
namlist	*b;
char	*c;
    b = a;
    while (b != NULL) {
	c = pop_name (&b);
	if (c != NULL) free (c);
    }
    return  b;
}


/*  ====  Prompting  ====  */
void  askme (char* buf, int blen, char* prompt)
/* Write prompt, then read buf from stdin */
{
int	l;
    printf ("%s", prompt);
    *buf = '\0';
    fgets (buf, blen, stdin);	/* Stores newline and nul */
    l = strlen (buf);
/* Kill any control characters at end of reply */
    while (l > 0) {
	if (buf[--l] < ' ') buf[l] = '\0';	/* ASCII specific */
    }
}

void  askmestr (char *buf, int blen, char *message)
/* Use message in prompt, obtain buf from user... */
{
int	l;
    printf (" >>%s ", message);
    *buf = '\0';
    fgets (buf, blen, stdin);	/* Stores newline and nul */
    l = strlen (buf);
/* Kill any control characters at end of reply */
    while (l > 0) {
	if (buf[--l] < ' ') buf[l] = '\0';	/* ASCII specific */
    }
}



/*  ====  I/O Routines  ====  */

FILE*	iopen (char* fs, char* tp)
/* Open file, provide error message if not possible */
/* Try opening without giving message, use fopen (library program) */
/* Value of tp is typically "r" or "w" or "rb" or "wb" */
{
FILE    *handle;
char    line[120];

    handle = fopen (fs, tp);

    if (handle == NULL)  {
        strcpy (line, "Opening \"");
        strcat (line, fs);
        strcat (line, "\"");
        perror (line);
/* The error message now prints on stderr. */
    } else strcpy (LastXopen, fs);
    return  handle;
}

FILE*	iopenx (char *fs, char *tp, namlist *exts)
/* Open file fs with I/O type tp.  If open fails, try adding extensions
  from all elements in queue exts in turn.  If still failed, generate
  the error message.  Return opened handle, or NULL if failed. */
{
FILE    *handle;
char    line[120];
char	*bigbuf;
namlist	*a;
int	l;

/* Try the bare thing first */
    handle = fopen (fs, tp);
    if (handle != NULL) {
	strcpy (LastXopen, fs);    // Save the open name
	return  handle;
    }

/* OK, we try any given extensions... */
    if (exts != NULL) {
	l = strlen (fs);
	bigbuf = (char *) malloc (2*l + 80);  /* Guessing that's good */
/* Now, loop around the exts list, see if adding ext to fs helps. */
	a = exts;
	while (a != NULL) {
	    strcpy (bigbuf, fs);
	    if (a->name != NULL) {	/* Added protection here */
		strcat (bigbuf, a->name);
		handle = fopen (bigbuf, tp);
		if (handle != NULL) break;
	    }
	    a = a->next;
	    if (a == exts) break;
	}	/* Loops over all extensions listed */
/* At this point, if handle is not null, we could print a line
  giving the name of the file found */
	if (handle != NULL) {
	    printf ("  Opened %s mode \"%s\".\n", bigbuf, tp);
/* May want to fprintf that to stderr... */
	    strcpy (LastXopen, bigbuf);    // Save actual open name
	}
/* Drop the dynamic buffer now that we're done */
	free (bigbuf);
    }  /* Checking non-null exts list */

    if (handle == NULL)  {
        strcpy (line, "Opening \"");
        strcat (line, fs);
        strcat (line, "\"");
        perror (line);
/* The error message now prints on stderr. */
    }
    return  handle;
}


/*  Finding of file.  Given a filename and a series of possible locations,
  including literals and environment names, find the proper file name and 
  return it as a dynamic string.  */

// /* The "my_access" and "checkfile" routines are local only.  Not to
//   be put into function prototypes section of ioutils.h  */

// int	my_access (char* path)
// /* Return bits 1 for exist, 2 for readable */
// {
// int	r=0;	// return
// int	k;
//     k = access (path, F_OK);
//     if (k == 0) {	// file exists
// 	r |= 1;
// 	k = access (path, R_OK);
// 	if (k == 0) r |= 2;
// 	return  r;
//     } else  return  0;
// }

// int	checkfile (char* spec, char* name, char* dir)
// /* Return 1 if exists, 2 if readable;
//   Put concatenated name into spec;
//   dir may be environment variable, check both ways.   */
// {
// int	k=0;	// return value
// char	*tev;	// translated environment variable
// 
// /* Try first directory/name */
//     sprintf (spec, "%s/%s", dir, name);
//     k = my_access (spec);
//     if (k & 1) return k;
// // printf (" -- Did not find \"%s\".\n", spec);		//debug
// 
// /* If not found, try directoryname without / */
//     sprintf (spec, "%s%s", dir, name);
//     k = my_access (spec);
//     if (k & 1) return k;
// 
// /* Now, try env translation of dir */
//     tev = getenv (dir);
//     if (tev != NULL) {
// /* Try it with / and again without */
// 	sprintf (spec, "%s/%s", tev, name);
// 	k = my_access (spec);
// 	if (k & 1) return k;
// // printf (" -- Did not find \"%s\".\n", spec);		//debug
// 	sprintf (spec, "%s%s", tev, name);
// 	k = my_access (spec);
// 	if (k & 1) return k;
//     }
// 
// /* If nothing found, we need to return a zero, and blank spec */
//     *spec = AN;
//     return  0;
// }


// char*	findrfile (char* fname, namlist* nl)
// /* Return dynamic string containing the readable filespec; nl is list
//   of locations to search.  NULL returned if not found. */
// /* NOTE - only returns files which exist AND are readable.  If writable
//   file is desired, write another similar program for that! */
// {
// char	spec[256];	// Constructed file spec
// int	k;		// Status of search
// char*	result;
// namlist	*p;
// 
// /* Check all those null cases... */
//     if (fname == NULL) return NULL;
// 
// /* Look at all the elements of nl, return any one found. */
//     for (p=nl; p != NULL; p = p->next) {
// 	if (checkfile (spec, fname, p->name) == 3) goto foundfile;
// 	if (p->next == nl) break;
//     }
// 
// /* Look also in the local directory and parent directory */
//     if (checkfile (spec, fname, "." ) == 3) goto foundfile;
//     if (checkfile (spec, fname, "..") == 3) goto foundfile;
// /* Check also bare file name, might have path already */
//     if (checkfile (spec, fname, ""  ) == 3) goto foundfile;
// 
// /* Fall through = not found; return a null value. */
//     return  NULL;
// 
// /*   foundfile  --  Here we return the dynamic string.  */
// foundfile:
//     result = dynamstr (spec);
//     return  result;
// }


/*  Finding of File, revised and generalized, Mark III */


void	findfile (char** spec, char* mode, char* name, char* dir)
/* Returns result in *spec as a dynamic string. */
/* Looks for file of name in dir, possibly translating these as
  environment variables.  The mode requires existance and if it
  contains r=>readable, w=>writable, e=>executable. */
/* Once spec is not NULL, subsequent calls give no action. */
{
char	buf[256];
char	*tev;
int	m;
int	nc;

/* Look for a current result, return if so. */
    if (*spec != NULL) return;

/* Without a name, we can't do a thing. */
    if (name == NULL) return;

/* We need to look.  Find the proper access argument... */
    m = 0;
    if (strchr (mode, 'r')) m |= R_OK;
    if (strchr (mode, 'w')) m |= W_OK;
    if (strchr (mode, 'e')) m |= X_OK;
    if (m == 0) m = F_OK;
/* When we call access(name, m) we use the "!" negation in the
  test, as access returns 0 for success... */

/* See if the dir shows we have a null case, special... */
    nc = (dir == NULL);
    if (!nc) nc = (strlen(dir) < 1);
    if (nc) {		// Null directory; look locally
	sprintf (buf, "./%s", name);		// ./name
	if (!access(buf,m)) goto foundfile;
	sprintf (buf, "%s", name);		// name
	if (!access(buf,m)) goto foundfile;
	sprintf (buf, "../%s", name);		// ../name
	if (!access(buf,m)) goto foundfile;
	sprintf (buf, "../../%s", name);	// ../../name
	if (!access(buf,m)) goto foundfile;
	sprintf (buf, "/%s", name);		// /name
	if (!access(buf,m)) goto foundfile;
	tev = getenv (name);
	if (tev != NULL) strcpy (buf, tev);	// getenv(name)
	if (!access(buf,m)) goto foundfile;
    } else {		// Good directory string
	sprintf (buf, "%s/%s", dir, name);	// dir/name
	if (!access(buf,m)) goto foundfile;
	sprintf (buf, "%s%s", dir, name);	// dirname
	if (!access(buf,m)) goto foundfile;
	tev = getenv (dir);
	if (tev != NULL) {
	    sprintf (buf, "%s/%s", tev, name);	// tev/name
	    if (!access(buf,m)) goto foundfile;
	    sprintf (buf, "%s%s", tev, name);	// tevname
	    if (!access(buf,m)) goto foundfile;
	}
    }
/* Well, if we didn't get to go to foundfile, we haven't found the
  file, so we just return... */
    return;

	/* - - - - - - *\
	|   foundfile   |
	\* - - - - - - */
foundfile:
/* The buffer contains a good file specification; create a dynamic
  string containing it, and return its pointer. */
    *spec = dynamstr (buf);
    return;
}



/*  ====  Parsing Routines  ====  */

char*   skipover (char* line, char* delim)
/* Return line pointer advanced over characters in the delim set. */
{
char*   c;
    if (line == NULL) return line;
    for (c=line; strchr (delim, *c) && *c != AN; c++) ;
    return  c;
}

/* Non-destructive parsing.  Given a pointer (within) a line, copy the
  delimited token to a (provided) buffer, and return the advanced
  pointer to the "rest" of the line.  Handle nulls and delimited
  strings intelegently. */
char*   parse (char* line, char* token, char* delim)
/* line points to string data; copy the first token delimited by
  characters from the delim string to the token buffer, and return
  a pointer to the remainder of the line.  */
{
char    *c;
int     j;

/* If token is not null, set its first character to nul */
    if (token != NULL) *token = AN;
/* If line is null, return a null */
    if (line == NULL) return line;
/* Advance line pointer while it points to delimiter characters. */
    line = skipover (line, delim);
/* If line points to nul, return a NULL pointer. */
    if (*line == AN) return NULL;
/* Find the next whitespace if any */
    c = strpbrk (line, delim);
/* If no whitespace left, we just copy to token, and return a
pointer to the terminating nul */
    if (c == NULL) {            // no whitespace in line
        j = strlen (line);
        c = line + j;
    } else {            // token and rest present
        j = c - line;   // length of token
        c = skipover (c, delim);
    }
/* Copy the stuff from line to token */
    strncpy (token, line, j);
    token[j] = AN;      // adds a nul at end
    return  c;
}


/* Find file name */
char*	cfilename (char* fspec)
/* Return a local pointer to a copy of a filename, extracted from
  the fspec text. */
/* Mostly used to extract the program name from execution command */
{
static	char	buf[32];
char	*c, *e;
int	n=0;

/* Take care of a few null cases */
//    if (fspec == NULL) return NULL;
    *buf = AN;
    if (fspec == NULL) return buf;

/* File name starts at start or following final / character */
    c = strrchr (fspec, '/');
    if (c == NULL) c = fspec;
    else c++;

/* File name ends at end or at first dot following start */
    e = strchr (c, '.');
    n = (e == NULL) ? strlen(c) : (int)(e-c);

/* Copy up to 32 characters to buf... */
    if (n > 31) n = 31;
    strncpy (buf, c, n);
    buf[n] = AN;	/* Supply the terminating null */
    return  buf;
}



/*  ====  Parameter File Reading  ====  */

// typedef  void*  user_prog (char* , void* );
// above may be necessary to prototype the user programs.

// Alternate below is to use... char* filespec, user_prog* userprog, ...
void*  get_data (char* filespec, void* (*userprog)(char*, void*),
	char* delim,	char* leadcom,
	char* anycom,	char* dxcom,
	char* xdcom,	char* dxdcom )
/* Filespec is the data file to read;
  userprog is a user program which will store the data,
  there are 6 comment character sets based on use:
  1) The delimiter list.  If null, strings 4,5,6 are ignored.  When
the data is read, leading and trailing delimiters are removed.
  2) The characters which, if present in the first position of a line
declare the entire line to be a comment and be ignored.
  3) The characters which introduce a comment if they appear anywhere
in the line, with or without delimiters.
  4) The characters which will start a comment only if they are preceeded
by any character from the delimiter set.
  5) The characters which will start a comment only if they are followed
by any character from the delimiter set.
  6) The characters which will start a comment only if both preceeded
and followed by any character from the delimiter set.
.. */
{
char	line[256];
FILE	*infile;
char	*p;
// int	i;
int	j;
void	*V = NULL;
int	kl = 0;		// debug

/* Open, read the file.  For each valid line, call the user
  program to decode it into its own structure.  Return the
  pointer last given by the user program. */

/* Initial stuff.  See if we have a valid file, and open it.
  Return NULL if anything fails. */
    infile = iopen (filespec, "r");
    if (infile == NULL) return NULL;

//    infile = fopen (filespec, "r");
//// NOTE this is redundant with "iopen" in ioutils
//    if (infile == NULL) {
//	strcpy (line, "Opening \"");
//	strcat (line, filespec);
//	strcat (line, "\"");
//	perror (line);
//	return NULL;
//    }
//    strcpy (LastXopen, filespec);

/* Loop on reading lines from the file. */
    while (fgets (line, 256, infile) ) {

/* Remove commentary -- this is the main thing this program does. */
	if (*line == AN) continue;	// Kill blank lines
	if (leadcom != NULL)	// Look for whole line comments
	    if (strchr(leadcom, (int)*line) ) continue;

	if ((p = strpbrk (line, anycom))) {	// unconditional comment
	    *p = AN;
	}


/* Comments preceeded by delimiter */
	for (p=line; p != NULL && *p != AN; p++ ) {
	    if ((p = strpbrk (p, dxcom))) {
		j = p - line;
		if (j < 1) continue;
		if (strchr(delim, line[j-1])) { *p = AN; break; }
	    } else break;
	}

/* Comments followed by delimiter */
	for (p=line; p != NULL && *p != AN; p++ ) {
	    if ((p = strpbrk (p, xdcom))) {
		j = p - line;
		if (line[j+1] == AN) continue;
		if (strchr(delim, line[j+1])) { *p = AN; break; }
	    } else break;
	}

/* Comments between delimiters */
	for (p=line; p != NULL && *p != AN; p++ ) {
	    if ((p = strpbrk (p, dxdcom))) {
		j = p - line;
		if (j < 1) continue;
		if (line[j+1] == AN) continue;
		if (strchr(delim, line[j-1]) &&
		    strchr(delim, line[j+1])) { *p = AN; break; }
	    } else break;
	}

/* Remove trailing cr, lf, delimiter characters if any. */
	j = strlen (line);
	while (j > 0) {
	    if (strchr ("\r\n", line[j-1]) ||
		strchr (delim,  line[j-1]) ) {
		line[j-1] = AN;
		j--;
	    }  else  break;
	}

/* Remove any leading delimiters */
	p = line;
	while ( *p != AN && strchr (delim, *p) ) p++;

	if (*p == AN) continue;	// Skip whole thing if nothing left

/* Call the user program with the remaining non-comment line */
// pointer p points to de-commented line in buffer...
	V = (*userprog) (p, V);
	kl++;	// debug
    }	// End of main loop
    fclose (infile);	// Close up input

//	printf (" ( Read %d lines.)\n", kl);	//debug

/* Final call to user program with NULL line for its final pointer. */
    V = (*userprog) (NULL, V);

/* Return the final pointer from user program */
    return  V;
}


/*  ====  Sexigisimal decode/encode  ====  */

double	posvalue (char*  field)
/* Recursively decode a position value field into degrees (hours)
  and optionally : minutes ... */
{
double	a;
int	k;
// int	j;
char	*fc;

/* Take any null cases out quickly at start */
    if (field == NULL) return 0.0;
    if (*field == '\0') return 0.0;

/* Note that this algorithm can be fooled by such errors as an embeded
 minus sign, or extra decimal points.  Also, non-numeric values
 get treated as zero.  All these should be obvious errors. */

/* Recursively evaluate any minus sign by returning the negative of
  what we return for all characters following that minus sign. */
    if (*field == '-')	/* Recursive minus sign feature */
	return -posvalue (field+1);

/* Evaluate floating point number before any colon, and add 1/60 of
  any recursively obtained value after a colon. */

    fc = strchr (field, ':');	/* Look for colon */
    if (fc == NULL) {		/* No colon, decode directly */
	k = sscanf (field, "%lf", &a);
	if (k != 1)  a = 0.0;
    } else {			/* Colon, split into two values */
	k = sscanf (field, "%lf:", &a);
	if (k != 1) a = 0.0;
	a += posvalue (fc+1) / 60.0;
    }
    return  a;
}

int	intvalue (char*  field)
/* Decode an integer value field */
{
int	k;
int	d;
    if (field == NULL) return 0;
    k = sscanf (field, "%d", &d);
    if (k != 1) d = 0;
    return  d;
}

double	floatvalue (char*  field)
/* Decode a floating value field */
{
int	k;
double	d;
    if (field == NULL) return 0.0;
    k = sscanf (field, "%lf", &d);
    if (k != 1) d = 0.0;
    return  d;
}


/*  ===  Encoding data  ===  */


char*	decomp (int i, char s, int n, int r)
/* Decompose i using radix r for n places, using s as a separator */
/* Return is a dynamic string */
{
char	*result;
char	*top;
int	m;

    if (n > 1) {	/* Multi-period, recursion needed */
	m = i % r;
	top = decomp (i/r, s, n-1, r);
	result = (char*)malloc(strlen(top)+12);
	sprintf (result, "%s%c%02d", top, s, m);
	free (top);
    } else {		/* Single period, return value */
	result = (char*)malloc(12);
	sprintf (result, "%02d%c", i, AN);
    }
    return  result;
}


char*	sexig (double val, int np, int nd)
/* Translate val into sexigisimal quantity with np periods and nd
  decimal points in fraction.  Usually np = 2 and nd is small */
/* Result is a static string here */
/* Note:  Since result is a static string, it needs to be fully used
BEFORE this program is called again.  When printing, if this is called
to format a string, only one such call may be done within a single
printf argument string. */
{
static	char	result[32];
char	*ip;
char	*fp;
double	e;
int	i;
int	neg;

/* Find the sign, store it.  Find the fraction as needed and
store that too.  Get the integer parts in periods using subroutine.
Finally, combine sign, integer and fraction.  */

/* Get the sign and absolute value of val */
    if (val < 0.0) {
	neg = 1;
	val = -val;
    } else  neg = 0;

/* Change radix by np periods of 60 */
    for (i=1; i<np; i++) val *= 60.0;

/* Get the fraction string, and integer part. */
    if (nd > 0) {
	if (nd > 12) nd = 12;
	e = 0.5 * exp (-nd * log(10.0));
	val += e;
	i = (int) val;
	val -= (double) i;
	if (val > e) val -= e;
	fp = (char*)malloc(nd+4);
	sprintf (fp, "%*.*f%c", nd+1, nd, val, AN);
    } else {
	fp = (char*)malloc(4);
	*fp = '\0';
	i = (int) (val + 0.5);
    }

/* Get the integer part from recursive subroutine */
    ip = decomp (i, ':', np, 60);

/* Find the decimal point, if any (often, fp will have the format
of 0.xx rather than .xx as desired.) */
    for (i=0; i<nd; i++) if (fp[i] == '.' || fp[i] == '\0') break;
/* Should be i < 3 above??? */

/* Combine integer parts and fraction parts */
    sprintf (result, "-%s%s%c", ip, fp+i, AN);
    free (ip);
    free (fp);
    if (neg) return result;
    return result+1;
}

void	sexigw (char* buf, double val, int np, int nd)
/* Similar to sexig, but avoids the use of a local buffer by
  writing the result into a user supplied character buffer;
  which must be of sufficient length... */
/* Translate val into sexigisimal quantity with np periods and nd
  decimal points in fraction.  Usually np = 3 and nd is small */
/* New feature!  If np is negative, write with a + sign in front
  when the value is positive.  */
/* Note:  Since result is a static string, it needs to be fully used
BEFORE this program is called again.  When printing, if this is called
to format a string, only one such call may be done within a single
printf argument string. */
{
//  static	char	result[32];
char	*ip;
char	*fp;
double	e;
int	i;
int	neg;
int	ps;

/* Find the sign, store it.  Find the fraction as needed and
store that too.  Get the integer parts in periods using subroutine.
Finally, combine sign, integer and fraction.  */

/* Get the sign and absolute value of val */
    if (val < 0.0) {
	neg = 1;
	val = -val;
    } else  neg = 0;
    if (np < 0) {
	ps = 1;
	np = -np;
    } else ps = 0;

/* Change radix by np periods of 60 */
    for (i=1; i<np; i++) val *= 60.0;

/* Get the fraction string, and integer part. */
    if (nd > 0) {
	if (nd > 12) nd = 12;
	e = 0.5 * exp (-nd * log(10.0));
	val += e;
	i = (int) val;
	val -= (double) i;
	if (val > e) val -= e;
	fp = (char*)malloc(nd+4);
	sprintf (fp, "%*.*f%c", nd+1, nd, val, AN);
    } else {
	fp = (char*)malloc(4);
	*fp = '\0';
	i = (int) (val + 0.5);
    }

/* Get the integer part from recursive subroutine */
    ip = decomp (i, ':', np, 60);

/* Find the decimal point, if any (often, fp will have the format
of 0.xx rather than .xx as desired.) */
    for (i=0; i<nd; i++) if (fp[i] == '.' || fp[i] == '\0') break;
/* Should be i < 3 above??? */

/* Combine integer parts and fraction parts */
//    sprintf (result, "-%s%s%c", ip, fp+i, AN);
//    free (ip);
//    free (fp);
//    if (neg) return result;
//    return result+1;

/* Write the result into the user buffer */
    if     (neg) sprintf (buf, "-%s%s", ip, fp+i);
    else if (ps) sprintf (buf, "+%s%s", ip, fp+i);
    else         sprintf (buf,  "%s%s", ip, fp+i);
    free (ip);
    free (fp);

}


/*  ====  Date decode and encode  ====  */

/* These programs encode and decode an integer date for contemporary
  Gregorian calendar dates.  These may be related to an actual J.D.
  with additional time zone and time of day information, but are
  used here only for calendar date integer arithmetic.  */

int     mjdate (int y, int m, int d)
{
int     i, j, c, jy, jm;

    i = m + 9;
    jy = y + (i/12) - 1;
    jm = i % 12;
    c = jy / 100;
    jy -= 100 * c;
    j = (146097*c)/4 + (1461*jy)/4 + (153*jm+2)/5 + d - 678882;
    return  j;
}


void    date4mjd (int mjd, int* year, int* month, int* day)
/* Convert mjd to Gregorian year, month, day */
{
long int        j, y, d, m;
 
//    j = (long)(jd-0.5) - 1721118;     /* De-bias JD */
 
    j = mjd + 678882;   // Convert mjd
/* j is days since 29 Feb. year"0" (1 BC) */
 
    y = (4*j-1)/146097;         /* Century */
    j = 4*j - 1 - 146097*y;     /* Quarter days in century */
    d = j/4;                    /* Days in century */
    j = (4*d+3)/1461;           /* Year */
    d = 4*d + 3 - 1461*j;       /* Quarter days in year */
    d = (d+4)/4;                /* Days */
    m = (5*d-3)/153;            /* Month */
    d = 5*d - 3 - 153*m;        /* 1/5 days in month */
    d = (d+5)/5;                /* Day in month */
    y = 100*y + j;              /* Combine year, century */
 
    *year = y + (m+2)/12;       /* Correct year */
    *month = 1 + ((m+2) % 12);  /* Correct month */
    *day = d;                   /* Day of month */
    return;
 
}


int	weekofyear (int mj)
{
int	yr, mo, dy;
int	nyd, dow;
int	mjs, wn;

/* Find week of year for given mjd */
/* mjd is 0 mod 7 on Wednesday */
/* Week 1 contains the first Thursday of the year */

// Find the mjd of jan 1 of same year
    date4mjd (mj, &yr, &mo, &dy);	// Get year
    nyd = mjdate (yr, 1, 1);		// mjd of Jan 1
    dow = (nyd+3) % 7;	// Sunday = 0
/* dow is day of week, sunday = 0, for Jan 1. */

// find mjd of start of week 0 of year
    mjs = nyd - dow;
    if (dow < 5) mjs -= 7;	// 5 is Friday
/* mjs is mjd of Sunday starting week 0 */

// find week number
    wn = (mj - mjs) / 7;
// return that result.
    return wn;
}




/*  ====  */

