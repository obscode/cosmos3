/* ---------------------------------------------------------------------
 *
 * utils.h
 *
 * Project: GMAP  (OCIW, Pasadena, CA)
 *
 * 1999-05-06  Christoph C. Birk, birk@ociw.edu
 *
 * --------------------------------------------------------------------- */

#ifndef INCLUDE_UTILS_H
#define INCLUDE_UTILS_H 

/* --------------------------------------------------------------------- */ 

#include <sys/time.h>                         /* time_t */

/* --------------------------------------------------------------------- */ 

char*    extend            (char*,char,int);
char*    uppercase         (char*);
char*    lowercase         (char*);
char*    cut_spaces        (char*);
 
char*    alpha_str         (char*,double);
char*    delta_str         (char*,double);
char*    time_str          (char*,double);
 
double   d_alpha           (const char*);
double   d_delta           (const char*);

void     precess           (double,double,double,double*,double*,double);
double   get_siderial      (time_t,double);
double   get_epoch         (time_t);
int      get_uts           (time_t);
char*    get_night         (char*,time_t);
double   get_airmass       (double,double,double);
double   get_ha            (time_t,double,double);

char*    dtime             (void);
void     adebug            (const char*,const char*);
void     tdebug            (const char*,const char*);
#if 0
void     append_logfile    (const char*,const char*,const char*);
#endif
time_t   cor_time          (int);
 
char*    inc_filename      (char*);
char*    dec_filename      (char*);
char*    extract_filename  (char*,const char*);
char*    extract_pathname  (char*,const char*);

void     msleep            (int); 

void     q_sort            (int*,int,int);
void     d_sort            (double*,int,int);

char*    genv2             (char*,char*);
char*    genv3             (char*,char*,char*);

int      imin              (int,int);
int      imax              (int,int);
double   dmin              (double,double);
double   dmax              (double,double);
 
char*    itos              (char*,char*,int);
char*    ltos              (char*,char*,long);
char*    ftos              (char*,char*,double);

int      get_list_index    (char*,char**);
int      str2int           (char*,char**);

char*    get_local_datestr (char*,time_t);
char*    get_ut_timestr    (char*,time_t);

void     put_string        (char*,char*,char*);
char*    get_string        (char*,char*,char*,char*);
long     get_long          (char*,char*,long);
double   get_double        (char*,char*,double);

/* --------------------------------------------------------------------- */

#endif  /* INCLUDE_UTILS_H */

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

