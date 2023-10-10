/*  Basic Statistics Package utility header */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
**                                                                   **
**                  Copyright (C) 2004 by                            **
**     Ken Clardy,  Carnegie Observatories,  Pasadena California.    **
**                                                                   **
**      This software is proprietary, and may be used or copied      **
**          only with the written consent of the author.             **
**              All copies must retain this notice.                  **
**                                                                   **
** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ---------------------------------------------------------------------- */

#ifndef	INCLUDE_STATPACK_H
#define	INCLUDE_STATPACK_H

/* INCLUDEs ------------------------------------------------------------- */


/* DEFINEs -------------------------------------------------------------- */

#define  SP_MM
/*  SP_MM is defined to support the computation of minimum and maximum
  values in the statistics structure.  For some applications, it may
  be necessary to save time by omitting this feature.  
  Also, the substat and difstat features cannot correctly keep the
  minimum and maximum values. */

/* TYPEDEFs ------------------------------------------------------------- */

struct	statpak {
	double	sum;
	double	ssq;
	double	bias;
#ifdef  SP_MM
	double	min;
	double	max;
#endif  /* - SP_MM */
	double	avg;
	double	std;
	int	n;	// could be unsigned except for difstat
	int	dummy;	// For alignment (needed in some Linux)
	};
typedef	struct statpak	stats;

struct	statres {
	double	avg;
	double	std;
	};
typedef struct statres	statr;

/* GLOBALs -------------------------------------------------------------- */

/* FUNCTION PROTOTYPEs -------------------------------------------------- */

/* The normal statistics package entry points */
stats	clrstat (double b);		/* Clear statpak structure */
void	addstat (stats* a, double v);	/* Enter a value in stats */
void	cmpstat (stats* a);		/* Compute avg. & s.d. */
stats	sumstat (stats a, stats b);	/* Combine 2 stats structures */

statr	asdstat (stats* s);		/* Obtain avg/sd separately */

/* Providing partial support, not recommended for use */
void	substat (stats* a, double v);	/* Partly remove value from stats */
stats	difstat (stats a, stats b);	/* Attempt to difference structs */
stats	difstatA (stats a, stats b);

/* ---------------------------------------------------------------------- */
/* End of statpack.h package */

#endif	/* INCLUDE_STATPACK_H */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
