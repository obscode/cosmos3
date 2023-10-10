/* Utility programs for performance timing - Header file */

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

#ifndef INCLUDE_TIMEPAK_H
#define INCLUDE_TIMEPAK_H

/* INCLUDEs ------------------------------------------------------------- */

#include	<time.h>
#include	<string.h>

/* DEFINEs -------------------------------------------------------------- */

#define  T_ReStart  1

/* TYPEDEFs ------------------------------------------------------------- */

struct	TS_timer {
	char	*title;		// Report identifier
	double	et;		// Total seconds found
	clock_t	start;		// Latest start time
	int	n;		// Interval total
	int	state;		// state = start count
	int	flag;		// flag bits
	int	er;		// count of errors
	};
typedef	struct	TS_timer	T_timer;

/* GLOBALs -------------------------------------------------------------- */

int	N_TX = 0;
T_timer	*T_TX = NULL;
 
/* FUNCTION PROTOTYPEs -------------------------------------------------- */

T_timer*  T_new (char* title);
T_timer*  T_kil (T_timer* b);

void	T_clr (T_timer* b);
void	T_beg (T_timer* b);
void	T_end (T_timer* b);
char*	T_rep (T_timer* b, int r);

void	TX_con (int n);
void	TX_beg (int n);
void	TX_end (int n);
void	TX_nam (int n, char* title);
void	TX_enm (int n, char* title);
void	TX_rep (int n);


/* ---------------------------------------------------------------------- */
/* End of timepak.h package */
 
#endif  /* INCLUDE_TIMEPAK_H */
 
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
