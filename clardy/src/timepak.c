/* Utility programs for performance timing */

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

#include	<stdlib.h>	// malloc, free
#include	<stdio.h>	// sprintf
#include	"timepak.h"

// Define subroutines called from here...
extern	char*	dynamstr (char* a);

/* ---------------------------------------------------------------------- */

// T_timer* T_new (title)
// void T_clr (T_timer a)
// beg end rep (begin, end, report)
//   twhole = T_new ("Overall");
//   tscan = T_new ("Scanning");
//  T_clr (twhole);
//  T_clr (tscan);
//  T_beg (tscan);
//  ...
//  T_end (tscan);


T_timer*  T_new (char* title) {
/* Create a timer block, clear it and return pointer */
T_timer	*b;
    b = (T_timer*) malloc (sizeof (T_timer) );
    b->title = dynamstr (title);
    T_clr (b);
    return  b;
}

// need dynamic string call

T_timer*  T_kil (T_timer* b) {
/* Remove title, return storage, return null */
    if (b == NULL) return b;
    if (b->title != NULL) free (b->title);
    b->title = NULL;
    T_clr (b);
    free (b);
    return NULL;
}

void	T_clr (T_timer* b) {
/* Clear count and status fields in block */
    if (b == NULL) return;
    b->et    = 0.0;
    b->start = 0;
    b->n     = 0;
    b->state = 0;
    b->flag  = 0;
    b->er    = 0;
}

void	T_beg (T_timer* b) {
/* Start timing */

// Check for null case
    if (b == NULL) return;

// See if already running (an error)
    if (b->state > 0) {
	(b->er)++;
	if (b->flag & T_ReStart) {
	    T_end (b);
	} else {
	    b->state = 0;
	}
    }

// Start the timer and set state
    b->start = clock();
    (b->state)++;

}

void	T_end (T_timer* b) {
/* Stop timing */
clock_t	e;

    if (b == NULL) return;

// See if running (error if not)
    if (b->state <= 0) {
	(b->er)++;
	b->state = 0;
	return;
    }

// Find the elapsed time and total it
    e = clock() - b->start;
    b->et += (double)e / (double)CLOCKS_PER_SEC;
    (b->n)++;
    (b->state)--;
}

char*	T_rep (T_timer* b, int r) {
/* Return dynamic buffer with formatted report */
char	*buf;
char	xbuf[32];
char	*tx;

    if (b == NULL) return  NULL;

    buf = malloc (128);
    buf[0] = '\0';
    if (r > 0) {
	tx = (b->title == NULL) ? " " : b->title;
	sprintf (buf, "%22s %.3f Sec. (%4d @ %6.3f ms)",
	    tx,  b->et,
	    b->n,  1.e3 * b->et / ( (b->n > 0) ? b->n : 1)  );
	if (b->er > 0) {
	    sprintf (xbuf, " (%2d errors)", b->er);
	    strcat (buf, xbuf);
	}
    }
    return  buf;
}


/* ---------------------------------------------------------------------- */

/*   ----    The Mark II system    ----   */

// Usage:
//	TX_con (20)
//	TX_beg (12);
//	TX_enm (12, "subject detail");
//	TX_rep (-1);

void	TX_con (int n) {
/* Control.  If n > 0, establish a set of n timers.
  If n == 0, remove the current timers.  */
T_timer	*t;
int	i;

/* Remove any existing timers before making any more */
    if (T_TX != NULL) {
	for (i=0; i<N_TX; i++) {
	    t = &(T_TX[i]);
	    if (t->title != NULL) free (t->title);
	    t->title = NULL;
	}
	free (T_TX);
    }
    T_TX = NULL;
    N_TX = 0;

/* If we have n > 0, establish n timers, all cleared. */
    if (n > 0) {		// Establish timers
	T_TX = (T_timer*) malloc (n * sizeof (T_timer) );
	N_TX = n;
	for (i=0; i<N_TX; i++) {
	    t = &(T_TX[i]);
	    t->title = NULL;
	    T_clr (t);
	}
    }
}

static	int	TX_ndx (int n) {
/* Common null case and index code */
    if (T_TX == NULL) return -1;
    if (N_TX < 1) return -1;
    if (n > N_TX) return -1;
    if (n < 0) return -1;
    return (n-1);
}

void	TX_beg (int n) {
/* Start timing for timer n.  */
T_timer	*t;
int	i;

    if ( (i = TX_ndx(n)) < 0) return;	// Null cases and index
    t = &(T_TX[i]);
    T_beg (t);
}

void	TX_end (int n) {
/* Stop timing for timer n.  */
T_timer	*t;
int	i;

    if ( (i = TX_ndx(n)) < 0) return;	// Null cases and index
    t = &(T_TX[i]);
    T_end (t);
}

void	TX_nam (int n, char* title) {
/* Apply given title to timer n. */
T_timer	*t;
int	i;

    if ( (i = TX_ndx(n)) < 0) return;	// Null cases and index
    if (title == NULL) return;
    t = &(T_TX[i]);
    if (t->title == NULL)   t->title = dynamstr (title);
}

void	TX_enm (int n, char* title) {
/* Both end and name - useful */
    TX_end (n);
    TX_nam (n, title);
}

void	TX_clr (int n) {
/* Clear timer n */
T_timer	*t;
int	i;

    if ( (i = TX_ndx(n)) < 0) return;	// Null cases and index

/* Clear the specified timer */
    t = &(T_TX[i]);
    T_clr (t);
}

void	TX_rep (int n) {
/* Report.  If n > 0, report timer n.
  If n == 0, report all timers.
  If n < 0, report all timers and remove them.  */
T_timer	*t;
char	*c;
int	i;
int	j, k;

/* Bypass any request for null data report */
    if (T_TX == NULL) return;
    if (N_TX < 1) return;

/* Establish report range, j to k */
    if (n > 0) {
	k = n;
	j = k - 1;
    } else {
	j = 0;
	k = N_TX;
    }

/* Report range of nodes */
    for (i=j; i < k; i++) {
	t = &(T_TX[i]);
	if (t->title != NULL || t->n > 0) {
	    c = T_rep (t, 1);
	    printf (" [%2d] %s\n", i+1, c);
	    free (c);
	}
    }

/* If requested, kill the data set */
    if (n < 0) TX_con (0);

}



/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
