#ifdef __cplusplus
extern "C" {
#endif

#ifndef COMBO_HEADER
#define COMBO_HEADER

	/* ======================================================================
	combo.h,    S.Martello, D.Pisinger, P.Toth     feb 1997
	====================================================================== */

	/* This is the header file of the COMBO algorithm described in 
	* S.Martello, D.Pisinger, P.Toth: "Dynamic Programming and Strong
	* Bounds for the 0-1 Knapsack Problem".
	*/

	/* ======================================================================
	definitions
	====================================================================== */

#define MINRUDI      1000    /* parameter M1 from paper: 1000 normally */
#define MINSET       2000    /* parameter M2 from paper: 2000 normally */
#define MINHEUR     10000    /* parameter M3 from paper: 10000 normally */
#define MAXSTATES 1500000

#undef HASCHANCE             /* should strong upper bounds be used? */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
///#include <malloc.h>


	/* ======================================================================
	macros
	====================================================================== */

#define srand(x)     srand48(x)
#define randm(x)    (lrand48() % (long) (x))

#define SYNC            5      /* when to switch to linear scan in bins */
#define SORTSTACK     200      /* stack for saving discarded intervals */
#define MINMED       1000      /* find exact median if larger size */
#define MAXV  (8*sizeof(long)) /* number of bits in a long integer */

#define TRUE  1                /* boolean value */
#define FALSE 0

#define LEFT  1                /* expansion of core in given direction */
#define RIGHT 2

#define PARTITION 1            /* should sort routine partition or sort */
#define SORTALL   2

#define DET(a1, a2, b1, b2)    ((a1) * (prod)(b2) - (a2) * (prod)(b1))
#define SWAP(a, b)   { register item q; q = *(a); *(a) = *(b); *(b) = q; }
#define NO(a,p)                ((int) ((p) - (a)->fitem + 1))
#define DIFF(a,b)              ((int) (((b)+1) - (a)))
#define TIME(t)                ((double) t / 1000)
#define MIN(a,b)               ((a) < (b) ? (a) : (b))
#define MAX(a,b)               ((a) > (b) ? (a) : (b))


	/* ======================================================================
	type declarations
	====================================================================== */

	typedef int           boolean; /* logical variable         */
	typedef int           ntype;   /* number of states/items   */
	typedef long          itype;   /* item profits and weights */
	typedef long          stype;   /* sum of profit or weight  */
	typedef unsigned long btype;   /* binary solution vector   */
	typedef double        prod;    /* product of state, item   */

	typedef int (*funcptr) (const void *, const void *);

	/* item record */
	//typedef struct {
	//  itype   p;              /* profit                  */
	//  itype   w;              /* weight                  */
	//  boolean x;              /* solution variable       */
	//} item;

	/* item record  modified*/
	typedef struct {
		itype   p;              /* profit                  */
		itype   w;              /* weight                  */
		boolean x;              /* solution variable       */
		int posizione;          /* posizione iniziale      */
	} item;

	/* interval record */
	typedef struct {
		item  *f;               /* first item in interval  */
		item  *l;               /* last item in interval   */
	} interval;

	/* state */
	typedef struct {
		stype psum;             /* profit sum of state     */
		stype wsum;             /* weight sum of state     */
		btype vect;             /* corresponding (partial) solution vector */
	} state;

	/* set of partial vectors */
	typedef struct pset {
		ntype size;             /* set size                */
		state *fset;            /* first element in set    */
		state *lset;            /* last element in set     */
		state *set1;            /* first element in array  */
		state *setm;            /* last element in array   */

		btype    vno;           /* current vector number   */
		item     *vitem[MAXV];  /* current last MAXV items */
		item     *ovitem[MAXV]; /* optimal set of items    */
		btype    ovect;         /* optimal solution vector */
	} partset;


	typedef struct { /* all info for solving separated problem */
		item     *fitem;        /* first item in problem */
		item     *litem;        /* last item in problem */
		item     *s;            /* core is [s,t] */
		item     *t;
		item     *b;            /* break item */
		item     *fpart;        /* first element in sorted core */
		item     *lpart;        /* last element in sorted core */
		stype    wfpart;        /* weight sum up to sorted core */
		item     *fsort;
		item     *lsort;
		stype    wfsort;
		stype    c;             /* capacity of problem */
		stype    z;             /* incumbent solution */
		stype    zwsum;         /* weight sum of incumbent solution */
		stype    lb;            /* lower bound */

		/* solutions may be represented in one of two ways: either a complete */
		/* array of items (fullsol = TRUE), or as last changes in dynamic     */
		/* programming enumeration (fullsol = FALSE) See description of partset */

		boolean  fullsol;       /* which representation of solution */
		item     *fsol;         /* first item in opt solution (fullsol=FALSE) */
		item     *lsol;         /* last item in opt solution (fullsol=FALSE) */
		item     *ffull;        /* start of item array (fullsol=TRUE) */
		item     *lfull;        /* end of item array (fullsol=TRUE) */
		partset  d;             /* set of states, including solution */

		stype    dantzig;       /* dantzig upper bound */
		stype    ub;            /* global upper bound */
		stype    psumb;         /* profit sum up to break item */
		stype    wsumb;         /* weight sum up to break item */

		stype    ps, ws, pt, wt;

		interval *intv1, *intv2;
		interval *intv1b, *intv2b;

		boolean relx;
		boolean master;
		int coresize;
	} allinfo;


	/* ======================================================================
	debug variables 
	====================================================================== */

	extern long simpreduced;
	extern long iterates;
	extern long maxstates;
	extern long coresize;
	extern long optsur;
	extern long relaxations;
	extern long relfeasible;
	extern long reltime;
	extern long pitested;
	extern long pireduced;
	extern long dynheur;


	/* ======================================================================
	forward declarations
	====================================================================== */

	extern stype combo(item *f, item *l, stype c, stype lb, stype ub,
		boolean def, boolean relx);
	/* f,l : first, last item                                               */
	/* c   : capacity of knapsack                                           */
	/* lb  : lower bound. Solution vector is only updated if better z found */
	/* ub  : upper bound. When upper bound is reached terminate immediately */
	/* def : should solution vector be defined or just find objective value */
	/* relx: relaxed problem is solved (no more relaxations will be made)   */
	/* returns the objective value of the problem                           */


#endif

#ifdef __cplusplus
}
#endif
