#ifndef COMBO_HEADER
#define COMBO_HEADER

#include <stdbool.h>
#include <stddef.h>

/* ======================================================================
combo.h,    S.Martello, D.Pisinger, P.Toth     feb 1997
====================================================================== */

/* This is the header file of the COMBO algorithm described in
* S.Martello, D.Pisinger, P.Toth: "Dynamic Programming and Strong
* Bounds for the 0-1 Knapsack Problem".
*/

/* ======================================================================
type declarations
====================================================================== */

typedef int           cmb_ntype;   /* number of states/items   */
typedef long          cmb_itype;   /* item profits and weights */
typedef long          cmb_stype;   /* sum of profit or weight  */
typedef unsigned long cmb_btype;   /* binary solution vector   */
typedef double        cmb_prod;    /* product of state, item   */

/* item record  modified */
typedef struct {
    cmb_itype   p;          /* profit                  */
    cmb_itype   w;          /* weight                  */
    _Bool   x;              /* solution variable       */
    size_t pos;             /* position in array       */
} cmb_item;

/* ======================================================================
forward declarations
====================================================================== */

extern cmb_stype combo(cmb_item *f, cmb_item *l, cmb_stype c, cmb_stype lb,
                       cmb_stype ub, _Bool def, _Bool relx);
/* f,l : first, last item                                               */
/* c   : capacity of knapsack                                           */
/* lb  : lower bound. Solution vector is only updated if better z found */
/* ub  : upper bound. When upper bound is reached terminate immediately */
/* def : should solution vector be defined or just find objective value */
/* relx: relaxed problem is solved (no more relaxations will be made)   */
/* returns the objective value of the problem                           */

#endif