/*
 * dwqp.h
 *
 *  Created on: May 26, 2011
 *      Author: zkk
 */

#ifndef DWQP_H_
#define DWQP_H_

typedef long int integer;
typedef double  real_T;

int dantzg(real_T *a, int *ndim, int *n, int *nuc, real_T *bv, integer *ib, integer *il, integer *maxiter);
int dantzg_zhang(real_T *a, int *ndim, int *n, int *nuc, real_T *bv, integer *ib, integer *il, integer *maxiter, integer *c_B, integer *c_N);


/* Subroutine */
int trsimp_(real_T *a, int *ndim, integer *m, int *n, real_T *bv, integer *ir, integer *ic);

#endif /* DWQP_H_ */




/*    DANTZIG QUADRATIC PROGRAMMING ALGORITHM. */

/*    N.L. RICKER    6/83 */

/*    ASSUMES THAT THE INPUT VARIABLES REPRESENT A FEASIBLE INITIAL */
/*    BASIS SET. */

/*    N       NUMBER OF CONSTRAINED VARIABLES (INCLUDING SLACK VARIABLES).*/

/*    NUC     NUMBER OF UNCONSTRAINED VARIABLES, IF ANY */

/*    BV      VECTOR OF VALUES OF THE BASIS VARIABLES.  THE LAST NUC */
/*            ELEMENTS WILL ALWAYS BE KEPT IN THE BASIS AND WILL NOT */
/*            BE CHECKED FOR FEASIBILITY. */

/*    IB      INDEX VECTOR, N ELEMENTS CORRESPONDING TO THE N VARIABLES. */
/*            IF IB(I) IS POSITIVE, THE ITH */
/*            VARIABLE IS BASIC AND BV(IB(I)) IS ITS CURRENT VALUE. */
/*            IF IB(I) IS NEGATIVE, THE ITH VARIABLE IS NON-BASIC */
/*            AND -IB(I) IS ITS COLUMN NUMBER IN THE TABLEAU. */

/*    IL      VECTOR DEFINED AS FOR IB BUT FOR THE N LAGRANGE MULTIPLIERS.*/

/*    A       THE TABLEAU -- SEE TRSIMP DESCRIPTION. */

/*    IRET    IF SUCCESSFUL, CONTAINS NUMBER OF ITERATIONS REQUIRED. */
/*            OTHER POSSIBLE VALUES ARE: */
/*            - I     NON-FEASIBLE BV(I) */
/*            -2N     NO WAY TO ADD A VARIABLE TO BASIS */
/*            -3N     NO WAY TO DELETE A VARIABLE FROM BASIS */
/*            NOTE:  THE LAST TWO SHOULD NOT OCCUR AND INDICATE BAD INPUT*/
/*            OR A BUG IN THE PROGRAM. */
