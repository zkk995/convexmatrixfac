/* DANTZGMP QP optimizer */

#include <stdlib.h>
#include<math.h>
#define EPS 0.f
#include<stdio.h>

//#define DEBUG
#include"dwqp.h"

int dantzg(real_T *a, int *ndim, int *n, int *nuc, real_T *bv, integer *ib, integer *il, integer *maxiter) {
	/* System generated locals */
	integer a_dim1, a_offset, i__1;

	/* Local variables */
	integer ichk, iter;
	real_T rmin, test;
	integer iout, i, ichki, ic, ir, nt, istand, irtest;
	extern int trsimp_(real_T *, int *, integer *, int *, real_T *, integer *, integer *);
	integer iad;
	real_T val, rat;
	int iret=-1;


	/* 	******************************************* */

	/*    VERSION MODIFIED 1/88 BY NL RICKER */
	/* 	  Modified 12/98 by NL Ricker for use as MATLAB MEX file */
	/* 	  Modified 03/01 by A Bemporad to introduce MAXITER */

	/* 	******************************************* */

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



	/*       CHECK FEASIBILITY OF THE INITIAL BASIS. */

	/* Parameter adjustments */
	--il;
	--ib;
	--bv;
	a_dim1 = *ndim;
	a_offset = a_dim1 + 1;
	a -= a_offset;

	/* Function Body */
	iter = 1;
	nt = *n + *nuc;
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
		if (ib[i] < 0 || bv[ib[i]] >= 0.f) {
			goto L50;
		}
		iret = -i;
		goto L900;
		L50:
		;
	}
	istand = 0;
	L100:

	/*       SEE IF WE ARE AT THE SOLUTION. */

	if (istand != 0) {
		goto L120;
	}
	val = 0.f;
	iret = iter;

	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
		if (il[i] < 0) {
			goto L110;
		}

		/*       PICK OUT LARGEST NEGATIVE LAGRANGE MULTIPLIER. */

		test = bv[il[i]];
		if (test >= val) {
			goto L110;
		}
		val = test;
		iad = i;
		ichk = il[i];
		ichki = i + *n;
		L110:
		;
	}

	/*       IF ALL LAGRANGE MULTIPLIERS WERE NON-NEGATIVE, ALL DONE. */
	/*       ELSE, SKIP TO MODIFICATION OF BASIS */

	if (val >= 0.f) {
		iret=iter;
		goto L900;
	}
	ic = -ib[iad];
	goto L130;

	/*       PREVIOUS BASIS WAS NON-STANDARD.  MUST MOVE LAGRANGE */
	/*       MULTIPLIER ISTAND INTO BASIS. */

	L120:
	iad = istand;
	ic = -il[istand - *n];

	/*       CHECK TO SEE WHAT VARIABLE SHOULD BE REMOVED FROM BASIS. */

	L130:
	ir = 0;

	/*       FIND SMALLEST POSITIVE RATIO OF ELIGIBLE BASIS VARIABLE TO */
	/*       POTENTIAL PIVOT ELEMENT.  FIRST TYPE OF ELIGIBLE BASIS VARIABLE
	 */
	/*       ARE THE REGULAR N VARIABLES AND SLACK VARIABLES IN THE BASIS. */

	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
		irtest = ib[i];

		/*       NO GOOD IF THIS VARIABLE ISN'T IN BASIS OR RESULTING PIVOT WOULD */
		/*       BE ZERO. */

		if (irtest < 0 || a[irtest + ic * a_dim1] == 0.f) {
			goto L150;
		}
		rat = bv[irtest] / a[irtest + ic * a_dim1];

		/*          THE FOLLOWING IF STATEMENT WAS MODIFIED 7/88 BY NL RICKER */
		/*          TO CORRECT A BUG IN CASES WHERE RAT=0. */

		if (rat < 0.f || (rat == 0.f && a[irtest + ic * a_dim1] < 0.f)   ) {
			goto L150;
		}
		if (ir == 0) {
			goto L140;
		}
		if (rat > rmin) {
			goto L150;
		}
		L140:
		rmin = rat;
		ir = irtest;
		iout = i;
		L150:
		;
	}

	/*      SECOND POSSIBLITY IS THE LAGRANGE MULTIPLIER OF THE VARIABLE ADDED*/
	/*       TO THE MOST RECENT STANDARD BASIS. */

	if (a[ichk + ic * a_dim1] == 0.f) {
		goto L170;
	}
	rat = bv[ichk] / a[ichk + ic * a_dim1];
	if (rat < 0.f) {
		goto L170;
	}
	if (ir == 0) {
		goto L160;
	}
	if (rat > rmin) {
		goto L170;
	}
	L160:
	ir = ichk;
	iout = ichki;

	L170:
	if (ir != 0) {
		goto L200;
	}
	iret = *n * -3;
	/* printf("** Fatal error in QP solver!\n"); */
	goto L900;

	L200:

	/*       SET INDICES AND POINTERS */

	if (iout > *n) {
		goto L220;
	}
	ib[iout] = -ic;
	goto L230;
	L220:
	il[iout - *n] = -ic;
	L230:
	if (iad > *n) {
		goto L240;
	}
	ib[iad] = ir;
	goto L250;
	L240:
	il[iad - *n] = ir;
	L250:

	/*       TRANSFORM THE TABLEAU */

	trsimp_(&a[a_offset], ndim, &nt, n, &bv[1], &ir, &ic);
	++iter;
	if (iter > *maxiter) {  /* No solution found within MAXITER iterations */
		iret = iter;
		goto L900;
	}

	/*       WILL NEXT TABLEAU BE STANDARD? */

	istand = 0;
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
		/* L260: */
		if (ib[i] > 0 && il[i] > 0) {
			goto L270;
		}
	}
	goto L280;
	L270:
	istand = iout + *n;
	L280:
	goto L100;

	L900:
	return iret;
} /* dantzg_ */


int dantzg_zhang(real_T *a, int *ndim, int *n, int *nuc, real_T *bv, integer *ib, integer *il, integer *maxiter, integer *c_B, integer *c_N) {
	/* System generated locals */
	integer a_dim1, a_offset, i__1;

	/* Local variables */
	integer iter;
	integer i, ic, ir, nt, b;
	extern int trsimp_(real_T *, int *, integer *, int *, real_T *, integer *, integer *);
	real_T val;
	int iret=-1;

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

	{/* Parameter adjustments */
		--il;--ib;
		--bv;
		a_dim1 = *ndim;
		a_offset = a_dim1 + 1;
		a -= a_offset;
		--c_B;--c_N;
	}


	/* Function Body */

	nt = *n + *nuc;
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
		if (ib[i] < 0 || bv[ib[i]] >= 0.f) {continue;}
		iret = -i;
		return iret;
	}

	for (i = 1; i <= *n; ++i) {
		if(ib[i]>0)c_B[ib[i]] = - i;else c_N[-ib[i]] = - i;
		if(il[i]>0)c_B[il[i]] = i;else c_N[-il[i]] = i;
	}
	iter = 0;
	while (iter++ < *maxiter){

		val = 0.f;
		for (i = 1; i <= *n; ++i) {
			if ( bv[i]< val && a[i+i*a_dim1]< -1e-14) {
				val = bv[i];
				ic = i;
			}
		}/* i */
		if(!(val<0)){break;}
		ir = ic;  /* j */

		/*       TRANSFORM THE TABLEAU */
		trsimp_(&a[a_offset], ndim, &nt, n, &bv[1], &ir, &ic);

		/*       Update c_B and c_N  */
		b = c_B[ic];c_B[ic] = c_N[ic];c_N[ic] = b;

	}/*end of while*/

	for (i = 1; i <= *n; ++i) {
		if(c_B[i]>0) il[c_B[i]] = i;else ib[-c_B[i]] = i;
		if(c_N[i]>0) il[c_N[i]] = -i;else ib[-c_N[i]] = -i;
	}

	return iter;
} /* dantzg_zhang */


int dantzg_(real_T *a, int *ndim, int *n, int *nuc, real_T *bv, integer *ib, integer *il, integer *maxiter) {
	/* System generated locals */
	integer a_dim1, a_offset, i__1;

	/* Local variables */
	integer ichk, iter;
	real_T rmin, test;
	integer iout, i, ichki, ic, ir, nt, istand, irtest;
	extern int trsimp_(real_T *, int *, integer *, int *, real_T *, integer *, integer *);
	integer iad;
	real_T val, rat;
	int iret=-1;

	int debug_k, debug_kk;/*//Debug*/

	/* 	******************************************* */

	/*    VERSION MODIFIED 1/88 BY NL RICKER */
	/* 	  Modified 12/98 by NL Ricker for use as MATLAB MEX file */
	/* 	  Modified 03/01 by A Bemporad to introduce MAXITER */

	/* 	******************************************* */

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



	/*       CHECK FEASIBILITY OF THE INITIAL BASIS. */

	/* Parameter adjustments */
	--il;
	--ib;
	--bv;
	a_dim1 = *ndim;
	a_offset = a_dim1 + 1;
	a -= a_offset;

	/* Function Body */
	iter = 1;
	nt = *n + *nuc;
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
		if (ib[i] < 0 || bv[ib[i]] >= 0.f) {continue;}
		iret = -i;
		return iret;
	}
	istand = 0;
	while (1){

		/*       SEE IF WE ARE AT THE SOLUTION. */

		if (istand != 0) {
			/*       PREVIOUS BASIS WAS NON-STANDARD.  MUST MOVE LAGRANGE */
			/*       MULTIPLIER ISTAND INTO BASIS. */
			iad = istand;
			ic = -il[istand - *n];
		}else{
			val = 0.f;
			iret = iter;
			/*       PICK OUT LARGEST NEGATIVE LAGRANGE MULTIPLIER. */
			for (i = 1; i <= *n; ++i) {
				if (il[i] >= 0) {
					test = bv[il[i]];
					if (test < val) {
						val = test;
						iad = i;
						ichk = il[i];
						ichki = i + *n;
					}
				}
			}

			/*       IF ALL LAGRANGE MULTIPLIERS WERE NON-NEGATIVE, ALL DONE. */
			/*       ELSE, SKIP TO MODIFICATION OF BASIS */
			if (val >= 0.f) {return iter;}
			ic = -ib[iad];
		}

#ifdef DEBUG
		printf("\n\n\n\n\n ");
		if(istand != 0)
			printf("ic =%6d  w_%d hat_i = %d \n", ic, (iad<*n?iad:(iad-*n)), ichk );
#endif
		ir = 0;
		/*
		 * //       FIND SMALLEST POSITIVE RATIO OF ELIGIBLE BASIS VARIABLE TO
		 * //       POTENTIAL PIVOT ELEMENT.  FIRST TYPE OF ELIGIBLE BASIS VARIABLE
		 * //       ARE THE REGULAR N VARIABLES AND SLACK VARIABLES IN THE BASIS.
		 */
		for (i = 1; i <= *n; ++i) {
			irtest = ib[i];
			/*	//       NO GOOD IF THIS VARIABLE ISN'T IN BASIS OR RESULTING PIVOT WOULD
			 * //       BE ZERO.//
			 */
			if (irtest > 0 && fabs(a[irtest + ic * a_dim1])>=EPS ) {
				rat = bv[irtest] / a[irtest + ic * a_dim1];

				if (rat > EPS && (ir == 0 || rat < rmin)) {
					rmin = rat;
					ir = irtest;
					iout = i;

#ifdef DEBUG
					printf("ratio =%6.2f, a_%d%d =%6.2f \n", rat, irtest, ic, a[irtest + ic * a_dim1]);
#endif
				}

			}
		}
#ifdef DEBUG
		if(istand != 0) printf("nonstand \n");
		for(debug_k=1;debug_k<=*ndim;debug_k++){
			printf("%6d\t", ib[debug_k]);
			for(debug_kk=1;debug_kk<=*ndim;debug_kk++)
				printf("%6.2f\t", a[debug_k + debug_kk * a_dim1]);
			printf("%6.2f\t", bv[debug_k]);
			printf("\n");
		}
		printf("%6d\t", 0);
		for(debug_kk=1;debug_kk<=*ndim;debug_kk++)
			printf("%6d\t", il[debug_kk]);
		printf("\n");
#endif
		/*      SECOND POSSIBLITY IS THE LAGRANGE MULTIPLIER OF THE VARIABLE ADDED*/
		/*       TO THE MOST RECENT STANDARD BASIS. */
		if (fabs(a[ichk + ic * a_dim1]) >=EPS) {
			rat = bv[ichk] / a[ichk + ic * a_dim1];
			if (rat >= EPS && (ir == 0 || rat < rmin)) {
				ir = ichk;
				iout = ichki;
#ifdef DEBUG
				printf(".ratio =%6.2f, a_%d%d =%6.2f \n", rat, ichk, ic, a[ichk + ic * a_dim1]);
#endif
			}
		}
		if (ir == 0) {
			iret = *n * -3;
			printf("** Fatal error in QP solver!\n");/**/
			return iret;
		}




		/*       SET INDICES AND POINTERS */
		if (iout > *n) {	il[iout - *n] = -ic;	}
		else{				ib[iout] = -ic;			}

		if (iad > *n) {		il[iad - *n] = ir;		}
		else{				ib[iad] = ir;			}
		/*       TRANSFORM THE TABLEAU */
		trsimp_(&a[a_offset], ndim, &nt, n, &bv[1], &ir, &ic);
		if (++iter > *maxiter) {return iter;} /* No solution found within MAXITER iterations */
		/*       WILL NEXT TABLEAU BE STANDARD? */
		istand = 0;
		for (i = 1; i <= *n; ++i) {
			if (ib[i] > 0 && il[i] > 0){istand = iout + *n;
#ifdef DEBUG
			printf(".w_%d and z_%d  w_%d \n", i, i, iout);
			for(debug_k=1;debug_k<=*ndim;debug_k++){printf("%6d\t", ib[debug_k]);}
			printf(" z \n");
			for(debug_k=1;debug_k<=*ndim;debug_k++){printf("%6d\t", il[debug_k]);}
			printf(" w\n");
#endif
			break;}
		}

	}/*end of while*/
} /* dantzg_ */

/* Subroutine */
int trsimp_(real_T *a, int *ndim, integer *m, int *n, real_T *bv, integer *ir, integer *ic) {
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2;

	/* Local variables */
	integer i, j;
	real_T ap;


	--bv;
	a_dim1 = *ndim;
	a_offset = a_dim1 + 1;
	a -= a_offset;

	/* Function Body */
	i__1 = *m;
	for (i = 1; i <= i__1; ++i) {
		if (i != *ir) {
			ap = a[i + *ic * a_dim1] / a[*ir + *ic * a_dim1];
			bv[i] -= bv[*ir] * ap;
			i__2 = *n;
			for (j = 1; j <= i__2; ++j) {
				if (j != *ic) {
					a[i + j * a_dim1] -= a[*ir + j * a_dim1] * ap;
				}
			}
		}
	}

	/*       NOW TRANSFORM THE PIVOT ROW AND PIVOT COLUMN. */

	ap = a[*ir + *ic * a_dim1];
	i__1 = *m;
	for (i = 1; i <= i__1; ++i) {
		a[i + *ic * a_dim1] = -a[i + *ic * a_dim1] / ap;
	}

	bv[*ir] /= ap;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		a[*ir + j * a_dim1] /= ap;
	}
	a[*ir + *ic * a_dim1] = 1.f / ap;

	return 0;
} /* trsimp_ */
