/*
 * solveqp.h
 *
 *  Created on: May 26, 2011
 *      Author: zkk
 */

#ifndef SOLVEQP_H_
#define SOLVEQP_H_

#include"dwqp.h"
#include"m_dw.h"
#ifdef __cplusplus
extern "C" {
#endif
// min_x 1/2 x^TGx -g^Tx s.t. x>=0 and x^T e=1.
int m_qp_convex(double *G, double*g, double* x, int n,int ng,
		integer * ib, integer *il);

// min_x 1/2 x^TGx -g^Tx s.t. x>=0.
int m_qp_nonnegative(double *G, double*g, double* x, int n, int ng,
		integer * ib, integer *il);

#ifdef __cplusplus
}
#endif

#endif /* SOLVEQP_H_ */
