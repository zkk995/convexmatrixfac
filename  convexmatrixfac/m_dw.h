/*
 * m_dw.h
 *
 *  Created on: May 30, 2011
 *      Author: zkk
 */

#ifndef M_DW_H_
#define M_DW_H_


#include"dwqp.h"

#ifdef __cplusplus
extern "C" {
#endif
// dw algorithm with multiple right hand sides.
int m_dw(double *M, double* q, int n,int nrhs,integer * ib, integer *il);

#ifdef __cplusplus
}
#endif

#endif /* M_DW_H_ */
