/*
 * solveqp.c
 *
 *  Created on: May 26, 2011
 *      Author: zkk
 */

const double lambda=1e-18;

#include"dwqp.h"
#include"m_dw.h"

#if defined(_WIN32)  || defined(_WIN64)

#else
#define dpotrf dpotrf_
#define dpotrs dpotrs_
#define dpotri dpotri_
#endif
//#define dpotrf FORTRAN_WRAPPER(dpotrf)
extern void dpotrf(char*uplo, integer *n, double *a, integer *lda, integer *info);
extern void dpotrs(char*uplo, integer *n, integer *nrhs, double *a, integer *lda,
		double *b, integer *ldb, integer *info);
extern void dpotri(char*uplo, integer *n, double *a, integer *lda, integer *info);

#include<string.h>

// min_x 1/2 x^TGx -g^T x s.t. x>=0 and x^T e=1. with multiple right hand sides
int m_qp_convex(double *G, double*g, double* x, int n, int ng,
		integer * ib, integer *il){

	integer n_=n-1,N=n;
	register int i, j,i_,j_;
	double *a=G+N*(N-1);
	double aa=G[N*N-1],bb;
	double *pG,*pg;
	integer nrhs = ng, nrhs1 = 1;
	char uplo = 'L';
	integer info=0;

	double sum_,*q=x;

	for(j=0;j<n_;j++){
		for(pG=G+N*j, i=0;i<n_;i++){
			pG[i]+=aa-a[i]-a[j];
		}

	}
	for(pg=g,i=0;i<nrhs;i++,pg+=N){
		bb=pg[N-1];
		for(j=0;j<n_;j++){
			pg[j]+=aa-a[j]-bb;
		}
	}

	for(j=0;j<n_;j++){
		G[j+N*j]+=lambda;
	} // pd matrix

	dpotrf(&uplo, &n_, G, &N, &info);
	dpotrs(&uplo, &n_, &nrhs, G, &N, g, &N, &info);

	for(i=0;i<N;i++){
		ib[i]=-(i+1);
		il[i]=(i+1);
		q[i]=1;
	}
	dpotrs(&uplo, &n_, &nrhs1, G, &N, q, &N, &info);

	dpotri(&uplo, &n_, G, &N, &info);


	for(j=1;j<n_;j++){
		for(i=0,i_=j*N,j_=j;i<j;i++){
			G[i_++]=G[j_];j_+=N;
		}
	} // symmetric


	for(sum_=0, i=0;i<n_;i++){
		sum_+=q[i];
	}
	G[N*N-1]=-sum_;
	memcpy(G+N*n_, q, n_*sizeof(double));

	for(j=0;j<n_;j++){
		for(pG=G+j*N, i=0;i<n_;i++){
			*pG = -*pG;pG++;
		}
		G[N*j+n_]=q[j];
	}  // now table is ready!

	memcpy(x, g, ng*N*sizeof(double));
	for(q=x,pg=g,j=0;j<ng;j++,q+=N,pg+=N){
		for(sum_=0, i=0;i<n_;i++){
			sum_+=pg[i];
		}
		q[n_]=1-sum_;
	}

//	return 0;
	return m_dw(G,x,n,ng,ib,il);
}


// min_x 1/2 x^TGx -g^Tx s.t. x>=0.
int m_qp_nonnegative(double *G, double*g, double* x, int n, int ng,
		integer * ib, integer *il){


	integer N = n;
	char uplo = 'L';
	integer nrhs = ng;
	integer info;


	// Temporary variables
	register int i, j,i_,j_;
	double sum_=0;


	for(j=0;j<n;j++){
		G[j+N*j]+=lambda;
	} // pd matrix

	dpotrf(&uplo, &N, G, &N, &info);
	dpotrs(&uplo, &N, &nrhs, G, &N, g, &N, &info);
	dpotri(&uplo, &N, G, &N, &info);

	for(j=1;j<n;j++){
		for(i=0,i_=j*N,j_=j;i<j;i++){
			G[i_++]=G[j_];j_+=N;
		}
	}

	memcpy(x, g, ng*n*sizeof(double));

	for(i=0;i<n*n;i++){
		G[i] = -G[i];
	}

	for(i=0;i<n;i++){
		ib[i]=-(i+1);
		il[i]=(i+1);
	}

//	return 1;
	return m_dw(G,x,n,ng,ib,il);
}

