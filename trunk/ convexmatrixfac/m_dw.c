/*
 * m_dw.c
 *
 *  Created on: May 30, 2011
 *      Author: zkk
 */

#include"m_dw.h"

#define NumThreads 4
#include<string.h>

int m_dw(double *M, double* q, int n,int ng,integer * ib, integer *il);

/*******************************/

#ifdef NumThreads

integer* ib;
integer*il;
integer *c_B;
integer *c_N;
double **A;

#include<pthread.h>
typedef struct{
	double *M;
	double *q;
	integer *ib;
	integer *il;
	int n;
	int nrhs;
	int *flag;
	int pid;
} myarg;

#include<malloc.h>

void mmalloc(int r){
	int i;
	il = (integer*)malloc(NumThreads*r*sizeof(integer));
	ib = (integer*)malloc(NumThreads*r*sizeof(integer));
	c_N = (integer*)malloc(NumThreads*r*sizeof(integer));
	c_B = (integer*)malloc(NumThreads*r*sizeof(integer));
	A = (double**)malloc(NumThreads*sizeof(double*));
	for(i=0;i<NumThreads;i++){
		A[i]=(double*)malloc(r*r*sizeof(double));
	}
}
void mfree(void){
	int i;
	for(i=0;i<NumThreads;i++){
		free(A[i]);
	}
	free(A);

	free(il);free(ib);
	free(c_N);free(c_B);
}


void* update_thread(void *argin){

	register int j, j_, i, i_, k, k_;

	myarg arg=*((myarg *)argin);
	double *M =arg.M;
	double *q_,*q =arg.q;
	const integer *cib=arg.ib;
	const integer *cil=arg.il;
	integer *ib_,*il_;
	int ng =arg.nrhs;
	int n=arg.n;
	int thread =arg.pid;
	integer maxiterint=3*n;
	int ndim=n,nuc =0;
	int iret;
	int* flag=arg.flag;*flag=1;

	for(j=thread;j<ng;j+=NumThreads){

		memcpy(A[thread],M,sizeof(double)*n*n);
		memcpy(ib+thread*n,cib,sizeof(integer)*n);
		memcpy(il+thread*n,cil,sizeof(integer)*n);

		ib_=ib+thread*n; il_=il+thread*n;q_=q+j*n;

#ifndef ZHANG
		iret= dantzg(A[thread],&ndim, &ndim, &nuc, q_, ib_, il_, &maxiterint);
#else
		iret = dantzg_zhang(A[thread],&ndim, &ndim, &nuc,q_, ib_, il_, &maxiterint, c_B+n*thread, c_N+n*thread);
#endif
		if(iret==n*-3) {
			*flag=0;
			break;
		}

		for(i=0;i<n;i++){
			if(il_[i]<0) {q_[i]=0;}
		}
	}
}

int m_dw(double*M, double*q, int n, int nrhs, integer*ib_in,integer*il_in){
	myarg thread_args[NumThreads];
	pthread_t pid[NumThreads];
	int i,flag=1;
	int *flags=(int*)malloc(sizeof(int)*NumThreads);
	mmalloc(n);

	for(i=0; i<NumThreads;i++){
		thread_args[i].M=M;
		thread_args[i].q=q;
		thread_args[i].ib=ib_in;
		thread_args[i].il=il_in;
		thread_args[i].n=n;
		thread_args[i].nrhs=nrhs;
		thread_args[i].flag=flags+i;
		thread_args[i].pid=i;

		pthread_create(&(pid[i]), NULL, update_thread, &thread_args[i]);
	}
	for(i=0;i<NumThreads;i++){
		pthread_join(pid[i], NULL);
	}

	for(i=0;i<NumThreads;i++){
		if(!flags[i]) {
			flag=0; break;
		}
	}
	mfree();
	free(flags);
	return flag;
}

#else
#include<malloc.h>

int m_dw(double*M, double*q, int n, int nrhs, integer*ib,integer*il){
	int i,j;
	double*q_=q;
	integer maxiterint=3*n;
	int ndim=n,nuc =0;
	int iret,flag=1;
	double *A= (double*)malloc(sizeof(double)*n*n);
	integer*il_ = (integer*)malloc(n*sizeof(integer));
	integer*ib_ = (integer*)malloc(n*sizeof(integer));
	integer*c_N = (integer*)malloc(n*sizeof(integer));
	integer*c_B = (integer*)malloc(n*sizeof(integer));

	for(j=0;j<nrhs;j++){

		memcpy(A,M,sizeof(double)*n*n);
		memcpy(ib_,ib,sizeof(integer)*n);
		memcpy(il_,il,sizeof(integer)*n);
		q_ = q+j*n;

#ifndef ZHANG
		iret= dantzg(A, &ndim, &ndim, &nuc, q_, ib_, il_, &maxiterint);
#else
		iret = dantzg_zhang(A, &ndim, &ndim, &nuc,q_, ib_, il_, &maxiterint,c_B,c_N);
#endif
		for(i=0;i<n;i++){
			if(il_[i]<0) {q_[i]=0;}
		}
		if(iret==n*-3) {
			flag=0;
		break;
		}
	}

	free(A);free(c_B);free(c_N);free(il_);free(ib_);
	return flag;
}


#endif
