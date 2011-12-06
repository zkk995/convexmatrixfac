#include<cstring>
#include<cmath>
#include<mex.h>
#include"solveqp.h"
#include "spMat.h"
const mxArray *X;
double nrm2=0;

#pragma comment(lib, "D:/Program Files/MATLAB/R2010b/extern/lib/win32/microsoft/libmwlapack.lib")

#define dgemm dgemm_
#define dgesvd dgesvd_
#define dsyev dsyev_ 
/**/
extern "C"{
	extern void dgemm(char*TRANSA,char*TRANSB,long int*M,long int*N,long int*K,double*ALPHA,
		double*A,long int*LDA,double*B,long int*LDB,double*BETA,double*C,long int*LDC);
	// C := alpha*op( A )*op( B ) + beta*C,
	extern void dgesvd( char*JOBU, char*JOBVT, long int*M, long int*N, double*A, long int*LDA, double*S, 
		double*U, long int*LDU, double*VT, long int*LDVT,double*WORK, long int*LWORK, long int*INFO );
	//JOBU = 'N', JOBVT='N';LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)) 
	extern void dsyev( char*JOBZ, char*UPLO, long int*N, double*A, long int*LDA, 
		double*W, double*WORK, long int*LWORK, long int*INFO );
	//JOBZ='N', UPLO='U' // LWORK >= max(1,3*N-1).For optimal efficiency, LWORK >= (NB+2)*N
}

/* call subroutine  b=b+U'*X; b=b+Y'*X';*/
void sub_calb(double* b,const int bn,int n,double*U,int m,bool isU){
	if(mxIsSparse(X)){
		spMat R;
		R.Ir =(long int*)mxGetIr(X);
		R.Jc =(long int*)mxGetJc(X);
		R.sr =mxGetPr(X);
		R.m  =mxGetM(X);
		R.n  =mxGetN(X);
		R.nnz=R.Jc[R.n];

		if (isU){
			register int j,i_,i,k;
			for (j=0;j<R.n;j++){
				for (i_=R.Jc[j];i_<R.Jc[j+1];i_++){
					i=R.Ir[i_];
					for(k=0;k<bn;k++){
						b[k+j*bn]+=U[i+k*m]*R.sr[i_];
					}
				}
			}
		}else{
			register int j,i_,i,k;
			for (j=0;j<R.n;j++){
				for (i_=R.Jc[j];i_<R.Jc[j+1];i_++){
					i=R.Ir[i_];
					for(k=0;k<bn;k++){
						b[k+i*bn]+=U[j+k*m]*R.sr[i_];
					}
				}
			}
		}

	}else{		
		char TRANSA='T',TRANSB='N';
		long int M = bn, K=m, N = n,LDX=mxGetM(X);
		if(!isU){TRANSB='T';}		
		double ALPHA=1,BETA=1;
		dgemm(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,
			U,&K,mxGetPr(X),&LDX,&BETA,b,&M);
	}
}

void calA(double*C, double *A, int m, int r){ //C = A'*A  A: m x r
	char TRANSA='T',TRANSB='N';
	long int M = r, N=r,K=m;
	double ALPHA=1.0,BETA=0.0;

	dgemm(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,
		A,&K,A,&K,&BETA,C,&M);
}

void calb(double*b, double*U,const int m, const int r,const int j,
		  const int bn,double*Y,const int n,double*work,bool isU){// b=U(I{i},:)*X-(U(I{i},:)*U(~I{i},:)')*Y(~I{i},:);
			  char TRANSA='T',TRANSB='N';
			  long int M,N,K;
			  double ALPHA,BETA;

			  memset(b,0,sizeof(double)*bn*n);
			  M = bn;K=m; N = j;ALPHA=1;BETA=0;
			  if (N>0){
				  dgemm(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,
					  U+j*m,&K,U,&K,&BETA,work,&M);
				  M = bn; K=N; N = n;ALPHA=-1;BETA=1;	   
				  dgemm(&TRANSB,&TRANSA,&M,&N,&K,&ALPHA,
					  work,&M,Y,&N,&BETA,b,&M);
			  }

			  M = bn; K=m; N = r-j-bn;ALPHA=1;BETA=0;
			  if(N>0){
				  dgemm(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,
					  U+j*m,&K,U+(j+bn)*m,&K,&BETA,work,&M);
				  M = bn; K=N; N = n;ALPHA=-1;BETA=1;	   
				  dgemm(&TRANSB,&TRANSA,&M,&N,&K,&ALPHA,
					  work,&M,Y+(j+bn)*n,&N,&BETA,b,&M);
			  }
			  // call subroutine  b=b+U(I{i},:)*X
			  sub_calb(b,bn,n,U+j*m,m,isU);

}




double cal_condition_A_lam1(double*A,int n,double& lam1){
	char JOBZ='N', UPLO='U';
	long int N=n, INFO,LWORK=n*20;
	double *WORK=new double [n*20];
	double *W=new double [n];
	dsyev( &JOBZ, &UPLO, &N, A, &N, 
		W, WORK, &LWORK,&INFO );

	double wmin=W[0],wmax=wmin;
	for(int i=n%2;i<n;i+=2){
		if(W[i]<W[i+1]){
			if(W[i+1]>wmax)wmax=W[i+1];
			if(W[i]  <wmin)wmin=W[i];
		}
		else{ 
			if(W[i]  >wmax)wmax=W[i];
			if(W[i+1]<wmin)wmin=W[i+1];
		}
	}
	double condA=wmin/wmax;

	lam1 = wmax;
	delete [] WORK;
	delete [] W;
	return condA;
}



/*   work: bn*(r-bn+1), y,b: bn*max(m,n), A: bn*bn, bn=(r+.5)/blk */
void block_nmf_update(Mat & U, Mat & Y, int r,int blk, double*A,double*A_,double* b, double tau,double*y,double *work){
	double lam1,ss;
	int  bn=r/blk,flag;		
	long int *il =new long int [bn+1];
	long int *ib =new long int [bn+1];

	int j=0;
	for(int i=0;i<blk;i++){
		calA( A,U.sr+U.m*j,U.m,bn);  
		calb(b, U.sr,U.m, r,j,bn,Y.sr,Y.m,work,true);
		memcpy(A_,A,sizeof(double)*bn*bn);
		ss=(cal_condition_A_lam1(A_,bn,lam1)<tau)?tau*lam1:0;

		//{
		//	mexPrintf("lam1=%8.6f\n",lam1);
		//	for (int k=0;k<bn;k++){
		//		mexPrintf("j=%d, %6.4e,b[%d]=%f\n",j,A[k+bn*k],k,b[k]);
		//	}
		//}

		if (ss>0){
			for (int k=0;k<bn;k++){
				A[bn*k+k]+=ss;			
				double*Yk=Y.sr+(j+k)*Y.m;
				for (int l=0;l<Y.m;l++){
					b[k+l*bn]+=Yk[l]*ss;
				}			
			}
		}


		flag=m_qp_nonnegative(A,b,y, bn,Y.m,ib,il);

		//for (int k=0;k<bn;k++){
		//	mexPrintf("j=%d, %6.5e \n",j,y[k]);
		//}

		for (int k=0;k<bn;k++){
			double*Yk=Y.sr+(k+j)*Y.m;
			for (int l=0;l<Y.m;l++){
				Yk[l] = y[k+l*bn];
			}			
		}

		//for (int k=0;k<bn;k++){
		//	mexPrintf("j=%d, %6.5e \n",j,Y.sr[k*Y.m]);
		//}

		/****************update U******************/
		calA( A,Y.sr+Y.m*j,Y.m,bn); 
		//	for (int k=0;k<bn;k++){
		//		mexPrintf("j=%d, %6.5e\n",j,A[k+bn*k]);
		//	}
		calb(b, Y.sr,Y.m, r,j,bn,U.sr,U.m,work,false);
		////calb_updateU(b, U,Y,j,bn,work);
		//for (int k=0;k<bn;k++){
		//	mexPrintf("j=%d, %6.5e\n",j,b[k]);
		//}


		memcpy(A_,A,sizeof(double)*bn*bn);
		ss=(cal_condition_A_lam1(A_,bn,lam1)<tau)?tau*lam1:0;  

		//{
		//	mexPrintf("lam1=%8.6f\n",lam1);
		//	for (int k=0;k<bn;k++){
		//		mexPrintf("j=%d, %6.5e,b[%d]=%f\n",j,A[k+bn*k],k,b[k]);
		//	}
		//}

		if (ss>0){
			for (int k=0;k<bn;k++){
				A[bn*k+k]+=ss;			
				double*Uk=U.sr+(j+k)*U.m;
				for (int l=0;l<U.m;l++){
					b[k+l*bn]+=Uk[l]*ss;
				}			
			}
		}
		flag=m_qp_nonnegative(A,b,y, bn,U.m,ib,il);
		for (int k=0;k<bn;k++){
			double*Yk=U.sr+(k+j)*U.m;
			for (int l=0;l<U.m;l++){
				Yk[l] = y[k+l*bn];
			}			
		}
		/**********************************/		
		j+=bn;
		if((i+1)==(blk-r%blk))bn=bn+1;
	} 	
	delete [] il;
	delete [] ib;
}

double CalculateObj(Mat&U,Mat&Y,double*work){ //obj = (nrm2+sum(sum((U'*U).*(Y'*Y))))/2- sum(sum(Y'.*(U'*X)));
	double obj=nrm2;
	int r= U.n, n=Y.m;
	double*U1 = new double [r*r];
	double*U2 = new double [r*r];
	calA(U1, U.sr, U.m, r);
	calA(U2, Y.sr, Y.m, r);
	for (int i=0;i<r*r;i++){
		obj+=U1[i]*U2[i];
	}
	obj/=2;

	memset(work,0,sizeof(double)*Y.m*r);
	sub_calb(work,r,Y.m,U.sr,U.m,true);
	for (int j=0;j<n;j++)
	for (int i=0;i<r;i++){
		obj-=Y.sr[j+i*n]*work[i+j*r];
	}

	delete [] U1;
	delete [] U2;
	return obj;
}
int max(int a,int b){return a>b?a:b;}
void block_nmf(Mat & U, Mat & Y, const int blk,double tol,double tau,const int maxiter,double*err,int disp){
	int r=U.n; 
	int bn=r/blk;
    double  obj=0,obj_;
	
	double *sr = mxGetPr(X);
	int n_;
	if (mxIsSparse(X)) {n_ = mxGetJc(X)[mxGetN(X)];}else{ n_ =mxGetN(X)*mxGetM(X); }
	nrm2=0;for (int i=0;i<n_;i++){nrm2+=sr[i]*sr[i];}
	
	double*work=new double[bn*(r-bn+1)];
	double*A=new double[(bn+1)*(bn+1)];
	double*A_=new double[(bn+1)*(bn+1)];
	double*b=new double[max((bn+1)*max(U.m,Y.m), r*Y.m)];
	double*y=new double[(bn+1)*max(U.m,Y.m)];

	for (int iter=0;iter<maxiter;iter++){
		block_nmf_update(U, Y, r,blk, A,A_,b, tau,y,work);

		//stop condition
		obj_=obj;
		obj=CalculateObj(U,Y,b);
		err[iter]=obj;
		if(disp==1){ mexPrintf("iter: %d obj:%6.4e \n",iter,obj);}
		if(fabs(obj_-obj)<tol*obj){break;}
	}

	delete[] work;
	delete[] A;
	delete[] A_;
	delete[] b;
	delete[] y;
}
void usageprint(void){
	mexPrintf(" [U,Y,err] = block_nmfmex(X,U,Y,blk,tol,maxiter,tau,display);"
		"\n  min |X - UY'|,  X : m x n "
		"\n call subroutine  obj=CalculateObj(X,U,Y); "
		"\n  subroutine [y,tabA] = qpmex(A,b,1) \n");
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	if (nrhs!=8||(nlhs!=3)){
		usageprint();
		mexErrMsgTxt("error: 8 inputs and 3 outputs are required\n");
	}
	X=prhs[0];
	plhs[0] = mxDuplicateArray(prhs[1]);
	plhs[1] = mxDuplicateArray(prhs[2]);

	Mat U,Y;
	U.m=mxGetM(plhs[0]);
	U.n=mxGetN(plhs[0]); 
	U.sr=mxGetPr(plhs[0]);
	Y.m=mxGetM(plhs[1]);
	Y.n=mxGetN(plhs[1]); 
	Y.sr=mxGetPr(plhs[1]);

	int blk = (int)mxGetScalar(prhs[3]);
	double tol = mxGetScalar(prhs[4]);
	int maxiter = (int)mxGetScalar(prhs[5]);
	double tau = mxGetScalar(prhs[6]);
	int disp =(int)mxGetScalar(prhs[7]);
	
	plhs[2] =mxCreateDoubleMatrix(maxiter,1,mxREAL);;
	double *err=mxGetPr(plhs[2]);
	block_nmf(U, Y, blk,tol,tau,maxiter,err,disp);
}
