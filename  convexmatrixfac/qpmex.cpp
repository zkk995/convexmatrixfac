#include"solveqp.h"
#include<cstring>
#include<mex.h>
void usageprint(void){
    mexPrintf(" [y,tabA] = qpmex(A,b,type);\n  min .5 x'Ax-bi'x s.t. ...;\n multiple right hand sides, b=[...bi...] \n nonnegative if type=1 else convex.\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs!=3||(nlhs!=2)){
        usageprint();
        mexErrMsgTxt("error: 3 inputs and 2 outputs are required\n");
    }
    plhs[0] = mxDuplicateArray(prhs[1]);
    double *y = mxGetPr(plhs[0]);
    int type = (int)mxGetScalar(prhs[2]);
    
    const int n= mxGetN(prhs[0]);
    const int n_= mxGetN(prhs[1]);
    
    plhs[1] = mxDuplicateArray(prhs[0]);
    
    long int *c_N=new long int [n];
    long int *c_B=new long int [n];
    long int *il =new long int [n];
    long int *ib =new long int [n];
    double *A = mxGetPr(plhs[1]);
    
    double *b = new double [n*n_]; 
    memcpy(b,mxGetPr(prhs[1]),sizeof(double)*n*n_);

    int flag=0;
    if(type==1){
    	flag=m_qp_nonnegative(A,b,y, n,n_,ib,il);
    }
    else{
    	flag=m_qp_convex(A,b,y, n,n_,ib,il);
    }
    if(!flag)mexErrMsgTxt("error: qp algorithm not convergent\n");
    
    
    delete [] b;
    delete [] il;
    delete [] ib;
    delete [] c_B;
    delete [] c_N;
    
    
}


