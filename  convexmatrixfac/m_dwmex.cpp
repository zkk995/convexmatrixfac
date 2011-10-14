#include"solveqp.h"
#include<cstring>
#include<mex.h>
void usageprint(void){
    mexPrintf(" x = m_dwmex(tab,b);\n multiple right hand sides, b=[...bi...] \n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs!=2||(nlhs!=1)){
        usageprint();
        mexErrMsgTxt("error: 2 inputs and 1 outputs are required\n");
    }
    plhs[0] = mxDuplicateArray(prhs[1]);
    double *x = mxGetPr(plhs[0]); 
    const int n= mxGetN(prhs[0]);
    const int ng= mxGetN(prhs[1]);   
    long int *il =new long int [n];
    long int *ib =new long int [n];
    double *G = mxGetPr(prhs[0]);
    
    for(int i=0;i<n;i++){
        ib[i]=-(i+1);
        il[i]=(i+1);
    }
   
    int flag=m_dw(G,x,n,ng,ib,il);
    delete [] il;
    delete [] ib;
    if(!flag)mexErrMsgTxt("error: qp algorithm not convergent\n");
    
}
