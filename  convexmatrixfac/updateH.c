/*
 *  updateH.c
 *
 *
 *  Created by Zhao on 11-10-1.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include<stdlib.h>
#include<string.h>
#include<math.h>

/* min f(H) = .5|V - W'H|_F   (W:r x m;) H:r x n
 * GH = nabla_H f  r x n
 * WW = WW' :r x r
 * WV = WV  :r x n
 */

#define WV(i, j) WV[i+j*r]
#define WW(i, j) WW[i+j*r]
#define GH(i, j) GH[i+j*r]
#define H(i, j)  H[i+j*r]
#define Hnew(i, j) Hnew[i+j*r]
#define SH(i, j) SH[i+j*r]

double tol =0.001;

int goHiter(double *WV, double *WW, double *GH, double *H, double*Hnew, int r, int n, int maxiter) {
    int i, j, p, q;
    double *maxh_v;
    int *maxh_ind;
    int hinner, flop = 0;
    double init = 0, s, ss, diffobj, pv;
    double*SH;
    
    maxh_v=(double*)malloc(n*sizeof(double));
    maxh_ind=(int*)malloc(n*sizeof(int));
    SH= (double*)malloc(r*n*sizeof(double));
    
    memset(Hnew, 0, sizeof(double)*r*n); /*set to zeros*/
    
    for ( j=0 ; j<n ; j++ ) {
        maxh_v[j] = 0;
        maxh_ind[j] = -1;
        for ( i=0 ; i<r ; i++ ) {
            s =fabs(WW(i , i))>1e-14 ? GH(i, j)/WW(i, i):0;/*avoid to divide 0 */
            s = H(i, j)-s;
            if (s < 0)s=0;
            
            s = s-H(i, j);
            SH(i, j) = s;
            diffobj = (-1)*s*GH(i, j)-0.5*WW(i, i)*s*s;
            if (diffobj > maxh_v[i]) {
                maxh_v[j] = diffobj;
                maxh_ind[j] = i;
            }
        }
        if ( maxh_v[j] > init )
            init = maxh_v[j];
    }
    flop += r*n *7;
    
    for (q=0 ; q<n; q++ ) {
        for ( hinner = 0 ; hinner <maxiter; hinner++) {
            pv=maxh_v[q];
            p = maxh_ind[q];
            s = SH(p, q);
            
            if (pv < init*tol)
                break;
            
            Hnew(p, q) += s;
            H(p, q) +=  s;
            maxh_v[q] = 0;
            maxh_ind[q]=-1;
            for ( i=0 ; i<r;  i++ ) {
                GH(i, q) = GH(i, q) + s*WW(i, p);
                ss = fabs(WW(i, i))>1e-14? GH(i, q)/WW(i, i):0; /*add*/
                ss = H(i, q)-ss;
                if (ss < 0)ss=0;
                ss = ss-H(i, q);
                SH(i, q) = ss;
                diffobj = (-1)*ss*GH(i, q)-0.5*WW(i, i)*ss*ss;
                if ( diffobj > maxh_v[q]) {
                    maxh_v[q] = diffobj;
                    maxh_ind[q] = i;
                }
            }
        }
        flop += hinner*3+7*r*hinner;
    }
    free(maxh_v);
    free(maxh_ind);
    free(SH);
    return flop;
}


#include<mex.h>
void usageprint(){
    mexPrintf("usage:\n"
            "[Hnew flopadd] = updateH(WV,WW,GH,H,tol);\n"
            "flop = flop + flopadd;\n"
            "H = H + Hnew;\n"
            "/* min f(H) = .5|V - W'H|_F   (W:r x m;) H:r x n \n"
            " * GH = nabla_H f  r x n\n"
            " * WW = WW' :r x r \n"
            " * WV = WV  :r x n \n"
            " */"
            "\n");
}

/* mex -largeArrayDims updateH.c */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double*WV, *WW, *GH, *H, *Hnew;
    int r, n;
    int maxiter;
    int flop=0;
    mxArray*mxH, *mxGH;
    
    if (nrhs!=5||nlhs!=2){
        usageprint();
        mexErrMsgTxt("error: 5 inputs and 2 outputs are required\n");
    }
    
    
    WV = mxGetPr(prhs[0]);
    r = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    WW = mxGetPr(prhs[1]);
    mxGH = mxDuplicateArray(prhs[2]);
    GH = mxGetPr(mxGH);
    mxH = mxDuplicateArray(prhs[3]);
    H = mxGetPr(mxH);
    tol =mxGetScalar(prhs[4]);
    
    
    plhs[0] = mxDuplicateArray(prhs[3]);
    Hnew = mxGetPr(plhs[0]);
    plhs[1] = mxDuplicateArray(prhs[4]);
    
    maxiter=r*2;
    flop=goHiter(WV, WW, GH, H, Hnew, r, n, maxiter);
    mxGetPr(plhs[1])[0]=flop;
    
    mxDestroyArray(mxH);/*free memeroy*/
    mxDestroyArray(mxGH);
}

