code for our paper which is about the topic of nonnegative matrix factorization
```
% convex matrix factorization
% [U,Y,err] = cmf(X,r,options)
% min |X - U'Y|^2 + t|U|^2 s.t. U>=0, Y>=0 and y'e=1;


% nonnegative matrix factorization
% [U,Y,err] = qp_nmf(X,r,options)
% min |X - U'Y| s.t. U>=0, Y>=0
% [tol,maxiter,U,Y,display] = getoption(options);
```