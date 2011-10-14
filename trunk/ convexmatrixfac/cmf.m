function [U,Y,err] = cmf(X,r,options)
% convex matrix factorization
% [U,Y,err] = cmf(X,r,options)
% min |X - U'Y| s.t. U>=0, Y>=0 and y'e=1;
% [tol,maxiter,U,Y,display] = getoption(options);

etol=1e-4*r;
if ~exist('options','var'),options=struct;end;
[tol,maxiter,U,Y,display]=getoption(options);
nrm2=norm(X,'fro')^2;
err=zeros(maxiter,1);obj=0;%Xt=X';
for iter =1:maxiter
    A = U*U'; b=U*X;
    d=svd(A);ss=0;
    if d(end)/d(1)<etol,ss=d(1)*etol;end
    [Y,tab] = qpmex(A+eye(size(A))*ss,b+Y*ss,2); % convex
    
    A = Y*Y';b= (X*Y')';%b= Y*Xt;
    d=svd(A);ss=0;
    if d(end)/d(1)<etol,ss=d(1)*etol;end
    [U,tab] = qpmex(A+eye(size(A))*ss,b,1); % nonnegative

    obj_=obj;
    obj=CalculateObj(U, Y);err(iter)=obj;
    if display,
        fprintf('iter: %d obj:%6.4e \n',iter,obj);
    end
    if abs(obj_-obj)<tol*obj,break;end
end % end of main
err=err(1:iter);

    function obj = CalculateObj(U, Y)
        obj = (nrm2+sum(sum((U*U').*(Y*Y'))))/2- sum(sum(Y.*(U*X)));
    end
    function [tol,maxiter,U,Y,display]=getoption(options)
        maxiter = 20;
        tol = 1e-4;
        display = 1;

        [m,n]=size(X);
        if isfield(options,'tol'),tol = options.tol;end
        if isfield(options,'maxiter'),maxiter = options.maxiter;end
        if isfield(options,'U'),U = options.U;else U=ones(r,m)+eye(r,m);end
        if isfield(options,'Y'),Y = options.Y;else Y=zeros(r,n);end
        if isfield(options,'display'),display=options.display;end

    end
end
