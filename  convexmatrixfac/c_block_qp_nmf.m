function [U,Y,err] = c_block_qp_nmf(X,r,blk,options)
% nonnegative matrix factorization
% [U,Y,err] = block_qp_nmf(X,r,blk,options)
% min |X - U'Y| s.t. U>=0, Y>=0
% [tol,maxiter,U,Y,display]=getoption(options);
etol=2e-4*r/blk;
if ~exist('options','var'),options=struct;end;
[tol,maxiter,U,Y,display]=getoption(options);

[U,Y,err] = block_nmfmex(X,U',Y',blk,tol,maxiter,etol,display);
U=U';Y=Y';
% normlize
U = diag(1./max(eps,max(U,[],2)))*U;
Y = diag(max(eps,max(U,[],2)))*Y;
err(err==0)=[];
    function [tol,maxiter,U,Y,display]=getoption(options)
        maxiter = 20;
        tol = 1e-4;
        display = 1;        
        [m,n]=size(X);
        if isfield(options,'tol'),tol = options.tol;end
        if isfield(options,'maxiter'),maxiter = options.maxiter;end
        if isfield(options,'U'),U = options.U;else U=ones(r,m)+eye(r,m);end
        if isfield(options,'Y'),Y = options.Y;else Y=ones(r,n);end
        if isfield(options,'display'),display=options.display;end        
    end
end
