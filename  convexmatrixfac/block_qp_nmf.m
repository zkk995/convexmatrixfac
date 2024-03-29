function [U,Y,err] = block_qp_nmf(X,r,blk,options)
% nonnegative matrix factorization
% [U,Y,err] = block_qp_nmf(X,r,blk,options)
% min |X - U'Y| s.t. U>=0, Y>=0
% [tol,maxiter,U,Y,display]=getoption(options);
etol=2e-4*r/blk;
if ~exist('options','var'),options=struct;end;
[tol,maxiter,U,Y,display]=getoption(options);
nrm2=norm(X,'fro')^2;
err=zeros(maxiter,1);obj=0;% Xt=X';

blk=fix(blk);
if blk<1||r/blk<10||r/blk>25,warning('nmf:blk',' choose a bad block number');end
rem_=rem(r,blk);
I=cell(blk,1);
bn=fix(r/blk);j=1;
for i=1:blk
    I_=false(r,1);I_(j:j+bn-1)=true;
    I{i}  = I_;j=j+bn;
    if i==blk-rem_,bn=bn+1;end
end

for iter =1:maxiter
    if 1
        for i =1:blk
            A = U(I{i},:)*U(I{i},:)'; b=U(I{i},:)*X-(U(I{i},:)*U(~I{i},:)')*Y(~I{i},:);
            d=svd(A);ss=0;
            if d(end)/d(1)<etol,ss=d(1)*etol;end
            [Y(I{i},:),~] = qpmex(A+eye(size(A))*ss,b+Y(I{i},:)*ss,1); %
            A = Y(I{i},:)*Y(I{i},:)';b=(X*Y(I{i},:)')'-Y(I{i},:)*Y(~I{i},:)'*U(~I{i},:);%b= Y*Xt;
            d=svd(A);ss=0;
            if d(end)/d(1)<etol,ss=d(1)*etol;end
            [U(I{i},:),~] = qpmex(A+eye(size(A))*ss,b+U(I{i},:)*ss,1); % nonnegative
        end;
    else
        for i =1:blk
            A = U(I{i},:)*U(I{i},:)'; b=U(I{i},:)*X-U(I{i},:)*U(~I{i},:)'*Y(~I{i},:);
            d=svd(A);ss=0;
            if d(end)/d(1)<etol,ss=d(1)*etol;end
            [Y(I{i},:),~] = qpmex(A+eye(size(A))*ss,b+Y(I{i},:)*ss,1); %
        end
        for i =1:blk
            A = Y(I{i},:)*Y(I{i},:)';b=(X*Y(I{i},:)')'-Y(I{i},:)*Y(~I{i},:)'*U(~I{i},:);%b= Y*Xt;
            d=svd(A);ss=0;
            if d(end)/d(1)<etol,ss=d(1)*etol;end
            [U(I{i},:),~] = qpmex(A+eye(size(A))*ss,b+U(I{i},:)*ss,1); % nonnegative
        end
        
    end
    % check stop condition  
    obj_=obj;
    obj=CalculateObj(U, Y);err(iter)=obj;
    if display,
        fprintf('iter: %d obj:%6.4e \n',iter,obj);
    end
    % normlize
    U = diag(1./max(eps,max(U,[],2)))*U;
    Y = diag(max(eps,max(U,[],2)))*Y;
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
        if isfield(options,'Y'),Y = options.Y;else Y=ones(r,n);end
        if isfield(options,'display'),display=options.display;end
        
    end
end

