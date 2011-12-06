function [W,H,difflist,timelist,floplist] = mGCD(V, maxiter, W, H)
% min |V - W'H| s.t.  W, H >=0
tol = 0.001;
if size(W,1)~=size(H,1),error('error:mgcd','incorrect size of W and H');end

total = 0;flop = 0;obj=0;
difflist = zeros(1,maxiter);timelist = difflist;floplist = difflist;
Hnew=zeros(size(H));Vt=V';nrmV2=norm(V,'fro')^2;
for iter = 1:maxiter
    begin = cputime;
    
    [W,Wnew,flop]=fupdateH(Vt,H,W,Hnew,flop,tol);
    [H,Hnew,flop]=fupdateH(V,W,H,Wnew,flop,tol);
    
    total = total + cputime - begin;
    obj_=obj;
    obj = CalculateObj(V,W, H,nrmV2);
    
    difflist(iter) = obj;
    timelist(iter) = total;
    floplist(iter) = flop;
    if abs(obj-obj_)<1e-16*obj;break;end
    fprintf('iter %d, gcd obj %4.3e \n',iter,obj);
end
% difflist=difflist(1:iter);timelist=timelist(1:iter);
% floplist=floplist(1:iter);
return
function obj = CalculateObj(X,U, Y,nrm2)
     obj = (nrm2+sum(sum((U*U').*(Y*Y'))))/2- sum(sum(Y.*(U*X)));
return

function [H,Hnew,flop]=fupdateH(V,W,H,Wnew,flop,tol)
% update variables of H
WV = W*V;
WW = W*W';
GH = -(WV-WW*H);
flop = flop + FlopMul(V,Wnew')+FlopMul(W,W)+FlopMul(WW,H);
flop = flop + numel(H);

% 	whos WV WW GH H tol
[Hnew flopadd] = updateH(WV,WW,GH,H,tol);
flop = flop + flopadd;
H = H + Hnew;
return

function f = FlopMul(A,B)

if issparse(A)
    f = nnz(A)*size(B,2)*2;
else
    f = size(A,1)*size(A,2)*size(B,2)*2;
end
return