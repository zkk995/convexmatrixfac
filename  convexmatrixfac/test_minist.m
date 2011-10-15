% function test_minist
X =combinedata('F:\ImageDatabase\MINIST_data\zip\Data');
Y =X; 
%%
% Y = 1-X;
 r= 20;
[m,n] = size(Y); maxiter =50;
opt.maxiter =maxiter;opt.tol = 1e-26;
U = ones(m,r)+120*eye(m,r)+2*rand(m,r);V = ones(r,n)+50*eye(r,n);
U0 = 1.0*Y(:,(1:r)+243)+0+0.*eye(m,r)+0*rand(m,r);% 23

opts=opt;opts.U=rand(r,m);opts.bound=1;
opts.lambda=1e-5;
tic;[U1,V1,err] = cmf(Y,r,opts);toc;
U1=U1';
%%
figure(2),subplot(3,1,1),showbaseimg(U0)
figure(2),subplot(3,1,2),showbaseimg(U1)
return
%%
width=.95; height=.45;
left=0.025;  bottom=0.5;
subplot('Position',[left bottom width height])
showbaseimg(U),title('some samples')
% width=1; %height=0.35;
% left=0; 
bottom=0.0;
subplot('Position',[left bottom width height])
showbaseimg(U1), title('recovered basis')

return
%%

[m,n] = size(Y); maxiter =30;
opt.maxiter =maxiter;opt.tol = 1e-26;
U = 1*ones(m,r)+5*eye(m,r);V = ones(r,n)+50*eye(r,n);
tic;[U1,V1,res]=qpas_nmf(Y,r,opt,U,V);toc;

%%
tic;[U1,V1,res]=cqpas_nmf(Y,r,opt,U1,V);toc;

%%
opt.maxiter =maxiter;opt.tol = 1e-26;
U = X(1:r,:);V = ones(r,n)+50*eye(r,n);
tic;[U1,V1,res]=cqpas_nmf(Y,r,opt,U,V);toc;

%%

