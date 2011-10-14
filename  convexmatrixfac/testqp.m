n=41;
A1=rand(n);
A =A1*A1'+eye(n);
yy=rand(n,20000)/n;

b= A*yy ;

tic;[y,tab] = qpmex(A,b,2);toc;