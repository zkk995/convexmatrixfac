function showbaseimg(U1)
basenum = size(U1,2);
sdim = sqrt(size(U1,1));
if (sdim^2)~= size(U1,1),error('error in dimension');end
T = zeros(sdim,basenum*sdim);Ind =1:sdim;
for i=1:basenum
    Ti=U1(:,i);
    T(:,Ind) = reshape(Ti,sdim,sdim)';Ind =Ind+sdim;
end
if mod(basenum,2)==1
    T=[T,zeros(sdim)];
end

TT=[T(:,1:end/2);T(:,1+end/2:end)];
imshow(TT)
