function [X,sn]=combinedata(baseDir)
% mixture of the minist digits  images
if nargin<1, baseDir='.';end;
sn=1000;
file=[baseDir '/mixtureminist.mat'];
try
    load(file);
catch ME1
    if strcmp(ME1.identifier,'MATLAB:load:couldNotReadFile')   
        files = cell(10,1);
        Ind = 1:sn;
        Xt = zeros(sn*10,784);
        for i=1:10,
            files{i} = [baseDir,'/digit' num2str(i-1) '.mat'];
            load(files{i},'D');
            Xt(Ind,:) = D(1:sn,:);Ind=Ind+sn;
        end;
        P = randperm(sn*10);
        X = Xt(P(1:end/2),:)+Xt(P(end/2+1:end),:);
        X = X/2;
        X = double(X')/255;
        save([baseDir '/mixtureminist.mat'], 'X', 'sn');
    else
        rethrow(ME1);
    end
end
return
%%
[X,sn]=combinedata;
for i=1:50,
    Xi=reshape(X(i,:),28,28);
    imshow(1-Xi');
    pause;
end

