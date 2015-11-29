clear; clc;
size = [2,2];
n = 100;
for i = 1:n
    M{i} = (i-1)*2/n*ones(size) + .1*rand(size);
end

t = linspace(0,1,n);
ti = linspace(0,1,10000);

tic
Mi = MatrixInterp(t,M,ti);
toc
tic
Mii = MatrixInterp(t,M,ti);
toc