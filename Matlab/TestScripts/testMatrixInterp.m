clear; clc;

M{2} = 2*ones(4,1)';
M{1} = zeros(4,1)';
t = [0,1];
ti = linspace(0,1,10000);

tic
Mi = MatrixInterp(t,M,ti);
toc
tic
Mii = MatrixInterp(t,M,ti);
toc