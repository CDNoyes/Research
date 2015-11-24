%% Test the constrained QP algorithm
clear; clc;

n = 40;
H = rand(n);
H = Regularize(0.5*(H+H'),.1); %Symmetric, pos-def hessian
g = -2+4*rand(n,1);
x0 = 10*rand(n,1);
% x0 = [0,5,-3,0]';
xu = 5*ones(n,1);
xl = -3*ones(n,1);

opt = optimset('Algorithm','interior-point-convex');
tic
[xMatlab,fMatlab] = quadprog(H,g,[],[],[],[],xl,xu,[],opt);
toc
tic
[x,f] = ProjectedNewtonQP(H,g,x0,xl,xu,1e-3);
toc