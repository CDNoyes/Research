%test out Submatrix on Jacobians and Hessians

%VDP Oscillator
clear; clc;
[F,J,H] = VDP();
[dim,ind] = Dimension(2,1,1);

x0 = [3;5];
u0 = 0;
mu = 0.25;
X0 = [x0;u0;mu];

%Analytical Results
f = F(0,x0,u0,mu);
j = J(x0,u0,mu);
h = H(x0,u0,mu);

jac = Submatrix(j,dim,ind);
hess = Submatrix(h,dim,ind);
% hess.xx == hxx?
fx = @(X) F(0,X,u0,mu);
hxx = Hessian(fx,x0);