%test out Submatrix on Jacobians and Hessians

%% VDP Oscillator
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

%% Test the scalar functionality
%Crazy made-up function of 3 states, 2 controls, and 1 parameter
clc; clear;
f = @(x,u,p) sum(x)*x(1)*u(1)/(u(2)+p^2);
[dim,ind] = Dimension(3,2,1);

F = @(X) f(X(ind.state),X(ind.control),X(ind.parameter));
X = [1;1;1;1;1;-1];
J = ComplexDiff(F,X);
H = Hessian(F,X);
j = Submatrix(J,dim,ind,true);
h = Submatrix(H,dim,ind,true);