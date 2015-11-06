function [x,fval,funcCount] = Penalty(f,c,x0)
%Implements a penalty method to solve constrained optimization. This
%implementation uses MATLAB's fminsearch as the unconstrained optimizer.

sigma = 0.01;
x = x0;
xold = 1e6;
tol = 1e-5;
funcCount = 0;
while abs(x-xold) > tol
phi = @(x) f(x) + sigma*min(0,c(x))^2;

xold = x;
[x,fval,~,output] = fminsearch(phi,xold);
sigma = sigma*10;
funcCount = output.funcCount + funcCount;
end

end