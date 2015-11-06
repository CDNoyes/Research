function [x,fval,funcCount,S] = Penalty(f,c,x0)
%Implements a penalty method to solve constrained optimization. This
%implementation uses MATLAB's fminsearch as the unconstrained optimizer.

sigma = f(x0)/norm(min(0,c(x0)).^2);
S.initial = sigma;
x = x0;
xold = 1e6;
tol = 1e-3;
funcCount = 0;
while abs(x-xold) > tol
    phi = @(x) f(x) + sum(sigma*min(0,c(x)).^2);
    xold = x;
    [x,fval,~,output] = fminsearch(phi,xold);
    sigma = sigma*10;
    funcCount = output.funcCount + funcCount;
end
S.final = sigma/10; %Since it gets multiplied by 10 before actually being used again
end