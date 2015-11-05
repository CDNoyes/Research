function H = Hessian( fun, x, stepSize )
%HESSIAN Numerically computes the hessian of a scalar function. 
%   The approximation is performed by using a complex step method to
%   estimate the gradient and using a central difference scheme to compute
%   the hessian. If the gradient of the function is known analytically,
%   ComplexDiff should be used on the gradient to obtain a more accurate
%   solution.

df = @(x) ComplexDiff(fun,x);
if nargin < 3
    H = CentralDiff(df,x);
else
    H = CentralDiff(df,x,stepSize);
end
end

