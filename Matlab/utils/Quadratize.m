% QUADRATIZE Computes a quadratic expansion.
%   QUADRATIZE(F,X,J,H) Returns a function handle that gives the second
%   order expansion of the function F around the point X. J and H are the
%   jacobian and hessian of F. They are optional, if not given they will be
%   estimated numerically.

function quadFunction = Quadratize(fun,point,firstDerivative,secondDerivative)

if nargin >= 3 && ~isempty(firstDerivative)
    df = firstDerivative;
    if nargin == 4 && ~isempty(secondDerivative)
        H = secondDerivative;
    else
        H = @(x) Hessian(fun,x);
    end
else
    df = @(x) ComplexDiff(fun,x);
    H = @(x) Hessian(fun,x);
end
f = fun(point);
n = length(f);
quadFunction = @(x) f+df(point)*(x-point)+0.5*HessianMultiply(H(point),x-point,n);

end

% function xHx = hessianMultiply(H,x,n)
% 
% %X' = full(kron(sparse(eye(n)),x'));
% %X = repmat(X,n,1)
% 
% xHx = full(kron(sparse(eye(n)),x'))*H*repmat(x,n,1);
% 
% end