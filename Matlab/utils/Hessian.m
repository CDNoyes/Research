function [H,h] = Hessian( fun, x, stepSize )
%HESSIAN Numerically computes the hessian of a function.
%   HESSIAN(FUN,X,STEPSIZE) approximates the hessian of function FUN at the
%   point X using optionally the stepsize STEPSIZE for finite differencing.
%   The approximation is performed by using a complex step method to
%   estimate the gradient/jacobian and using a central difference scheme to
%   compute the hessian. If the gradient/jacobian of the function is known
%   analytically, ComplexDiff should be used on the gradient to obtain a
%   more accurate solution.
%
%   The Hessian is output in two formats. The first is a special 2D matrix
%   that facilitates a simple matrix multiplication form (for e.g. x'*H*x)
%   in the event that fun is a vector-valued function; if the function is
%   scalar then this representation is equivalent to the standard form of
%   the hessian. The second output has each "page" of the Hessian in its
%   own cell.

df = @(x) ComplexDiff(fun,x);
fx = fun(x);
if isscalar(fx)
    if nargin < 3
        H = CentralDiff(df,x);
        h{1} = H;
    else
        H = CentralDiff(df,x,stepSize);
        h{1} = H;
    end
else %Vector valued function
    n = length(fx);
    h = cell(1,n);
    for i = 1:n
        h{i} = CentralDiff(@(X)wrapper(df,X,i),x);
        h{i} = 0.5*(h{i}+h{i}.'); %Remove small asymmetries due to numerical error
    end
    H = sparse(blkdiag(h{:})); %Convert from an array of cells to a single 2d matrix

end
end

function dfi = wrapper(df,x,index)
dvalue = df(x); %A matrix
dfi = dvalue(index,:); %The gradient of the scalar "index" function wrt the the states

end