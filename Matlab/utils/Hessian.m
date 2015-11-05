function H = Hessian( fun, x, stepSize )
%HESSIAN Numerically computes the hessian of a function.
%   The approximation is performed by using a complex step method to
%   estimate the gradient/jacobian and using a central difference scheme to
%   compute the hessian. If the gradient/jacobian of the function is known
%   analytically, ComplexDiff should be used on the gradient to obtain a
%   more accurate solution.

df = @(x) ComplexDiff(fun,x);
fx = fun(x);
if isscalar(fx)
    if nargin < 3
        H = CentralDiff(df,x);
    else
        H = CentralDiff(df,x,stepSize);
    end
else %Vector valued function
    n = length(fx);
    m = length(x);
    H = nan(m,m,n);
    for i = 1:n
        H(:,:,i) = CentralDiff(@(X)wrapper(df,X,i),x);
    end
    
end
end

function dfi = wrapper(df,x,index)
dvalue = df(x); %A matrix
dfi = dvalue(index,:); %The gradient of the scalar "index" function wrt the the states

end