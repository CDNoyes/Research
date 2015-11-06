%LINEARIZE Linearize a function.
%   The function FUN is linearized around the point POINT. The output is a
%   function handle that can be evaluated at points near POINT.
%   Linearization is the first order Taylor expansion of the function
%   around the input point. If the first derivative of the function is
%   known, it can be provided; otherwise, it will be approximated
%   numerically.

function linearizedFunction = Linearize(fun,point,firstDerivative)
if nargin == 3 && ~isempty(firstDerivative)
    df = firstDerivative;
else
    df = @(x) ComplexDiff(fun,x);
end
f = fun(point);

linearizedFunction = @(x) f+df(point)*(x-point);

end