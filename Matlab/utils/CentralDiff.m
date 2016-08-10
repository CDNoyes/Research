%CENTRALDIFF Computes the central difference approximation to the
%derivative of a vector-valued function F at the point X. The step size H
%is an optional input.
%The function f should return a column vector.
function dfdx = CentralDiff(f,x,h)

if nargin < 3
    h = 1e-8;
end
n = length(x);
x = x(:);
xhf = repmat(x,1,n)+0.5*h*eye(n);
xhb = repmat(x,1,n)-0.5*h*eye(n);

for i = 1:n
    dfdx(:,i) = (f(xhf(:,i))-f(xhb(:,i)))/h;
end
end