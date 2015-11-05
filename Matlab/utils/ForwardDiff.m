function dfdx = ForwardDiff(f,x,h)

if nargin < 3
    h = 1e-8;
end
n = length(x);
fx = f(x);
xh = repmat(x(:),1,n)+h*eye(n);
for i = 1:n
    dfdx(:,i) = (f(xh(:,i))-fx)/h;
end
end