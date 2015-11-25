%From Control-Limited DDP by Tass, Mansard, Todorov
%XUPPER and XLOWER represent the bounds on the solution, and each is a vector with the same length as x0.
function [x,fOpt,hff] = ProjectedNewtonQP(hess,grad,x0,xlower,xupper,tol)

%Preliminary problem info:
n = length(x0);
if nargin == 5 || isempty(tol)
    tol = 1e-6;
end
iter = 0;
iterMax = Saturate(1e-4/tol,50,500);

%Input checking
if isscalar(xlower) 
    xlower = xlower*ones(n,1);
elseif length(xlower) ~= n
    error('Lower bound must be scalar or the same length as x.')
end
if isscalar(xupper) 
    xupper = xpper*ones(n,1);
elseif length(xupper) ~= n
    error('Upper bound must be scalar or the same length as x.')
end

%Make the initial point feasible:
x = Saturate(x0(:),xlower(:),xupper(:));

while iter < iterMax
    % Gradient of the objective function:
    g = grad + hess*x;
    
    %Determine the complimentary sets based on constraints:
    idu = (x-xupper) == 0;  % The elements on the upper boundary
    idl = (xlower-x) == 0;  % The elements on the lower boundary
    gu = g>0;               % Elements where the objective function is increasing
    gl = g<0;               % Elements where the objective function is decreasing
    %Find the elements that are at max and growing, or at min and decreasing.
    ci = (gl+idu==2)+(gu+idl==2);
    c = find(ci,n);
    f = find(~ci,n);
    
    gf = grad(f) + hess(f,f)*x(f) + hess(f,c)*x(c);
    if norm(gf) < tol
        disp(['ProjectedNewton QP terminated successfully in ',num2str(iter),' iterations.'])
        break
    end
    dx = zeros(n,1);
    dx(f,1) = -inv(hess(f,f))*gf;

    hff = hess(f,f);
    alpha = armijo(@(X) fQuad(hess,grad,X),x,dx,g,xlower,xupper);
    x = Saturate(x+alpha*dx,xlower,xupper);
    iter = iter+1;
end
if iter == iterMax
    disp('Max number of iterations reached.')
end
fOpt = fQuad(hess,grad,x);

end

function alpha = armijo(f,x,dx,g,xl,xu)
gamma = 0.1;
c = 0.5;
alpha = 2*max(xu-xl);
r = gamma/2; % Initialize
while r < gamma
    alpha = c*alpha;
    xa = Saturate(x+alpha*dx,xl,xu);
    r = (f(x)-f(xa))/(g'*(x-xa));
end

end

function f = fQuad(H,q,x)
f = 0.5*x'*H*x + q'*x;
end
