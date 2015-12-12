function sol = SDRE(x0,tf,A,B,C,Q,R,F,z)

nPoints = 100;
t = linspace(0,tf,nPoints);
n = length(x0);

%Check the input matrices
if ~isa(A,'function_handle') %This should probably never get used
    warning('Constant A matrix implies the problem is linear and can be solved optimally by a different method')
    A = @(x) A;
end
if ~isa(R,'function_handle')
    R = @(x) R;
end
m = size(R(x0),1);
if ~isa(B,'function_handle')
    B = @(x) B;
elseif isa(B,'function_handle') && (nargin(B) == 2)
    warning('Nonaffine B matrix will be affinized by setting u = 0')
    B = @(x) B(x,zeros(m,1));
else
    error('Input B must be constant, function of state B(x), or state/control B(x,u)')
end
if isempty(C)
    C = @(x) eye(n);
elseif ~isa(C,'function_handle')
    C = @(x) C;
end
if ~isa(Q,'function_handle')
    Q = @(x) Q;
end

if ~isa(F,'function_handle')
    F = @(x) F;
end
%Determine if we're tracking or regulating
tracking = false;
if nargin < 9 || isempty(z)
    % Have to check nargin but don't need to actually do anything here
elseif ~isa(z,'function_handle')
    if any(z ~= 0)
        tracking = true;
        z = @(t) z(:);
    end
else
    tracking = true;  
end
if tracking
    for i = 1:nPoints
        Z(:,i) = z(t(i)); % Loop in case z isn't vectorized.
    end
else
    Z = zeros(n,1);
end

x = [x0(:),zeros(n,nPoints-1)];
for i = 1:nPoints-1
    xc = x(:,i);
    Pss = -care(-A(xc),B(xc),C(xc)'*Q(xc)*C(xc),R(xc)); %We want negative def. Pss
    Acl = A(xc)-B(xc)/R(xc)*B(xc)'*Pss;
    D = lyap(Acl,-B(xc)/R(xc)*B(xc)');
    Kf = inv(C(Z(:,end))'*F(xc)*C(Z(:,end))-Pss);
    K = expm(Acl*(t(i)-tf))*(Kf-D)*expm(Acl'*(t(i)-tf)) + D;
    P{i} = inv(K) + Pss;
    u = @(X) -R(X)\B(X)'*P{i}*X;
    [~,xi] = ode45(@dynamics,[t(i),t(i+1)],xc,[],A,B,u);
    x(:,i+1) = xi(end,:)';
end

sol.state = x;
sol.P = P;
sol.time = t;

end

function dx = dynamics(t,x,A,B,U)
a = A(x);
u = U(x);
b = B(x);
dx = a*x + b*u;
end