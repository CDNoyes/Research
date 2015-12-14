%ASRE Solves a nonlinear regulation/tracking problem via an Approximating
%Sequence of Riccati Equations methodology.
%   ASRE(X0,TF,M,A,B,C,Q,R,F,z)
%   X0 is the initial state
%   TF is the fixed final time
%   M is the dimension of the control vector
%   A(x) is the state-dependent coefficient matrix
%   B(x,u) is the state-dependent control matrix
%   C(x) is the state-dependent output matrix, optional
%   Q(x) is the state/output/error weight matrix in the Lagrange cost
%   R(x) is the control weight matrix in the Lagrange cost
%   F(x) is the state/output weight matrix of the Mayer cost
%   z(t) is a function defining a reference trajectory to be followed and
%        is optional. It should return a column vector with (1 <= l <= n)
%        elements. It z is a constant vector the problem is interpreted as
%        regulating the state with fixed final conditions Cx = z(tf). If
%        instead a constant reference trajectory is desired, z(t) should be
%        a function handle returning a constant value or vector.
%
% References:
% State Regulation:
% Global optimal feedback control for general nonlinear systems with
% nonquadratic performance criteria, Cimen and Banks 2004
%
% Tracking Control:
% Nonlinear optimal tracking control with application to super-tankers for
% autopilot design, Cimen and Banks, 2004

function sol = ASRE(x0,tf,A,B,C,Q,R,F,z)

% tf = OCP.bounds.upper.finalTime;
% x0 = OCP.bounds.upper.initialState;
% x0 = x0(:);
% n = OCP.dimension.state;
% m = OCP.dimension.control;

interpType = 'nearest';
nPoints = 500; % Very little impact on solution quality
t = linspace(0,tf,nPoints);
tb = tf - t;

n = length(x0);

%Check the input matrices
if ~isa(A,'function_handle') %This should probably never get used
    warning('Constant A matrix implies the problem is linear and can be solved optimally by a different method')
    A = @(x) A;
end
if ~isa(B,'function_handle')
    B = @(x,u) B;
elseif isa(B,'function_handle') && (nargin(B) == 1)
    B = @(x,u) B(x);
elseif isa(B,'function_handle') && (nargin(B) == 2)
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
if ~isa(R,'function_handle')
    R = @(x) R;
end
if ~isa(F,'function_handle')
    F = @(x) F;
end
%Determine if we're tracking or regulating
tracking = false;
fixedFS = false; % Any endpoint constraints?
if nargin < 9 || isempty(z)
    disp('Problem Type: Regulation')
elseif ~isa(z,'function_handle')
    fixedFS = true;
    disp('Problem Type: Regulation with state constraints at final time.')
    l = length(z);
else
    disp('Problem Type: Reference Tracking')
    tracking = true;  
    for i = 1:nPoints
        Z(:,i) = z(t(i)); % Loop in case z isn't vectorized.
    end
end

m = size(R(x0),1);
iter = 0;
iterMax = 25;
tol = 0.01; 
diff = tol+1;

while iter < iterMax && diff > tol
    if ~iter % First iteration is LTI
        
        % Integrate the finite-time riccati equation backward
        if fixedFS
            Pf = reshape(F(x0),[],1);
            Cfun = @(X) eye(n);
        else
            Pf = reshape(C(x0)'*F(x0)*C(x0),[],1);
            Cfun = @(X) C(x0);
        end
        [~,Pvec] = ode45(@Riccati,tb,Pf,[],...
            @(X) A(x0),...
            @(X,U) B(x0,0),...
            Cfun,...
            @(X) Q(x0),...
            @(X) R(x0),...
            @(X)x0,...
            @(T)zeros(m,1));
        
        % Compute the feedforward term if we're tracking
        if tracking
            sf = C(x0)'*F(x0)*z(tf);
            [~,s] = ode45(@feedforward,tb,sf,[],...
                @(X) A(x0),...
                @(X,U) B(x0,zeros(m,1)),...
                @(X) C(x0),...
                @(X) Q(x0),...
                @(X) R(x0),...
                @(T) interp1(tb,Pvec,T,interpType),...
                @(X)x0,...
                @(T)zeros(m,1),...
                z);
        elseif fixedFS
            Vf = C(x0)';
            Gf = zeros(l);
            [~,VG] = ode45(@fixedGains,tb,[Vf(:);Gf(:)],[],...
                @(X) A(x0),...
                @(X,U) B(x0,zeros(m,1)),...
                @(X) R(x0),...
                @(T) interp1(tb,Pvec,T,interpType),...
                @(X) x0,...
                @(T) zeros(m,1),l);
        else
            s = zeros(nPoints,n);
        end      
        
        if fixedFS
            % Compute the states:
            
            % Compute the controls
            computeFixedControl(B,R,Pvec,x,u,Vvec,Gvec,r)
        else
            % Compute the states
            [~,x{iter+1}] = ode45(@dynamics,t,x0,[],...
                @(X) A(x0),...
                @(X,U) B(x0,0),...
                @(X) R(x0),...
                @(T) interp1(tb,Pvec,T,interpType),...
                @(T) 0,...
                @(T) interp1(tb,s,T,interpType)');
            
            % Compute the controls
            u{iter+1} = computeControl(@(X,U) B(x0,0), @(X) R(x0),Pvec,...
                x{iter+1}',zeros(m,nPoints),s');
        end
    else % Recursive LTV systems
        
        %Integrate the finite-time riccati equation backward
        Pf = reshape(C(x{iter}(end,:))'*F(x{iter}(end,:))*C(x{iter}(end,:)),[],1);
        [~,Pvec] = ode45(@Riccati,tb,Pf,[],A,B,C,Q,R,...
            @(T) interp1(t,x{iter},T,interpType),@(T) interp1(t,u{iter},T,interpType));
        
        % Compute the feedforward term if we're tracking
        if tracking
            sf = C(x{iter}(end,:))'*F(x{iter}(end,:))*z(tf);
            [~,s] = ode45(@feedforward,tb,sf,[],A,B,C,Q,R,...
                @(T) interp1(tb,Pvec,T,interpType),...
                @(T) interp1(t,x{iter},T,interpType),...
                @(T) interp1(t,u{iter},T,interpType),...
                z);
        else
            s = zeros(nPoints,n);
        end      
        
        % Integrate the states
        [~,x{iter+1}] = ode45(@dynamics,t,x0,[], A,B,R,...
            @(T) interp1(tb,Pvec,T,interpType),...
            @(T) interp1(t,u{iter},T,interpType),...
            @(T) interp1(tb,s,T,interpType)');

        % Compute the controls
        u{iter+1} = computeControl(B,R,Pvec,x{iter}',u{iter},s');
        
        % Compute the convergence criteria
        diff = norm(x{iter+1}-x{iter});
    end
    
    iter = iter + 1;
    disp(['Iteration ',num2str(iter),' complete.'])
end

if (iter == iterMax) && diff > tol
    fprintf(['ASRE failed to converge in ',num2str(iterMax),' iterations.\n',...
        'Consider raising the maximum number of iterations or\ndecreasing',...
        ' the convergence tolerance (currently set to ',num2str(tol),').\n'])
end

for i = nPoints:-1:1 
    P{i} = reshape(Pvec(i,:),n,n);
end

sol.state = x{end};
sol.control = u{end}; %Open loop control
sol.P = P;
sol.s = s';
sol.time = t;
sol.history.state = x;
sol.history.control = u;

end

function U = computeControl(B,R,Pvec,x,u,s)
l = sqrt(size(Pvec,2));
for i = 1:size(x,2)
    P = reshape(Pvec(end-i+1,:),l,l);
    U(:,i) = -R(x(:,i))\B(x(:,i),u(:,i))'*(P*x(:,i)-s(:,end-i+1));
end
end

function dP = Riccati(t,P,A,B,C,Q,R,X,U)
n = sqrt(length(P));
p = reshape(P,n,n);
x = X(t);
u = U(t);
a = A(x);
b = B(x,u);
c = C(x);
q = Q(x);
r = R(x);
s = (b/r)*b';
dP = reshape(-c'*q*c - p*a - a'*p + p*s*p,[],1);

end

function dx = dynamics(t,x,A,B,R,P,U,s)
Pt = P(t);
n = sqrt(length(Pt));
p = reshape(Pt,n,n);
a = A(x);
u = U(t);
b = B(x,u);
r = R(x);
S = (b/r)*b';
dx = (a - S*p)*x + S*s(t);

end

function ds = feedforward(t,s,A,B,C,Q,R,P,X,U,z)
Pt = P(t);
n = sqrt(length(Pt));
p = reshape(Pt,n,n);
x = X(t);
u = U(t);
a = A(x);
b = B(x,u);
c = C(x);
r = R(x);
q = Q(x);
S = (b/r)*b';
ds = -(a - S*p)'*s - c'*q*z(t);

end

% Additional Functions for the Fixed final state version of the problem
function dVG = fixedGains(t,VG,A,B,R,P,X,U,l)
Pt = P(t);
n = sqrt(length(Pt));
p = reshape(Pt,n,n);
V = reshape(VG(1:l*n),l,n); %Turn vector into matrix if needed
x = X(t);
u = U(t);
a = A(x);
b = B(x,u);
r = R(x);
K = r\b'*p;
dV = -(a-b*K)'*V;
dG = V'*b/r*b'*V;
dVG = [dV(:);dG(:)];
end

function U = computeFixedControl(B,R,Pvec,x,u,Vvec,Gvec,r)
n = sqrt(size(Pvec,2));
l = sqrt(size(Gvec,2));
for i = 1:size(x,2)
    P = reshape(Pvec(end-i+1,:),n,n);
    V = reshape(Vvec(end-i+1,:),n,l);
    G = reshape(Gvec(end-i+1,:),l,l);
    U(:,i) = -R(x(:,i))\B(x(:,i),u(:,i))'*((P + V/G*V')*x(:,i) + V/G*r);
end
end