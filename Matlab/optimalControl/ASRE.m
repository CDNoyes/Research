%ASRE Solves a nonlinear regulation/tracking problem via an Approximating
%Sequence of Riccati Equations methodology.
%   ASRE(X0,TF,M,A,B,C,Q,R,F)
%   A(x) is the state-dependent coefficient matrix
%   B(x,u) is the 
%
% References:
% State Regulation:
% Global optimal feedback control for general nonlinear systems with
% nonquadratic performance criteria, Cimen and Banks 2004
%
% Tracking Control:
% Nonlinear optimal tracking control with application to super-tankers for
% autopilot design, Cimen and Banks, 2004

function sol = ASRE(x0,tf,m,A,B,C,Q,R,F)

%Handle constant matrices
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
    C = @(x) eye(length(x0));
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


iter = 0;
iterMax = 25;
tol = 0.01; %Seems large
diff = tol+1;

% tf = OCP.bounds.upper.finalTime;
% x0 = OCP.bounds.upper.initialState;
% x0 = x0(:);
% n = OCP.dimension.state;
% m = OCP.dimension.control;

nPoints = 10000;
t = linspace(0,tf,nPoints);

while iter < iterMax && diff > tol
    if ~iter % First iteration is LTI
        %Integrate the finite-time riccati equation backward
        [tb,Pvec] = ode45(@Riccati,tf-t,reshape(C(x0)'*F(x0)*C(x0),[],1),[],...
            @(X) A(x0),...
            @(X,U) B(x0,0),...
            @(X) C(x0),...
            @(X) Q(x0),...
            @(X) R(x0),...
            @(X)x0,...
            @(T)zeros(m,1));
        
        %Compute the states
        [~,x{iter+1}] = ode45(@dynamics,t,x0,[],...
            @(X) A(x0),...
            @(X,U) B(x0,0),...
            @(X) R(x0),...
            @(T) interp1(tb,Pvec,T),...
            @(T) 0);
        
        %Compute the controls
        u{iter+1} = computeControl(@(X,U) B(x0,0), @(X) R(x0),Pvec,x{iter+1}',zeros(m,nPoints));
        
    else % Recursive LTV systems
        %Integrate the finite-time riccati equation backward
        Pf = reshape(C(x{iter}(end,:))'*F(x{iter}(end,:))*C(x{iter}(end,:)),[],1);
        [tb,Pvec] = ode45(@Riccati,tf-t,Pf,[],A,B,C,Q,R,...
            @(T) interp1(t,x{iter},T),@(T) interp1(t,u{iter},T));
        
        %Compute the states
        [~,x{iter+1}] = ode45(@dynamics,t,x0,[],...
            A,...
            B,...
            R,...
            @(T) interp1(tb,Pvec,T),...
            @(T) interp1(t,u{iter},T));

        %Compute the controls
        u{iter+1} = computeControl(B,R,Pvec,x{iter}',u{iter});
        diff = norm(x{iter+1}-x{iter});
    end
    
    iter = iter + 1;
    disp(['Iteration ',num2str(iter),' complete.'])
end

sol.state = x;
sol.control = u;
sol.time = t;

end

function U = computeControl(B,R,Pvec,x,u)
l = sqrt(size(Pvec,2));
for i = 1:size(x,2)
    P = reshape(Pvec(end-i+1,:),l,l);
    U(:,i) = -R(x(:,i))\B(x(:,i),u(:,i))'*P*x(:,i);
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

function dx = dynamics(t,x,A,B,R,P,U)
Pt = P(t);
n = sqrt(length(Pt));
p = reshape(Pt,n,n);
a = A(x);
u = U(t);
b = B(x,u);
r = R(x);
s = (b/r)*b';
dx = (a - s*p)*x;

end
