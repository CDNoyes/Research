%CTDDP Continous-time differential dynamic programming implementation.
%   CTDDP(OptimalControlProblem) Solves an OCP via an implementation of the
%   algorithm described in "Continous-time DDP with Terminal Constraints"
%   by Sun et al.

function sol = CTDDP(OCP,order)

if ~OCP.fixedFinalTime
    error('This continous time DDP formulation currently requires a fixed final time.')
end
if nargin < 2 || isempty(order)
    OCP.order = 0; % Order [0,1]
else
    OCP.order = order;
end

%Check gradients/jacobians/hessians, if empty compute them numerically


n = OCP.dimension.state;
m = OCP.dimension.control;
p = OCP.dimension.adjoint;

[~,OCP.ind] = Dimension(n,m,0);

nPoints = 250;
t = linspace(0,OCP.bounds.upper.finalTime,nPoints);
u = zeros(m,nPoints); %Initial control guess of all zeros.
x0 = OCP.bounds.upper.initialState;
fixed = ~OCP.free.finalState; %Logical indices of the state vector elements that are fixed at tfinal, length p
if p == 0 && any(fixed) %This means we treat constraints softly
    p = sum(fixed);
    OCP.dimension.adjoint = p;
    softConstraints = true; 
end
lambda = zeros(p,1);
constraintGrad = Psi_X(fixed,p,n);

if OCP.order == 0
    eta = .3; % A parameter to be set
    gamma = 0.4;
else
    eta = 1;
    gamma = 1;
end

tol = 1e-4;
iter = 0;
iterMax = 15;
constraintViolation = tol+1;
opt = odeset('AbsTol',1e-8,'RelTol',1e-8);

% The DDP loop
while constraintViolation > tol && iter < iterMax
    
    %Integrate dynamics forward with current control
    [~,x] = ode45(OCP.dynamics,t,x0,[], @(T) interp1(t,u,T));
    constraintViolation = max(abs(x(end,fixed)-OCP.bounds.upper.finalState(fixed)'));
    
    %Compute V and its partials, eq. 28, at tFinal
    V = OCP.cost.mayer(t(end),x(end,:)') +...
        (x(end,fixed)-OCP.bounds.upper.finalState(fixed)')*lambda;        % Scalar
    Vx = ComplexDiff(@(X)OCP.cost.mayer(t(end),X),x(end,:))...
        + lambda'*constraintGrad;                                   % Row vector of lenth n
    Vlambda = (x(end,fixed)-OCP.bounds.upper.finalState(fixed)');         % Row vector of length p
    Vxx = OCP.hessian.mayer(t(end),x(end,:)'); % Hessian(@(X)OCP.cost.mayer(t(end),X),x(end,:)');           % nxn matrix
    Vxlambda = constraintGrad;                                      % pxn matrix
    Vlambda2 = zeros(p);                                            % pxp matrix
    
    %Integrate the Value functions backward
    Vf = [V;Vx(:);Vlambda(:);Vxx(:);Vlambda2(:);Vxlambda(:)];
    [tb,Vstate] = ode45(@valueFunction,linspace(t(end),0,nPoints),Vf,opt,OCP,@(T) interp1(t,x,T),@(T) interp1(t,u,T));
    
    %Parse Vstate
    [V,Vx,Vl,Vxx,Vll,Vxl] = parseValueStates(Vstate,OCP);
    
    %Compute the delta lambda quantity
    if ~softConstraints
        dlambda = -eta*inv(Vll{1})*Vlambda';
    else
        dlambda = zeros(p,1);
    end
    
    %Compute gains for the whole trajectory
    for i = nPoints:-1:1
        [H,L,F] = computePartials(OCP,x(i,:)',u(:,i));
        [l(1:m,i),Kx{i},Kl{i}] = computeGains(H,L,F,Vx{i},Vxx{i},Vxl{i});
    end
    
    %Integrate the linearize dynamics forward using new delta-control
    [~,deltaX] = ode45(@linearDynamics,t,zeros(n,1),opt,...
        dlambda,...
        @(X) ComplexDiff(@(Y)OCP.dynamics(0,Y(OCP.ind.state),Y(OCP.ind.control)),X),...
        @(T) interp1(t,x,T),...
        @(T) interp1(t,u,T),...
        @(T) interp1(t,l,T),...
        @(T) MatrixInterp(t,Kx,T),...
        @(T) MatrixInterp(t,Kl,T));
    
    %Compute delta-control
    for i = 1:nPoints
        du(1:m,i) = l(:,i)+Kx{i}*deltaX(i,:)' + Kl{i}*dlambda;
    end
    
    %Update control and adjoint
    u = u+gamma*du;
    lambda = lambda + dlambda;
    
    iter = iter+1; %Increment the counter and continue the loop!
    hist(iter).state = x;
    hist(iter).control = u;
    disp(['Iteration ', num2str(iter),' complete'])
end

sol.state = x;
sol.time = t;
sol.control = u;
sol.l = l;
sol.Kx = Kx;
sol.Kl = Kl;

end

function Vdot = valueFunction(t,State,ocp,x,u)

len = ocp.dimension;

%Parse states:
% V = State(1);
Vx = State(1+(1:len.state));
% Vl = State(1+len.state+(1:len.adjoint));
Vxx = reshape(State(1+len.state+len.adjoint+(1:len.state^2)),len.state,len.state);
% Vll = reshape(State(1+len.state+len.adjoint+len.state^2+(1:len.adjoint^2)),len.adjoint,len.adjoint);
Vxl = reshape(State(1+len.state+len.adjoint+len.state^2+len.adjoint^2+1:end),len.state,len.adjoint);

%Compute the Jacobian and Hessian terms:
F = Submatrix(ComplexDiff(@(X)ocp.dynamics(t,X(ocp.ind.state),X(ocp.ind.control)),[x(t)';u(t)]),ocp.dimension,ocp.ind);
JL = ComplexDiff(@(X)ocp.cost.lagrange(X(ocp.ind.state),X(ocp.ind.control)),[x(t)';u(t)]);
% L = Submatrix(JL,ocp.dimension,ocp.ind);
L.x = JL(ocp.ind.state);
L.u = JL(ocp.ind.control);
HL = ocp.hessian.lagrange(x,u); %Hessian(@(X)ocp.cost.lagrange(X(ocp.ind.state),X(ocp.ind.control)),[x(t)';u(t)]); %Hessian of the lagrange cost
% H = Submatrix(HL); %+ ocp.order*(SOMEHESSIANCRAP); %Only 0th order for now
H.xx = HL(ocp.ind.state,ocp.ind.state);
H.xu = HL(ocp.ind.state,ocp.ind.control);
H.ux = H.xu';
H.uu = HL(ocp.ind.control,ocp.ind.control);

%Compute the Gain terms:
Huui = inv(H.uu);
l = -Huui*(L.u+F.control'*Vx);
Kx = -Huui*(0.5*H.ux + 0.5*H.xu' + F.control'*Vxx); %This can be simplified
Kl = -Huui*F.control'*Vxl;

%The value function itself:
dV = 0.5*l'*H.uu*l - ocp.cost.lagrange(x(t)',u(t));
%First partial derivatives:
dVx = -L.x - (F.state*Vx - Kx'*H.uu*l)';
dVl = -Kl'*L.u;
%Second partial derivatives:
dVxx = -H.xx + Kx'*H.uu*Kx - Vxx*F.state' - F.state*Vxx;
dVll = Kl'*H.uu*Kl;
dVxl = -H.xu*Kl - F.state*Vxl - Vxx*F.control*Kl;

Vdot = [dV;dVx(:);dVl(:);reshape(dVxx,[],1);reshape(dVll,[],1);reshape(dVxl,[],1)];

end

function constraintGrad = Psi_X(fixed,p,n)
constraintGrad = zeros(p,n);
j = 0;
for i=1:p
    j = j+1;
    constraintGrad(j,fixed(i)) = 1;
end

end

function [l,Kx,Kl] = computeGains(H,L,F,Vx,Vxx,Vxl)
Huui = inv(H.uu);
l = -Huui*(L.u+Vx*F.control);
Kx = -Huui*(0.5*H.ux + 0.5*H.xu' + F.control'*Vxx); %This can be simplified
Kl = -Huui*F.control'*Vxl;

end

function [H,L,F] = computePartials(ocp,x,u)

F = Submatrix(ComplexDiff(@(X)ocp.dynamics(0,X(ocp.ind.state),X(ocp.ind.control)),[x;u]),ocp.dimension,ocp.ind);
JL = ComplexDiff(@(X)ocp.cost.lagrange(X(ocp.ind.state),X(ocp.ind.control)),[x;u]);
L.x = JL(ocp.ind.state);
L.u = JL(ocp.ind.control);
HL = ocp.hessian.lagrange(x,u); %Hessian(@(X)ocp.cost.lagrange(X(ocp.ind.state),X(ocp.ind.control)),[x;u]); %Hessian of the lagrange cost
% H = Submatrix(HL); %+ order*(SOMEHESSIANCRAP); %Only 0th order for now
H.xx = HL(ocp.ind.state,ocp.ind.state);
H.xu = HL(ocp.ind.state,ocp.ind.control);
H.ux = H.xu';
H.uu = HL(ocp.ind.control,ocp.ind.control);
end

function [V,Vx,Vl,Vxx,Vll,Vxl] = parseValueStates(VState,ocp)
len = ocp.dimension;
for i = 1:size(VState,1)
    State = VState(i,:);
    V(i) = State(1);
    Vx{i} = State(1+(1:len.state));
    Vl{i} = State(1+len.state+(1:len.adjoint));
    Vxx{i} = reshape(State(1+len.state+len.adjoint+(1:len.state^2)),len.state,len.state);
    Vll{i} = reshape(State(1+len.state+len.adjoint+len.state^2+(1:len.adjoint^2)),len.adjoint,len.adjoint);
    Vxl{i} = reshape(State(1+len.state+len.adjoint+len.state^2+len.adjoint^2+1:end),len.state,len.adjoint);
end
end

function dx = linearDynamics(t,x,lambda,jac,state,control,l,Kx,Kl)
J = jac([state(t)';control(t)]);
du = l(t) + Kx(t)*x + Kl(t)*lambda;
dx = J*[x;du]; %+HESSIAN TERM

end