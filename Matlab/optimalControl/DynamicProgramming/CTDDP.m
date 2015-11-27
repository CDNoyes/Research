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
n = OCP.dimension.state;
m = OCP.dimension.control;
p = OCP.dimension.adjoint;

nPoints = 100;
t = linspace(0,tf,nPoints);
u = zeros(m,nPoints); %Initial control guess of all zeros.
x0 = OCP.state.initial;
lambda = zeros(p,1); %Initial adjoint guess of zeros
fixed = ~OCP.free.finalState; %Logical indices of the state vector elements that are fixed at tfinal, length p

constraintGrad = Psi_X(fixed,p,n);
if OCP.order == 0
    eta = .3; % A parameter to be set
    gamma = 0.1;
else
    eta = 1;
    gamma = 1;
end

tol = 1e-4;
iter = 0;
iterMax = 40;
constraintViolation = tol+1;
% The DDP loop

while constraintViolation > tol && iter < iterMax
    
    %Integrate dynamics forward with current control
    [t,x] = ode45(OCP.dynamics,t,x0,[], @(T) interp1(t,u,T));
    constraintViolation = max(abs(x(end,fixed)-OCP.bounds.finalState(fixed)'));
    
    %Compute V and its partials, eq. 28, at tFinal
    V = OCP.cost.mayer(t(end),x(end,:)) +...
        (x(end,fixed)-OCP.bounds.finalState(fixed)')*lambda;        % Scalar
    Vx = ComplexDiff(@(X)OCP.cost.mayer(t(end),X),x(end,:))...
        + constraintGrad'*lambda;                                   % Row vector of lenth n
    Vlambda = (x(end,fixed)-OCP.bounds.finalState(fixed)');         % Row vector of length p
    Vxx = Hessian(@(X)OCP.cost.mayer(t(end),X),x(end,:));           % nxn matrix
    Vxlambda = constraintGrad;                                      % pxn matrix
    Vlambda2 = zeros(p);                                            % pxp matrix
    
    %Integrate the Value functions backward
    Vf = [V;Vx(:);Vlambda(:);Vxx(:);Vxlambda(:);Vlambda2(:)];
    [tb,Vstate] = ode45(valueFunction,linspace(tf,0,nPoints),Vf,[]);
    
    %Parse Vstate
    [V,Vx,Vl,Vxx,Vll,Vxl] = parseValueStates(Vstate,ocp);
    
    %Compute the delta lambda quantity
    dlambda = -eta*inv(Vll{1})*Vlambda;
    
    %Compute gains for the whole trajectory
    for i = nPoints:-1:1
        [l{i},Kx{i},Kl{i}] = computeGains(H,L,F,Vx,Vxx,Vxl);
    end
    %Integrate the linearize dynamics forward using new delta-control
    
    %Computer delta-control
    
    %Update control and adjioint
    u = u+lambda*du;
    lambda = lambda + dlambda;
    iter = iter+1; %Increment the counter and continue the loop!
end

sol.state = x;
sol.time = t;
sol.control = u;
sol.l = l;
sol.Kx = Kx;
sol.Kl = Kl;

end

function Vdot = valueFunction(t,State,ocp,x)

len = ocp.dimension;

%Parse states:
% V = State(1);
Vx = State(1+(1:len.state));
% Vl = State(1+len.state+(1:len.adjoint));
Vxx = reshape(State(1+len.state+len.adjoint+(1:len.state^2)),len.state,len.state);
% Vll = reshape(State(1+len.state+len.adjoint+len.state^2+(1:len.adjoint^2)),len.adjoint,len.adjoint);
Vxl = reshape(State(1+len.state+len.adjoint+len.state^2+len.adjoint^2:end),len.state,len.adjoint);

%Compute the Jacobian and Hessian terms:
F = Submatrix(ComplexDiff(ocp.dynamics,x(t)));
JL = ComplexDiff(ocp.cost.lagrange,x(t));
L = Submatrix(JL);
HL = Hessian(ocp.cost.lagrange,x(t)); %Hessian of the lagrange cost
H = Submatrix(HL); %+ order*(SOMEHESSIANCRAP); %Only 0th order for now

%Compute the Gain terms:
Huui = inv(H.uu);
l = -Huui*(L.u+F.u*Vx);
Kx = -Huui*(0.5*H.ux + 0.5*H.xu' + F.u*Vxx); %This can be simplified
Kl = -Huui*Fu*Vxl;

%The value function itself:
dV = 0.5*l'*H.uu*l - L.matrix;
%First partial derivatives:
dVx = -L.x - F.x*Vx + Kx'*H.uu*l;
dVl = -Kl'*L.u;
%Second partial derivatives:
dVxx = -H.xx + Kx'*H.uu*Kx - Vxx*F.x' - F.x*Vxx;
dVll = Kl'*H.uu*Kl;
dVxl = -H.xu*Kl - Fx*Vxl - Vxx*Fu'*Kl;

Vdot = [dV;dVx;dVl;reshape(dVxx,[],1);reshape(dVll,[],1);reshape(dVxl,[],1)];

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
l = -Huui*(L.u+F.u*Vx);
Kx = -Huui*(0.5*H.ux + 0.5*H.xu' + F.u*Vxx); %This can be simplified
Kl = -Huui*Fu*Vxl;

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
    Vxl{i} = reshape(State(1+len.state+len.adjoint+len.state^2+len.adjoint^2:end),len.state,len.adjoint);
end
end