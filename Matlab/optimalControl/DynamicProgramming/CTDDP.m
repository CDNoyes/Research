%CTDDP Continous-time differential dynamic programming implementation.
%   CTDDP(OptimalControlProblem) Solves an OCP via an implementation of the
%   algorithm described in "Continous-time DDP with Terminal Constraints"
%   by Sun et al.

function sol = CTDDP(OCP)

if ~OCP.fixedFinalTime
    error('This continous time DDP formulation currently requires a fixed final time.')
end

n = OCP.dimension.state;
m = OCP.dimension.control;
p = OCP.dimension.adjoint;

nPoints = 100;
u = @(t) zeros(m,1); %Initial control guess of all zeros.
x0 = OCP.state.initial;
lambda = zeros(p,1); %Initial adjoint guess of zeros
fixed = ~OCP.free.finalState; %Logical indices of the state vector elements that are fixed at tfinal, length p

constraintGrad = Psi_X(fixed,p,n);

% The DDP loop
while max(abs(something)) > tol
    
    %Integrate dynamics forward with current control
    [t,x] = ode45(OCP.dynamics,linspace(0,tf,nPoints),x0,[], u);
    
    %Compute V and its partials, eq. 28, at tFinal
    V = OCP.cost.mayer(t(end),x(end,:)) +...
        (x(end,fixed)-OCP.bounds.finalState(fixed)')*lambda;    % Scalar
    Vx = ComplexDiff(@(TX)OCP.cost.mayer(TX(1),TX(2:end)),[t(end),x(end,:)])...
        + constraintGrad*lambda;                                 % Row vector of lenth n
    Vlambda = (x(end,fixed)-OCP.bounds.finalState(fixed)');     % Row vector of length p
    Vxx = Hessian(@(X)OCP.cost.mayer(t(end),X),x(end,:));       % nxn matrix
    Vxlambda = constraintGrad;                                    % pxn matrix
    Vlambda2 = zeros(p);                                        % pxp matrix
    
    %Integrate the Value functions backward
    Vf = [V;Vx(:);Vlambda(:);Vxx(:);Vxlambda(:);Vlambda2(:)];
    [t,Vstate] = ode45(valueFunction,linspace(tf,0,nPoints),Vf,[]);
    
    
    
end

end

function Vdot = valueFunction(t,State,len)

%Parse states:
% V = State(1);
Vx = State(1+(1:len.state));
% Vl = State(1+len.state+(1:len.adjoint));
Vxx = reshape(State(1+len.state+len.adjoint+(1:len.state^2)),len.state,len.state);
% Vll = reshape(State(1+len.state+len.adjoint+len.state^2+(1:len.adjoint^2)),len.adjoint,len.adjoint);
Vxl = reshape(State(1+len.state+len.adjoint+len.state^2+len.adjoint^2:end),len.state,len.adjoint);

%The value function itself:
dV = 0.5*l'*Huu*l - L;
%First partial derivatives:
dVx = -Lx - Fx*Vx + Kx'*Huu*l;
dVl = -Kl'*Lu;
%Second partial derivatives:
dVxx = -Hxx + Kx'*Huu*Kx - Vxx*Fx' - Fx*Vxx;
dVll = Kl'*Huu*Kl;
dVxl = -Hxu*Kl - Fx*Vxl - Vxx*Fu'*Kl;

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