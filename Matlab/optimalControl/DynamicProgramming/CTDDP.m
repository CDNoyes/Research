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
u = @(t) 0; %Initial control guess of all zeros. Implicitly assumes m = 1
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
    
    
    
    
end

end

function Vdot = valueFunction(t,V,T,X)

end

function constraintGrad = Psi_X(fixed,p,n)
constraintGrad = zeros(p,n);
j = 0;
for i=1:p
    j = j+1;
    constraintGrad(j,fixed(i)) = 1;
end

end