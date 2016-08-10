% We have to be very particular about the way the inputs are treated. Since
% the function may have state, control, and parameters for which we wish to
% compute sensitivities via STM, the function F should take only one input
% vector X = [x,u,p] and should have trivial dynamics for the constant
% parameters P. Perhaps we can write a wrapper that handles the input
% conversion and trivial dynamics for us. Some thought is required. 

%After thinking about it, the incoming f should take only one argument, the
%user should augment the dynamics with trivial dynamics as needed for
%parameters. This should also be reflecteed in the jac and hess inputs, if
%the user provides them. This way if the algorithm computes them
%numerically everything will work out fine.

%Essentially, I am requiring that the jacobian of the system be square.
function [fAugmented,xAugmented] = AppendSTM(f,x0,order,ind,jac,hess)

%Order 1 or 2, the order of the state transition matrix stuff to include
% nInputs = nargin(f);

n = length(x0);
xAugmented = [x0(:);reshape(eye(n),[],1)];
if nargin < 5 || isempty(jac)
    jac = @(x,u) ComplexDiff(@(X)f(0,X(ind.state),u(ind.control)),[x;u]);
end

fAugmented = @(t,x,u) [f(t,x,u); stmDynamics(x,u,length(ind.state)+length(ind.control),jac)];

end

function dPhi = stmDynamics(x,u,n,J)
Phi = reshape(x(n+1:end),n,n);
dPhi = reshape(J(x,u)*Phi,[],1);
end