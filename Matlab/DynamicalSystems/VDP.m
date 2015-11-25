%VDP Creates the system dynamics of a Van Der Pol oscillator.
%   VDP(MU) returns the system dynamics of an oscillator with nonlinear
%   parameter MU. Also returns the system's Jacobian matrix and Hessian
%   tensor. The dynamics output is a function handle that can be integrated
%   via Matlab's builtin integrators.
%
% dynamics: f(time,state,control)
% jacobian: j(state,control)

function [fVDP,jVDP,hVDP] = VDP()

fVDP = @vdpDynamics;
jVDP = @vdpJacobian;
hVDP = @vdpHessian;

end

function dx = vdpDynamics(t,x,u,mu)

dx = [ x(2)
        -x(1) + mu*(1-x(1).^2)*x(2) + ControlWrapper(t,u) ];

end

function J = vdpJacobian(x,u,mu)

J = [0, 1, 0, 0
    -1-2*mu*x(1)*x(2), mu*(1-x(1).^2), 1,(1-x(1).^2)*x(2)];

end

function [H,h] = vdpHessian(x,u,mu)
h{2} = [-2*mu*x(2),-2*mu*x(1),0,-2*x(1)*x(2);
        -2*mu*x(1),0,0,(1-x(1).^2)
        0,0,0,0
        -2*x(1)*x(2),(1-x(1).^2),0,0];
h{1} = zeros(4,4);
H = sparse(blkdiag(h{:}));
end