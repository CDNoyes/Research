% ADAPTIVESMCDYNAMICS Implements the dynamics from 
% "Adaptive Sliding-Mode Control for Nonlinear Systems With Uncertain
% Parameters" (IEEE Transactions on System, Man, and Cybernetics, April
% 2008)

function dX = AdaptiveSMCDynamics(t,X)
w = 0; %Noise term

xd = [sin(t);cos(t);-sin(t)]; %The reference trajectory to be tracked

x = X(1:3);     % Actual, corrupted dynamics
x0 = X(4:6);    % Our system model of the dynamics
gain = X(7);    % An adaptive gain

% Define the sliding surface
c = [12,7,1]';  % The sliding mode coefficients
e = x-xd;
sig = c'*e;

% Nominal system dynamics: f0(x0) + b0(x0)*u
b0 = @(x) [0;0;3];
f0 = @(y) [y(2:3);-1*y(1)^2-(1.5)*y(2)-y(3)];

% Compute the controls
u0 = -(c'*b0(x0))\(c'*f0(x0)-c'*[xd(2:3);-cos(t)]);
uas = -(c'*b0(x0))\gain*tanh(10*sig);
u = u0+uas;

% Compute the dynamics
dx = [x(2:3);
    -(1+0.3*sin(t))*x(1)^2-(1.5+0.2*cos(t))*x(2)-(1+0.4*sin(t))*x(3)+(3+cos(t))*u+w];
dx0 = f0(x0) + b0(x0)*u;
dgain = 5*abs(sig);

dX = [dx;dx0;dgain];
end