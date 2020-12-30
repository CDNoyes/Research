function [A,B] = EntryJacobians(x, u, discrete)
%%% Computes linearization matrices around a trajectory (x,u)

% Assumes x = [v, h, fpa, s] and u = cos(sigma)

% state and control indices
ix = 1:4;
iu = 5;

dv = x(1,2)-x(1,1);

% dynamics first derivatives
xu_dyn  = @(xu) entry_dynamics(xu(ix,:),xu(iu,:),dv,discrete);
J       = complex_difference(xu_dyn, [x; u]);
if discrete
    A      = J(2:end,2:ix(end),:);
    B      = J(2:end,iu,:);
else
    A      = J(1:end,1:ix(end),:);
    B      = J(1:end,iu,:);
end

function y = entry_dynamics(x,u,dv,discrete)

% === states and controls:
% x = [h dv gamma]'
% u = [u]'     = [cos(bank)]

% constants
rp = 3396.2e3;

% states
v = x(1,:);
h = x(2,:);
fpa = x(3,:);
% s = x(4,:);

[g,L,D] = entry_accels(x);

% Derivs
sdot = v.*cos(fpa)/1000; 
hdot = v.*sin(fpa);
vdot = -D-g.*sin(fpa);
fpadot = L./v.*u + (v./(rp+h) - g./v).*cos(fpa);

xdot = [vdot; hdot; fpadot; sdot];            % change in state

dt = dv./(vdot);   % just for estimate
if discrete
    y  = x + xdot.*dt;  % new state
else
    y = xdot;
end

function [g,L,D] = entry_accels(x)
% constants
cd  = 1.46;
cl  = 0.35;
rp = 3396.2e3;
m = 7200;
S = 15.8;
mu = 4.2830e13;

% states
h = x(2,:);
v = x(1,:);

% Accels
rho = 0.0158*exp(-h/9354.5);
f = 0.5*rho.*v.^2 * S/m;
D = f*cd;
L = f*cl;
g = mu./(rp+h).^2;

function J = complex_difference(fun, x)
% simple finite-difference derivatives
% assumes the function fun() is vectorized

h = 0.00001j;

[n, K]  = size(x);
H       = h*eye(n);
H       = permute(H, [1 3 2]);
X       = pp(x, H);
X       = reshape(X, n, K*n);
Y       = fun(X);
m       = numel(Y)/(K*n);
Y       = reshape(Y, m, K, n);
J       = imag(Y)/imag(h);
J       = permute(J, [1 3 2]);

% utility functions: singleton-expanded addition and multiplication
function c = pp(a,b)
c = bsxfun(@plus,a,b);