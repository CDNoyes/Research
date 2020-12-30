function [V,X,U] = UnscentedEntry(V, X0, u, K, W, Xr)
%%% Xr should be interpolant of [D, fpa, s]
global n_samples weights
n_samples = length(X0)/3;
weights = W;

%%% Integrates longitudinal equations of motion wrt velocity using
%%% feedforward control u and feedback gains K
fun = @(v,x) entry_dynamics(v,x,u,K,Xr);
[V,X] = ode45(fun, [V(1), V(end)], X0);

for i = 1:length(V)
    [~, U(:,i)] = entry_dynamics(V(i), X(i,:)', u, K, Xr);
end

function [h,fpa,s,g,L,D] = entry_accels(v,x)
% constants
cd  = 1.408;
cl  = 0.357;
rp = 3396.2e3;
m = 5000;
S = 15.8;
mu = 4.2830e13;

% states
[h,fpa,s] = get_states(x);

% Accels
rho = 0.0158*exp(-h/9354.5);
f = 0.5*rho.*v.^2 * S/m;
D = f*cd;
L = f*cl;
g = mu./(rp+h).^2;

function [dx, uout] = entry_dynamics(v, x, u, K, Xr)

% === states and controls:
% x = [h gamma s]'
% u = [u]'     = [cos(bank)]

% constants
rp = 3396.2e3;

% states + current accelerations
[h, fpa, s, g, L, D] = entry_accels(v,x);

% feedback terms
k = K(v);
kd = k(1);
ks = k(2);
kf = k(3);

% Ref, aka mean from the optimization routine 
xr = Xr(v);
eD = D - xr(1);
ef = fpa-xr(2);
es = s - xr(3);

du = kd*eD + ks*es + kf*ef;
% u_cl = Saturate(u(v)+du, 0, 1);
u_cl = smooth_sat(u(v)+du);


% Derivs
sdot = v.*cos(fpa)/1000; % in km
hdot = v.*sin(fpa);
vdot = -D-g.*sin(fpa);
fpadot = L./v.*u_cl + (v./(rp+h) - g./v).*cos(fpa);

xdot = [hdot; fpadot; sdot];            % change in state

dt = 1./vdot;   % just for estimate
dx =  xdot.*repmat(dt, 3, 1);  % new state
if nargout == 2
    uout = u_cl;
end


function [h,fpa,s] = get_states(x)
global n_samples


ih = 1:n_samples;
ig = n_samples + ih;
is = n_samples + ig;

h = x(ih,:); % m
fpa = x(ig,:);
s = x(is,:); % km

% function er = get_error(state, mean)
% er = state - get_stats(state);


function y = smooth_sat(x)
K = 5;
q = 2*x - 1;
y =  0.5/K * log(cosh(K*(q+1))./cosh(K*(q-1))); % saturates between [-1,1]
y = 0.5 + 0.5*y; % map back to [0,1]

function [m,v] = get_stats(x)
global weights
y = x.*weights;
m = sum(y, 1);
if nargout > 1
    v = sum(weights.*(x-m).^2, 1);
end