function [V,X,U] = UnscentedEntry(V, X0, u, K, W, Xr)
%%% Xr should be interpolant of [D, fpa, s]
% W is the UT sigma point weights to compute mean/var

global n_samples weights scale
n_samples = length(X0)/6;
weights = W;
hscale = 120e3;
vscale = 5500;
fpascale = 0.2;
rangescale = 500e3;
scale = [hscale, vscale, fpascale, rangescale];

%%% Integrates longitudinal equations of motion wrt velocity using
%%% feedforward control u and feedback gains K
fun = @(v,x) entry_dynamics(v,x,u,K,Xr);
[V,X] = ode45(fun, [V(1), V(end)], X0);

for i = 1:length(V)
    [~, U(:,i)] = entry_dynamics(V(i), X(i,:)', u, K, Xr);
end

% Unscale
[h, fpa, s, ~,~,~] = get_states(X');
X = [h; fpa; s/1000]';

function [h,fpa,s,g,L,D] = entry_accels(v,x)
% constants
rp = 3396.2e3;
[m, S, cl, cd] = aero_const();
mu = 4.2830e13;

% states
[h, fpa, s, fcl, fcd, frho] = get_states(x);

% Accels
rho = frho.*0.0158.*exp(-h./9354.5);
f = 0.5*rho.*v.^2*S/m;
D = f.*cd.*fcd;
L = f.*cl.*fcl;
g = mu./(rp+h).^2;

function [dx, uout] = entry_dynamics(v, x, u, K, Xr)
global scale

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
es = s/1000 - xr(3);

du = kd*eD + ks*es + kf*ef;
% u_cl = Saturate(u(v)+du, 0, 1);
LoD = L./D;
u_cl = smooth_sat(u(v)./LoD + du);


% Derivs
sdot = v.*cos(fpa);
hdot = v.*sin(fpa);
vdot = -D-g.*sin(fpa);
fpadot = L./v.*u_cl + (v./(rp+h) - g./v).*cos(fpa);


xdot = [hdot/scale(1); fpadot/scale(3); sdot/scale(4)];            % change in state

dt = 1./vdot;   % just for estimate
dx =  xdot.*repmat(dt, 3, 1);  % new state
dx = [dx; dx*0]; % append the zero dynamics for parameters 
if nargout == 2
    uout = u_cl;
end


function [h,fpa,s, fcl, fcd,frho] = get_states(x)
global scale n_samples

ih = 1:n_samples;
ig = n_samples + ih;
is = n_samples + ig;

h = x(ih,:)*scale(1); % m
fpa = x(ig,:)*scale(3);
s = x(is,:)*scale(4); % m
icl = n_samples + is;
icd = n_samples + icl;
irho = n_samples + icd;
fcl = x(icl,:);
fcd = x(icd,:);
frho = x(irho,:);



function y = smooth_sat(x)
K = 20;
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