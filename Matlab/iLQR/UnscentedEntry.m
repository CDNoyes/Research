function [V,X] = UnscentedEntry(X0, u, K)

global n_samples
n_samples = length(X0)/3;

%%% Integrates longitudinal equations of motion wrt velocity using
%%% feedforward control u and feedback gains K
fun = @(v,x) entry_dynamics(v,x,u,K);
[V,X] = ode45(fun, [5461.4, 550], X0);


function [h,fpa,s,g,L,D] = entry_accels(v,x)
% constants
cd  = 1.46;
cl  = 0.35;
rp = 3396.2e3;
m = 7200;
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

function dx = entry_dynamics(v, x, u, K)

% === states and controls:
% x = [h gamma s]'
% u = [u]'     = [cos(bank)]

% constants
rp = 3396.2e3;

% states + current accelerations
[h, fpa, s, g, L, D] = entry_accels(v,x);

% feedback terms
if 1
    % Tuned for decent performance in the guess trajectory
    kd = 0.1; %more lift up when too much drag
    ks = -0.1;  % less lift up when too close
    kf = -50.0;% left lift up when to shallow

    ef = get_error(fpa);
    es = get_error(s);
    eD = get_error(D);
    du = kd*eD + ks*es + kf*ef;
    u_cl = Saturate(u(v)+du, 0, 1);
else
    u_cl = u(v);
end

% Derivs
sdot = v.*cos(fpa)/1000; % in km
hdot = v.*sin(fpa);
vdot = -D-g.*sin(fpa);
fpadot = L./v.*u_cl + (v./(rp+h) - g./v).*cos(fpa);

xdot = [hdot; fpadot; sdot];            % change in state

dt = 1./vdot;   % just for estimate
dx =  xdot.*repmat(dt, 3, 1);  % new state


function [h,fpa,s] = get_states(x)
global n_samples


ih = 1:n_samples;
ig = n_samples + ih;
is = n_samples + ig;

h = x(ih,:); % m
fpa = x(ig,:);
s = x(is,:); % km

function er = get_error(state)
er = [state(1,:)*0; state(2:end,:)-state(1,:)]; % subtracts the 'prime' point from the rest
% in reality this should be the mean subtracted from all the points 