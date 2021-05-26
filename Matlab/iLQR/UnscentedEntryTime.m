function UnscentedEntryTime()
% X = [r,s,v,gamma]

global n_samples weights

x0 = [120e3, 0, 5505, -15.5*pi/180, 1, 1, 1]';
P0 = diag([0.001, 8000, 1, 0.2*pi/180, 0.05,0.05,0.07]).^2;
[X0,weights] = UnscentedTransform(x0,P0,15);
n_samples = length(weights);

%%% Integrates longitudinal equations of motion wrt velocity using
%%% feedforward control u and feedback gains K
Xvec = X0';
Xvec = Xvec(:);
u = 0; % bank = 90
fun = @(t,x) entry_dynamics(t,x,u);
nsteps = 500;
[t,X] = ode45(fun, linspace(0,100,nsteps), Xvec);

X3d = reshape(X, nsteps, 15,[]);
% X3d(1,:,1) % should be altitudes

X3d(:,:,5:7) = []; % Remove the parameters
X3d(:,:,3) = []; % Remove velocity state

[h,v,fpa,s,g,L,D] = entry_accels(X');
vdot = -D - g.*sin(fpa);

% Find a velocity for which vdot < threshold < 0 for all of the
% samples.
    % First, cut off the increasing velocity portion of the trajectory 
decrease = all(vdot <= -0.2, 1);
min(t(decrease))
% td = t(all(decrease, 1)); % Finds all times at which all points are decreasing
% td = td(1); % the first time at which all the points are decreasing
% Now take the lowest velocity
vd = v(:,decrease); % velocities at the first time that all are sufficiently decreasing
v0 = min(vd(:,1)) - 1;
v0 = 5500;% hardcoded works as long as each traj reached it and its sufficiently low

% Then find the state of each sigma point at that velocity
for i = 1:n_samples
    data = squeeze(X3d(decrease,i,:));
    X_v0(i,:) = interp1(v(i,decrease), data, v0);
end
% Finally, compute the UT stats for those states
[hm, hs] = get_stats_(X_v0(:,1)/1000);
[fm, fs] = get_stats_(X_v0(:,3)*180/pi);
[sm, ss] = get_stats_(X_v0(:,2)/1000);

disp(['v0 = ', num2str(v0)])
disp(['h0 = ', num2str(hm), ' km, sigma = ', num2str(hs)])
disp(['s0 = ', num2str(sm), ' km, sigma = ', num2str(ss)])
disp(['fpa0 = ', num2str(fm), ' km, sigma = ', num2str(fs)])

function [h,v,fpa,s,g,L,D] = entry_accels(x)
% constants
rp = 3396.2e3;
[m, S, cl, cd] = aero_const();
mu = 4.2830e13;

% states
[h, v, fpa, s, fcl, fcd, frho] = get_states(x);

% Accels
rho = frho.*0.0158.*exp(-h./9354.5);
f = 0.5*rho.*v.^2*S/m;
D = f.*cd.*fcd;
L = f.*cl.*fcl;
g = mu./(rp+h).^2;

function xdot = entry_dynamics(t, x, u)

% === states and controls:
% x = [h, s, v, gamma]'
% u = [u]'     = [cos(bank)]

% constants
rp = 3396.2e3;

% states + current accelerations
[h, v, fpa, s, g, L, D] = entry_accels(x);


% Derivs
sdot = v.*cos(fpa);
hdot = v.*sin(fpa);
vdot = -D-g.*sin(fpa);
fpadot = L./v.*u + (v./(rp+h) - g./v).*cos(fpa);


xdot = [hdot; sdot; vdot; fpadot; 0*h; 0*h; 0*h];  % change in state


function [h, v, fpa, s, fcl, fcd,frho] = get_states(x)
global n_samples

ih = 1:n_samples;
is = n_samples + ih;
iv = n_samples + is;
ig = n_samples + iv;

h = x(ih,:); % m
v = x(iv,:);
fpa = x(ig,:);
s = x(is,:); % m
icl = n_samples + ig;
icd = n_samples + icl;
irho = n_samples + icd;
fcl = x(icl,:);
fcd = x(icd,:);
frho = x(irho,:);


function [m,s] = get_stats_(x)
global weights
y = x.*weights;
m = sum(y, 1);
if nargout > 1
    s = sum(weights.*(x-m).^2, 1).^0.5;
end