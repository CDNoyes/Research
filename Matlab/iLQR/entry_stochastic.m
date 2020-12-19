function sol = entry_stochastic(closed_loop, W, bounds)
% A demo of iLQG/DDP with stochastic entry dynamics
% clc;
% close all

fprintf(['\nUse of the iLQG algorithm '...
    'with entry dynamics, velocity loss as independent variable.\n'])

% Set full_DDP=true to compute 2nd order derivatives of the
% dynamics. This will make iterations more expensive, but
% final convergence will be much faster (quadratic)

global ds v0 full_DDP n_samples n_states scale weights closed wh ws wu
% note: weights are the sigma point weights
% w_ are the cost weights

n_states = 3;
full_DDP = 1;
if nargin == 0
    closed = 1;
    wh = 2;
    ws = 1;
    wu = 0.2;
    bounds = [0 1];
else
    closed = closed_loop;
    wh = W(1);
    ws = W(2);
    wu = W(3);
    full_DDP = 1; %~closed;
end

% Scale states for numerical condition
hscale = 120e3;
vscale = 5500;
fpascale = 0.2;
rangescale = 500e3;
scale = [hscale, vscale, fpascale, rangescale];
% v0 = 5461.4;
v0 = 5525;

% set up the optimization problem
DYNCST  = @(x,u,i) entry_dyn_cst(x,u,full_DDP);
V       = 600; % terminal velocity, m/s
dV      = v0-V;
T       = 1000;              % horizon
ds      = dV/T;

% Initial States
% x0      = [0, 39.4497e3/hscale, (-10.604*pi/180)/fpascale, 0]';   % nominal initial state
x0 = [0 , 51.5356e3/hscale, (-11.2179428280964*pi/180)/fpascale, 0]';
% n_samples = 3;
% delta = zeros(3, n_samples);
% delta(2,:) = [0, -1, 1];
% X0 = [x0(2:end),x0(2:end),x0(2:end)]' + delta'*0.75*pi/180/fpascale;
% weights = ones(n_samples,1)/n_samples;
% X0_vector = [0;X0(:)];

[X0, weights] = UnscentedTransform(x0(2:end), diag([2500/hscale, 0.25*pi/180/fpascale, 5e3/rangescale]).^2, 4);
n_samples = length(weights);
X0 = [zeros(n_samples,1), X0'];
X0_vector = X0(n_samples:end)';

fprintf(['Samples = ',num2str(n_samples),'\nTotal State Dimension = ',num2str(1+n_samples*n_states),'\n'])

% Initial Control
% u0      = ones(1,floor(T))*0.4;
% u0 =     [zeros(1, floor(1600/ds)), ones(1,T-floor((1600)/ds))]; % good guess
u0 = linspace(0.35, 1, T);
% u0 = linspace(1, 0, T);


Op.lims  = bounds;
Op.plot = 1;               % plot the derivatives as well
Op.maxIter = 50;
Op.parallel = 0;


% === run the optimization!
[x, u, ~, ~, ~, cost] = iLQG(DYNCST, X0_vector, u0, Op);

% Compute stats and setup the output structure 
[h,v,fpa,s,~,L,D] = entry_accels(x);
v=v';
s = s/1000;
[hm, hv] = get_stats(h);
[sm, sv] = get_stats(s);
[fm, fv] = get_stats(fpa);
[Dm, Dv] = get_stats(D);

sol.u = u;
sol.v = v;
sol.mean = [hm; fm; sm];
sol.var = [hv; fv; sv];
sol.h = h;
sol.fpa = fpa;
sol.s = s;
sol.weights = [wh,ws,wu];
sol.cost = sum(cost);
sol.bounds = bounds;
sol.X0 = [h(:,1)',fpa(:,1)',s(:,1)'];
sol.ddp = full_DDP;
sol.sigma_weights = weights;
sol.L = L;
sol.D = D;
sol.Dm = Dm;
sol.Dv = Dv;


print_stats(x(:,end));
lw = 2;


figure
x2 = [v', fliplr(v')];
inBetween = [(sm-3*sv.^0.5), fliplr((sm+3*sv.^0.5))];
fill(x2, inBetween, 'c');
hold all
plot(v, s', 'k--','linewidth', 1)
plot(v, sm, 'm', 'linewidth', lw)

xlabel('Velocity, m/s')
ylabel('Downrange, km')
grid on

figure
x2 = [v', fliplr(v')];
inBetween = [(hm-3*hv.^0.5)/1000, fliplr((hm+3*hv.^0.5)/1000)];
fill(x2, inBetween, 'c');
hold all
plot(v, h'/1000, 'k--','linewidth', 1)
plot(v, hm/1000, 'm', 'linewidth', lw)

xlabel('Velocity m/s')
ylabel('Altitude km')
grid on

figure
inBetween = [(Dm-3*Dv.^0.5), fliplr((Dm+3*Dv.^0.5))];
fill(x2, inBetween, 'c');
hold all
plot(v, D', 'k--','linewidth', 1)
plot(v, Dm, 'm', 'linewidth', lw)
xlabel('Velocity, m/s')
ylabel('Drag, m/s^2')
grid on

function [h,v,fpa,s,g,L,D] = entry_accels(x)
% constants
cd  = 1.46;
cl  = 0.35;
rp = 3396.2e3;
m = 7200;
S = 15.8;
mu = 4.2830e13;
rho0 = 0.0158;

% states
[h,v,fpa,s] = get_states(x);

% Accels
rho = rho0*exp(-h/9354.5);
f = 0.5*rho.*v.^2 * S/m;
D = f*cd;
L = f*cl;
g = mu./(rp+h).^2;

function y = entry_dynamics(x,u)
global ds full_DDP scale n_states closed

% === states and controls:
% x = [h dv gamma]'
% u = [u]'     = [cos(bank)]

% constants
rp = 3396.2e3;

% states + current accelerations
[h,v,fpa,s,g,L,D] = entry_accels(x);

% feedback terms
if closed
    % Tuned for decent performance in the guess trajectory
    kd = 1*0.1; %more lift up when too much drag
    ks = 1*-0.1;  % less lift up when too close
    kf = 1*-50.0;% left lift up when to shallow
    
    ef = get_error(fpa);
    es = get_error(s/1000);
    eD = get_error(D);
    du = kd*eD + ks*es + kf*ef;
            u_cl = smooth_sat(u + du); % cosh
%     u_cl = smooth_sat2(u + du); % sqrt
    
else
    u_cl = u;
end

% Derivs
sdot = v.*cos(fpa); % in km
hdot = v.*sin(fpa);
vdot = -D-g.*sin(fpa);
% fpadot = L./v.*u_cl.^(1+full_DDP) + (v./(rp+h) - g./v).*cos(fpa);
fpadot = L./v.*u_cl + (v./(rp+h) - g./v).*cos(fpa);

xdot = [hdot/scale(1); fpadot/scale(3); sdot/scale(4)];            % change in state

dt = -ds./vdot;   % just for estimate
y  = x + [ds/scale(2)+v*0; xdot.*repmat(dt, n_states, 1)];  % new state

function y = smooth_sat(x)
K = 20;
q = 2*x - 1;
y =  0.5/K * log(cosh(K*(q+1))./cosh(K*(q-1))); % saturates between [-1,1]
y = 0.5 + 0.5*y; % map back to [0,1]

function y = smooth_sat2(x)
k = 1e-3;
q = x-0.5;
c = 0.5;
y = c + c*0.5*(sqrt(k + (q/c+1).^2) - sqrt(k + (q/c-1).^2));

function er = get_error(state)
% er = [state(1,:)*0; state(2:end,:)-state(1,:)]; % subtracts the 'prime' point from the rest
er = state - get_stats(state);
% in reality this should be the mean subtracted from all the points



function c = entry_cost(x, u)
global ds wh ws wu

cost_scale = 10000;

% states
[h,v,fpa,s,g,L,D] = entry_accels(x);

% cost function
% sum of 3 terms:
% lu: quadratic cost on controls
% lf: final cost
% lx: running cost on altitude loss

final = isnan(u(1,:));
u(:,final)  = 0;


% running cost terms
hdot = v.*sin(fpa);
vdot = -D - g.*sin(fpa);
hprime = hdot./vdot * ds; %it's actually -hdot over -vdot
hprime = get_stats(hprime); % expected value of running cost

sdot = v.*cos(fpa); % meters/second
[smean, svar] = get_stats(s);
[hmean, hvar] = get_stats(h);

vs_dot = -2*get_stats(s.*sdot./vdot) + 2*smean.*get_stats(sdot./vdot);
vh_dot = -2*get_stats(h.*hdot./vdot) + 2*hmean.*get_stats(hdot./vdot);

%convert variance rates to std rates
vs_dot = 0.5*vs_dot./(svar.^0.5);
vh_dot = 0.5*vh_dot./(hvar.^0.5);

lx = hprime + ws*vs_dot*ds + wh*vh_dot*ds;

% control cost
lu    = wu*(u-0.5).^2 * ds.* v/2500; % this term is for smoothness

% final cost
if any(final)
    llf = 0*s;
    %        llf      = ws*1000*svar.^0.5 + wh*hvar.^0.5;
    lf = zeros(size(u));
    lf(1,final)= llf(1,final);
else
    lf    = 0;
end

% total cost
c     = (lu + lx + lf)/cost_scale;


function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = entry_dyn_cst(x,u,full_DDP)
global n_states n_samples
% combine dynamics and cost
% use helper function finite_difference() to compute derivatives

if nargout == 2
    f = entry_dynamics(x,u);
    c = entry_cost(x,u);
else
    % state and control indices
    ix = 1:(1+n_states*n_samples);
    iu = ix(end)+1;
    
    % dynamics first derivatives
    xu_dyn  = @(xu) entry_dynamics(xu(ix,:),xu(iu,:));
    J       = complex_difference(xu_dyn, [x; u]);
    fx      = J(:,ix,:);
    fu      = J(:,iu,:);
    
%     J2 = EntryJacobian(x, u, 0);
    
    % dynamics second derivatives
    if full_DDP
        N_J = size(J);
        xu_Jcst = @(xu) complex_difference(xu_dyn, xu);
        JJ      = finite_difference(xu_Jcst, [x; u]);
        %         JJ      = reshape(JJ, [4 6 size(J)]); % Original code
        if N_J <= 2
            JJ = reshape(JJ,[ix(end) iu N_J(2)]);
        else
            JJ = reshape(JJ, [ix(end) iu N_J(2) N_J(3)]);
        end
        JJ      = 0.5*(JJ + permute(JJ,[1 3 2 4])); %symmetrize
        fxx     = JJ(:,ix,ix,:);
        fxu     = JJ(:,ix,iu,:);
        fuu     = JJ(:,iu,iu,:);
    else
        [fxx,fxu,fuu] = deal([]);
    end
    
    % cost first derivatives
    xu_cost = @(xu) entry_cost(xu(ix,:),xu(iu,:));
    J       = squeeze(complex_difference(xu_cost, [x; u]));
    cx      = J(ix,:);
    cu      = J(iu,:);
    
    %     cost second derivatives
    xu_Jcst = @(xu) squeeze(complex_difference(xu_cost, xu));
    JJ      = finite_difference(xu_Jcst, [x; u]);
    JJ      = 0.5*(JJ + permute(JJ,[2 1 3])); %symmetrize
    cxx     = JJ(ix,ix,:);
    cxu     = JJ(ix,iu,:);
    cuu     = JJ(iu,iu,:);
    
    [f,c] = deal([]);
end


function J = EntryJacobian(x, u, dv)
%%% Jacobian of longitudinal entry dynamics [h,fpa,s] written wrt to velocity
%%% as the independent variable. Also includes control jacobian. 

% constants
hs = 9354.5;
cd  = 1.46;
cl  = 0.35;
rp = 3396.2e3;
m = 7200;
S = 15.8;
mu = 4.2830e13;
rho0 = 0.0158;

[h,v,fpa,s,g,L,D] = entry_accels(x);
r = h + rp;

% velocity rate  = 1 wrt velocity 
J0 = 0*h;

% altitude rate terms
J11 = -(v.*sin(fpa).*((2*mu*sin(fpa))./r.^3 + (S*cd*rho0*v.^2.*exp(-h./hs))./(2*hs*m)))./((mu.*sin(fpa))./(h + rp).^2 + (S*cd*rho0.*v.^2.*exp(-h./hs))./(2*m)).^2;
J12 = (mu.*v.*cos(fpa).*sin(fpa))./(r.^2.*((mu.*sin(fpa))./(r).^2 + (S*cd*rho0*v.^2.*exp(-h./hs))./(2*m)).^2) - (v.*cos(fpa))./((mu.*sin(fpa))./(r).^2 + (S*cd*rho0.*v.^2.*exp(-h./hs))./(2*m));
J13 = 0*h;
J14 = 0*h;

J = [J11; J12; J13; J14];


function J = finite_difference(fun, x, h)
% simple finite-difference derivatives
% assumes the function fun() is vectorized

if nargin < 3
    h = 2^-17;
end

[n, K]  = size(x);
H       = [zeros(n,1) h*eye(n)];
H       = permute(H, [1 3 2]);
X       = pp(x, H);
X       = reshape(X, n, K*(n+1));
Y       = fun(X);
m       = numel(Y)/(K*(n+1));
Y       = reshape(Y, m, K, n+1);
J       = pp(Y(:,:,2:end), -Y(:,:,1)) / h;
J       = permute(J, [1 3 2]);

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


function [h,v,fpa,s] = get_states(x)
global scale n_samples v0

iv = 1;
ih = iv + 1:n_samples+1;
ig = n_samples + ih;
is = n_samples + ig;

h = x(ih,:)*scale(1); % m
v = v0 - x(iv,:)*scale(2); % m/s
fpa = x(ig,:)*scale(3); % rad
s = x(is,:)*scale(4); % m


function [m,v] = get_stats(x)
global weights
y = x.*weights;
m = sum(y, 1);
if nargout > 1
    v = sum(weights.*(x-m).^2, 1);
end


function print_stats(x)
[h,~,~,s] = get_states(x);

[hm, hv] = get_stats(h/1000);
[sm, sv] = get_stats(s/1000);

disp(['hf = ',num2str(hm),' +/- 3*',num2str(hv^0.5),' km'])
disp(['sf = ',num2str(sm),' +/- 3*', num2str(sv^0.5),' km'])