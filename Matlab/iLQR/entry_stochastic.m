function entry_stochastic
% A demo of iLQG/DDP with stochastic entry dynamics
clc; clear;
close all

fprintf(['\nUse of the iLQG algorithm '...
'with entry dynamics, velocity loss as independent variable.\n'...
'for details see\nTassa, Mansard & Todorov, ICRA 2014\n'...
'\"Control-Limited Differential Dynamic Programming\"\n'])

% Set full_DDP=true to compute 2nd order derivatives of the 
% dynamics. This will make iterations more expensive, but 
% final convergence will be much faster (quadratic)
 
global ds v0 full_DDP n_samples n_states scale weights
n_states = 3;
full_DDP = 0;

hscale = 120e3;
vscale = 5500;
fpascale = 0.2;
v0 = 5461.4;
rangescale = 500;
scale = [hscale, vscale, fpascale, rangescale];

% set up the optimization problem
DYNCST  = @(x,u,i) entry_dyn_cst(x,u,full_DDP);
V       = 550; % terminal velocity, m/s
dV      = v0-V; 
T       = 1000;              % horizon
ds      = dV/T;

% Initial States 
x0      = [0, 39.4497e3/hscale, (-10.604*pi/180)/fpascale, 0]';   % nominal initial state
% n_samples = 5;
% delta = zeros(4,n_samples);
% delta(3,:) = [0, -1, 1, -2, 2];
% X0 = [x0,x0,x0,x0,x0]' + delta'*0.75*pi/180/fpascale;
% weights = ones(n_samples,1)/n_samples;

[X0, weights] = UnscentedTransform(x0(2:end), diag([2500/hscale, 0.25*pi/180/fpascale, 1/rangescale]).^2, 3);
n_samples = length(weights);
X0 = [zeros(n_samples,1), X0'];
X0_vector = X0(n_samples:end)';

fprintf(['Samples = ',num2str(n_samples),'\nTotal State Dimension = ',num2str(1+n_samples*n_states),'\n'])
% Initial Control 
% u0      = ones(1,floor(T));    
% u0 =     [zeros(1, floor(1600/ds)), ones(1,T-floor((1600)/ds))]; % good guess
u0 = linspace(0, 1, T);

Op.lims  = [0 1];        
Op.plot = 1;               % plot the derivatives as well
Op.maxIter = 50;
Op.parallel = 0;


% === run the optimization!
[x, u, L, Vx, Vxx, cost, trace, stop] = iLQG(DYNCST, X0_vector, u0, Op);

[h,v,fpa,s] = get_states(x);

hm = get_stats(h);
sm = get_stats(s);
print_stats(x(:,end));

figure
plot(s', h'/1000)
hold on
plot(sm, hm/1000, 'k--')
xlabel('Downrange km')
ylabel('Altitude km')
grid on

figure
plot(v', h'/1000)
hold on
plot(v(1,:), hm/1000, 'k--')
xlabel('Velocity m/s')
ylabel('Altitude km')
grid on

function [h,v,fpa,s,g,L,D] = entry_accels(x)
% constants
cd  = 1.46;     
cl  = 0.35;    
rp = 3396.2e3;
m = 7200;
S = 15.8;
mu = 4.2830e13;

% states
[h,v,fpa,s] = get_states(x);

% Accels
rho = 0.0158*exp(-h/9354.5);
f = 0.5*rho.*v.^2 * S/m;
D = f*cd;
L = f*cl;
g = mu./(rp+h).^2;

function y = entry_dynamics(x,u)
global ds full_DDP hscale vscale fpascale rangescale n_states

% === states and controls:
% x = [h dv gamma]' 
% u = [u]'     = [cos(bank)]

% constants
rp = 3396.2e3;

% states + current accelerations
[h,v,fpa,s,g,L,D] = entry_accels(x);

% Derivs
sdot = v.*cos(fpa)/1000; % in km 
hdot = v.*sin(fpa);
vdot = -D-g.*sin(fpa);
fpadot = L./v.*u.^(1+full_DDP) + (v./(rp+h) - g./v).*cos(fpa);

xdot = [hdot/hscale; fpadot/fpascale; sdot/rangescale];            % change in state

dt = -ds./vdot;   % just for estimate 
y  = x + [ds/vscale+v*0;xdot.*repmat(dt, n_states, 1)];  % new state


function c = entry_cost(x, u)
global ds 

cost_scale = 10000;
% targets
% hf = 7.5;
sf = 334;

% weights
wu = 0.0;
wh = 0;
ws = 0;

% states
[h,v,fpa,s,g,L,D] = entry_accels(x);

% h_mean = get_stats(h);

% cost function 
% sum of 3 terms:
% lu: quadratic cost on controls
% lf: final cost 
% lx: running cost on altitude loss

final = isnan(u(1,:));
u(:,final)  = 0;

% control cost
lu    = wu*u.^2;

% running cost
hdot = v.*sin(fpa); 
vdot = -D - g.*sin(fpa);
lx = hdot./vdot * ds;
lx = get_stats(lx); % expected value of running cost 


% final cost
if any(final)
   llf      = ws*(sf-s).^2; % in real coordinates
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
    
    % dynamics second derivatives
    if full_DDP
        N_J = size(J);
        xu_Jcst = @(xu) finite_difference(xu_dyn, xu);
        JJ      = finite_difference(xu_Jcst, [x; u]);
%         JJ      = reshape(JJ, [4 6 size(J)]); % Original code
        if N_J <= 2 
            JJ = reshape(JJ,[3 4 N_J(2)]); 
        else 
            JJ = reshape(JJ, [4 5 N_J(2) N_J(3)]); 
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
H       = [h*eye(n)];
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

function c = tt(a,b)
c = bsxfun(@times,a,b);


function [h,v,fpa,s] = get_states(x)
global scale n_samples v0

iv = 1;
ih = iv + 1:n_samples+1;
ig = n_samples + ih;
is = n_samples + ig;

h = x(ih,:)*scale(1); % m
v = v0 - x(iv,:)*scale(2); % m/s
fpa = x(ig,:)*scale(3);
s = x(is,:)*scale(4); % km

function [m,v] = get_stats(x)
global weights
y = x.*weights;
m = sum(y, 1);
if nargout > 1
v = sum(weights.*(x-m).^2, 1);
end

function print_stats(x)
[h,v,fpa,s] = get_states(x);

[hm, hv] = get_stats(h/1000);
% [vm, vv] = get_stats(v);
% [fm, fv] = get_stats(fpa);
[sm, sv] = get_stats(s);


disp(['hf = ',num2str(hm),' +/- 3*',num2str(hv^0.5),' km'])
% disp(['fpa = ',num2str(vm),'=/- 3*',num2str(hv),' m/s'])
disp(['sf = ',num2str(sm),' +/- 3*', num2str(sv^0.5),' km'])