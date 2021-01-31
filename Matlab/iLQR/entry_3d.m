function sol = entry_3d(inp, lagrange, penalty)
% Nominal 3D Entry Dynamics (downrange instead of longitude is used)
% Allows for iterations with different values of the lagrange multiplier
% and penalty value 


% Set full_DDP=true to compute 2nd order derivatives of the
% dynamics. This will make iterations more expensive, but
% final convergence will be much faster (quadratic)

if nargin == 0
    inp = DDPInput([0,0,0]);
    lagrange = 0;
    penalty = 0.1;
end

% state vector: [h, theta, phi, v, gamma, psi]


global scale ds
full_DDP = 1;

hscale = 120e3;
vscale = 5500;
fpascale = 0.2;
v0 = 5525;
rangescale = 500e3;
scale = [hscale, rangescale, 1, vscale, fpascale, 1];

% set up the optimization problem
DYNCST  = @(x,u,i) entry_dyn_cst(x,u,full_DDP,lagrange,penalty);
V       = 460; % terminal velocity, m/s
dV      = V-v0;
T       = inp.horizon;              % horizon
ds      = dV/T;

x0      = [54.5e3/hscale, 0, 0, v0/vscale, (-11.5*pi/180)/fpascale, 0]';   % initial state
if isfield(inp, 'guess') && ~isempty(inp.guess)
    u0 = inp.guess;
else
%     u0 = acos(ones(1,floor(T))*0.7).*sin(linspace(0,2*pi,T))*0;
    u0 = ones(1,T);
end

Op.lims  = [-90*pi/180 90*pi/180];         % angle limits (radians)
Op.plot = inp.running_plots;               % plot the derivatives as well
Op.maxIter = inp.max_iterations;
Op.parallel = 1;


% === run the optimization!

    [x, u, L, Vx, Vxx, cost, trace, stop] = iLQG(DYNCST, x0, u0, Op);
    
[g,L,D] = entry_accels(x);
h = x(1,:)*scale(1);
s = x(2,:)*scale(2)/1000; 
phi = x(3,:)*scale(3);
v = x(4,:)*scale(4);
fpa = x(5,:)*scale(5);
psi = x(6,:)*scale(6);
% s = 3396.2*theta; % TODO compute the actual traj length

constraint = entry_constraints(x,u);

disp(['hf = ',num2str(h(end)/1000),' km'])
disp(['Vf = ',num2str(v(end)),' m/s'])
disp(['sf = ',num2str(s(end)),' km'])
disp(['final crossrange = ',num2str(-phi(end)*3396.2),' km'])

if inp.terminal_plots
    figure
    plot(s, h/1000)
    xlabel('Downrange km')
    ylabel('Altitude km')
    grid on

    figure
    plot(v, h/1000)
    xlabel('Velocity m/s')
    ylabel('Altitude km')
    grid on
end

sol.u = u;
sol.v = v;
sol.state = [h; s; phi; v; fpa; psi];
% sol.h = h;
% sol.fpa = fpa;
% sol.s = s;
% sol.weights = [wh,ws,wu];
sol.cost = sum(cost);
sol.constraint = constraint(:,end); 
sol.inp = inp;
sol.L = L;
sol.D = D;


function [g,L,D] = entry_accels(x)
global scale
% constants
[m,S,cl,cd] = aero_const();

rp = 3396.2e3;
mu = 4.2830e13;

% states
h = x(1,:)*scale(1);
v = x(4,:)*scale(4);

% Accels
rho = 0.0158*exp(-h/9354.5);
f = 0.5*rho.*v.^2 * S/m;
D = f*cd;
L = f*cl;
g = mu./(rp+h).^2;

function y = entry_dynamics(x,u)
global ds scale

% === states and controls:
% x = [h dv gamma]'
% u = [u]'     = [cos(bank)]

% constants
rp = 3396.2e3;

% states
h = x(1,:)*scale(1);
% theta = x(2,:); 
phi = x(3,:)*scale(3);
v = x(4,:)*scale(4);
fpa = x(5,:)*scale(5);
psi = x(6,:)*scale(6);

[g,L,D] = entry_accels(x);

% Derivs
% thetadot = v./(rp+h).*cos(fpa).*cos(psi)./cos(phi); 
sdot = v.*cos(fpa);
phidot = v./(rp+h).*cos(fpa).*sin(psi);
hdot = v.*sin(fpa);
vdot = -D-g.*sin(fpa);
fpadot = L./v.*cos(u) + (v./(rp+h) - g./v).*cos(fpa);
psidot = -1./v./cos(fpa).*(L.*sin(u) + v.^2./(rp+h).*cos(fpa).^2.*cos(psi).*tan(phi));
xdot = [hdot; sdot; phidot; vdot; fpadot; psidot]./scale';            % change in state

dt = ds./(vdot);   % just for estimate
y  = x + xdot.*dt;  % new state


function c = entry_cost(x, u, lagrange, penalty)
global scale ds

cost_scale = 10000;
% targets


% weights
wu = 0.0;
wh = 0.0;
ws = 0.0;

% states
% h = x(1,:)*scale(1);
% s = x(2,:)*scale(2); 
phi = x(3,:)*scale(3);
v = x(4,:)*scale(4);
fpa = x(5,:)*scale(5);
% psi = x(6,:)*scale(6);

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
[g,~,D] = entry_accels(x);

hdot = v.*sin(fpa);
vdot = -D - g.*sin(fpa);
lx = -hdot./vdot * ds;


% final cost
if any(final)
    lf = 0;
%     llf      = 10000*3397*phi.^2; % in real coordinates
%     lf       = zeros(size(u));
%     lf(1,final)= llf(1,final);
else
    lf    = 0;
end

% Constraint terms
cons = entry_constraints(x,u);
if any(final)
    con = zeros(size(cons,1),size(u,2));
    con(:,final) = cons(:,final);
else
    con = 0*lagrange;
end

%Lagrange should be a vector same length as con
%Penalty is generally a scalar even when con is not 

% total cost
c     = (lu + lx + lf + sum(lagrange.*con, 1) + 0.5*penalty*sum(con.^2))/cost_scale ;


function c = entry_constraints(x,u)
global scale 
% h = x(1,:)*scale(1);
s = x(2,:)*scale(2); 
phi = x(3,:)*scale(3);
% v = x(4,:)*scale(4);
% fpa = x(5,:)*scale(5);
% psi = x(6,:)*scale(6);
c = [phi*3397; s/1000-350];% final latitude = 0


function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = entry_dyn_cst(x,u,full_DDP,lagrange,penalty)
% combine dynamics and cost
% use helper function finite_difference() to compute derivatives

if nargout == 2
    f = entry_dynamics(x,u);
    c = entry_cost(x,u,lagrange,penalty);
else
    % state and control indices
    ix = 1:6;
    iu = 7;
    
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
    xu_cost = @(xu) entry_cost(xu(ix,:),xu(iu,:),lagrange,penalty);
    J       = squeeze(complex_difference(xu_cost, [x; u]));
    cx      = J(ix,:);
    cu      = J(iu,:);
    
    %     cost second derivatives
        xu_Jcst = @(xu) squeeze(complex_difference(xu_cost, xu));    
    
        JJ      = finite_difference(xu_Jcst, [x; u]);
        JJ      = 0.5*(JJ + permute(JJ,[2 1 3])); %symmetrize
        cxx     = JJ(ix,ix,:);
        cxu     = JJ(ix,iu,:); % all zeros for alt obj
        cuu     = JJ(iu,iu,:); % all zeros for alt obj

    
    
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
