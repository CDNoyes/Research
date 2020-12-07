function entry_vel
% A demo of iLQG/DDP with car-parking dynamics
clc;
close all

fprintf(['\nA demonstration of the iLQG algorithm '...
'with entry dynamics, velocity loss as independent variable.\n'...
'for details see\nTassa, Mansard & Todorov, ICRA 2014\n'...
'\"Control-Limited Differential Dynamic Programming\"\n'])

% Set full_DDP=true to compute 2nd order derivatives of the 
% dynamics. This will make iterations more expensive, but 
% final convergence will be much faster (quadratic)
 
global hscale vscale fpascale ds v0 rangescale full_DDP
full_DDP = 1;

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

x0      = [39.4497e3/hscale, 0, (-10.604*pi/180)/fpascale, 0]';   % initial state
% u0      = ones(1,floor(T));    % initial controls
% u0 =     [zeros(1, floor(1600/ds)), ones(1,T-floor((1600)/ds))]; % good guess
u0 = linspace(0, 1, T);

Op.lims  = [0 1];         % wheel angle limits (radians)
Op.plot = 1;               % plot the derivatives as well
Op.maxIter = 28;
Op.parallel = 0;


% === run the optimization!
[x, u, L, Vx, Vxx, cost, trace, stop] = iLQG(DYNCST, x0, u0, Op);

h = x(1,:)*hscale/1000;
v = v0 - x(2,:)*vscale;
s = x(4,:)*rangescale;

disp(['hf = ',num2str(h(end)),' km'])
disp(['Vf = ',num2str(v(end)),' m/s'])
disp(['sf = ',num2str(s(end)),' km'])


figure
plot(s, h)
xlabel('Downrange km')
ylabel('Altitude km')
grid on

figure
plot(v, h)
xlabel('Velocity m/s')
ylabel('Altitude km')
grid on

function [g,L,D] = entry_accels(x)
global hscale vscale v0 
% constants
cd  = 1.46;     
cl  = 0.35;    
rp = 3396.2e3;
m = 7200;
S = 15.8;
mu = 4.2830e13;

% states
h = x(1,:)*hscale;
v = v0 - x(2,:)*vscale;

% Accels
rho = 0.0158*exp(-h/9354.5);
f = 0.5*rho.*v.^2 * S/m;
D = f*cd;
L = f*cl;
g = mu./(rp+h).^2;

function y = entry_dynamics(x,u)
global hscale vscale ds fpascale v0 rangescale full_DDP

% === states and controls:
% x = [h dv gamma]' 
% u = [u]'     = [cos(bank)]

% constants
rp = 3396.2e3;

% states
h = x(1,:)*hscale;
v = v0 - x(2,:)*vscale;
fpa = x(3,:)*fpascale;

[g,L,D] = entry_accels(x);

% Derivs
sdot = v.*cos(fpa)/1000; % in km 
hdot = v.*sin(fpa);
vdot = -D-g.*sin(fpa);
fpadot = L./v.*u.^(1+full_DDP) + (v./(rp+h) - g./v).*cos(fpa);

xdot = [hdot/hscale; -vdot/vscale; fpadot/fpascale; sdot/rangescale];            % change in state

dt = ds./(-vdot);   % just for estimate 
y  = x + xdot.*dt;  % new state


function c = entry_cost(x, u)
global hscale vscale wu hf cost_scale wh fpascale ds v0 rangescale

cost_scale = 10000;
% targets
hf = 7.5;
sf = 334;

% weights
wu = 0.0;
wh = 0;
ws = 0;

% states
h = x(1,:)*hscale;
v = v0 - x(2,:)*vscale;
fpa = x(3,:)*fpascale;
s = x(4,:)*rangescale;

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
lx = hdot./vdot * ds;


% final cost
if any(final)
   llf      = ws*(sf-s).^2; % in real coordinates
   lf       = zeros(size(u));
   lf(1,final)= llf(1,final);
else
   lf    = 0;
end

% total cost
c     = (lu + lx + lf)/cost_scale;

function [cx, cu, cxx, cxu, cuu] = cost_derivs(x, u)
global hscale vscale wu vf hf cost_scale wv wh fpascale ds v0

v = v0 - x(2,:)*vscale;
fpa = x(3,:)*fpascale;

cu = 2*wu*u/cost_scale;
cu(:,end) = 0;
cuu = repmat(2*wu/cost_scale, 1,1,length(u));
cuu(:,:,end) = 0;

% For running altitude cost and terminal quadratic cost alt
% cx = x*0;
% cx = []
cx(:,end) = [2*wh*hscale/1000*(x(1,end)*hscale/1000 - hf),0, 0]/cost_scale; % gradient of terminal cost 

cxu = zeros(3,1,length(u));
cxx = zeros(3,3,length(u));
% cxx(3,3,:) = -2*fpascale.^2*sec(fpa).^2.*tan(fpa)/cost_scale*ds/1000;

cxx(:,:,end) = [2*wh*(hscale/1000)^2, 0, 0;0, 0, 0; 0, 0, 0]/cost_scale;% hessian of terminal cost 



function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = entry_dyn_cst(x,u,full_DDP)
% combine dynamics and cost
% use helper function finite_difference() to compute derivatives

if nargout == 2
    f = entry_dynamics(x,u);
    c = entry_cost(x,u);
else
    % state and control indices
    ix = 1:4;
    iu = 5;
    
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
%     cx_ = gradient(x);
    
%     cost second derivatives
    xu_Jcst = @(xu) squeeze(complex_difference(xu_cost, xu));
    JJ      = finite_difference(xu_Jcst, [x; u]);
    JJ      = 0.5*(JJ + permute(JJ,[2 1 3])); %symmetrize
    cxx     = JJ(ix,ix,:);
    cxu     = JJ(ix,iu,:);
    cuu     = JJ(iu,iu,:);
%     [cx_, cu_, cxx_, cxu_, cuu_] = cost_derivs(x, u);

%     [cx, cu, cxx, cxu, cuu] = cost_derivs(x, u);

    
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

function lx = gradient(x)
% gradient of running cost wrt to [h,v,fpa]
global hscale vscale fpascale v0 ds
h = x(1,:)*hscale;
v = v0 - x(2,:)*vscale;
fpa = x(3,:)*fpascale;
rp = 3396.2e3;
mu = 4.2830e13;
cd  = 1.46;     
% cl  = 0.35;    
m = 7200;
S = 15.8;
rho0 = 0.0158;
hs = 9354.5;
lx = [-(v.*sin(fpa).*((2*mu*sin(fpa))./(h + rp).^3 + (S*cd*rho0*v.^2.*exp(-h/hs))./(2*hs*m)))./((mu*sin(fpa))./(h + rp).^2 + (S*cd*rho0*v.^2.*exp(-h/hs))./(2*m)).^2 
    (S*cd*rho0*v.^2.*exp(-h/hs).*sin(fpa))./(m*((mu*sin(fpa))./(h + rp).^2 + (S*cd*rho0*v.^2.*exp(-h/hs))./(2*m)).^2) - sin(fpa)./((mu.*sin(fpa))./(h + rp).^2 + (S*cd*rho0*v.^2.*exp(-h/hs))./(2*m))
    (mu*v.*cos(fpa).*sin(fpa))./((h + rp).^2.*((mu*sin(fpa))./(h + rp).^2 + (S*cd*rho0*v.^2.*exp(-h/hs))./(2*m)).^2) - (v.*cos(fpa))./((mu.*sin(fpa))./(h + rp).^2 + (S*cd*rho0*v.^2.*exp(-h/hs))./(2*m))]*ds;
 lx = lx./[hscale, vscale, fpascale]';