function entry
% A demo of iLQG/DDP with car-parking dynamics
clc;
close all

fprintf(['\nA demonstration of the iLQG algorithm '...
'with entry dynamics.\n'...
'for details see\nTassa, Mansard & Todorov, ICRA 2014\n'...
'\"Control-Limited Differential Dynamic Programming\"\n'])

% Set full_DDP=true to compute 2nd order derivatives of the 
% dynamics. This will make iterations more expensive, but 
% final convergence will be much faster (quadratic)
full_DDP = 1;
 
global hscale vscale fpascale ds
hscale = 120e3;
vscale = 5500;
fpascale = 0.2;
% set up the optimization problem
DYNCST  = @(x,u,i) entry_dyn_cst(x,u,full_DDP);
S       = 720; % terminal downrange, km
T       = 1000;              % horizon
ds      = S/(T) * 1000; % meters
x0      = [127e3/hscale, 5505/vscale, (-14.5*pi/180)/fpascale]';   % initial state
% u0      = ones(1,floor(T));    % initial controls
u0 =     [ones(1, floor(75000/ds)), zeros(1,floor((550000-75000)/ds)), ones(1,T-floor((550000-75000)/ds)-floor(75000/ds))]*0.9;
Op.lims  = [0 1];         % wheel angle limits (radians)
Op.plot = 1;               % plot the derivatives as well
Op.maxIter = 150;
Op.parallel = 0;


% === run the optimization!
[x, u, L, Vx, Vxx, cost, trace, stop] = iLQG(DYNCST, x0, u0, Op);

Jf = entry_cost(x(:,end), nan);
h = x(1,:)*hscale/1000;
v = x(2,:)*vscale;
s = linspace(0, S, T+1);
% cost
disp(['Jf = ',num2str(Jf)])

disp(['hf = ',num2str(h(end)),' km'])
disp(['Vf = ',num2str(v(end)),' m/s'])

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

function y = entry_dynamics(x,u)
global hscale vscale ds fpascale

% === states and controls:
% x = [h v gamma]' = [x; y; car_angle; front_wheel_velocity]
% u = [u]'     = [cos(bank)]

% constants
cd  = 1.46;      % d = distance between back and front axles
cl  = 0.35;     % h = timestep (seconds)
rp = 3396.2e3;
m = 7200;
S = 15.8;

% states
h = x(1,:)*hscale;
v = x(2,:)*vscale;
fpa = x(3,:)*fpascale;
% 
stop = v <= 2;
v(stop) = 2;
h(stop) = 0;

% Accels
rho = 0.0158*exp(-h/9354.5);
f = 0.5*rho.*v.^2 * S/m;
D = f*cd;
L = f*cl;
g = 4.2830e13./(rp+h).^2;

% Derivs
sdot = v.*cos(fpa);
hdot = v.*sin(fpa);
vdot = -D-g.*sin(fpa);
fpadot = L./v.*u.^2 + (v./(rp+h) - g./v).*cos(fpa);

xdot = [hdot/hscale; vdot/vscale; fpadot/fpascale];            % change in state

dt = ds./sdot; % just for estimate 
y  = x + xdot.*dt;                % new state

% function [fx, fu, fxx, fxu, fuu] = entry_derivs(x, u)
% h = x(1,:)*hscale;
% v = x(2,:)*vscale;
% fpa = x(3,:);
% 
% z = zeros(1,length(u))
% fx = [sec(fpa).^2


function c = entry_cost(x, u)
global hscale vscale wu vf hf cost_scale wh wv

cost_scale = 1;
hf = 7.5;
vf = 550;
wu = 0.01;
wv = 1;
% cost function 
% sum of 3 terms:
% lu: quadratic cost on controls
% lf: final cost on distance from target parking configuration
% lx: running cost on distance from origin to encourage tight turns

final = isnan(u(1,:));
u(:,final)  = 0;

% control cost
lu    = wu*u.^2;

% final cost
if any(final)
   llf      = (hf-x(1,end)*hscale/1000).^2 + wv*(x(2,end)*vscale-vf).^2;
%    llf = x(2,end)*vscale - x(1,end)*hscale/1000;
%    lf       = double(final);
   lf = zeros(size(u));
   lf(1,end)= llf;
else
   lf    = 0;
end

% running cost
lx = 0; %cx*sabs(x(1:2,:),px);

% total cost
c     = (lu + lx + lf)/cost_scale;

function [cx, cu, cxx, cxu, cuu] = cost_derivs(x, u)
global hscale vscale wu vf hf cost_scale wv

cu = 2*wu*u/cost_scale;
cu(:,end) = 0;
cuu = repmat(wu/cost_scale, 1,1,length(u));
cuu(:,:,end) = 0;

% For quadratic cost on alt and vel
cx = 0*x;
cx(:,end) = [2*hscale/1000*(hscale/1000*x(1,end)-hf), wv*2*vscale*(x(2,end)*vscale-vf), 0]/cost_scale; % gradient of terminal cost 
cxu = zeros(3,1,length(u));
cxx = zeros(3,3,length(u));
cxx(:,:,end) = [2*(hscale/1000)^2, 0, 0;0, wv*2*vscale^2, 0; 0, 0, 0]/cost_scale;% hessian of terminal cost 

% For linear cost on both 
% cx = 0*x;
% cx(:,end) = [-hscale/1000, vscale, 0]/cost_scale; % gradient of terminal cost 
% cxu = zeros(3,1,length(u));
% cxx = zeros(3,3,length(u));

function y = sabs(x,p)
% smooth absolute-value function (a.k.a pseudo-Huber)
y = pp( sqrt(pp(x.^2,p.^2)), -p);


function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = entry_dyn_cst(x,u,full_DDP)
% combine dynamics and cost
% use helper function finite_difference() to compute derivatives

if nargout == 2
    f = entry_dynamics(x,u);
    c = entry_cost(x,u);
else
    % state and control indices
    ix = 1:3;
    iu = 4;
    
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
            JJ = reshape(JJ, [3 4 N_J(2) N_J(3)]); 
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
%     
%     % cost second derivatives
%     xu_Jcst = @(xu) squeeze(finite_difference(xu_cost, xu));
%     JJ      = finite_difference(xu_Jcst, [x; u]);
%     JJ      = 0.5*(JJ + permute(JJ,[2 1 3])); %symmetrize
%     cxx     = JJ(ix,ix,:);
%     cxu     = JJ(ix,iu,:);
%     cuu     = JJ(iu,iu,:);
    
    [cx, cu, cxx, cxu, cuu] = cost_derivs(x, u);
    
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

h = 0.001j;

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