function sol = entry_stochastic_gains(input)
% A demo of iLQG/DDP with stochastic entry dynamics and joint gain
% optimization
% clc;
% close all

fprintf(['\nUse of the iLQG algorithm '...
    'with entry dynamics, velocity loss as independent variable.\n'...
    'Performing joint gain optimization.\n'])

% Set full_DDP=true to compute 2nd order derivatives of the
% dynamics. This will make iterations more expensive, but
% final convergence will be much faster (quadratic)

global ds v0 full_DDP n_samples n_states n_controls scale weights wh ws wu
% note: weights are the sigma point weights
% w_ are the cost weights

n_states = 3;
n_controls = 4; % open loop + 3 feedback terms 

if nargin == 0
    input = DDPInput([2, 1, 0.2]);
    input.ut_scale = 5;

end
full_DDP = input.ddp;
wh = input.weights(1);
ws = input.weights(2);
wu = input.weights(3);

% Scale states for numerical condition
hscale = 120e3;
vscale = 5500;
fpascale = 0.2;
rangescale = 500e3;
scale = [hscale, vscale, fpascale, rangescale];

v0 = input.v0;

% set up the optimization problem
DYNCST  = @(x,u,i) entry_dyn_cst(x, u, full_DDP);
V       = input.vf; % terminal velocity, m/s
dV      = v0-V;
T       = input.horizon;              % horizon
ds      = dV/T;

% Initial States
x0 = input.x0./[1;hscale;fpascale;1];
Pscale = scale([1,3,4])'*scale([1,3,4]);
P0 = input.P0./Pscale; % scaled initial covariance 

[X0, weights] = UnscentedTransform(x0(2:end), P0, input.ut_scale);
n_samples = length(weights);
X0 = [zeros(n_samples,1), X0'];
X0_vector = X0(n_samples:end)';

fprintf(['Samples = ',num2str(n_samples),'\nTotal State Dimension = ',num2str(1+n_samples*n_states),'\n'])

% Initial Control
% u0      = ones(1,floor(T))*0.4;
% u0 =     [zeros(1, floor(1600/ds)), ones(1,T-floor((1600)/ds))]; % good guess
% u0 = [linspace(0.35, 1, T-200), ones(1,200)];
u0 = linspace(0.35, 1, T);

kd = input.gains(1); 
ks = input.gains(2);  
kf_scale = 200; % dont touch this
kf = input.gains(3)/kf_scale;

L = ones(T,1);
% L = linspace(1, 2, T)';
K0 = [L*kd, L*ks, L*kf]';
% K0 = 0*K0;

% Op.lims  = [input.bounds; [0,1]; [-1,0];[-200, 0]/kf_scale]; % final gain is scaled in dynamics
Op.lims  = [input.bounds; [-0.25,1]; [-0.5, 0.1];[-100, 50]/kf_scale]; % allow 'bad' sign gains

Op.plot = input.running_plots;               % plot the derivatives as well
Op.maxIter = input.max_iterations;
Op.parallel = 1;


% === run the optimization!
[x, u, ~, ~, ~, cost] = iLQG(DYNCST, X0_vector, [u0;K0], Op);


% Compute stats and setup the output structure 
[h,v,fpa,s,~,L,D] = entry_accels(x);

v=v';
s = s/1000;
u(4,:) = u(4,:)*kf_scale; % allows the gains to be used outside this program

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
sol.X0 = [h(:,1)',fpa(:,1)',s(:,1)'];
sol.P0 = P0;
sol.input = input;
sol.sigma_weights = weights;
sol.L = L;
sol.D = D;
sol.Dm = Dm;
sol.Dv = Dv;

print_stats(x(:,end));

% Plots
lw = 2;

if input.terminal_plots
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

    figure
    plot(v(1:end-2), u(1,1:end-1), 'k', 'linewidth', lw)
    hold all
    us = u;
    us(end,:) = us(end,:)/kf_scale;
    plot(v(1:end-2), us(2:end,1:end-1), 'linewidth', 1)
    xlabel('Velocity, m/s')
    ylabel('Controls')
    legend('Open Loop','K_D','K_s',['K_f / ',num2str(kf_scale)])
    grid on
end

function [h,v,fpa,s,g,L,D] = entry_accels(x)
% constants
cd  = 1.408;
cl  = 0.357;
rp = 3396.2e3;
m = 5000;
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

function y = entry_dynamics(x,U)
global ds scale n_states

% === states and controls:
% x = [h dv gamma]'
% u = [cos(bank), K1, K2, K3]'

% constants
rp = 3396.2e3;

% states + current accelerations
[h,v,fpa,s,g,L,D] = entry_accels(x);

% u = U(1,:);
% K = U(2:4,:);

% feedback terms
ef = get_error(fpa);
es = get_error(s/1000);
eD = get_error(D);
du = U(2,:).*eD + U(3,:).*es + U(4,:).*ef*200;
u_cl = smooth_sat(U(1,:) + du); % cosh


% Derivs
sdot = v.*cos(fpa); 
hdot = v.*sin(fpa);
vdot = -D-g.*sin(fpa);
fpadot = L./v.*u_cl + (v./(rp+h) - g./v).*cos(fpa);

xdot = [hdot/scale(1); fpadot/scale(3); sdot/scale(4)];            % change in state

dt = -ds./vdot;   % just for estimate
y  = x + [ds/scale(2)+v*0; xdot.*repmat(dt, n_states, 1)];  % new state

function y = smooth_sat(x)
K = 20;
q = 2*x - 1;
y =  0.5/K * log(cosh(K*(q+1))./cosh(K*(q-1))); % saturates between [-1,1]
y = 0.5 + 0.5*y; % map back to [0,1]


function er = get_error(state)
er = state - get_stats(state);
% the mean subtracted from all the points

function c = entry_cost(x, u)
global ds wh ws wu

cost_scale = 10000;

% states
[h,v,fpa,s,g,~,D] = entry_accels(x);

% cost function
% sum of terms:
% lu: quadratic cost on controls
% lx: running cost on altitude loss


% Don't delete these two lines even if unused 
final = isnan(u(1,:));
u(:,final)  = 0;

% running cost terms
hdot = v.*sin(fpa);
vdot = -D - g.*sin(fpa);
hprime = hdot./vdot * ds; %it's actually -hdot over -vdot
hprime = get_stats(hprime); % expected value of running cost

sdot = v.*cos(fpa); % meters/second
[smean, svar] = get_stats(s); % m
[hmean, hvar] = get_stats(h); % m

vs_dot = -2*get_stats(s.*sdot./vdot) + 2*smean.*get_stats(sdot./vdot);
vh_dot = -2*get_stats(h.*hdot./vdot) + 2*hmean.*get_stats(hdot./vdot);

%convert variance rates to std rates
vs_dot = 0.5*vs_dot./(svar.^0.5);
vh_dot = 0.5*vh_dot./(hvar.^0.5);

lx = hprime + ws*vs_dot*ds + wh*vh_dot*ds;

% control cost
lu    = wu*(u(1,:)-0.5).^2 * ds.* v/2500; % this term is for smoothness
% lu    = wu*(sum(u(2:end,:).^2,1)+ (u(1,:)-0.5).^2) * ds.* v/2500; % terms
% on gains don't seem to help much 

% total cost
c     = (lu + lx)/cost_scale;


function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = entry_dyn_cst(x,u,full_DDP)
global n_states n_samples n_controls
% combine dynamics and cost
% use helper function finite_difference() to compute derivatives

if nargout == 2
    f = entry_dynamics(x,u);
    c = entry_cost(x,u);
else
    % state and control indices
    ix = 1:(1+n_states*n_samples);
    iu = ix(end) + (1:n_controls);
    
    % dynamics first derivatives
    xu_dyn  = @(xu) entry_dynamics(xu(ix,:),xu(iu,:));
    J       = complex_difference(xu_dyn, [x; u]);
    fx      = J(:,ix,:);
    fu      = J(:,iu,:);
    
    % dynamics second derivatives
    if full_DDP
        if 1
            N_J = size(J);
            xu_Jcst = @(xu) complex_difference(xu_dyn, xu);
            JJ      = finite_difference(xu_Jcst, [x; u]);
            %         JJ      = reshape(JJ, [4 6 size(J)]); % Original code
            if N_J <= 2
                JJ = reshape(JJ,[ix(end) iu(end) N_J(2)]);
            else
                JJ = reshape(JJ, [ix(end) iu(end) N_J(2) N_J(3)]);
            end
            JJ      = 0.5*(JJ + permute(JJ,[1 3 2 4])); %symmetrize
            fxx     = JJ(:,ix,ix,:);
            fxu     = JJ(:,ix,iu,:);
            fuu     = JJ(:,iu,iu,:);
        else
            u_dyn  = @(u) entry_dynamics(repmat(x,1,size(u,2)/size(x,2)), u);
            u_Jcst = @(u) complex_difference(u_dyn, u);
            fuu = finite_difference(u_Jcst, u);
            fuu = reshape(fuu, [ix(end) 4 4 size(J,3)]);
            fuu = 0.5*(fuu + permute(fuu,[1 3 2 4]));
            fxx = [];
            fxu = [];  
            
        end
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
fpa = x(ig,:)*scale(3);
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