function sol = entry_stochastic(input)
fprintf(['\nUse of the iLQG algorithm '...
    'with entry dynamics, velocity loss as independent variable.\n'])

% Set full_DDP=true to compute 2nd order derivatives of the
% dynamics. This will make iterations more expensive, but
% final convergence will be much faster (quadratic)

global ds v0 full_DDP n_samples n_states scale weights closed wh ws wu gains QN
% note: weights are the sigma point weights
% w_ are the cost weights

n_states = 3;
if nargin == 0
    input = DDPInput([3, 3, 0.2]);
    input.ut_scale = 17; % pretty good for a lot
%     input.ut_scale = 6;
%     input.gains = [0,0,0];
end
closed = input.closed_loop;
wh = input.weights(1);
ws = input.weights(2);
wu = input.weights(3);
full_DDP = input.ddp;
gains = @(v) input.gains;
QN = input.qn;

% Scale states for numerical condition
hscale = 120e3;
vscale = 5500;
fpascale = 0.2;
rangescale = 500e3;
scale = [hscale, vscale, fpascale, rangescale];
v0 = input.v0;

% set up the optimization problem
DYNCST  = @(x,u,i) entry_dyn_cst(x, u, full_DDP,i);
V       = input.vf; % terminal velocity, m/s
dV      = v0-V;
T       = input.horizon;              % horizon
ds      = dV/T;

% Initial States
x0 = input.x0./[1;hscale;fpascale;1];
Pscale = scale([1,3,4])'*scale([1,3,4]);
P0 = input.P0./Pscale; % scaled initial covariance

if 1
    correlate = 0;
    if correlate
        d = diag(P0).^0.5; % std devs
        rho = [0, d(1)*d(2), d(1)*d(3)
               d(1)*d(2), 0, d(2)*d(3)
               d(1)*d(3), d(2)*d(3),0]; % maximum correlations
       P0 = P0 - 0.3*rho;
    end
    [X0, weights] = UnscentedTransform(x0(2:end), P0, input.ut_scale);
    n_samples = length(weights);
    X0 = [zeros(n_samples,1), X0'];
    X0_vector = X0(n_samples:end)';
    
else % Random sampling
    n_samples = 10;
    
    X0 = repmat(x0(2:end)', n_samples, 1) + randn(n_samples, 3)*chol(P0);
    X0 = [zeros(n_samples, 1), X0];
    X0_vector = X0(n_samples:end)';
    weights = ones(n_samples,1)/n_samples;
end

fprintf(['Samples = ',num2str(n_samples),'\nTotal State Dimension = ',num2str(1+n_samples*n_states),'\n'])

% Initial Control
% u0      = ones(1,floor(T))*0.7;
% u0 =     [zeros(1, floor(1600/ds)), ones(1,T-floor((1600)/ds))]; % good guess
u0 = linspace(0.35, 1, T);
% u0 = linspace(1, 0, T);


Op.lims  = input.bounds;
Op.plot = input.running_plots;               % plot the derivatives as well
Op.maxIter = input.max_iterations;
Op.parallel = 1;
Op.regType = 1; % default 1

% Kopt = optimize_gain(X0_vector, u0, [0,0,0]);
LQR = 0;
if LQR
    initialize_gains(X0_vector, u0)
    x = forward_pass(X0_vector, u0, [], T);
    u = u0;
    print_stats(x(:,end));

    cost = 0;
    ddp_gains = [];
end


% === run the optimization!
if QN
    [x, u, ddp_gains, ~, ~, cost] = QNDDP(DYNCST, X0_vector, u0, Op);
else
    [x, u, ddp_gains, ~, ~, cost] = iLQG(DYNCST, X0_vector, u0, Op);
end

% Kopt = optimize_gain(X0_vector, u);
% sol.gains = Kopt;

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
sol.X0 = [h(:,1)',fpa(:,1)',s(:,1)'];
sol.input = input;
sol.sigma_weights = weights;
sol.L = L;
sol.D = D;
sol.Dm = Dm;
sol.Dv = Dv;
sol.stats = [hm(end)/1000, (hm(end)-3*hv(end)^0.5)/1000, 3*sv(end).^0.5]; % 3-sigma low altitude, range std

print_stats(x(:,end));
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
    
    G = squeeze(ddp_gains(1,2:end,1:end-1));
    G = reshape(G, n_samples, n_states, input.horizon-1);
    G = squeeze(get_stats(G));
    
    %     G0 = squeeze(ddp_gains(1,2:n_samples:end,1:end-1));
    
    figure
    plot(v(1:end-2), u(1,1:end-1), 'k', 'linewidth', lw)
    hold all
    %     plot(v(1:end-2), G/100, 'linewidth', 1)
    %     plot(v(1:end-2), G0/100, '--','linewidth', 1)
    
    xlabel('Velocity, m/s')
    ylabel('Controls')
    %     legend('Open Loop',['K / ',num2str(100)])
    grid on
    
end

function [h,v,fpa,s,g,L,D] = entry_accels(x)
% constants
[m,S,cl,cd] = aero_const();

rp = 3396.2e3;

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
global ds full_DDP scale n_states closed gains

% === states and controls:
% x = [h dv gamma]'
% u = [u]'     = [cos(bank)]

% constants
rp = 3396.2e3;

% states + current accelerations
[h,v,fpa,s,g,L,D] = entry_accels(x);

% feedback terms
% Tuned for decent performance in the guess trajectory
K = gains(real(v));
%     kd = gains(1); %more lift up when too much drag
%     ks = gains(2);  % less lift up when too close
%     kf = gains(3);% left lift up when to shallow

ef = get_error(fpa);
% ef = get_error(v.*cos(fpa));
es = get_error(s/1000);
eD = get_error(D);
N = size(ef,1);
du = eD.*repmat(K(:,1)',N,1) + es.*repmat(K(:,2)',N,1) + ef.*repmat(K(:,3)',N,1); % NOTE TO SELF: Removed a mult by 200 on ef when I added LQR gains
if size(ef,2)  == 1 % Use hard saturation to propagate?
    u_cl = Saturate(u+du,0,1);
else
    %     du = Saturate(du,-2,2);
    u_cl = smooth_sat(u + du); % cosh
end


% Derivs
sdot = v.*cos(fpa); % in km
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


function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = entry_dyn_cst(x,u,full_DDP,iter)
global n_states n_samples QN
% combine dynamics and cost
% use helper function finite_difference() to compute derivatives

if nargout == 2
    f = entry_dynamics(x,u);
    c = entry_cost(x,u);
else
    
%     update_gains(x,u)
    
    % state and control indices
    ix = 1:(1+n_states*n_samples);
    iu = ix(end)+1;
    
    % dynamics first derivatives
    xu_dyn  = @(xu) entry_dynamics(xu(ix,:),xu(iu,:));
    J       = complex_difference(xu_dyn, [x; u]);
    fx      = J(:,ix,:);
    fu      = J(:,iu,:);
    
    % dynamics second derivatives
    if full_DDP || (QN && iter <= 2)
        if 1
            N_J = size(J);
            xu_Jcst = @(xu) complex_difference(xu_dyn, xu);
            JJ      = finite_difference(xu_Jcst, [x; u]);
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
            u_dyn  = @(u) entry_dynamics(repmat(x,1,size(u,2)/size(x,2)), u);
            u_Jcst = @(u) complex_difference(u_dyn, u);
            fuu = finite_difference(u_Jcst, u);
            fuu = reshape(fuu, [ix(end) 1 1 size(J,3)]);
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
    
    if QN && iter > 2
        % if not QuasiNewton, otherwise, these can also be approximated, right?
        % we can either compute first iteration (probably wise)
        % or set to zero or identity for first iteration
        [cxx,cxu,cuu] = deal([]);
    else
        %     cost second derivatives
        xu_Jcst = @(xu) squeeze(complex_difference(xu_cost, xu));
        JJ      = finite_difference(xu_Jcst, [x; u]);
        JJ      = 0.5*(JJ + permute(JJ,[2 1 3])); %symmetrize
        cxx     = JJ(ix,ix,:);
        cxu     = JJ(ix,iu,:);
        cuu     = JJ(iu,iu,:);
    end
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

function std = gain_opt_obj(x0, u, K, N)

x = forward_pass(x0, u, K, N);
[h,v,fpa,s,g,L,D] = entry_accels(x(:,end));
[~,std] = get_stats(s(:,end));
[~,hstd] = get_stats(h(:,end));
std = std + 0.1*hstd;

function x = forward_pass(x0, u, K, N)
global gains
if ~isempty(K)
    gains = @(V) K;
end
x = x0(:);

for i = 1:N-1
    x(:,i+1) = entry_dynamics(x(:,i), u(i));
end

function initialize_gains(x0, u0)
global gains
x = forward_pass(x0, u0, gains(0), length(u0));
update_gains(x,u0);

function update_gains(x,u)
global gains

[h,v,fpa,s,g,L,D] = entry_accels(x);

% LQR
if 0
    % This expects UNSCALED states
    [A,B] = EntryJacobians([v;get_stats(h);get_stats(fpa);get_stats(s)], u, 1);
    
    Qf = diag([0, 0, 100]).^2;
    Q  = diag([0.01, 1, 1]).^2;
    R  = 1e3;
    
    K = LQRGains(A,B,Qf,Q,R); % gains in h, fpa, s
    
else % Apollo Gains
    t = cumtrapz(v, -1./get_stats(D+g.*sin(fpa)));
%     [A,B] = EntryJacobians([v;get_stats(h);get_stats(fpa);get_stats(s)], u, 0);
    
%     K = ApolloGains(t,v,[get_stats(h);get_stats(fpa);get_stats(s)],u,get_stats(L),get_stats(D),A,B); % gains in h, fpa, s
    K = ApolloGains2(t,v,[get_stats(h);get_stats(fpa);get_stats(s)],u,get_stats(L),get_stats(D));
    
end

% K(1,:) = -K(1,:)*9354.5./get_stats(D); % gain in Drag
% K(3,:) = K(3,:)*1000; % to account for using km in place of m

gain_hold = 600;
Kf = K(:,v<gain_hold);
K(:, v<gain_hold) = repmat(Kf(:,1),1,sum(v<gain_hold));
if 0
    figure
    plot(v, K(1:3,:)')
    legend('Drag','Range', 'FPA')
    grid on
end
k3 = 1;
gains = @(vn) interp1(v, k3*K', vn);


function Kopt = optimize_gain(x0, u, guess)
options = optimoptions('fminunc');
% options.Algorithm = 'trust-region';
N = length(u);
obj = @(K) gain_opt_obj(x0, u, K, N);

Kopt = fminunc(obj, guess, options);
