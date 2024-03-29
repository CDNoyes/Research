function sol = FixedGainOpt(input)

global ds v0 full_DDP n_samples n_states n_controls scale weights wh ws wu n_params
n_states = 3;
n_controls = 1;
n_params = 3;

if nargin == 0
    wh = 1;
    ws = 1;
    wu = 0.3;
    input = DDPInput([wh,ws,wu]);
    input.ut_scale = 16;
    input.horizon = 250;
    load solutions_cl_ddp_max_guess
    input.guess = sols{6}.u;
end
W = input.weights;
wh = W(1);
ws = W(2);
wu = W(3);
full_DDP = input.ddp;

% Scale states for numerical condition
hscale = 120e3;
vscale = 5500;
fpascale = 0.2;
rangescale = 500e3;
scale = [hscale, vscale, fpascale, rangescale];

v0 = input.v0;

% set up the optimization problem
V       = input.vf; % terminal velocity, m/s
dV      = v0-V;
T       = input.horizon;
ds      = dV/T;

% Initial States
x0      = [0, 39.4497e3/hscale, (-10.604*pi/180)/fpascale, 0, ones(1, n_params)]';   % nominal initial state

% Note that horizon must be set lower because of the increased state
% dimension
[X0, weights] = UnscentedTransform(x0(2:end), 1*diag([2500/hscale, 0.25*pi/180/fpascale, 10000/rangescale, 5/100, 5/100, 7/100]).^2, input.ut_scale);
n_samples = length(weights);
X0 = [zeros(n_samples,1), X0'];
X0_vector = X0(n_samples:end)';

fprintf(['Samples = ',num2str(n_samples),'\nTotal State Dimension = ',num2str(length(X0_vector)),'\n'])

% Initial Control
if isempty(input.guess)
    % u0      = ones(1,floor(T))*0.4;
    % u0 =     [zeros(1, floor(1600/ds)), ones(1,T-floor((1600)/ds))]; % good guess
    % u0 = [linspace(0.35, 1, T-floor(T/4)), ones(1,floor(T/4))];
    u0 = linspace(0.35, 1, T);
else
    u0 = input.guess(1,:);
    disp('Loaded reference control from guess')
end

% Tuned for decent performance in the guess trajectory
kd = input.gains(1); %more lift up when too much drag
ks = input.gains(2);  % less lift up when too close
kf_scale = 200; % dont touch this
kf = input.gains(3)/kf_scale;% less lift up when too shallow

K0 = [kd, ks, kf];
% K0 = zeros(3,1);
Kopt = optimize_gain(X0_vector, u0, K0, W);
x = forward_pass(X0_vector, u0, Kopt, T+1);
u = u0;
Kopt(3) = Kopt(3)*kf_scale; % lol dont do this till after the forward pass
input.gains = Kopt;

% Compute stats 
[h,v,fpa,s,~,L,D] = entry_accels(x);
v=v';
s = s/1000;
[hm, hv] = get_stats(h);
[sm, sv] = get_stats(s);
[fm, fv] = get_stats(fpa);
[Dm, Dv] = get_stats(D);
[Lm, Lv] = get_stats(L);

sol.u = u;
sol.v = v;
sol.mean = [hm; fm; sm];
sol.var = [hv; fv; sv];
sol.h = h;
sol.fpa = fpa;
sol.s = s;
sol.weights = [wh,ws,wu];
sol.X0 = [h(:,1)',fpa(:,1)',s(:,1)'];
sol.input = input;
sol.sigma_weights = weights;
sol.L = L;
sol.D = D;
sol.Lm = Lm;
sol.Dm = Dm;
sol.Dv = Dv;
sol.Lv = Lv;
[m,S,cl,cd] = aero_const();
sol.mass = m;
sol.area = S;
sol.cl = cl;
sol.cd = cd;
sol.BC = m/(S*cd);

print_stats(x(:,end));

% Plots
lw = 2;
if input.terminal_plots
figure
x2 = [v', fliplr(v')];
inBetween = [(sm-3*sv.^0.5)/1000, fliplr((sm+3*sv.^0.5)/1000)];
fill(x2, inBetween, 'c');
hold all
plot(v, s'/1000, 'k--','linewidth', 1)
plot(v, sm/1000, 'm', 'linewidth', lw)

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
plot(v(1:end-1), u(1,1:end-1), 'k', 'linewidth', lw)
hold all

xlabel('Velocity, m/s')
ylabel('Controls')
grid on
end

function [h,v,fpa,s,g,L,D] = entry_accels(x)
% constants
[m,S,cl,cd] = aero_const();

rp = 3396.2e3;
mu = 4.2830e13;

% states
[h,v,fpa,s, fcl, fcd, frho] = get_states(x);

% Accels
rho = frho.*0.0158.*exp(-h./9354.5);
f = 0.5*rho.*v.^2*S/m;
D = f.*cd.*fcd;
L = f.*cl.*fcl;
g = mu./(rp+h).^2;

function x = entry_dynamics(x, U, K)
% Calls the discrete time Euler approx multiple times for better accuracy
global ds
n_int_steps = 4;
dsn = ds/n_int_steps;

k = size(U,2)/size(x,2);
x = repmat(x,1,k);

for i = 1:n_int_steps
    dx = entry_step(x, U, K);
    x = x + dsn*dx;
end

function x = entry_step(x, U, K)
global scale n_states n_params

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

du = K(1)*eD + K(2)*es + K(3)*ef*200; % [0.0725, -0.025, -4]


LoD = L./D;
LoDr = get_stats(LoD);
u_cl = smooth_sat(LoDr./LoD.*U(1,:) + du); % cosh

% Derivs
sdot = v.*cos(fpa);
hdot = v.*sin(fpa);
vdot = -D-g.*sin(fpa);
fpadot = L./v.*u_cl + (v./(rp+h) - g./v).*cos(fpa);
xdot = [hdot/scale(1); fpadot/scale(3); sdot/scale(4)];            % change in state

% dt = -ds./vdot;   % just for estimate
% x  = x + [ds/scale(2)+v*0; xdot.*repmat(dt, n_states, 1); 0*x(2+n_states*size(dt,1):end, :)];  % new state
x = [1/scale(2)+v*0; xdot.*repmat(-1./vdot, n_states, 1); zeros(size(L).*[n_params,1])]; % really dx/dv

function y = smooth_sat(x)
K = 20;
q = 2*x - 1;
y =  0.5/K * log(cosh(K*(q+1))./cosh(K*(q-1))); % saturates between [-1,1]
y = 0.5 + 0.5*y; % map back to [0,1]


function er = get_error(state)
er = state - get_stats(state);
% the mean subtracted from all the points

% utility functions: singleton-expanded addition and multiplication
function c = pp(a,b)
c = bsxfun(@plus,a,b);


function [h,v,fpa,s, fcl, fcd,frho,fhs] = get_states(x)
global scale n_samples v0 n_params

iv = 1;
ih = iv + 1:n_samples+1;
ig = n_samples + ih;
is = n_samples + ig;

h = x(ih,:)*scale(1); % m
v = v0 - x(iv,:)*scale(2); % m/s
fpa = x(ig,:)*scale(3);
s = x(is,:)*scale(4); % m
if nargout > 4
    icl = n_samples + is;
    icd = n_samples + icl;
    irho = n_samples + icd;
    fcl = x(icl,:);
    fcd = x(icd,:);
    if n_params == 2
        frho = 1;
    else
        frho = x(irho,:);
    end
end


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

disp(['hf = ',num2str(hm),' +/- 3*',num2str(hv^0.5),' km (3s low = ',num2str(hm-3*hv^0.5) ,')'])
disp(['sf = ',num2str(sm),' +/- 3*', num2str(sv^0.5),' km'])

function Kopt = optimize_gain(x0, u, guess, W)
options = optimoptions('fminunc');
% options = optimoptions('fmincon');

% options.Algorithm = 'trust-region';
options.MaxFunctionEvaluations = 1000;

N = length(u);
obj = @(K) gain_opt_obj(x0, u, K, N, W);

Kopt = fminunc(obj, guess, options);

function out = gain_opt_obj(x0, u, K, N,W)

x = forward_pass(x0, u, K, N);
[h,v,fpa,s,g,L,D] = entry_accels(x(:,end));
[~,std] = get_stats(s(:,end));
[hm,hstd] = get_stats(h(:,end));
out = -hm + W(2)*std^0.5 + W(1)*hstd^0.5; % match the actual cost fun
% out = (std^0.5 + 0.02*hstd^0.5)/1000;

function x = forward_pass(x0, u, K, N)
x = x0(:);

for i = 1:N-1
    x(:,i+1) = entry_dynamics(x(:,i), u(i), K);
end