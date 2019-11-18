function sol = pd1d()

global k
k = 5.5; % this needs to match what's used in the dynamics
% k = 0;

% close all;

output = SolveOCP();

sol = output.result.solution.phase(1);


mass = sol.state(:, 3);
H = Hamiltonian(sol.state, sol.costate, sol.control, k);
x = sol.state;
sol.state(:,1:2) = sol.state(:,1:2)*4000/1050;


disp(['Prop used: ',num2str(mass(1)-mass(end)), ' kg.'])
disp(['ToF: ',num2str(sol.time(end)), ' s.'])

% Plots
figure(1)
plot(sol.state(:,1), sol.state(:,2))
hold all
set(gcf,'name','Phase Portrait','numbertitle','off')
set(gcf,'WindowStyle','docked')

c1 = (1-sol.costate(:,3)*k) + sol.costate(:,2)*1050./sol.state(:,3);
c2 = (1-sol.costate(:,3)*k) - sol.costate(:,2)*1050./sol.state(:,3);
    t1 = (sol.costate(1,2) + (1-k*sol.costate(1,3))*mass(1)/1050)/sol.costate(1,1);
    disp(['u=1 at t=', num2str(t1), ' s']); % only true is u(0) = 0, not -1
c3 =  sol.costate(:,2)*1050./sol.state(:,3);
figure(3)
plot(sol.time,sol.control)
hold all
plot(sol.time, c1)
plot(sol.time, c2)
plot(sol.time, c3)

title('Control')
set(gcf,'name','Control dirs vs Time','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(5)
plot(sol.time, mass)
title('Mass vs Time')
set(gcf,'name','Mass','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(6)
plot(sol.time,sol.control(:,1)./mass)
title('Thrust Acceleration (m/s^2)')
set(gcf,'name','Thrust Acceleration','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(7)
plot(sol.time,sol.costate)
hold on
plot(sol.time, sol.costate(:,1).*x(:,2)/(1.62*1050/4000),'o')
title('Costates')
set(gcf,'name','Costates','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(8)
plot(sol.time, H)
title('Hamiltonian')
set(gcf,'name','H','numbertitle','off')
set(gcf,'WindowStyle','docked')

% p = sol.costate;
% u = sol.control;
% u(abs(u)<1e-3) = 0;
% s = sign(u);
% S = (1 - p(:,3)*k).*s + p(:,2)*1050./sol.state(:,3);
% s(s==0)=-1;
% Sneg = s.* ((1 - p(:,3)*k).*s + p(:,2)*1050./sol.state(:,3));
% figure(9)
% plot(sol.time, S)
% hold all
% plot(sol.time, Sneg, '--')
% legend("True","If negative thrust were used")
% set(gcf,'name','Switch','numbertitle','off')
% set(gcf,'WindowStyle','docked')
end

function output = SolveOCP()
t0 = 0;
tmin = 5;
tmax = 140;

p = 4000;

m0 = 1050;
z0 = 4000*m0/p;
v0 = -70*m0/p;
% v0 = 120*m0/p;


zmin = 0;
zmax = 100e3;

vmin = -1000;
vmax = 100;

mmin = 200; % dry mass, should be set to >0
mmax = m0;


zf = 0;
vf = 0;


% Thrust acceleration limits
Tmax = 1;
Tmin = -1;

%% Bounds Definition
iphase = 1;

bounds.phase(iphase).initialtime.lower    = t0;
bounds.phase(iphase).initialtime.upper    = t0;
bounds.phase(iphase).finaltime.lower      = tmin;
bounds.phase(iphase).finaltime.upper      = tmax;

bounds.phase(iphase).initialstate.lower          = [z0 v0 m0];
bounds.phase(iphase).initialstate.upper          = [z0 v0 m0];

bounds.phase(iphase).state.lower          = [zmin vmin mmin];
bounds.phase(iphase).state.upper          = [zmax vmax mmax];

bounds.phase(iphase).finalstate.lower          = [zf vf mmin];
bounds.phase(iphase).finalstate.upper          = [zf vf mmax];

bounds.phase(iphase).control.lower = [Tmin];
bounds.phase(iphase).control.upper = [Tmax];

bounds.phase(iphase).integral.lower = 0;
bounds.phase(iphase).integral.upper = tmax*Tmax;

%% Guess
guess.phase(iphase).time      = [0; 70];
guess.phase(iphase).state     = [bounds.phase(iphase).initialstate.lower; bounds.phase(iphase).finalstate.lower];
guess.phase(iphase).control   = [Tmax; Tmax];
guess.phase(1).integral = tmax*Tmax;

%% NLP Parameters and Optimization Call
setup.name = 'SRP Optimization';
setup.functions.continuous = @Dynamics;
setup.functions.endpoint = @Cost;
setup.nlp.solver = 'snopt';
if strcmp(setup.nlp.solver,'snopt')
    setup.nlp.snoptoptions.tolerance = 1e-6; % Default 1e-6
    setup.nlp.snoptoptions.maxiterations = 2000; % Default 2000
elseif strcmp(setup.nlp.solver,'ipopt')
    setup.nlp.ipoptoptions.linear_solver = 'mumps'; % Default 'mumps' ('ma57')
    setup.nlp.ipoptoptions.tolerance = 1e-7; % Default 1e-7
    setup.nlp.ipoptoptions.maxiterations = 2000; % Default 2000
end
setup.bounds = bounds;
setup.guess = guess;
setup.derivatives.supplier = 'sparseCD'; % 'sparseFD' or 'sparseBD' or 'sparseCD' or 'adigator'
setup.derivatives.derivativelevel = 'first'; % 'first' or 'second'
setup.derivatives.dependencies = 'sparseNaN'; % 'full', 'sparse' or 'sparseNaN'
setup.scales.method = 'automatic-guessUpdate'; % 'none' or 'automatic-bounds' or 'automatic-guess' or 'automatic-guessUpdate' or 'automatic-hybrid' or 'automatic-hybridUpdate' or 'defined'
setup.mesh.method = 'hp-PattersonRao'; % 'hp-PattersonRao' or 'hp-DarbyRao' or 'hp-LiuRao'
setup.mesh.tolerance = 1e-5; % Default 1e-3
setup.mesh.maxiterations = 30; % Default 10
setup.method = 'RPM-Differentiation'; % 'RPM-Differentiation' or 'RPM-Integration'
setup.displaylevel = 1; % 0 = no output. 1 = only mesh refinement. 2 = NLP solver iteration output and mesh refinement

output = gpops2(setup);

end

function output = Cost(input)

% output.objective = -input.phase(1).finalstate(3);
output.objective = input.phase(1).integral;
% output.objective = input.phase(1).finaltime; % Time optimal control 
end

function output = Dynamics(input)
global k 

%---------------------%
% Dynamics in Phase 1 %
%        SRP          %
%---------------------%
t = input.phase(1).time;
s1 = input.phase(1).state;
control1 = input.phase(1).control;

% Variables
z     = s1(:,1);
v     = s1(:,2);
m     = s1(:,3);

u     = control1(:,1);

if k == 0
    aT = u;
else
    aT    = 1050.*u./m;
end
gravity1 = 1.62*1050/4000;


% EOM
zdot = v;
vdot = aT - gravity1;
mdot = -k*dabs(u);

% Output
output(1).dynamics = [zdot vdot mdot];
output(1).integrand = dabs(u);

end

function H = Hamiltonian(x, p, u, k)

H = abs(u) + p(:,1).*x(:,2) + p(:,2).*(1050./x(:,3).*u-1.62*1050/4000) - p(:,3).*k.*abs(u);
end

function xa = dabs(x)
% xa = sqrt(x.^2 + 1e-6);
xa = abs(x);
end