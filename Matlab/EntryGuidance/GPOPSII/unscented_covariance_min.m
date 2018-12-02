function output = unscented_covariance_min()
% Computes the optimal solution to variance minimization of scalar
% nonlinear system
if 1
    close all;
    
end
clc;

global P0 CLOSED_LOOP
P0 = 0.0069; % 0.0069 = +/-0.25 3-sigma
CLOSED_LOOP = false;

output = SolveOCP();

sol = output.result.solution.phase(1);

D = Dynamics(output.result.solution);
std = D.std;
P = std.^2;

% H = sol.costate(:,1).*dx + sol.costate(:,2).*dstm;

% u_theory = -3*sign(sol.costate(:,1));
% u_theory(abs(sol.costate(:,1))< 0.1) = 0;

% Plots
figure(1)
hold all
plot(sol.time,sol.state(:,1))
plot(sol.time, D.mean)
plot(sol.time,sol.state(:,2),'k--')
plot(sol.time,sol.state(:,3),'k--')
plot(sol.time, 1.8*ones(size(sol.time)))
legend('Nominal', 'Mean', 'sigma pts')
title('States')
set(gcf,'name','Position vs Time','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(2)
hold all
plot(sol.time,P)
% legend('P')
title('Variance')
set(gcf,'name','Variance vs Time','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(3)
hold all
plot(sol.time,sol.control)
% plot(sol.time, u_theory)
% legend('Numerical','Theory')
title('Control')
set(gcf,'name','Control vs Time','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(4)
hold all
plot(sol.time,sol.costate)
legend('x','stm')
title('Costates')
set(gcf,'name','Costates vs Time','numbertitle','off')
set(gcf,'WindowStyle','docked')

% figure(5)
% hold all
% plot(sol.time, H)
% title('Hamiltonian')
% set(gcf,'name','H vs Time','numbertitle','off')
% set(gcf,'WindowStyle','docked')
% if CLOSED_LOOP
%     figure(6)
%     hold all
%     plot(sol.time,D.gains)
%     title('Gain')
%     set(gcf,'name','Gain vs Time','numbertitle','off')
%     set(gcf,'WindowStyle','docked')
% end
output = output.result.solution;

end

function output = SolveOCP()
global CLOSED_LOOP P0

t0 = 0;
tmin = 1;
tmax = 1;

x0 = 1;
n = 3;
x01 = x0 +n*0.1175;
x02 = x0 -n*0.1175;

xmin = -10;
xmax = 1.8;

xf = 0;
umax =  3;

%% Bounds Definition
iphase = 1;

bounds.phase(iphase).initialtime.lower    = t0;
bounds.phase(iphase).initialtime.upper    = t0;
bounds.phase(iphase).finaltime.lower      = tmin;
bounds.phase(iphase).finaltime.upper      = tmax;

bounds.phase(iphase).initialstate.lower          = [x0 x01 x02];
bounds.phase(iphase).initialstate.upper          = [x0 x01 x02];

bounds.phase(iphase).state.lower          = [xmin xmin xmin];
bounds.phase(iphase).state.upper          = [xmax xmax xmax];

bounds.phase(iphase).finalstate.lower          = [xf xmin xmin];
bounds.phase(iphase).finalstate.upper          = [xf xmax xmax];

bounds.phase(iphase).control.lower = [-umax];
bounds.phase(iphase).control.upper = [umax];

    bounds.phase(iphase).path.lower = [-umax, -umax];
    bounds.phase(iphase).path.upper = [umax, umax];

% bounds.phase(iphase).integral.lower = 0;
% bounds.phase(iphase).integral.upper = tmax*amax^2;

%% Guess
guess.phase(iphase).time      = [0; tmax]; % 0.5*(tmin+tmax)
guess.phase(iphase).state     = [x0, x01 x02; xf, xf, xf];
guess.phase(iphase).control   = [3; -3];

%% NLP Parameters and Optimization Call
setup.name = 'VarMin';
setup.functions.continuous = @Dynamics;
setup.functions.endpoint = @Cost;
setup.nlp.solver = 'snopt';
if strcmp(setup.nlp.solver,'snopt')
    setup.nlp.snoptoptions.tolerance = 1e-6; % Default 1e-6
    setup.nlp.snoptoptions.maxiterations = 2000; % Default 2000
elseif strcmp(setup.nlp.solver,'ipopt')
    setup.nlp.ipoptoptions.linear_solver = 'mumps'; % Default 'mumps' ('ma57')
    setup.nlp.ipoptoptions.tolerance = 1e-6; % Default 1e-7
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
xf = input.phase(1).finalstate;
m = 0.5*xf(1) + 0.25*(xf(2)+xf(3));

delta = [0.5,0.25,0.25].*(xf-m).^2;
Pf = sum(delta);
output.objective = Pf;
% output.objective = input.phase(1).finaltime;

end

function output = Dynamics(input)
global CLOSED_LOOP

t = input.phase(1).time;
s = input.phase(1).state;
control1 = input.phase(1).control;

% Variables
x  = s(:,1);
x1 = s(:,2);
x2 = s(:,3);
m = 0.5*x + 0.25*(x1+x2);

delta1 = x1-m;
delta2 = x2-m;
P = 0.5*(delta1.^2 + delta2.^2); % Variance is relative to mean, not nominal 

c = 0.25;
dx = -c*dabs(x).*x + control1;


if CLOSED_LOOP
    K = -1;
else
    K = 0;
end
dx1 = -c*dabs(x1).*x1 + control1 + K*(x1-x); % Feedback is relative to nomina traj, not mean 
dx2 = -c*dabs(x2).*x2 + control1 + K*(x2-x);
u1 = control1 + K*(x1-x);
u2 = control1 + K*(x2-x);

output(1).path = [u1, u2];
output(1).std = sqrt(P);
output(1).mean = m;
output(1).dynamics = [dx dx1 dx2];

end

function y= dabs(x)
y = sqrt(x.^2 + 1e-9);
end
