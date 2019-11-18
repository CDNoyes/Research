function output = covariance_min()
% Computes the optimal solution to variance minimization of scalar
% nonlinear system
if 1
    close all;
    
end
clc;

global P0 CLOSED_LOOP Pc
P0 = 0.0069; % 0.0069= +/-0.25 3-sigma
CLOSED_LOOP = true;
Pc = (0.1/3)^2;

output = SolveOCP();

sol = output.result.solution.phase(1);
STM = sol.state(:,2);
STM2 = sol.state(:,3);
P = STM.^2 * P0 + STM2.^2*Pc;
std = sqrt(P);

D = Dynamics(output.result.solution);
dx = D.dynamics(:,1);
dstm = D.dynamics(:,2);

H = sol.costate(:,1).*dx + sol.costate(:,2).*dstm;

% u_theory = -3*sign(sol.costate(:,1));
% u_theory(abs(sol.costate(:,1))< 0.1) = 0;

% Plots
figure(1)
hold all
plot(sol.time,sol.state(:,1))
plot(sol.time,sol.state(:,1)+3*std,'k')
plot(sol.time, 1.8*(sol.time+1e-10)./(sol.time+1e-10), 'r--', 'linewidth', 3)
plot(sol.time,sol.state(:,1)-3*std,'k')
xlabel('Time')
ylabel('State')
legend('x', 'x\pm 3\sigma','constraint')
% title('States')
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
plot(sol.time, D.u)
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

figure(5)
hold all
plot(sol.time, H)
title('Hamiltonian')
set(gcf,'name','H vs Time','numbertitle','off')
set(gcf,'WindowStyle','docked')
if CLOSED_LOOP
    figure(6)
    hold all
    plot(sol.time,D.gains)
    title('Gain')
    set(gcf,'name','Gain vs Time','numbertitle','off')
    set(gcf,'WindowStyle','docked')
end
output = output.result.solution;

% save("closed_loop_constrained.mat", "sol")

end

function output = SolveOCP()
global CLOSED_LOOP

t0 = 0;
tmin = 2;
tmax = 2;

x0 = 1;
stm0 = 1;

xmin = -1.8;
xmax = 1.8;

xf = 0;
umax =  3;

%% Bounds Definition
iphase = 1;

bounds.phase(iphase).initialtime.lower    = t0;
bounds.phase(iphase).initialtime.upper    = t0;
bounds.phase(iphase).finaltime.lower      = tmin;
bounds.phase(iphase).finaltime.upper      = tmax;

bounds.phase(iphase).initialstate.lower          = [x0 stm0, 0];
bounds.phase(iphase).initialstate.upper          = [x0 stm0, 0];

bounds.phase(iphase).state.lower          = [xmin -1000, -1000];
bounds.phase(iphase).state.upper          = [xmax 1000, 1000];

bounds.phase(iphase).finalstate.lower          = [xf -1000, -1000];
bounds.phase(iphase).finalstate.upper          = [xf 1000, 1000];

bounds.phase(iphase).control.lower = [-umax];
bounds.phase(iphase).control.upper = [umax];

if CLOSED_LOOP
    bounds.phase(iphase).path.lower = [-1.8, -1.8, -umax, -umax];
    bounds.phase(iphase).path.upper = [1.8, 1.8, umax, umax];
end
% bounds.phase(iphase).integral.lower = 0;
% bounds.phase(iphase).integral.upper = tmax*amax^2;

%% Guess
if 1
guess.phase(iphase).time           = [0; tmax]; % 0.5*(tmin+tmax)
guess.phase(iphase).state     = [x0, 1, 0; xf, 1e-3, 0];
guess.phase(iphase).control   = [3; -3];
% guess.phase(iphase).integral = [ 0.5*(amin+amax)*15];
else
   load open_loop_constrained.mat 
   guess.phase = sol;
end

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
setup.mesh.tolerance = 1e-8; % Default 1e-3
setup.mesh.maxiterations = 30; % Default 10
setup.method = 'RPM-Differentiation'; % 'RPM-Differentiation' or 'RPM-Integration'
setup.displaylevel = 1; % 0 = no output. 1 = only mesh refinement. 2 = NLP solver iteration output and mesh refinement

output = gpops2(setup);

end

function output = Cost(input)
global P0;
stm = input.phase(1).finalstate(2);
Pf = stm*P0*stm;
output.objective = Pf;
% output.objective = input.phase(1).finaltime;

end

function output = Dynamics(input)
global P0 Pc CLOSED_LOOP

t = input.phase(1).time;
s = input.phase(1).state;
control1 = input.phase(1).control;

% Variables
x     = s(:,1);
stm   = s(:,2);
stm2   = s(:,3);

alpha = 0.1;
dx = -alpha*dabs(x).*x + control1;
A = -2*alpha*dabs(x);
B = 1;
dstm = A.*stm;
dstm2 = A.*stm2 -dabs(x).*x;

output.u = control1;

if CLOSED_LOOP
    
    Qf = 10;
    Q = 0;
    R = 0.1;
    K = t*0;
    S = Qf;
    dt = diff(t);
%     for i = length(t)-1:-1:1
%         if false % Finite horizon, discrete implementation 
%             Ad = 1 + dt(i)*A(i);
%             Bd = B*dt(i);
%             K(i) = -Bd*S*Ad/(R + Bd*S*Bd);
%             S = Ad^2*S + (Ad*S*Bd)*K(i) + Q;
%         else
%             if isnan(A(i)) || isnan(t(i))
%                 K(i) = nan;
%             else
%                 K(i) = -lqr(A(i), B, Qf, R);
%             end
%         end
%     end
%     K(end) = K(end-1);
            K = -10;
%         K = 0; % This allows for open loop covariance, still respecting state bounds
    dstm = (A+B*K).*stm;
    dstm2 = (A+B*K).*stm2 -dabs(x).*x;
    P = stm.^2 * P0 + stm2.^2*Pc;
    std = sqrt(P);
    delta = 3*std;
    u1 = control1 + K.*delta;
    u2 = control1 + -K.*delta;
    output(1).path = [x+delta, x-delta, u1, u2];
    output(1).gains = K;
    output(1).std = delta;
    output(1).u = [u1, u2];
end

output(1).dynamics = [dx dstm dstm2];

end

function y= dabs(x)
y = sqrt(x.^2 + 1e-9);
end
