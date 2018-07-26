function sol = optimize_srp(x0)
% Computes the optimal solution to SRP problems treating thrust
% magnitude and unit vector as the control variables
x0 = cell2mat(x0);
% load('E:\Documents\GitHub\Research\Matlab\data\srp_z3200.mat')
% guess.phase = sol;
% output = SolveOCP(x0, guess);
output = SolveOCP(x0);

sol = output.result.solution.phase(1);
% save(['./data/srp_z',num2str(int16(x0(3))),'.mat'], 'sol');

if 0
    mass = sol.state(:, 7); %getMass(sol);
    tSwitch = ThrustPlaneAnalysis(sol);
    H = hamiltonian(sol);
    [pitch,azimuth] = getAnglesFromUnitVector(sol.control(:,2:4));
    pitch = pitch*180/pi;
    azimuth = azimuth*180/pi;

    disp(['Prop used: ',num2str(mass(1)-mass(end)), ' kg.'])

    % Plots
    figure(1)
    plot(sol.time,sol.state(:,1:3))
    legend('x','y','z')
    title('Positions')
    set(gcf,'name','Position vs Time','numbertitle','off')
    set(gcf,'WindowStyle','docked')

    figure(2)
    plot(sol.time,sol.state(:,4:6))
    legend('u','v','w')
    title('Velocities')
    set(gcf,'name','Velocity vs Time','numbertitle','off')
    set(gcf,'WindowStyle','docked')

    figure(3)
    plot(sol.time,sol.control(:,2:4))
    hold all
    plot(sol.time, -sol.costate(:,4:6)./sum(sol.costate(:,4:6).^2, 2).^0.5,'--');
    legend('u_x','u_y','u_z', '\lambda_u','\lambda_v','\lambda_w')
    title('Control Directions')
    set(gcf,'name','Control dirs vs Time','numbertitle','off')
    set(gcf,'WindowStyle','docked')

    figure(4)
    legText = {'T'};
    plot(sol.time,sol.control(:,1))
    hold all
    for i = 1:length(tSwitch)
        if isreal(tSwitch(i))
            legText{end+1} = ['Theoretical Switch ',num2str(length(legText))];
            plot(tSwitch(i)*ones(size(sol.time)),linspace(min(sol.control(:,1)),max(sol.control(:,1)),length(sol.time)),'--')
        end
    end
    legend(legText{:})
    title('Thrust Magnitude (N)')
    set(gcf,'name','Thrust Arcs','numbertitle','off')
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
    legend('x','y','z','u','v','w','m','Location','best')
    title('Costates')
    set(gcf,'name','Costates','numbertitle','off')
    set(gcf,'WindowStyle','docked')

    figure(8)
    plot(sol.time,H)
    title('Hamiltonian')
    set(gcf,'name','Hamiltonian','numbertitle','off')
    set(gcf,'WindowStyle','docked')

    figure(9)
    plot(sol.time,pitch)
    hold on
    plot(sol.time,azimuth)
    legend('Pitch','Yaw','Location','best')
    title('Control Angles')
    set(gcf,'name','Control Angles','numbertitle','off')
    set(gcf,'WindowStyle','docked')


    figure(12)
    plot3(sol.state(:,1),sol.state(:,2),sol.state(:,3),'o-')
    hold on
    plot3(sol.state(:,1),sol.state(:,2),0*sol.state(:,3),'k--') % Ground track
    title('Trajectory')
    axis('equal')
    set(gcf,'name','3D Trajectory','numbertitle','off')
    set(gcf,'WindowStyle','docked')
end
end

function output = SolveOCP(X, guess)
t0 = 0;
tmin = 16;
tmax = 16;

x0 = X(1);
y0 = X(2);
z0 = X(3);
u0 = X(4);
v0 = X(5);
w0 = X(6);
m0 = X(7);

xmin = -15e3;
xmax = 15e3;
ymin = xmin;
ymax = xmax;
zmin = 0;
zmax = 100e3;

umin = -2e3;
umax = 2e3;
vmin = umin;
vmax = umax;
wmin = -3e3;
wmax = 100;

mmin = 1; % dry mass, should be set to >0
mmax = m0;

xf = 0;
yf = 0;
zf = 0;
uf = 0;
vf = 0;
wfmin = 0;
wfmax = 0;

% Thrust acceleration limits
amin = 40*0.001;
amax = 70;
Tmin = amin*m0;
Tmax = amax*m0;

%% Bounds Definition
iphase = 1;

bounds.phase(iphase).initialtime.lower    = t0;
bounds.phase(iphase).initialtime.upper    = t0;
bounds.phase(iphase).finaltime.lower      = tmin;
bounds.phase(iphase).finaltime.upper      = tmax;

bounds.phase(iphase).initialstate.lower          = [x0 y0 z0 u0 v0 w0 m0];
bounds.phase(iphase).initialstate.upper          = [x0 y0 z0 u0 v0 w0 m0];

bounds.phase(iphase).state.lower          = [xmin ymin zmin umin vmin wmin mmin];
bounds.phase(iphase).state.upper          = [xmax ymax zmax umax vmax wmax mmax];

bounds.phase(iphase).finalstate.lower          = [xf yf zf uf vf wfmin mmin];
bounds.phase(iphase).finalstate.upper          = [xf yf zf uf vf wfmax mmax];

bounds.phase(iphase).control.lower = [Tmin, -1,-1,-1];
bounds.phase(iphase).control.upper = [Tmax,  1, 1, 1];

bounds.phase(iphase).path.lower = 1;
bounds.phase(iphase).path.upper = 1;

% bounds.phase(iphase).integral.lower = 0;
% bounds.phase(iphase).integral.upper = tmax*Tmax^2;

%% Guess
if nargin ==1 || isempty(guess)
    guess.phase(iphase).time      = [0; 16];
    guess.phase(iphase).state     = [bounds.phase(iphase).initialstate.lower; bounds.phase(iphase).finalstate.lower];
    guess.phase(iphase).control   = [Tmax 1 0 0; Tmin 0 0 -1];
end
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
    setup.nlp.ipoptoptions.tolerance = 1e-6; % Default 1e-7
    setup.nlp.ipoptoptions.maxiterations = 2000; % Default 2000
end
setup.bounds = bounds;
setup.guess = guess;
setup.derivatives.supplier = 'sparseFD'; % 'sparseFD' or 'sparseBD' or 'sparseCD' or 'adigator'
setup.derivatives.derivativelevel = 'first'; % 'first' or 'second'
setup.derivatives.dependencies = 'sparseNaN'; % 'full', 'sparse' or 'sparseNaN'
setup.scales.method = 'automatic-guessUpdate'; % 'none' or 'automatic-bounds' or 'automatic-guess' or 'automatic-guessUpdate' or 'automatic-hybrid' or 'automatic-hybridUpdate' or 'defined'
setup.mesh.method = 'hp-PattersonRao'; % 'hp-PattersonRao' or 'hp-DarbyRao' or 'hp-LiuRao'
setup.mesh.tolerance = 1e-5; % Default 1e-3
setup.mesh.maxiterations = 20; % Default 10
setup.method = 'RPM-Differentiation'; % 'RPM-Differentiation' or 'RPM-Integration'
setup.displaylevel = 0; % 0 = no output. 1 = only mesh refinement. 2 = NLP solver iteration output and mesh refinement

output = gpops2(setup);

end

function output = Cost(input)

output.objective = -input.phase(1).finalstate(7);

end

function output = Dynamics(input)


%---------------------%
% Dynamics in Phase 1 %
%        SRP          %
%---------------------%
t1 = input.phase(1).time;
s1 = input.phase(1).state;
control1 = input.phase(1).control;

% Variables
x1     = s1(:,1);
y1     = s1(:,2);
z1     = s1(:,3);
u1     = s1(:,4);
v1     = s1(:,5);
w1     = s1(:,6);
m1     = s1(:,7);


T     = control1(:,1);
ux    = control1(:,2);
uy    = control1(:,3);
uz    = control1(:,4);
aT    = T./m1;

gravity1 = 3.71;
Isp = 290;
g0 = 9.81;
ve = Isp*g0;

% EOM
xdot1 = u1;
ydot1 = v1;
zdot1 = w1;
udot1 = aT.*ux;
vdot1 = aT.*uy;
wdot1 = aT.*uz - gravity1;
mdot1 = -T/ve;

% Output
% output(1).integrand = T;
output(1).dynamics = [xdot1 ydot1 zdot1 udot1 vdot1 wdot1 mdot1];
output(1).path = (sum(control1(:,2:4).^2,2)); %Norm of the control directions

end


function tSwitch = ThrustPlaneAnalysis(sol)
Isp = 290;
g0 = 9.81;
ve = Isp*g0;

lv = sol.costate(:,4:6);
% lr = sol.costate(1,1:3);
lm = sol.costate(:,7);
m = sol.state(:,7);
u = sol.control(:,2:4);

S = sign(-lm/ve + dot(lv, u, 2)./m);
s1 = S(1:end-1);
s2 = S(2:end);
iSwitch = s1.*s2 < 0;
iSwitch = [iSwitch; false];
tSwitch = sol.time(iSwitch);

end

function H = hamiltonian(sol)
T = sol.control(:,1);
CR = sol.costate(1,1:3);
gvec = [0;0;-3.7].';
aT = T./sol.state(:,7);

Isp = 290;
g0 = 9.81;
ve = Isp*g0;

for i = 1:length(T)
    H(i) = -T(i)*sol.costate(i,7)/ve + dot(CR,sol.state(i,4:6)) + dot(sol.costate(i,4:6),aT(i)*sol.control(i,2:4)+gvec);
    
end
end

function [pitch,azimuth] = getAnglesFromUnitVector(u)

for i = 1:size(u,1)
    pitch(i) = atan2(u(i,3),norm(u(i,1:2)));
    azimuth(i) = atan(u(i,2)/u(i,1));
end
end