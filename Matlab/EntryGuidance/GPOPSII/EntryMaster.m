%% GPOPSII, 3DOF Optimal Solution
clear; clc; close all;

dtr = pi/180;

% Reference Trajectory
% load('EntryGuidance/HighElevationPlanner/Trajectories/ReferenceTrajectory_DR780_CR0.mat');




% System Models:
mars = Mars();
vm = VehicleModel();

% Initial states:
%[radius long lat velocity fpa heading]
x0 = [3540e3; -90.07*dtr; -43.90*dtr; 5505; -14.15*dtr; 4.99*dtr]';

% Target info:
target.DR = 780;
target.CR = 0;
[target.lon,target.lat] = FinalLatLon(x0(2),x0(3),x0(6),target.DR,target.CR);

% auxdata.ref = ref;
auxdata.target = target;
auxdata.planet = mars;
auxdata.vehicle = vm;
auxdata.dimension.state = 7;
auxdata.dimension.control = 1;
auxdata.dimension.parameter = 0;
auxdata.delta.CD = -.14;


% Bounds
t0 = 0;
tfl = 0;
tfu = 400;

bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfl;
bounds.phase.finaltime.upper = tfu;

lb = [mars.radiusEquatorial+6e3
   -pi
    -pi/2
    300
    -25*pi/180 %-25 degrees
    -pi/2
    -pi/2]';
ub = [x0(1)+100
    0
    0
    x0(4)+500
    5*dtr %5 degrees
    pi/2
    pi/2]';
tol = 0*dtr;
xfl = [lb(1), target.lon, target.lat, lb(4), lb(5), lb(6), lb(7)];
xfu = [ub(1), target.lon, target.lat, 475, ub(5), ub(6), ub(7)];
bounds.phase.initialstate.lower =[x0, lb(7)];
bounds.phase.initialstate.upper = [x0, ub(7)];
bounds.phase.state.lower = [lb];
bounds.phase.state.upper = [ub];
bounds.phase.finalstate.lower = [xfl];
bounds.phase.finalstate.upper = [xfu];

bounds.phase.control.lower = -20*dtr;
bounds.phase.control.upper = 20*dtr;
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 500;

% Initial Guess
tGuess = [t0;0.5*(tfu+tfl)];
guess.phase.time = tGuess;
guess.phase.state = [[x0,0];0.5*(xfu+xfl)];
guess.phase.control = ones(2,1);
guess.phase.integral = 0;

% Initial Mesh
meshphase.colpoints = 4*ones(1,10);
meshphase.fraction = .1*ones(1,10);
setup.name = 'OptimalDragTracking-GPOPS2';
setup.functions.continuous = @EntryConstraints;
setup.functions.endpoint = @EntryCost;
setup.auxdata = auxdata;
setup.mesh.phase = meshphase;
setup.bounds = bounds;
setup.guess = guess;
setup.nlp.solver = 'snopt';
% setup.nlp.solver = 'ipopt';
% setup.nlp.options.ipopt.linear_solver = 'ma57';
%setup.nlp.options.ipopt.linear_solver = 'mumps';
setup.derivatives.supplier = 'sparseFD';
setup.derivatives.derivativelevel = 'first';
setup.scales.method = 'automatic-bounds';
setup.method = 'RPMintegration';
setup.mesh.method = 'hp1';
setup.mesh.tolerance = 1e-4; % default 1e-3
setup.mesh.colpointsmin = 3;
setup.mesh.colpointsmax = 10;
setup.mesh.maxiterations = 20;

sol = gpops2(setup);
Sol = sol.result.solution.phase;

%% Plot the solution
traj = TrajectorySummary(Sol.time,Sol.state,Sol.state(:,7),target.DR,target.CR);
EntryPlots(traj)

figure
plot(Sol.time, Sol.control/dtr)