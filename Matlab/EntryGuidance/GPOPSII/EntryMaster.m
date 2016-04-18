%% GPOPSII, 3DOF Optimal Drag Tracking With Reversals
clear; clc;

% Reference Trajectory
load('EntryGuidance/HighElevationPlanner/Trajectories/ReferenceTrajectory_DR780_CR0.mat');

% System Models:
mars = Mars();
vm = VehicleModel();

% Initial states:
x0 = ref.state(1,1:6);

auxdata.ref = ref;
auxdata.planet = mars;
auxdata.vehicle = vm;
auxdata.dimension.state = 7;
auxdata.dimension.control = 1;
auxdata.dimension.parameter = 0;

dtr = pi/180;

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
    0 %0 degrees
    pi/2
    pi/2]';
tol = 0*dtr;
xfl = [lb(1), lb(2), lb(3), lb(4), lb(5), lb(6), lb(7)];
xfu = [ub(1), ub(2), ub(3), 485, ub(5), ub(6), ub(7)];
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
% setup.nlp.options.ipopt.linear_solver = 'ma57';
%setup.nlp.options.ipopt.linear_solver = 'mumps';
setup.derivatives.supplier = 'sparseFD';
setup.derivatives.derivativelevel = 'first';
setup.scales.method = 'automatic-bounds';
setup.method = 'RPMintegration';
setup.mesh.method = 'hp1';
setup.mesh.tolerance = 1e-2; % default 1e-3
setup.mesh.colpointsmin = 3;
setup.mesh.colpointsmax = 10;
setup.mesh.maxiterations = 20;

sol = gpops2(setup);

%% Plot the solution

