function traj = optimize_edl()

dtr = pi/180;

DR = 900;

% System Models:
mars = Mars();
vm = VehicleModel();

% Initial states:
%[radius long lat velocity fpa heading]
x0 = [3540e3; 0*dtr; 0*dtr; 5505; -14.15*dtr; 0*dtr]';
m0 = vm.mass;

% Target info:
target.DR = DR;
target.CR = 0.0;
[target.lon,target.lat] = FinalLatLon(x0(2),x0(3),x0(6),target.DR,target.CR);

% auxdata.ref = ref;
auxdata.target = target;
auxdata.planet = mars;
auxdata.vehicle = vm;
auxdata.dimension.state = 7;
auxdata.dimension.control = 1;
auxdata.dimension.parameter = 0;
auxdata.delta.CD = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase 1 - Entry 
% Bank angle rate as control variable 
phase = 1;

% Bounds
t0 = 0;
tfl = 0;
tfu = 700;

bounds.phase(1).initialtime.lower = t0;
bounds.phase(1).initialtime.upper = t0;
bounds.phase(1).finaltime.lower = tfl;
bounds.phase(1).finaltime.upper = tfu;

lb = [mars.radiusEquatorial+0e3
    -pi/2
    -pi/2
    100
    -45*pi/180 %-25 degrees, FPA
    -pi/2      % Heading
    -pi/2]';   % Bank angle

ub = [x0(1)+100
    pi
    pi/2
    x0(4)+500
    15*dtr %5 degrees
    pi
    pi/2]';

    
% Free lat/lon
xfl = lb;
xfu = [ub(1), ub(2:3), 470, ub(5), ub(6), ub(7)]; 

bounds.phase(phase).initialstate.lower =[x0, lb(7)];
bounds.phase(phase).initialstate.upper = [x0, ub(7)];
bounds.phase(phase).state.lower = lb;
bounds.phase(phase).state.upper = ub;
bounds.phase(phase).finalstate.lower = lb;
bounds.phase(phase).finalstate.upper = ub;

bank_rate = 20*dtr;
bounds.phase(phase).control.lower = -bank_rate;
bounds.phase(phase).control.upper = bank_rate;

bounds.eventgroup(phase).lower = zeros(1,7);
bounds.eventgroup(phase).upper = zeros(1,7);


% Initial Guess
guess.phase(phase).time = [0; 240];
guess.phase(phase).state = [x0,-pi/4; xfl(1)+3000, target.lon-.02, target.lat, 500, -12*pi/180, 0*pi/180, m0];
guess.phase(phase).control = zeros(2,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase 2 - SRP 
phase=2;


% Bounds
t0 = 0;
tfl = 0;
tfu = 700;

bounds.phase(phase).initialtime.lower = t0;
bounds.phase(phase).initialtime.upper = tfu;
bounds.phase(phase).finaltime.lower = tfl;
bounds.phase(phase).finaltime.upper = tfu;

lb = [mars.radiusEquatorial+0e3
    -pi/2
    -pi/2
    0
    -90*pi/180 %-25 degrees, FPA
    -pi/2      % Heading
    0]';       % Mass

ub = [x0(1)+100
    pi
    pi/2
    900
    15*dtr      %5 degrees
    pi
    m0]';

    
bounds.phase(phase).initialstate.lower =lb;
bounds.phase(phase).initialstate.upper = ub;
bounds.phase(phase).state.lower = lb;
bounds.phase(phase).state.upper = ub;
bounds.phase(phase).finalstate.lower = [lb(1), target.lon-0.1, target.lat-0.02, 0, -pi/2, -pi/2, 0];
bounds.phase(phase).finalstate.upper = [lb(1), target.lon+0.1, target.lat+0.002, 1, 0, pi/2, m0];

bounds.phase(phase).control.lower = [vm.thrust*0.1, 0];
bounds.phase(phase).control.upper = [vm.thrust, 2*pi];

guess.phase(phase).time = [260; 280];
guess.phase(phase).state = [xfl(1)+3000, target.lon-.05, target.lat, 500, -12*pi/180, 0*pi/180, m0
                            xfl(1), target.lon, target.lat, 1, -90*pi/180, 0*pi/180, m0*0.7];
guess.phase(phase).control = [vm.thrust, pi; vm.thrust, pi];

% End Phase Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial Mesh
% meshphase.colpoints = 4*ones(1,10);
% meshphase.fraction = .1*ones(1,10);
setup.name = 'OptimalEDL';
setup.functions.continuous = @EDLConstraints;
setup.functions.endpoint = @EDLPhaseLink;
setup.auxdata = auxdata;
% setup.mesh.phase(1) = meshphase;
% setup.mesh.phase(2) = meshphase;

setup.bounds = bounds;

setup.guess = guess;
setup.nlp.solver = 'snopt';
setup.derivatives.supplier = 'sparseFD';
setup.derivatives.derivativelevel = 'first';
setup.scales.method = 'automatic-guessUpdate';
setup.method = 'RPM-Differentiation';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-2; % default 1e-3
setup.mesh.colpointsmin = 3;
setup.mesh.colpointsmax = 10;
setup.mesh.maxiterations = 50;
setup.displaylevel = 1; % default 2
sol = gpops2(setup);
Sol = sol.result.solution.phase(1);

traj = TrajectorySummary(Sol.time,Sol.state,Sol.state(:,7),target.DR,target.CR);
% traj.sigma_dot = Sol.control; % left in radians 
EntryPlots(traj);
end