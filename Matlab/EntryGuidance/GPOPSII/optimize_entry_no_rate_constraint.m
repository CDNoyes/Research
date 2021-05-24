function traj = optimize_entry_no_rate_constraint(DR, CR)
if nargin == 0
    DR = 350;
    CR = 0;
end
dtr = pi/180;

% System Models:
mars = Mars();
vm = VehicleModel(5000);


% Initial states:
%[radius long lat velocity fpa heading]
% x0 = [3540e3; 0*dtr; 0*dtr; 5505; -14.5*dtr; 0*dtr]';
x0 = [54.5e3 + 3396.2e3; 0*dtr; 0*dtr; 5525; -11.5*dtr; 0*dtr]';

% Target info:
target.DR = DR;
target.CR = CR;
[target.lon,target.lat] = FinalLatLon(x0(2), x0(3), x0(6),target.DR,target.CR);

% auxdata.ref = ref;
auxdata.target = target;
auxdata.planet = mars;
auxdata.vehicle = vm;
auxdata.dimension.state = 6;
auxdata.dimension.control = 1;
auxdata.dimension.parameter = 0;
auxdata.delta.CD = 0;


% Bounds
t0 = 0;
tfl = 100;
tfu = 500;

bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfl;
bounds.phase.finaltime.upper = tfu;

lb = [mars.radiusEquatorial + -20e3
    -pi/2
    -pi/2
    300
    -45*dtr %-25 degrees, FPA
    -pi/2      % Heading
    ]';   %

ub = [x0(1) + 20e3
    0.4  % like 1500 km on Mars 
    pi/2
    x0(4)+500
    30*dtr % FPA constraint degrees
    pi/2
    ]';

heading_max = 30. * dtr;
fpa_min = -40 * dtr;
vel_max = 460;
if 1 % Fixed lat/lon
    xfl = [lb(1), target.lon, target.lat, lb(4), fpa_min, -heading_max];
    xfu = [ub(1), target.lon, target.lat, vel_max, ub(5), heading_max];
    
elseif 0 % Fixed lon i.e. downrange
    xfl = [lb(1), target.lon, lb(3), lb(4), lb(5), lb(6)];
    xfu = [ub(1), target.lon, ub(3), vel_max, ub(5), ub(6)];
    
elseif 0 % Fixed lat i.e. crossrange
    xfl = [lb(1:2), target.lat, lb(4), lb(5), lb(6)];
    xfu = [ub(1:2), target.lat, vel_max, ub(5), ub(6)];
    
else % Free lat/lon
    xfl = lb;
    xfu = [ub(1), ub(2:3), vel_max, ub(5), ub(6)];
end


bounds.phase.initialstate.lower =[x0 ];
bounds.phase.initialstate.upper = [x0];
bounds.phase.state.lower = [lb];
bounds.phase.state.upper = [ub];
bounds.phase.finalstate.lower = [xfl];
bounds.phase.finalstate.upper = [xfu];

bank_limit = 90*dtr;
bounds.phase.control.lower = -bank_limit;
bounds.phase.control.upper = bank_limit;
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 5000;


% Initial Guess
tGuess = [t0;0.5*(tfu+tfl)];
guess.phase.time = tGuess;
guess.phase.state = [x0; [3396200, 0.27, 0, vel_max, -0.26, 0]];
guess.phase.control = zeros(2,1);
guess.phase.integral = 0;

% Initial Mesh
% meshphase.colpoints = 4*ones(1,10);
% meshphase.fraction = .1*ones(1,10);
setup.name = 'OptimalEntry-GPOPS2';
setup.functions.continuous = @EntryConstraints;
setup.functions.endpoint = @EntryCost;
setup.auxdata = auxdata;
% setup.mesh.phase = meshphase;
setup.bounds = bounds;

setup.guess = guess;
setup.nlp.solver = 'snopt';
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'first';
setup.scales.method = 'automatic-guessUpdate';
setup.method = 'RPM-Differentiation';
setup.mesh.method = 'hp-PattersonRao'; % 'hp-PattersonRao'
setup.mesh.tolerance = 1e-5; % default 1e-3
% setup.mesh.colpointsmin = 3;
% setup.mesh.colpointsmax = 10;
setup.mesh.maxiterations = 30;
setup.displaylevel = 1; % default 2
sol = gpops2(setup);
Sol = sol.result.solution.phase;

traj = TrajectorySummary(Sol.time,Sol.state,Sol.control,target.DR,target.CR);
EntryPlots(traj);

end