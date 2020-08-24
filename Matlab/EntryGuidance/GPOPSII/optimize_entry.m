function traj = optimize_entry()

dtr = pi/180;

glimit = 20; % Earth gs, set to [] for no limit 
minimum_altitude = 1000; % meters relative to MOLA, can be negative,  set to [] for no limit 

mass = 7200;

% System Models:
mars = Mars();
vm = VehicleModel(mass);


% Initial states:
%[radius long lat velocity fpa heading]
x0 = [3524e3; 0*dtr; 0*dtr; 5500; -16*dtr; 0*dtr]';


DR = 0;
CR = 0;

% Target info:
target.DR = DR;
target.CR = CR;
[target.lon,target.lat] = FinalLatLon(x0(2), x0(3), x0(6),target.DR,target.CR);

% auxdata.ref = ref;
auxdata.target = target;
auxdata.planet = mars;
auxdata.vehicle = vm;
auxdata.dimension.state = 7;
auxdata.dimension.control = 1;
auxdata.dimension.parameter = 0;
auxdata.delta.CD = 0;


% Bounds
t0 = 0;
tfl = 40;
tfu = 500;

bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfl;
bounds.phase.finaltime.upper = tfu;

lb = [mars.radiusEquatorial + -20e3
    -pi/2
    -pi/2
    300
    -45*pi/180 %-25 degrees, FPA
    -pi/2      % Heading
    -90*pi/180]';   % Bank angle

ub = [x0(1) + 20e3
    0.4  % like 1500 km on Mars
    pi/2
    x0(4)+500
    10*dtr % FPA constraint degrees
    pi/2
    90*pi/180]';

% heading_max = 2. * dtr;
% fpa_min = -14.5 * dtr;
vel_max = 500;
update_range = 0;
if 0 % Fixed lat/lon
    xfl = [lb(1), target.lon, target.lat, lb(4), fpa_min, -heading_max, lb(7)];
    xfu = [ub(1), target.lon, target.lat, vel_max, ub(5), heading_max, ub(7)];
    
elseif 0 % Fixed lon i.e. downrange
    xfl = [lb(1), target.lon, lb(3), lb(4), lb(5), lb(6), lb(7)];
    xfu = [ub(1), target.lon, ub(3), vel_max, ub(5), ub(6), ub(7)];
    
elseif 0 % Fixed lat i.e. crossrange
    xfl = [lb(1:2), target.lat, lb(4), lb(5), lb(6), lb(7)];
    xfu = [ub(1:2), target.lat, vel_max, ub(5), ub(6), ub(7)];
    update_range = 1;
    
else % Free lat/lon
    xfl = lb;
    xfu = [ub(1), ub(2:3), vel_max, ub(5), ub(6), ub(7)];
    update_range = 1;
end

x0l = x0;
x0u = x0;
if 1    % This allows the optimizer to choose the initial FPA
    x0l(5) = -25*dtr;
    x0u(5) = -10*dtr;
end
bounds.phase.initialstate.lower =[x0l, lb(7)];
bounds.phase.initialstate.upper = [x0u, ub(7)];
bounds.phase.state.lower = [lb];
bounds.phase.state.upper = [ub];
bounds.phase.finalstate.lower = [xfl];
bounds.phase.finalstate.upper = [xfu];

bank_rate = 20*dtr;
bounds.phase.control.lower = -bank_rate;
bounds.phase.control.upper = bank_rate;
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 5000;
if glimit
    bounds.phase.path.lower = 0;
    bounds.phase.path.upper = glimit;
end
% Initial Guess
tGuess = [t0;0.5*(tfu+tfl)];
guess.phase.time = tGuess;
guess.phase.state = [[x0,-pi/4]; [3396200, 0.27, 0, vel_max, -0.25, 0, 0]];
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
setup.mesh.tolerance = 1e-6; % default 1e-3
% setup.mesh.colpointsmin = 3;
% setup.mesh.colpointsmax = 10;
setup.mesh.maxiterations = 20;
setup.displaylevel = 1; % default 2
sol = gpops2(setup);
Sol = sol.result.solution.phase;

if update_range
    [target.DR, target.CR] = Range(x0(2), x0(3), x0(6),Sol.state(end,2),Sol.state(end,3));
    target.lon = Sol.state(end,2);
end

traj = TrajectorySummary(Sol.time,Sol.state,Sol.state(:,7),target.DR,target.CR);
traj.sigma_dot = Sol.control; % left in radians

end