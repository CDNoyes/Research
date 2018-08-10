function traj = optimal_entry_guidance(x0, bank, lon_final, lat_final)

dtr = pi/180;

% System Models:
mars = Mars();
vm = VehicleModel();

% Initial states:
%[radius long lat velocity fpa heading]
x0 = cell2mat(x0);
% disp(x0)

% Target info:
target.lon = lon_final;
target.lat = lat_final;
[target.DR,target.CR] = Range(x0(2), x0(3), x0(6), lon_final, lat_final);
% disp(target)

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
tfl = 1;
tfu = 500;

bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfl;
bounds.phase.finaltime.upper = tfu;

lb = [mars.radiusEquatorial+1e3
    -pi/2
    -pi/2
    300
    -45*pi/180 %-25 degrees, FPA
    -pi/2      % Heading
    -pi/2]';   % Bank angle

ub = [x0(1)+100
    pi
    pi/2
    x0(4)+500
    10*dtr % FPA constraint degrees
    pi
    pi/2]';

heading_max = 5. * dtr;
fpa_min = -45 * dtr;
vel_max = 475;

xfl = [lb(1), target.lon, target.lat, lb(4), fpa_min, -heading_max, lb(7)];
xfu = [lb(1)+20e3, target.lon, target.lat, vel_max, ub(5), heading_max, ub(7)];
    

bounds.phase.initialstate.lower =[x0, bank];
bounds.phase.initialstate.upper = [x0, bank];
bounds.phase.state.lower = [lb];
bounds.phase.state.upper = [ub];
bounds.phase.finalstate.lower = [xfl];
bounds.phase.finalstate.upper = [xfu];

bank_rate = 20*dtr;
bounds.phase.control.lower = -bank_rate;
bounds.phase.control.upper = bank_rate;
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 5000;

% Initial Guess
if 1
    load('Guess.mat')

    for i = 1:length(guess.phase.time)
       if guess.phase.state(i,4) <= x0(4)
           break
       end
    end
    disp(['Starting from index = ',num2str(i)])
    guess.phase.state = guess.phase.state(i:end,:);
    guess.phase.time = guess.phase.time(i:end);
    guess.phase.time = guess.phase.time - guess.phase.time(1);
    guess.phase.control = guess.phase.control(i:end);

else
    tGuess = [t0;0.5*(tfu+tfl)];
    guess.phase.time = tGuess;
    guess.phase.state = [[x0,bank];0.5*(xfu+xfl)];
    guess.phase.control = ones(2,1);
    guess.phase.integral = 0;
end

% Initial Mesh
meshphase.colpoints = 4*ones(1,10);
meshphase.fraction = .1*ones(1,10);
setup.name = 'OptimalEntry-GPOPS2';
setup.functions.continuous = @EntryConstraints;
setup.functions.endpoint = @EntryCost;
setup.auxdata = auxdata;
setup.mesh.phase = meshphase;
setup.bounds = bounds;

setup.guess = guess;
setup.nlp.solver = 'snopt';
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'first';
setup.scales.method = 'automatic-guessUpdate';
setup.method = 'RPM-Differentiation';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-3; % default 1e-3
setup.mesh.colpointsmin = 3;
setup.mesh.colpointsmax = 10;
setup.mesh.maxiterations = 20;
setup.displaylevel = 0; % default 2
sol = gpops2(setup);
% guess = sol.result.solution;
% save("Guess.mat", 'guess');
Sol = sol.result.solution.phase;
% traj = {};
traj = TrajectorySummary(Sol.time,Sol.state,Sol.state(:,7),target.DR,target.CR);
traj.sigma_dot = Sol.control; % left in radians

end