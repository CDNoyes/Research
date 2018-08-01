%% GPOPSII, 3DOF Optimal Solution
clear; clc; close all;
[lineSpecs,textSpecs,figSpecs] = PlotSpecs();

dtr = pi/180;

% DR = linspace(740,1030,20);
% DR = 885;
DR = 815;

for i = 1:length(DR)

% System Models:
mars = Mars();
vm = VehicleModel();

% Initial states:
%[radius long lat velocity fpa heading]
x0 = [3540e3; 0*dtr; 0*dtr; 5505; -14.15*dtr; 0*dtr]';

% Target info:
target.DR = DR(i);
target.CR = 0;
[target.lon,target.lat] = FinalLatLon(x0(2),x0(3),x0(6),target.DR,target.CR);

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
tfl = 0;
tfu = 700;

bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfl;
bounds.phase.finaltime.upper = tfu;

lb = [mars.radiusEquatorial+0e3
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
    15*dtr %5 degrees
    pi
    pi/2]';

tol = 0*dtr;
heading_max = 5.5 * dtr;
fpa_min = -25 * dtr;
vel_max = 470;
if 1 % Fixed lat/lon
    xfl = [lb(1), target.lon, target.lat, lb(4), fpa_min, -heading_max, lb(7)];
    xfu = [ub(1), target.lon, target.lat, vel_max, ub(5), heading_max, ub(7)];
elseif 0 % Fixed lon i.e. downrange
    xfl = [lb(1), target.lon, lb(3), lb(4), lb(5), lb(6), lb(7)];
    xfu = [ub(1), target.lon, ub(3), 470, ub(5), ub(6), ub(7)];
elseif 0 % Fixed lat i.e. crossrange
    xfl = [lb(1:2), target.lat, lb(4), lb(5), lb(6), lb(7)];
    xfu = [ub(1:2), target.lat, vel_max, ub(5), ub(6), ub(7)];    
else % Free lat/lon
    xfl = lb;
    xfu = [ub(1), ub(2:3), 470, ub(5), ub(6), ub(7)]; 
end
bounds.phase.initialstate.lower =[x0, lb(7)];
bounds.phase.initialstate.upper = [x0, ub(7)];
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
tGuess = [t0;0.5*(tfu+tfl)];
guess.phase.time = tGuess;
guess.phase.state = [[x0,-pi/4];0.5*(xfu+xfl)];
guess.phase.control = ones(2,1);
guess.phase.integral = 0;

% Initial Mesh
meshphase.colpoints = 4*ones(1,10);
meshphase.fraction = .1*ones(1,10);
setup.name = 'OptimalEntry-GPOPS2';
setup.functions.continuous = @EntryConstraints;
setup.functions.endpoint = @EntryCost;
setup.auxdata = auxdata;
setup.mesh.phase = meshphase;
setup.bounds = bounds;
if i>1
    setup.guess = sol.result.solution;
else
    setup.guess = guess;
end
setup.nlp.solver = 'snopt';
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'first';
setup.scales.method = 'automatic-bounds';
setup.method = 'RPM-Differentiation';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-6; % default 1e-3
setup.mesh.colpointsmin = 3;
setup.mesh.colpointsmax = 10;
setup.mesh.maxiterations = 20;
setup.displaylevel = 1; % default 2
disp(['Solving DR = ',num2str(DR(i))]);
sol = gpops2(setup);
Sol = sol.result.solution.phase;

solutions{i} = Sol;
end
traj = TrajectorySummary(Sol.time,Sol.state,Sol.state(:,7),target.DR,target.CR);
EntryPlots(traj)
figure
grid on
box on
set(gcf,'name','Bank Rate', figSpecs{:})
plot(Sol.time, Sol.control/dtr)
ylabel('Bank Rate (deg/s)', textSpecs{:})
xlabel('Time (s)', textSpecs{:})
% p = polyfit(traj.energy_norm',traj.D,5);
% D_ = polyval(p,linspace(0,1));
% figure
% plot(traj.energy_norm,traj.D)
% hold on
% plot(linspace(0,1),D_);
% 
% vsplit = 2700;
% vi = traj.state(:,4) < 5400;
% vii = traj.state(:,4) > vsplit;
% vi = vi & vii;
% vj = traj.state(:,4) <= vsplit;
% sf = 500;
% pi = polyfit(traj.state(vi,4)'/sf,traj.D(vi),4);
% pj = polyfit(traj.state(vj,4)'/sf,traj.D(vj),7);
% Di = polyval(pi,linspace(vsplit,5400)/sf);
% Dj = polyval(pj,linspace(470,vsplit)/sf);
% figure
% plot(traj.state(:,4),traj.D)
% hold on
% plot(linspace(vsplit,5400),Di);
% plot(linspace(470,vsplit),Dj);

% vd = traj.D > 1;
% p = polyfit(traj.DR(vd)'/1000, traj.D(vd),2);
% D_ = polyval(p,linspace(min(traj.DR(vd)),max(traj.DR(vd)))/1000);
% figure
% plot(traj.DR(vd),traj.D(vd))
% hold on
% plot(linspace(min(traj.DR(vd)),max(traj.DR(vd))),D_)
return
%% Plot the reachable set solutions
load RCHB.mat;
DR = linspace(740,1030,20);
CR = 0;

for i = 1:1:(length(solutions)-4)
    Sol = solutions{i};
    target.DR = DR(i);
    target.CR = 0;
    traj = TrajectorySummary(Sol.time,Sol.state,Sol.state(:,7),target.DR,target.CR);
    EntryPlots(traj)
end

% figure
% plot(Sol.time, Sol.control/dtr)
return
%% 4 Phase Parametrized Version

clear; clc; close all;

dtr = pi/180;
nPhases = 4;
auxdata.nPhases = nPhases;
% DR = linspace(740,1030,20);
DR = 850;
for i = 1:length(DR)

% System Models:
mars = Mars();
vm = VehicleModel();

% Initial states:
%[radius long lat velocity fpa heading]
x0 = [3540e3; 0*dtr; 0*dtr; 5505; -14.15*dtr; 0*dtr]';

% Target info:
target.DR = DR(i);
target.CR = 0;
[target.lon,target.lat] = FinalLatLon(x0(2),x0(3),x0(6),target.DR,target.CR);

auxdata.target = target;
auxdata.planet = mars;
auxdata.vehicle = vm;
auxdata.dimension.state = 6;
auxdata.dimension.control = 0;
auxdata.dimension.parameter = 3;
auxdata.delta.CD = 0;


% Bounds
t0 = 0;
tfl = 0;
tfu = 400;



lb = [mars.radiusEquatorial+6e3
    -pi/2
    -pi/2
    0
    -45*pi/180 %-25 degrees, FPA
    -pi/2      % Heading
    ]';   
ub = [x0(1)+100
    pi
    pi/2
    x0(4)+500
    15*dtr %5 degrees
    pi/2
    ]';

if 0 % Fixed lat/lon
    xfl = [lb(1), target.lon, target.lat, lb(4), lb(5), lb(6)];
    xfu = [ub(1), target.lon, target.lat, 440, ub(5), ub(6)];
elseif 0 % Fixed lon i.e. downrange
    xfl = [lb(1), target.lon, lb(3), lb(4), lb(5), lb(6)];
    xfu = [ub(1), target.lon, ub(3), 440, ub(5), ub(6)];
elseif 1 % Fixed lat i.e. crossrange
    xfl = [lb(1:2), target.lat, lb(4), lb(5), lb(6)];
    xfu = [ub(1:2), target.lat, 440, ub(5), ub(6)];    
else % Free lat/lon
    xfl = lb;
    xfu = [ub(1), ub(2:3), 470, ub(5), ub(6)]; 
end
for p = 1:nPhases
    bounds.phase(p).initialtime.lower = t0;
    bounds.phase(p).initialtime.upper = t0;
    bounds.phase(p).finaltime.lower = tfl;
    bounds.phase(p).finaltime.upper = tfu;
    bounds.phase(p).initialstate.lower =[lb];
    bounds.phase(p).initialstate.upper = [ub];
    bounds.phase(p).state.lower = [lb];
    bounds.phase(p).state.upper = [ub];
    bounds.phase(p).finalstate.lower = [lb];
    bounds.phase(p).finalstate.upper = [ub];
    if p<nPhases
        bounds.eventgroup(p).lower = zeros(1,8);
        bounds.eventgroup(p).upper = zeros(1,8);
    end

end
bounds.phase(1).initialstate.lower =[x0];
bounds.phase(1).initialstate.upper = [x0];
bounds.phase(nPhases).finalstate.lower = [xfl];
bounds.phase(nPhases).finalstate.upper = [xfu];
bounds.parameter.lower = [0,0,0];
bounds.parameter.upper = [90,150,350];

% Initial Guess
guess.phase(1).time = [0;50];
guess.phase(2).time = [50;120];
guess.phase(3).time = [120;160];
guess.phase(4).time = [160;250];

a = 0.7;b = 1-a;
guess.phase(1).state = [x0;a*x0+b*0.5*(xfu+xfl)];
c = 0.45;d = 1-c;
guess.phase(2).state = [a*x0+b*0.5*(xfu+xfl);c*x0+d*0.5*(xfu+xfl)];
a = 0.2;b=1-a;
guess.phase(3).state = [c*x0+d*0.5*(xfu+xfl);a*x0+b*0.5*(xfu+xfl)];
guess.phase(4).state = [a*x0+b*0.5*(xfu+xfl);0.5*(xfu+xfl)];
xf_guess = [mars.radiusEquatorial+6e3, target.lon,target.lat,470,-0.26,0];

for j = 1:6
    sol_guess(1:5,j) = linspace(x0(j),xf_guess(j),5);
end
for j = 1:nPhases
    guess.phase(j).state = sol_guess(j:j+1,:);
end

guess.parameter = [50, 120, 160];

% Initial Mesh
meshphase.colpoints = 4*ones(1,10);
meshphase.fraction = .1*ones(1,10);
setup.name = 'OptimalEntry-GPOPS2';
setup.functions.continuous = @EntryConstraintsMultiphase;
setup.functions.endpoint = @EntryCostMultiphase;
setup.auxdata = auxdata;
setup.mesh.phase(1:nPhases) = meshphase;
setup.bounds = bounds;
if i>1
    setup.guess = sol.result.solution;
else
    setup.guess = guess;
end
setup.nlp.solver = 'snopt';
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.dependencies = 'sparse';
setup.derivatives.derivativelevel = 'first';
setup.scales.method = 'automatic-bounds';
setup.method = 'RPM-Differentiation';
setup.mesh.method = 'hp-LiuRao'; %'hp-PattersonRao'
setup.mesh.tolerance = 1e-3; % default 1e-3
setup.mesh.colpointsmin = 3;
setup.mesh.colpointsmax = 10;
setup.mesh.maxiterations = 20;
setup.displaylevel = 2; % default 2
disp(['Solving DR = ',num2str(DR(i))]);
sol = gpops2(setup);
Sol = sol.result.solution;
mysol.time = vertcat(Sol.phase(:).time);
mysol.state = vertcat(Sol.phase(:).state);
mysol.parameter = Sol.parameter;

solutions{i} = Sol;
end

%%
% Sol = solutions{i};
% DR = linspace(740,1030,20);
target.DR = DR(i);
traj = TrajectorySummary(mysol.time,mysol.state,zeros(size(mysol.state(:,5))),target.DR,target.CR);
EntryPlots(traj)