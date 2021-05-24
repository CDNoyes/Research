function sol = optimize_edl(efpa, min_alt, input_guess)

%%% efpa is in degrees
%   Set efpa is [] for free variable
%   Set min_alt to low value to essentially eliminate it
%
%

opt_efpa = isempty(efpa);
if opt_efpa
    efpa = -15; % just need any initial guess really
end


dtr = pi/180;

target_altitude = 0;
min_altitude_srp = min_alt;

signed_bank_angle = 1;

DR = 710;

% System Models:
mars = Mars();
m0 = 7200;
vm = VehicleModel(m0);

% Initial states:
%[radius long lat velocity fpa heading]
x0 = [3522e3; 0*dtr; 0*dtr; 5500; efpa*dtr; 0*dtr]';
% m0 = vm.mass;

% Target info:
target.DR = DR;
target.CR = 0.0;
[target.lon, target.lat] = FinalLatLon(x0(2), x0(3), x0(6), target.DR, target.CR);
target.altitude = target_altitude;

auxdata.target = target;
auxdata.planet = mars;
auxdata.vehicle = vm;

auxdata.dimension.state = 7;
auxdata.dimension.control = 1;
auxdata.dimension.parameter = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase 1 - Entry
% Bank angle rate as control variable
phase = 1;

% Bounds
t0 = 0;
tfl = 0;
tfu = 400;

bounds.phase(1).initialtime.lower = t0;
bounds.phase(1).initialtime.upper = t0;
bounds.phase(1).finaltime.lower = tfl;
bounds.phase(1).finaltime.upper = tfu;

lb = [mars.radiusEquatorial + target_altitude + min_altitude_srp
    -pi/2
    -pi/2
    100
    -35*pi/180 % FPA
    -pi/2      % Heading
    -pi/2*signed_bank_angle]';   % Bank angle

ub = [x0(1)+100
    pi
    pi/2
    x0(4)+500
    6*dtr %5 degrees
    pi
    pi/2]';


% Free lat/lon
xfl = lb;
xfu = [ub(1), ub(2:3), 700, ub(5), ub(6), ub(7)];

bounds.phase(phase).initialstate.lower =[x0, lb(7)];
bounds.phase(phase).initialstate.upper = [x0, ub(7)];

if 0 % fixed initial bank angle
    bank = 90*dtr;
    bounds.phase(phase).initialstate.lower(7) = bank;
    bounds.phase(phase).initialstate.upper(7) = bank;
end

if opt_efpa % free EFPA
    bounds.phase(phase).initialstate.lower(5) = -25*pi/180;
    bounds.phase(phase).initialstate.upper(5) = 0;
end
bounds.phase(phase).state.lower = lb;
bounds.phase(phase).state.upper = ub;
bounds.phase(phase).finalstate.lower = lb;
bounds.phase(phase).finalstate.upper = ub;
if 0  % enforce zero latitude at ignition
    bounds.phase(phase).finalstate.lower(3) = 0;
    bounds.phase(phase).finalstate.upper(3) = 0;
end

bank_angle_control = 1;

if ~bank_angle_control
    bank_rate = 100*dtr;
    bounds.phase(phase).control.lower = -bank_rate;
    bounds.phase(phase).control.upper = bank_rate;
else
%     bank_max = 90;
%     bounds.phase(phase).control.lower = -bank_max;
%     bounds.phase(phase).control.upper = bank_max;
    bounds.phase(phase).control.lower = 0;
    bounds.phase(phase).control.upper = 1;
end

bounds.eventgroup(phase).lower = zeros(1,6);
bounds.eventgroup(phase).upper = zeros(1,6);


% Initial Guess
guess.phase(phase).time = [0; 240];
guess.phase(phase).state = [x0,-pi/4; xfl(1)+3000, target.lon-.02, target.lat, 500, -12*pi/180, 0*pi/180, m0];
guess.phase(phase).control = zeros(2,1);

if bank_angle_control
  % strip seventh state in bounds and guess
  bounds.phase(1).state.lower = bounds.phase(1).state.lower(1:6);
  bounds.phase(1).state.upper = bounds.phase(1).state.upper(1:6);
  bounds.phase(1).initialstate.lower = bounds.phase(1).initialstate.lower(1:6);
  bounds.phase(1).initialstate.upper = bounds.phase(1).initialstate.upper(1:6);
  bounds.phase(1).finalstate.lower = bounds.phase(1).finalstate.lower(1:6);
  bounds.phase(1).finalstate.upper = bounds.phase(1).finalstate.upper(1:6);
  guess.phase(phase).state = guess.phase(phase).state(:,1:6);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase 2 - SRP in Cartesian Coords
phase=2;

amin = 8;
amax = 15;
Tmin = amin*m0;
Tmax = amax*m0;
disp(['T/W = ', num2str(amax/3.71) ])

t0 = 0;
tf = 0;
tmin = tf;

if tf == 0
    tmax = (m0-1)*9.81*300/Tmin; % At this point the vehicle is nearly massless if using the minimal thrust
else
    tmax = tf;
end

xmin = -45e3;
xmax = 45e3;
ymin = xmin;
ymax = xmax;
zmin = 0;
zmax = 100e3;

umin = -1e3;
umax = 1e3;
vmin = umin;
vmax = umax;
wmin = -1e3;
wmax = 100;

mmin = 10; % dry mass, should be set to >0
mmax = m0;

xf = 0;
yf = 0;
zf = 0;
uf = 0;
vf = 0;
wfmin = 0;
wfmax = 0;

pinpoint = 0; % Used to initially make the problem easier, just null the velocity optimally

bounds.phase(phase).initialtime.lower    = t0;
bounds.phase(phase).initialtime.upper    = t0;
bounds.phase(phase).finaltime.lower      = tmin;
bounds.phase(phase).finaltime.upper      = tmax;

bounds.phase(phase).initialstate.lower          = [xmin ymin zmin umin 0 wmin m0];
bounds.phase(phase).initialstate.upper          = [xmax ymax zmax umax 0 wmax m0];

bounds.phase(phase).state.lower          = [xmin ymin zmin umin vmin wmin mmin];
bounds.phase(phase).state.upper          = [xmax ymax zmax umax vmax wmax mmax];
if pinpoint
    bounds.phase(phase).finalstate.lower    = [xf yf zf uf vf wfmin mmin];
    bounds.phase(phase).finalstate.upper    = [xf yf zf uf vf wfmax mmax];
else
    bounds.phase(phase).finalstate.lower    = [xmin ymin zf uf vf wfmin mmin];
    bounds.phase(phase).finalstate.upper    = [xmax ymax zf uf vf wfmax mmax];
end

bounds.phase(phase).control.lower = [Tmin, -1,-1,-1];
bounds.phase(phase).control.upper = [Tmax,  1, 1, 1];

bounds.phase(phase).path.lower = 1;
bounds.phase(phase).path.upper = 1;

guess.phase(phase).time      = [0; tmax/3];
guess.phase(phase).state     = [0, 0, min_altitude_srp, -450, 0, -130, m0; -10000, 0, 0, 0, 0, 0, 0.75*m0];
guess.phase(phase).control   = [Tmax 1 0 0; Tmin 0 0 -1];

if ~isempty(input_guess)
    guess = input_guess;
end

if bank_angle_control
  guess.phase(1).state = guess.phase(1).state(:,1:6);
end

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
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'first';
setup.scales.method = 'automatic-guessUpdate';
setup.method = 'RPM-Differentiation';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-5; % default 1e-3
setup.mesh.colpointsmin = 3;
setup.mesh.colpointsmax = 10;
setup.mesh.maxiterations = 30;
setup.displaylevel = 1; % default 2
sol = gpops2(setup);
sol.result.solution.auxdata = auxdata;
Sol = sol.result.solution.phase(1);

% L = Tx'*[costate_srp(3); costate_srp(6)];

% traj = TrajectorySummary(Sol.time,Sol.state,Sol.state(:,7),target.DR,target.CR);
% traj.sigma_dot = Sol.control; % left in radians
% EntryPlots(traj);

% srp = sol.result.solution.phase(2);
% SRPPlots(srp, 20);

end

function SRPPlots(sol, n0)
mass = sol.state(:, 7); %getMass(sol);

disp(['Prop used: ',num2str(mass(1)-mass(end)), ' kg.'])

% Plots
figure(n0+1)
plot(sol.time,sol.state(:,1:3))
legend('x','y','z')
title('Positions')
set(gcf,'name','Position vs Time','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(n0+2)
plot(sol.time,sol.state(:,4:6))
legend('u','v','w')
title('Velocities')
set(gcf,'name','Velocity vs Time','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(n0+3)
plot(sol.time,sol.control(:,2:4))
hold all
plot(sol.time, -sol.costate(:,4:6)./sum(sol.costate(:,4:6).^2, 2).^0.5,'--');
legend('u_x','u_y','u_z', '\lambda_u','\lambda_v','\lambda_w')
title('Control Directions')
set(gcf,'name','Control dirs vs Time','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(n0+4)
legText = {'T'};
plot(sol.time,sol.control(:,1))
hold all
legend(legText{:})
title('Thrust Magnitude (N)')
set(gcf,'name','Thrust Arcs','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(n0+5)
plot(sol.time, mass)
title('Mass vs Time')
set(gcf,'name','Mass','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(n0+6)
plot(sol.time,sol.control(:,1)./mass)
title('Thrust Acceleration (m/s^2)')
set(gcf,'name','Thrust Acceleration','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(n0+7)
plot(sol.time,sol.costate)
legend('x','y','z','u','v','w','m','Location','best')
title('Costates')
set(gcf,'name','Costates','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(n0+8)
plot3(sol.state(:,1),sol.state(:,2),sol.state(:,3),'o-')
hold on
plot3(sol.state(:,1), sol.state(:,2), 0*sol.state(:,3),'k--') % Ground track
title('Trajectory')
axis('equal')
set(gcf,'name','3D Trajectory','numbertitle','off')
set(gcf,'WindowStyle','docked')

end