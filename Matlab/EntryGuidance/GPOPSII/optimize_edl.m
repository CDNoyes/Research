function sol = optimize_edl()

dtr = pi/180;

target_altitude = 0;
min_altitude_srp = 3000;

DR = 710;

% System Models:
mars = Mars();
m0 = 7200;
vm = VehicleModel(m0);

% Initial states:
%[radius long lat velocity fpa heading]
x0 = [3522e3; 0*dtr; 0*dtr; 5500; -14.5*dtr; 0*dtr]';
% m0 = vm.mass;

% Target info:
target.DR = DR;
target.CR = 0.0;
[target.lon,target.lat] = FinalLatLon(x0(2),x0(3),x0(6),target.DR,target.CR);
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
    400
    -35*pi/180 % FPA
    -pi/2      % Heading
    -pi/2]';   % Bank angle

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
bounds.phase(phase).state.lower = lb;
bounds.phase(phase).state.upper = ub;
bounds.phase(phase).finalstate.lower = lb;
bounds.phase(phase).finalstate.upper = ub;

bank_rate = 20*dtr;
bounds.phase(phase).control.lower = -bank_rate;
bounds.phase(phase).control.upper = bank_rate;

bounds.eventgroup(phase).lower = zeros(1,6);
bounds.eventgroup(phase).upper = zeros(1,6);


% Initial Guess
guess.phase(phase).time = [0; 240];
guess.phase(phase).state = [x0,-pi/4; xfl(1)+3000, target.lon-.02, target.lat, 500, -12*pi/180, 0*pi/180, m0];
guess.phase(phase).control = zeros(2,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Phase 2 - SRP in Entry Coords
% phase=2;
%
%
% % Bounds
% t0 = 0;
% tfl = 0;
% tfu = 700;
%
% bounds.phase(phase).initialtime.lower = t0;
% bounds.phase(phase).initialtime.upper = tfu;
% bounds.phase(phase).finaltime.lower = tfl;
% bounds.phase(phase).finaltime.upper = tfu;
%
% lb = [mars.radiusEquatorial+0e3
%     -pi/2
%     -pi/2
%     0
%     -90*pi/180 %-25 degrees, FPA
%     -pi/2      % Heading
%     1]';       % Mass
%
% ub = [x0(1)+100
%     pi
%     pi/2
%     900
%     15*dtr      %5 degrees
%     pi
%     m0]';
%
%
% bounds.phase(phase).initialstate.lower =lb;
% bounds.phase(phase).initialstate.upper = ub;
% bounds.phase(phase).state.lower = lb;
% bounds.phase(phase).state.upper = ub;
% bounds.phase(phase).finalstate.lower = [lb(1), target.lon-0.5, target.lat-0.2, 1, -pi/2, -pi/2, 1];
% bounds.phase(phase).finalstate.upper = [lb(1), target.lon+0.5, target.lat+0.2, 5, 0, pi/2, m0];
%
% bounds.phase(phase).control.lower = [vm.thrust*0.4, 0];
% bounds.phase(phase).control.upper = [vm.thrust, 2*pi];
%
% guess.phase(phase).time = [240; 280];
% guess.phase(phase).state = [xfl(1)+3000, target.lon-.05, target.lat, 500, -12*pi/180, 0*pi/180, m0
%                             xfl(1), target.lon, target.lat, 1, -90*pi/180, 0*pi/180, m0*0.7];
% guess.phase(phase).control = [vm.thrust, pi; vm.thrust, pi];
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
setup.mesh.tolerance = 1e-3; % default 1e-3
setup.mesh.colpointsmin = 3;
setup.mesh.colpointsmax = 10;
setup.mesh.maxiterations = 30;
setup.displaylevel = 1; % default 2
sol = gpops2(setup);
Sol = sol.result.solution.phase(1);

traj = TrajectorySummary(Sol.time,Sol.state,Sol.state(:,7),target.DR,target.CR);
% traj.sigma_dot = Sol.control; % left in radians
EntryPlots(traj);

srp = sol.result.solution.phase(2);
SRPPlots(srp, 20);

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