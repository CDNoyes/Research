function Sol = optimize_entry_experimental(c)
%%% Longitudinal dynamics only, but with experimental stuff 


dtr = pi/180;

% System Models:
mars = Mars();
vm = VehicleModel(7200);


% Target info:
DR = 600e3;
CR = 0;
target.DR = DR;
target.CR = CR;
[target.lon,target.lat] = FinalLatLon(0, 0, 0,target.DR,target.CR);


% Initial states:
%[altitude range_to_go velocity fpa]
x0 = [127.8e3; DR; 5505; -14.5*dtr]';

auxdata.target = target;
auxdata.planet = mars;
auxdata.vehicle = vm;
auxdata.dimension.state = 4;
auxdata.dimension.control = 1;
auxdata.dimension.parameter = 0;
auxdata.c = c;


% Bounds
t0 = 0;
tfl = 200;
tfu = 500;

bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfl;
bounds.phase.finaltime.upper = tfu;

lb = [-10e3
    0
    300
    -35*dtr %-25 degrees, FPA
    ]';   %

ub = [x0(1) + 1
    1000e3
    x0(3)+500
    10*dtr % FPA constraint degrees
    ]';

vel_max = 750;
if 0 % Fixed  downrange
    xfl = [lb(1), 0, lb(3), lb(4)];
    xfu = [ub(1), 0, vel_max, ub(4)];
    x0l = x0;
    x0l(2) = DR;
    x0u = x0;
    x0u(2) = DR; 
    
else % Free downrange
    x0l = x0;
    x0l(2) = 300e3;
    x0u = x0;
    x0u(2) = ub(2); % probably no longer than 1200 km eh?
    xfl = lb;
    xfu = [ub(1), 0, vel_max, ub(4)];
end

bounds.phase.initialstate.lower =[x0l];
bounds.phase.initialstate.upper = [x0u];
bounds.phase.state.lower = [lb];
bounds.phase.state.upper = [ub];
bounds.phase.finalstate.lower = [xfl];
bounds.phase.finalstate.upper = [xfu];

bank_limit = 90*dtr;
u_limit = cos(bank_limit);
bounds.phase.control.lower = u_limit;
bounds.phase.control.upper = 1;
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 5000;


% Initial Guess
tGuess = [t0; 0.5*(tfu+tfl)];
guess.phase.time = tGuess;
guess.phase.state = [x0; [0e3, 0, vel_max, -15*pi/180]];
guess.phase.control = 1*ones(2,1);
guess.phase.integral = 0;

% Initial Mesh
% meshphase.colpoints = 4*ones(1,10);
% meshphase.fraction = .1*ones(1,10);
setup.name = 'OptimalEntry-GPOPS2';
setup.functions.continuous = @EoM;
setup.functions.endpoint = @Cost;
setup.auxdata = auxdata;
% setup.mesh.phase = meshphase;
setup.bounds = bounds;

setup.guess = guess;
setup.nlp.solver = 'snopt';
setup.derivatives.supplier = 'sparseFD';
setup.derivatives.derivativelevel = 'first';
setup.scales.method = 'automatic-guessUpdate';
setup.method = 'RPM-Differentiation';
setup.mesh.method = 'hp-PattersonRao'; % 'hp-PattersonRao'
setup.mesh.tolerance = 1e-4; % default 1e-3
% setup.mesh.colpointsmin = 3;
% setup.mesh.colpointsmax = 10;
setup.mesh.maxiterations = 20;
setup.displaylevel = 1; % default 2
sol = gpops2(setup);
Sol = sol.result.solution.phase;

x = Sol.state';
Lon = cumtrapz(Sol.time, x(3,:).*cos(x(4,:))./(3396.2e3 + x(1,:)));
X = [Sol.state(:,1)+3396.2e3, Lon', Sol.state(:,2)*0, Sol.state(:,3:4), Sol.state(:,1)*0];

if Sol.state(1,2) ~= DR
    target.DR = Sol.state(1,2);
end

traj = TrajectorySummary(Sol.time, X, acos(Sol.control),target.DR/1000,target.CR);
EntryPlots(traj);
figure
plot(x(2,1)-x(2,:), Sol.control)
xlabel('DR')

end

function output = EoM( input )

S       = input.auxdata.vehicle.area;
m       = input.auxdata.vehicle.mass;
mu      = input.auxdata.planet.mu;
rp      = input.auxdata.planet.radiusEquatorial;

x = input.phase.state;
% nx = size(x,2);

u = input.phase.control;

h = x(:,1);
r = h + rp;
v = x(:,3);
gamma = x(:,4);

g = mu./(r.^2);
[rho,a] = MarsAtmosphericDensity(h);
M = v./a;
[C_D, C_L] = input.auxdata.vehicle.aerodynamics(M);
q = 0.5*rho.*v.^2*S/m;
drag = q.*C_D;
lift = q.*C_L;


dr = v.*sin(gamma);
ds = -v.*cos(gamma);
dv = -drag - g.*sin(gamma);
dgamma = lift.*u.^2./v - (g./v - v./r).*cos(gamma);
output.dynamics = [dr ds dv dgamma];
% output.integrand = zeros(size(r));
output.integrand = u.^2;



end

function output = Cost( input )
xf = input.phase.finalstate;
% tf = input.phase.finaltime;

% output.objective = (xf(1)-8000)^2 + input.auxdata.c*input.phase.integral;
output.objective = (7.5-xf(1)/1000)^2 + 1*(xf(3)-550)^2 + 0.1*input.phase.integral;
end