function EntryPredictorCorrector()
tic;
% Reference Data:
load('EntryGuidance/HighElevationPlanner/Trajectories/ReferenceTrajectory_DR780_CR0.mat');

% Flags
CONSTANTBANK = 1;
SCALE = 0;
DENSE = 1; % Whether to keep data at every guidance cycle or every integration step taken

% Compute the initial range to go:
S = TrajLength(ref);
S0 = S(end);

% System Models:
mars = Mars();
vm = VehicleModel();
sf = getScaleFactors(mars);

% Initial States:
x0 = [ref.state(1,1:6)';S0];
e0 = eFromState(x0(1)/sf.radius,x0(4)/sf.velocity);

% Final energy:
xf = ref.state(end,1:6)';
ef = eFromState(xf(1)/sf.radius,xf(4)/sf.velocity);

% Guidance System Specs:
f = 1; %Hz
dt = 1/f; % Guidance cycle length, seconds
tmax = 350; % while loop max iteration condition
tol = 1e-3; % Newton iteration tolerance in the predictor-corrector process

% Create a data structure for convenience
param.dt = dt;
param.planet = mars;
param.vehicle = vm;
param.sf = sf;
param.CONSTANTBANK = CONSTANTBANK;
param.SCALE = SCALE;
param.DENSE = DENSE;
param.e0 = e0;
param.ef = ef;
param.tol = tol;
param.sigmaf = 60*pi/180;

% Initialization
x = x0';
t = 0;
sigma0 = 20*pi/180;
bank = sigma0;
L = ref.L(1);
D = ref.D(1);
M = ref.M(1);
terminate = 0;
tCurrent = 0;

% The main loop occurs in the time domain so that the guidance cycle can be
% enforced simply, while the prediction/correction steps happen in the
% energy domain. EDIT: Actually the predict/correct shit can also happen in
% the time domain if we want, simply compute e and stop if e > ef.

while x(end,7) > 0 && ~terminate && tCurrent < tmax
    % Be careful if working with normalized states or not!
    
    if all(D < 0.1*mars.g0)                             % Pre-entry phase
        sigma = @(E) sigma0*ones(size(E));
        
    else                                                % Predictor-corrector phase
        
        sigma0 = abs(computeSigma0(x(end,:),sigma0,param)); %These calls should be combined
        sigma = @(E) getBankAngle(E,sigma0,param);
        % lateral control logic tbd later
    end

    [tnew,xnew,terminate,banknew,Lnew,Dnew,Mnew] = Step(x(end,:),sigma,param);
    
    % Update the solution vectors
    if param.DENSE
        if terminate
            idx = 2:length(tnew)-2;
        else
            idx = 2:length(tnew);
        end
    else
        if terminate
            idx = length(tnew)-2;
        else
            idx = length(tnew);
        end
    end
    x = [x;xnew(idx,:)]; %Augment the state
    t = [t;tnew(idx)+t(end)];
    L = [L;Lnew(idx)'];
    D = [D;Dnew(idx)'];
    M = [M;Mnew(idx)'];
    bank = [bank;banknew(idx)];
    
    tCurrent = t(end)
end
t_elapsed = toc;

% Display and analyze the results
disp(['Simulation Terminated. Time elapsed: ',num2str(t_elapsed),' s.'])
traj = TrajectorySummary(t,x,bank,ref.target.DR,ref.target.CR);
EntryPlots(traj);

end

%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = Predict(xc,sigmaFun,param)

e = @(x)  eFromState(x(:,1),x(:,4),param.sf);
[~,xn] = ode45(@(t,x) dynamics(t,x,sigmaFun(e(x')),param), [0,350],xc,[]);

z = xn(end,7);

end

function [sigma,znew,dznew] = Correct(xc,sigma0,param)
% e = @(x)  eFromState(x(:,1),x(:,4),param.sf);

%Algorithm Specs:
deltaSigma = 1e-9; % Finite difference stepsize
nMax = 5;
n = 0;

% Repeatedly call predict while changing "i" to compute the new bank angle
z1 = Predict(xc,@(e) getBankAngle(e,sigma0,param),param);
f1 = 0.5*z1^2; % Current cost
while n < nMax
    % Compute sensitivity to sigma0
    z2 = Predict(xc,@(e) getBankAngle(e,sigma0+deltaSigma,param),param);
    dz = (z2-z1)/deltaSigma;
    
    sigma = Saturate(sigma0 - (0.5^n)*z1/dz, -pi/2,pi/2); % Eqn 24 to update sigma
    
    znew = Predict(xc,@(e) getBankAngle(e,sigma,param),param);
    
    if 0.5*znew^2 < f1
        
        break;
    else
        n = n+1;
        %         sigma0 = sigma; % Don't do this
    end
    
end

% Compute the new sensitivity at the new sigma value, used for convergence check
znew_delta = Predict(xc,@(e) getBankAngle(e,sigma+deltaSigma,param),param);
dznew = (znew_delta-znew)/deltaSigma;


end

function [tn,xn,terminate,Sigma,L,D,M] = Step(xc, sigma, param,debug)
% Step method will need to output any data the user may need to switch
% phases - for example the measured Drag force.

% Note that if we scale t,r,v appropriately, and pass sf to EntryForces, we
% can perform all the integrations in dimensionless coordinates, if
% desired.

e = @(x)  eFromState(x(:,1),x(:,4),param.sf);
[tn,xn] = ode45(@(t,x) dynamics(t,x,sigma(e(x')),param), [0,param.dt],xc,[]);

Sigma = sigma(e(xn));

for i = 1:length(tn)
%     if i == 29 && debug
%         Debug = 1;
%     end
    [~,L(i),D(i),M(i),done(i)] = dynamics(tn(i),xn(i,:),Sigma(i),param);
    
end
terminate = any(done);

end

function [dX,L,D,M,done] = dynamics(t,x,sigma,param)

e = eFromState(x(1),x(4),param.sf);

if e > param.ef
    dX = zeros(size(x));
    done = true;
    L = nan;
    D = nan;
    M = nan;
else
    done = false;
    [g,L,D,~,M,a,rho] = EntryForces(x,param.planet,param.vehicle);
    dx = EntryDynamics(x,sigma,g,L,D);
    dX = [dx;-x(4)*cos(x(5))/x(1)]; % Augment the state with the range computation
end


end

function sf = getScaleFactors(planet)
if nargin == 1
    sf.time = sqrt(planet.radiusEquatorial/planet.gm);
    sf.velocity = sqrt(planet.gm*planet.radiusEquatorial);
    sf.radius = planet.radiusEquatorial;
else
    sf.time = 1;
    sf.velocity = 1;
    sf.radius = 1;
end
end

function xout = Scale(x,sf)

%x = [r,theta,phi,v,gamma,psi,s]
%xout = x dimensionless

xout = x./[sf.radius,1,1,sf.velocity,1,1,1];
end

function e = eFromState(r,v,sf)
% Computes the current energylike variable e from the dimensional radius
% and velocity
if nargin < 3 || isempty(sf)
    e = 1./r-0.5*v.^2;
else
    e =  1./(r/sf.radius)-0.5*(v/sf.velocity).^2;
end
end

function sigma0 = computeSigma0(xc,sigma0,param) % Add tons of inputs here
znew = 1;
dznew = param.tol*1.1;
maxIter = 10;
iter = 1;

while iter < maxIter && abs(znew*dznew) > param.tol
    [sigma0,znew,dznew] = Correct(xc,sigma0,param);
end

end

function sigma = getBankAngle(e,sigma0,param)


if param.CONSTANTBANK
    sigma = sigma0*ones(size(e));
else
    e0 = param.e0;
    ef = param.ef;
    sigmaf = param.sigmaf;
    sigma = sigma0 + (e-e0)/(ef-e0)*(sigmaf-sigma0);
end


end