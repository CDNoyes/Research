function EntryPredictorCorrector()
close all; clc;
tic;

dtr = pi/180;


% Reference Data:
load('EntryGuidance/HighElevationPlanner/Trajectories/ReferenceTrajectory_DR780_CR0.mat');

% Flags
CONSTANTBANK = 1; % Use a constant bank or the linear (in energy) parametrization
SCALE = 0; % Perform all calculations and integrations in dimensionless units?
DENSE = 1; % Whether to keep data at every guidance cycle or every integration step taken


% System Models:
mars = Mars();
vm = VehicleModel();
sf = getScaleFactors(mars);


% Initial States:
S0 = ref.target.DR*1e3/sf.radius; % Compute the initial range to go:
x0 = [ref.state(1,1:6)';S0];
e0 = eFromState(x0(1),x0(4),sf);

% Final energy:
% xf = ref.state(end,1:6)';
% ef = eFromState(xf(1),xf(4),sf);
ef = eFromState(10e3+sf.radius,400,sf);

% Guidance System Specs:
f = .2; % Hz
dt = 1/f; % Guidance cycle length, seconds
tmax = 350; % while loop max iteration condition
tol = 1e-5; % Newton iteration tolerance in the predictor-corrector process

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
param.sigmaf = 20*dtr;
param.phif = ref.target.lat;
param.thetaf = ref.target.lon;
param.nReversal = 0;
param.PCInit = 0;
param.CONSTANTBANKSTOP = 2500; % Stop using constant bank at this velocity
param.FDStep = .000175; % Finite Difference Step Size in Radians

% Bank Angle Dynamics, Limits, Gains, etc
param.ENFORCEBANK = 0;
sigmaMin = 0*dtr;
sigmaMax = 90*dtr;
K = 10*[0.5643, 1.2934,0]; %Bank Angle dynamics gains
param.bankGains = K;
rateMax = 20*dtr;
accMax = 10*dtr;
param.bankLimits.rate = rateMax;
param.bankLimits.acceleration = accMax;
param.bankLimits.angleMax = sigmaMax;
param.bankLimits.angleMin = sigmaMin;

% Initialization
t = 0;
sigma0 = -20*pi/180;
bank = sigma0;
if param.ENFORCEBANK
    x = [x0',bank,0];
else
    x = x0';
end
L = ref.L(1);
D = ref.D(1);
M = ref.M(1);
terminate = 0;
tCurrent = 0;
banksign = sign(sigma0); %Current sign of the bank angle
limit = 0.1; % Heading angle error limit
headingError = 0; 

% The main loop occurs in the time domain so that the guidance cycle can be
% enforced simply, while the prediction/correction steps happen in the
% energy domain. EDIT: Actually the predict/correct shit can also happen in
% the time domain if we want, simply compute e and stop if e > ef.

while ~terminate && tCurrent < tmax
    
    [~,CR] = Range(x(1,2),x(1,3),x(1,6),x(end,2),x(end,3));
    [sCurrent,lim,headError] = DeadbandLateralGuidance(x(end,:),banksign(end),sign(CR),param);
    if sCurrent ~= banksign(end) 
        param.nReversal = param.nReversal+1;
    end
    if ~any(sqrt(L.^2 + D.^2) > 0.15*mars.g0)                             % Pre-entry phase
        sigma = @(E,V) sigma0*ones(size(E));
        sCurrent = sign(sigma0);
        
    else                                                % Predictor-corrector phase
        if ~param.PCInit
            param.PCInit = 1;
            sCurrent = -sCurrent;
            param.nReversal = param.nReversal+1;
        end
        sigma0 = abs(computeSigma0(x(end,:),sigma0,param)); %These calls could be combined
        sigma = @(E,V) sCurrent*getBankAngle(E,V,sigma0,param);
    end
    
    [tnew,xnew,terminate,banknew,Lnew,Dnew,Mnew] = Step(x(end,:),sigma,param);
    
    % Update the solution vectors
    if param.DENSE
        idx = ~isnan(Dnew);
        idx(1) = false;
    else
        if terminate
            temp = ~isnan(Dnew);
            idx = find(temp,1,'last');
        else
            idx = length(tnew);
        end
    end
    x = [x;xnew(idx,:)]; %Augment the state
    t = [t;tnew(idx)+t(end)];
    L = [L;Lnew(idx)'];
    D = [D;Dnew(idx)'];
    M = [M;Mnew(idx)'];
    bank = [bank;banknew(idx)']; % Bank angle command by the PC and Lateral guidance
    banksign = [banksign;sCurrent*ones(length(tnew(idx)),1)];
    limit = [limit;lim*ones(length(tnew(idx)),1)];
    headingError = [headingError;headError*ones(length(tnew(idx)),1)];
    
    tCurrent = t(end)
end
t_elapsed = toc;

% Display and analyze the results
disp(['Simulation Terminated. Time elapsed: ',num2str(t_elapsed),' s.'])
traj = TrajectorySummary(t,x,bank,ref.target.DR,ref.target.CR);
EntryPlots(traj);
[lineSpecs,textSpecs,figSpecs] = PlotSpecs();

figure
plot(traj.state(:,4),headingError, traj.state(:,4),banksign,traj.state(:,4), traj.sigma,traj.state(:,4), sign(traj.CR), traj.state(:,4),limit,traj.state(:,4),-limit)
hold on
plot(ref.state(:,4),ref.headingError,'k--',lineSpecs{:})
xlabel('Velocity (m/s)',textSpecs{:})
ylabel('Heading Error (rad)',textSpecs{:})
grid on
box on
set(gcf,'name','Crossrange Deadband vs Velocity', figSpecs{:})
legend('Trajectory Heading Error','Commanded Sign','Bank Angle (rad)','CR','CR Deadband','CR Deadband')

if param.ENFORCEBANK
figure
plot(traj.time,bank/dtr,traj.time,Saturate(traj.state(:,8),-sigmaMax,sigmaMax)/dtr)
end

end

%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = Predict(xc,sigmaFun,param)

e = @(x)  eFromState(x(:,1),x(:,4),param.sf);
[~,xn] = ode45(@(t,x) dynamics(t,x,sigmaFun(e(x.'),x(4)),param), [0,350],xc,[]);

z = xn(end,7);

end

function [sigma,znew,dznew] = Correct(xc,sigma0,param)

%Algorithm Specs:
deltaSigma = param.FDStep; % Finite difference stepsize
nMax = 5;
n = 0;

% Repeatedly call predict while changing "i" to compute the new bank angle
z1 = Predict(xc,@(e,v) getBankAngle(e,v,sigma0,param),param);
f1 = 0.5*z1^2; % Current cost
while n < nMax
    % Compute sensitivity to sigma0
%     dz = ComplexDiff(@(SIGMA) Predict(xc,@(e,v) getBankAngle(e,v,SIGMA,param),param), sigma0);
    dz = ForwardDiff(@(SIGMA) Predict(xc,@(e,v) getBankAngle(e,v,SIGMA,param),param), sigma0,deltaSigma);

    sigma = Saturate(sigma0 - (0.5^n)*z1/dz, -pi/2,pi/2); % Eqn 24 to update sigma
    
    znew = Predict(xc,@(e,v) getBankAngle(e,v,sigma,param),param);
    
    if 0.5*znew^2 < f1
        
        break;
    else
        n = n+1;
    end
    
end

% Compute the new sensitivity at the new sigma value, used for convergence check
znew_delta = Predict(xc,@(e,v) getBankAngle(e,v,sigma+deltaSigma,param),param);
dznew = (znew_delta-znew)/deltaSigma;


end

function [tn,xn,terminate,Sigma,L,D,M] = Step(xc, sigma, param)
% Step method will need to output any data the user may need to switch
% phases - for example the measured Drag force.

% Note that if we scale t,r,v appropriately, and pass sf to EntryForces, we
% can perform all the integrations in dimensionless coordinates, if
% desired.

e = @(x)  eFromState(x(:,1),x(:,4),param.sf);
[tn,xn] = ode45(@(t,x) dynamics(t,x,sigma(e(x'),x(4)),param), [0,param.dt],xc,[]);


for i = 1:length(tn)
    Sigma(i) = sigma(e(xn(i,:)),xn(i,4));
    [~,L(i),D(i),M(i),done(i)] = dynamics(tn(i),xn(i,:),Sigma(i),param);
    
end
terminate = any(done);

end

function [dX,L,D,M,done] = dynamics(t,x,sigma,param)

% These stopping conditions should be in param
hmin = 6; %km
vmax = 480; %m/s
if (x(4) < vmax && x(7) <= 0) || (x(1)-param.planet.radiusEquatorial)/1000 < hmin %e > param.ef
    dX = zeros(size(x));
    done = true;
    L = nan;
    D = nan;
    M = nan;
else
    done = false;
    [g,L,D,~,M,a,rho] = EntryForces(x,param.planet,param.vehicle);
    if param.ENFORCEBANK
        sigma_ex = Saturate(x(8),-param.bankLimits.angleMax,param.bankLimits.angleMax);
        dx = EntryDynamics(x,sigma_ex,g,L,D);
        du = BankAngleDynamics(t,x(8:9),sigma,param.bankLimits,param.bankGains);

    else
        dx = EntryDynamics(x,sigma,g,L,D);
        du = [];
    end

    dX = [dx;-x(4)*cos(x(5))/x(1);du]; % Augment the state with the range computation and bank dynamics
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
    iter = iter + 1;
end
if iter == maxIter
    disp('Max iter reached in sigma0 computation.')
end
end

function sigma = getBankAngle(e,v,sigma0,param)


if param.CONSTANTBANK && v > param.CONSTANTBANKSTOP
    sigma = sigma0*ones(size(e));
else
    e0 = param.e0;
    ef = param.ef;
    sigmaf = param.sigmaf;
    sigma = sigma0 + (e-e0)/(ef-e0)*(sigmaf-sigma0);
end


end

function [snew,limit,e] = DeadbandLateralGuidance(state,sBank,sCR,param)

% Velocity-varying heading corridor
% lim0 = 0.08;
% limf = 0.07;
% v0 = 5505;
% vf = 2500;
lim0 = 0.1;
limf = 0.02;
v0 = 5505;
vf = 500;

M = (lim0-limf)/(v0-vf);
CR0 = lim0-M*v0;
limit = limf*(state(4)<vf && state(4)>0) + (CR0 + state(4)*M).*(state(4)>=vf);
psi_d = DesiredHeading(state(2),state(3),param.thetaf,param.phif);
e = psi_d-state(6);

if state(4) < 180 
    snew = -sCR;
else
    
    if abs(e) > limit &&  sBank == -sign(e)
        snew = -sBank;
    else
        snew = sBank;
    end
    
end

end