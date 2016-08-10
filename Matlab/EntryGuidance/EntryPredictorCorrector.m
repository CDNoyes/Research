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
sf = getScaleFactors(mars,SCALE);

% Guidance System Specs:
f = 1; % Hz
dt = 1/f./sf.time; % Guidance cycle length, seconds
tmax = 350/sf.time; % while loop max iteration condition
tol = 1e-5; % Newton iteration tolerance in the predictor-corrector process

% Create a data structure for convenience
param.dt = dt;

param.planet = mars;
param.vehicle = vm;
param.sf = sf;

param.CONSTANTBANK = CONSTANTBANK;
param.SCALE = SCALE;
param.DENSE = DENSE;

param.tol = tol;
param.sigmaf = 20*dtr;
param.phif = ref.target.lat;
param.thetaf = ref.target.lon;
param.nReversal = 0;
param.PCInit = 0;
param.CONSTANTBANKSTOP = 2500/sf.velocity; % Stop using constant bank at this velocity
param.FDStep = .000000175; % Finite Difference Step Size in Radians

% Bank Angle Dynamics, Limits, Gains, etc
param.ENFORCEBANK = 1;
sigmaMin = 10*dtr;
sigmaMax = 90*dtr;
K = 10*[0.5643, 1.2934,0]; %Bank Angle dynamics gains
param.bankGains = K;
rateMax = 20*dtr;
accMax = 10*dtr;
param.bankLimits.rate = rateMax;
param.bankLimits.acceleration = accMax;
param.bankLimits.angleMax = sigmaMax;
param.bankLimits.angleMin = sigmaMin;

% Initial States:
S0 = ref.target.DR*1e3/mars.radiusEquatorial; % Compute the initial range to go:
x0 = Scale([ref.state(1,1:6)';S0],param);

% Initial and Final Energy:
e0 = eFromState(x0(1),x0(4));
ef = eFromState((10e3+mars.radiusEquatorial)/sf.radius,400/sf.velocity); %Approximate final state, used in the linear parametrization of bank angle
param.e0 = e0;
param.ef = ef;

% Initialization
t = 0;
sigma0 = -16*pi/180;
bank = sigma0;
if param.ENFORCEBANK
    x = [x0',bank,0];
else
    x = x0';
end
L = ref.L(1)*sf.radius/sf.velocity^2;
D = ref.D(1)*sf.radius/sf.velocity^2;
M = ref.M(1);
terminate = 0;
tCurrent = 0;
banksign = sign(sigma0); %Current sign of the bank angle
limit = 0.1; % Heading angle error limit
headingError = 0;

% The main loop 

while ~terminate && tCurrent < tmax
    
    [~,CR] = Range(x(1,2),x(1,3),x(1,6),x(end,2),x(end,3));
    [sCurrent,lim,headError] = DeadbandLateralGuidance(x(end,:),banksign(end),sign(CR),param);
    if sCurrent ~= banksign(end)
        param.nReversal = param.nReversal+1;
    end
    if ~any(sqrt((L/sf.radius*sf.velocity^2).^2 + (D/sf.radius*sf.velocity^2).^2) > 0.15*mars.g0)                             % Pre-entry phase
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
    if param.DENSE                      % Take all non-nan elements
        idx = ~isnan(Dnew);
        idx(1) = false;
    else
        if  terminate
            temp = ~isnan(Dnew);
            idx = find(temp,1,'last');  % Take final non-nan element
        else
            idx = length(tnew);         % Take final element
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
    
    tCurrent = t(end);
    if SCALE || ~mod(tCurrent*param.sf.time,5)
        disp(['Current sim time: ',num2str(tCurrent*param.sf.time),' s'])
    end
end

% Compute a better endpoint:
% xf = interp1(xnew(1:2,7),[xnew(1:2,1:6),tnew(1:2)],0);
% t = [t;t(end)+xf(7)];
% if param.ENFORCEBANK
%     xf = [xf(1:6),0,x(end,8:9)];
% else
%     xf = [xf(1:6),0];
% end
% x = [x;xf];
% 
% L(end+1) = L(end);
% D(end+1) = D(end); %Lazy, for now
% M(end+1) = M(end);
% bank(end+1) = bank(end);
% banksign(end+1) = banksign(end);
% limit(end+1) = limit(end);
% headingError(end+1) = headingError(end);

%Unscale if needed:
if SCALE
    t = t*param.sf.time;
    x = Unscale(x,param);
    L = L/sf.radius*sf.velocity^2;
    D = D/sf.radius*sf.velocity^2;
end
t_elapsed = toc;

% Display and analyze the results
disp(['Simulation Terminated. Time elapsed: ',num2str(t_elapsed),' s.'])
traj = TrajectorySummary(t,x,bank,ref.target.DR,ref.target.CR);
traj.DR = ref.target.DR-traj.state(:,7).*traj.state(:,1)/1000;
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
    set(gcf,'name','Commanded and Flown', figSpecs{:})

end

end

%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = Predict(xc,sigmaFun,param)

e = @(x)  eFromState(x(:,1),x(:,4));
[~,xn] = ode45(@(t,x) dynamics(t,x,sigmaFun(e(x.'),x(4)),param,false), [0,350],xc,[]);

z = xn(end,7);

end

function [sigma,znew,dznew] = Correct(xc,sigma0,param)

%Algorithm Specs:
deltaSigma = param.FDStep; % Finite difference stepsize
nMax = 4;
n = 0;

% Repeatedly call predict while changing "i" to compute the new bank angle
z1 = Predict(xc,@(e,v) getBankAngle(e,v,sigma0,param),param);
f1 = 0.5*z1^2; % Current cost

% Compute sensitivity to sigma0
dz = ForwardDiff(@(SIGMA) Predict(xc,@(e,v) getBankAngle(e,v,SIGMA,param),param), sigma0,deltaSigma);
while n < nMax
    
    
    sigma = Saturate(sigma0 - (0.5^n)*z1/dz, -param.bankLimits.angleMax, param.bankLimits.angleMax); % Eqn 24 to update sigma
    
    znew = Predict(xc,@(e,v) getBankAngle(e,v,sigma,param),param);
    
    if 0.5*znew^2 < f1
        
        break;
    else
        n = n+1;
%         if n == nMax
%             disp('Max minor iterations reached.')
%         end
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

e = @(x)  eFromState(x(:,1),x(:,4));
[tn,xn] = ode45(@(t,x) dynamics(t,x,sigma(e(x'),x(4)),param,true), [0,param.dt],xc,[]);

for i = 1:length(tn)
    Sigma(i) = sigma(e(xn(i,:)),xn(i,4)); % sigma command 
    [~,L(i),D(i),M(i),done(i)] = dynamics(tn(i),xn(i,:),Sigma(i),param,false);
    
end
terminate = any(done);

end

function [dX,L,D,M,done] = dynamics(t,x,sigma,param,step)
h = (x(1)*param.sf.radius-param.planet.radiusEquatorial)/1000;
% These stopping conditions should be in param
hmin = 6; %km
% vmax = 480/param.sf.velocity; %m/s
% vmin = 312/param.sf.velocity;
% if (x(4) < vmax && x(7) <= 0) || (h < hmin) || (x(4) < vmin)
if (SatisfiesParachute(h,x(4)) && x(7) <= 0) || ((x(1)*param.sf.radius-param.planet.radiusEquatorial)/1000 < hmin) 
    dX = zeros(size(x));
    done = true;
    L = nan;
    D = nan;
    M = nan;
else
    done = false;
    [g,L,D,~,M,a,rho] = EntryForces(x,param.planet,param.vehicle,param.sf);
%     L = L + 0.1;
%     D = D + 0.1;
    if param.ENFORCEBANK && step
        sigma_ex = Saturate(x(8),-param.bankLimits.angleMax,param.bankLimits.angleMax);
        dx = EntryDynamics(x,sigma_ex,g,L,D);
        du = BankAngleDynamics(t,x(8:9),sigma,param.bankLimits,param.bankGains);
        
    else
        dx = EntryDynamics(x,sigma,g,L,D);
        if param.ENFORCEBANK % Enforce bank but we aren't stepping so no bank angle computation is needed
            du = [0;0];
        else
            du = [];
        end
    end
    
    dX = [dx;-x(4)*cos(x(5))/x(1);du]; % Augment the state with the range computation and bank dynamics
end


end

function sf = getScaleFactors(planet,scale)
if scale
    sf.time = sqrt(planet.radiusEquatorial/planet.gm);
    sf.velocity = sqrt(planet.gm*planet.radiusEquatorial);
    sf.radius = planet.radiusEquatorial;
else
    sf.time = 1;
    sf.velocity = 1;
    sf.radius = 1;
end
end

function xout = Scale(x,param)
%x = [r,theta,phi,v,gamma,psi,s] (+[sigma, sigma_dot] if bank constraints)
%xout = x scaled appropriately based on flags

xout = x;
xout(1) = x(1)/param.sf.radius;
xout(4) = x(4)/param.sf.velocity;
end
function xout = Unscale(x,param)
%x = [r,theta,phi,v,gamma,psi,s] (+[sigma, sigma_dot] if bank constraints)
%xout = x scaled appropriately based on flags

xout = x;
xout(:,1) = x(:,1)*param.sf.radius;
xout(:,4) = x(:,4)*param.sf.velocity;
end
function e = eFromState(r,v)
% Computes the current energylike variable e from the radius and velocity
e = 1./r-0.5*v.^2;

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
v0 = 5505/param.sf.velocity;
vf = 500/param.sf.velocity;

M = (lim0-limf)/(v0-vf);
CR0 = lim0-M*v0;
limit = limf*(state(4)<vf && state(4)>0) + (CR0 + state(4)*M).*(state(4)>=vf);
psi_d = DesiredHeading(state(2),state(3),param.thetaf,param.phif);
e = psi_d-state(6);

if state(4) < 180/param.sf.velocity
    snew = -sCR;
else
    
    if abs(e) > limit &&  sBank == -sign(e)
        snew = -sBank;
    else
        snew = sBank;
    end
    
end

end