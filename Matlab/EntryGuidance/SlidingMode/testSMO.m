% Run script for open loop simulation to demonstrate observer performance

%% First we compute the reference trajectory and the drag profile
%High Elevation Planner
clear; clc; close all

%Initial Guess
p = [5,51,250]; % 12.5 km for dr 780

%Setup the optimization tolerances
opt = optimset('tolX',1e-3,'tolFun',1e-8);

%Create the standard models to be used
mars = Mars();
vm = VehicleModel();

%Define the target location
DR = 780;
CR = 0;

%Call the optimization routine
[p,fval,flag,output] = fminsearch(@(p) HighElevationCostFunction(p,mars,vm,DR,CR), p, opt);
%Compute the state trajectory
[cost,t,x] = HighElevationCostFunction(p, mars, vm,DR,CR);

dtr = pi/180;
sigmaMin = 18.19*dtr;
sigmaMax = 87.13*dtr;
rateMax = 20*dtr;
accMax = 5*dtr;
Sigma = @(T) BankAngleProfile(T,p(1),p(2),p(3), sigmaMin, sigmaMax);

ref = TrajectorySummary(t,x,DR,CR);
EntryPlots(ref)

%% Fly the open loop trajectory with the drag observer

x0 = ref.state(1,:)';
xhat0 = [ref.D(1);0;0];
[t,x] = ode45(@(T,X) EntrySMO(T,X,Sigma(T),mars,vm),[0,tf],[x0;xhat0]);

tf = EntryAnalysis(t,x,DR,CR);
obs = TrajectorySummary(t,x,DR,CR);
EntryPlots(obs)

%Compare observed quantities with those computed directly:
figure
plot(t,abs(obs.D'-obs.state(:,7)))
xlabel('Time (s)')
ylabel('Error in Drag [model-observer] (m/s^2)')
figure
plot(t,x(:,9),'k')