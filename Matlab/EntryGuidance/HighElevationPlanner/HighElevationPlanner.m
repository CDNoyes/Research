%High Elevation Planner
clear; clc; close all

% p = [15;130;160]; %Works well for 780

%Initial Guess
p = [5,51,250]; % 12km

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

%Show statistics and plots
tf = EntryAnalysis(t,x,DR,CR);
% EntryPlots(t,x)
% close all

%% Bank Angle Dynamics - Applying constraints
dtr = pi/180;
sigmaMin = 18.19*dtr;
sigmaMax = 87.13*dtr;
rateMax = 20*dtr;
accMax = 5*dtr;
Sigma = @(T) BankAngleProfile(T,p(1),p(2),p(3),sigmaMin, sigmaMax); %The ideal bank profile with no constraints
lim.rate = rateMax;
lim.acceleration = accMax;
lim.angleMax = sigmaMax;
lim.angleMin = sigmaMin;
% Sigma = @(t) sin(t/10);
[T,S] = ode45(@BankAngleDynamics,[0,tf],[Sigma(0);0],[],Sigma,lim);

figure
plot(T,Sigma(T),'r--')
hold all
plot(T,Saturate(S(:,1),-sigmaMax,sigmaMax))
legend('Ideal','Constrained')