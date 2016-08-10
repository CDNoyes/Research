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

ref = TrajectorySummary(t,x,Sigma(t),DR,CR);
EntryPlots(ref)


save(['EntryGuidance/HighElevationPlanner/Trajectories/ReferenceTrajectory_DR',num2str(DR),'_CR',num2str(CR),'.mat'],'ref','Sigma');
%% Fly the open loop trajectory with the drag observer
% Conclusions: The tracking is perfect. Without acceleration constraint,
% the error is 0.388 meters. With acceleration constraints but optimized
% gains and prediction horizon, the total range error is 225 meters.
clc; close all; clear;
load('EntryGuidance/HighElevationPlanner/Trajectories/ReferenceTrajectory_DR780_CR0.mat');
mars = Mars();
vm = VehicleModel();
[gains.alpha, gains.k] = computeSMOGains(10);

x0 = ref.state(1,:)';
xhat0 = [ref.D(1);ref.D_dot(1);0];
u0 = [Sigma(0);0];
tic;
[t,x] = ode45(@(T,X) EntrySMO(T,X,Sigma(T+2.42),Mars(),VehicleModel(),gains,ref),[0,max(ref.time)+20],[x0;xhat0;u0]);
t_obs = toc;
disp(['Simulation time: ',num2str(t_obs),' s.'])

for i = 1:length(t)
    [~,D_true(i)] = EntrySMO(t(i),x(i,:)',Sigma(t(i)),mars,vm,gains,ref);
end

obs = TrajectorySummary(t,x,x(:,10),ref.target.DR,ref.target.CR);
obs.observer = 1;
EntryPlots(obs)
D_true = D_true(1:length(obs.time));
save(['EntryGuidance/HighElevationPlanner/Trajectories/ObserverTrajectory_DR',num2str(ref.target.DR),'_CR',num2str(ref.target.CR),'.mat'],'obs','Sigma');
%Compare observed quantities with those computed directly:
lineWidth = 2;
markerSize = 10;
fontSize = 12;
fontWeight = 'bold';
fontColor = 'k';

figure
semilogy(obs.time,abs(D_true'-obs.state(:,7)),'LineWidth',lineWidth)
xlabel('Time (s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Error in Drag [true-observer] (m/s^2)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
title('Observer Performance')
set(gcf,'name','Observer Performance', 'numbertitle','off','WindowStyle','docked')
grid on




