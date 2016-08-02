%High Elevation Planner
clear; clc; close all

%Initial Guess
p = [5,51,250]; % 12.5 km for dr 780

%Setup the optimization tolerances
opt = optimset('tolX',1e-2,'tolFun',1e-4,'Display','Iter');

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
% tf = EntryAnalysis(t,x,DR,CR);
dtr = pi/180;
sigmaMin = 18.19*dtr;
sigmaMax = 87.13*dtr;
Sigma = @(T) BankAngleProfile(T,p(1),p(2),p(3), sigmaMin, sigmaMax); %The ideal bank profile with only rate constraint

ref = TrajectorySummary(t,x,Sigma(t),DR,CR);
EntryPlots(ref)
return
save('EntryGuidance/HighElevationPlanner/Trajectories/controlHEP.mat','ref','p');
%% Robust Version
ref.P0 = zeros(9);
ref.P0(7,7) = 0.04^2;
[CD_err,W] = SigmaPoints(0,0.04^2); % +-12% error in drag coefficient
CD_err = linspace(-.16,.16,9);
delta.CD = CD_err;
delta.rho = zeros(size(CD_err));
delta.p0 = normpdf(delta.CD,0,0.04);
x0 = repmat(ref.state(1,:)',1,length(delta.CD));
ref.sp.state = x0;
ref.sp.W = W;
ref.sp.delta = delta;
ref.planet = mars;
ref.vehicle = vm;
ref.method = 'PF';
p_rand = sort(rand(1,3)*175); %48.7372   95.7043  167.5637
f_sto = HighElevationCostFunctionStochastic(p,mars,vm,DR,CR,ref);
[p_rob,f_sto_rob,flag,output] = fminsearch(@(p) HighElevationCostFunctionStochastic(p,mars,vm,DR,CR,ref), p, opt);


% Plot nominal, robustified, and expected values on one plot
[~,t,x] = HighElevationCostFunction(p_rob, mars, vm, DR, CR);

%Show statistics and plots
dtr = pi/180;
sigmaMin = 18.19*dtr;
sigmaMax = 87.13*dtr;
Sigma = @(T) BankAngleProfile(T,p_rob(1),p_rob(2),p_rob(3), sigmaMin, sigmaMax); %The ideal bank profile with only rate constraint

rob = TrajectorySummary(t,x,Sigma(t),DR,CR);
rob.method = ref.method;
EntryPlots(rob)
return
ref = rob;
p = p_rob;
save(['EntryGuidance/HighElevationPlanner/Trajectories/controlHEP-',ref.method,'.mat'],'ref','p');

%% Bank Angle Dynamics - Applying acceleration constraints

rateMax = 20*dtr;
accMax = 5*dtr;
lim.rate = rateMax;
lim.acceleration = accMax;
lim.angleMax = sigmaMax;
lim.angleMin = sigmaMin;
KpKdH = 10*[0.5,1.4]; %Initial Guess
KpKdH = fminsearch(@(g)optimizeBankAngleDynamics(g,Sigma,lim,tf),KpKdH);
if length(KpKdH) == 2
    KpKdH(3) = 0;
end
[T,S] = ode45(@BankAngleDynamics,[0,tf],[Sigma(0);0],[],Sigma,lim,KpKdH);
err = (Sigma(T)'/dtr-Saturate(S(:,1),-sigmaMax,sigmaMax)/dtr).^2;
e = trapz(T,err);
disp(['Norm of error between commanded and executed bank angle: ',num2str(e)])
figure
subplot 211
plot(T,Sigma(T)/dtr,'r--')
hold all
plot(T,Saturate(S(:,1),-sigmaMax,sigmaMax)/dtr)
plot(T,Saturate(S(:,2),-lim.rate,lim.rate)/dtr)

legend('Commanded','Executed','Bank Rate')
subplot 212
plot(T,err)

 