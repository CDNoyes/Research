% Test bed for drag tracking control schemes. There is also a drag observer
% and disturbance estimator available for use in the controllers if
% desired.

function traj = EntryDragTracking()

% Reference Data:
load('EntryGuidance/HighElevationPlanner/Trajectories/ObserverTrajectory_DR780_CR0.mat');
% ref = obs;
load('EntryGuidance/HighElevationPlanner/Trajectories/ReferenceTrajectory_DR780_CR0.mat');
ref.state(:,9) = ReplaceNAN(interp1(obs.time,obs.state(:,9),ref.time),0);

% System Models:
mars = Mars();
vm = VehicleModel();

% e0 = [ref.D(1)-0;ref.D_dot(i)-0];

[S,controller,T] = DesignSM(0,1,'Velocity');
[gains.alpha, gains.k] = computeSMOGains(10);

x0 = ref.state(1,1:6)';
xhat0 = [ref.D(1);ref.D_dot(1);0];
tic
[t,x] = ode45(@Tracking, linspace(0, max(ref.time)*1.3,5000), [x0;xhat0],[],ref,Mars(),VehicleModel(),gains,S,controller);
t_elapsed = toc;
disp(['Simulation time: ',num2str(t_elapsed),' s.'])

for i=1:length(t)
    [~,sigma(i),D_true(i),~,s(i),d_hat(i)] = Tracking(t(i),x(i,:)',ref,mars,vm,gains,S,controller);
end

traj = TrajectorySummary(t,x,sigma,ref.target.DR,ref.target.CR);
EntryPlots(traj)
n = length(traj.time);
lineWidth = 2;
markerSize = 10;
fontSize = 12;
fontWeight = 'bold';
fontColor = 'k';

figure
plot(traj.time,abs(D_true(1:n)'-traj.state(:,7)),'LineWidth',lineWidth)
xlabel('Time (s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Error in Drag [true-observer] (m/s^2)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
title('Observer Performance','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)

figure
plot(ref.time,ref.D,'LineWidth',lineWidth)
hold all
plot(traj.time,D_true(1:n),'LineWidth',lineWidth)
plot(traj.time,traj.D,'LineWidth',lineWidth)
legend('Desired Drag','Actual Drag','Expected Based on Model')
xlabel('Time (s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Drag (m/s^2)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)

figure
plot(traj.energy_norm,(D_true(1:n)'-interp1(ref.energy_norm,ref.D,traj.energy_norm,'spline')),'LineWidth',lineWidth)


figure
plot(traj.time,s(1:n),'LineWidth',lineWidth)
title('Sliding Surface', 'FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)

figure
plot(traj.time,d_hat(1:n),'LineWidth',lineWidth)
title('Disturbance Estimate', 'FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
end

function [dX,sigma,D_perturbed,D_noisy,s,dhat] = Tracking(t,x,ref,planetModel,vehicleModel,observerGains,S,controller)
if ~mod(floor(t),20)
    t
end
% Easier variable names
[r,theta,phi,V,gamma,~] = ParseState(x');

% Current energy
E = 0.5*V^2 + planetModel.mu/planetModel.radiusEquatorial - planetModel.mu./r; 

[g,L,D,~,~,~,rho,rhodot,~,C_D,C_D_dot] = EntryForces(x,planetModel,vehicleModel);

% Model uncertainty, control uncertainty, and "noise" acting in the input channel are applied here:
delta.CD = .1; % Constant offset
% delta.rho = .1; % Constant offset
delta.rho = .1*rho; % Multiplicative factor
D_perturbed = D + delta.CD*D/C_D + D*delta.rho/rho + D/C_D*delta.rho/rho*delta.CD;

% Measurement noise
D_noisy = D_perturbed;%D_perturbed*(1 + .0001*sin(t)); % probably very small in truth


% Check parachute constraints - stop if slow enough or too low.
hmin = 6; %km
vmax = 480; %m/s
[DR,CR] = Range(ref.state(1,2),ref.state(1,3),ref.state(1,6),theta,phi);
if (V < vmax && DR >= ref.target.DR) || (r-planetModel.radiusEquatorial)/1000 < hmin
    dX = zeros(size(x));
    sigma = 0;
    s = 0;
    dhat = 0;
else
    % Nominal model drag dynamics
    % Could also used measured drag directly if we think its better than our
    % observer
    if t < 1
        D_dot = [];
        D_hat = D;
    else
        D_dot = x(8);
        D_hat = x(7);
    end
    [a,b] = DragFBL(g,L,D_hat,r,V,gamma,rho,rhodot,D_dot,C_D,C_D_dot); %g,L,D,r,V,gamma,rho,rho_dot,D_dot
    
    % Reference trajectory
    D_des = interp1(ref.energy,ref.D,E,'spline');
    D_des_dot = interp1(ref.energy,ref.D_dot,E,'spline');
    D_des_ddot = interp1(ref.energy,ref.D_ddot,E,'spline');
    
    % Error states
    e = [x(7)-D_des; x(8)-D_des_dot]; %Should we (always) use the observer states? Maybe wait for convergence like above?

    %Sliding Surface Computation
    s = S(e,t);

    % Control

%     dhat = x(9)-interp1(ref.energy,ref.state(:,9),E,'spline');
%     dmax = Saturate(abs(dhat),0,.01); %0.001;%
    [a_p,b_p] = DragFBL(g,L,D_hat,r,V,gamma,rho+delta.rho,rhodot,D_dot,C_D+delta.CD,C_D_dot);
    dhat = (a-a_p) + (b-b_p)*interp1(ref.energy,ref.control,E,'spline');
%     dmax = abs(dhat)+1e-3*(~dhat);
    dmax = 0.12;
    u0 = D_des_ddot - a -dmax*Saturate(s/.001,-1,1);
    u = controller(b,u0,e,t);
    U = 1;
    u = Saturate(u,-U,U);

    % Compute the bank angle
    sigma = acos(u)*sign(interp1(ref.energy,ref.sigma,E,'spline',ref.sigma(end)));
    
if D < 0.1
    if delta.rho > 0
        sigma = 0;
%         sigma = interp1(ref.energy,ref.sigma,E,'linear',ref.sigma(end));
    else
        sigma = pi/3;
    end
    u = cos(sigma);
end

    % Compute the dynamics
    dxhat = SMO(x(7:9),D_noisy,a,b,u,observerGains.k,observerGains.alpha); % Observer states
    dx = EntryDynamics(x(1:6),sigma,g,L,D_perturbed);

    dX = [dx;dxhat];
end

end