% Test bed for drag tracking control schemes. There is also a drag observer
% and disturbance estimator available for use in the controllers if
% desired.

function traj = EntryDragTracking()

% Reference Data:
load('EntryGuidance/HighElevationPlanner/Trajectories/ReferenceTrajectory_DR780_CR0.mat');

% System Models:
mars = Mars();
vm = VehicleModel();

% Initial states:
x0 = ref.state(1,1:6)';
xhat0 = [ref.D(1);ref.D_dot(1);0];
sig0 = [ref.sigma(1);0];

% Control Design:
e0 = [xhat0(1)-ref.D(1);xhat0(2)-ref.D_dot(1)];
[S,controller,T] = DesignSM(e0(1),1,'Velocity');
[gains.alpha, gains.k] = computeSMOGains(10);

% Enforce Bank angle constraints on rate and acceleration?
global BANK_CONSTRAINTS
BANK_CONSTRAINTS = 0;

tic
[t,x] = ode45(@Tracking, linspace(0, max(ref.time)*1.3,5000), [x0;xhat0;sig0],[],ref,Mars(),VehicleModel(),gains,S,controller);
t_elapsed = toc;
disp(['Simulation time: ',num2str(t_elapsed),' s.'])

for i=1:length(t)
    [~,sigmaCom(i),D_true(i),~,s(i),d_hat(i),phase(i)] = Tracking(t(i),x(i,:)',ref,mars,vm,gains,S,controller);
end

bank = sigmaCom'*(~BANK_CONSTRAINTS) + x(:,10)*(BANK_CONSTRAINTS);

traj = TrajectorySummary(t,x,bank,ref.target.DR,ref.target.CR);
for i=1:length(traj.time)
    if i == 1
        csgn = sign(sig0(1));
    else
        csgn = sgn(i-1);
    end
    [sgn(i),headingErrorLimit(i)] = LateralGuidance(traj.CR(i),traj.D(i),csgn); % Drag based
    
    %     [sgn(i),headingErrorLimit(i)] = LateralGuidance(traj.CR(i),traj.state(i,4),csgn); %Velocity based
end
traj.observer = 1;
save('EntryGuidance/HighElevationPlanner/Trajectories/tempDragTrackingSolution.mat'); 
%%
EntryPlots(traj)
n = length(traj.time);
lineWidth = 2;
markerSize = 10;
fontSize = 12;
fontWeight = 'bold';
fontColor = 'k';

figure
plot(traj.energy_norm,sgn(1:n).*sigmaCom(1:n)*180/pi,'LineWidth',lineWidth)
hold all
plot(traj.energy_norm,x(1:n,10)*180/pi,'LineWidth',lineWidth)
legend('Command','Actual')
xlabel('Normalized Energy (-)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Bank Angle (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
grid on
box on
set(gcf,'name','Command Bank Profile', 'numbertitle','off','WindowStyle','docked')

% figure
% semilogy(traj.time,abs(D_true(1:n)'-traj.state(:,7)),'LineWidth',lineWidth)
% xlabel('Time (s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
% ylabel('Error in Drag [|true-observer|] (m/s^2)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
% title('Observer Performance','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
% grid on
% box on
% set(gcf,'name','Observer Performance', 'numbertitle','off','WindowStyle','docked')

figure
plot(ref.energy_norm,ref.D,'LineWidth',lineWidth)
hold all
plot(traj.energy_norm,D_true(1:n),'LineWidth',lineWidth)
plot(traj.energy_norm,traj.D,'LineWidth',lineWidth)
legend('Desired Drag','Actual Drag','Expected Based on Model')
xlabel('Normalized Energy (-)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Drag (m/s^2)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
grid on
box on
set(gcf,'name','Drag vs Energy', 'numbertitle','off','WindowStyle','docked')

% figure
% plot(traj.energy_norm,s(1:n)/max(traj.state(:,7)),'LineWidth',lineWidth)
% title('Normalized Sliding Surface', 'FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
% grid on
% box on
% axis([0,1, -max(1,1.2*max(abs(s(1:n)/max(traj.state(:,7))))),max(1,1.2*max(abs(s(1:n)/max(traj.state(:,7))))) ])
% set(gcf,'name','Normalized Sliding Surfance vs Energy', 'numbertitle','off','WindowStyle','docked')
% xlabel('Normalized Energy (-)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
% ylabel('s/D_{max}','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)


figure
plot(traj.state(:,4),traj.CR, traj.state(:,4),sgn,traj.state(:,4), traj.sigma, traj.state(:,4),headingErrorLimit,traj.state(:,4),-headingErrorLimit)
xlabel('Velocity (m/s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Cross Range (km)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
grid on
box on
set(gcf,'name','Crossrange Deadband vs Velocity', 'numbertitle','off','WindowStyle','docked')
legend('Trajectory CR','Commanded Sign','Bank Angle (rad)','CR Deadband')

if 1
    temp = 1;
end
end

function [dX,sigma,D_perturbed,D_noisy,s,dhat,phase] = Tracking(t,x,ref,planetModel,vehicleModel,observerGains,S,controller)
global BANK_CONSTRAINTS;
persistent banksign;
persistent currentPhase;
phaseText = {'Phase 1: Pre-entry -  Constant bank until .1g in sensed drag.'
    'Phase 2: Drag Tracking - SMO+SMC with lateral guidance.'
    'Phase 3: Heading Alignment - Nulling horizontal error.'
    'Phase 4: Parachute Deployment - Terminating simulation.'};
if isempty(banksign)
    banksign = sign(x(10));
end
if isempty(currentPhase)
    currentPhase = 0;
end


% Easier variable names
[r,theta,phi,V,gamma,~] = ParseState(x');
dtr = pi/180;
sigmaMin = 0*dtr;
sigmaMax = 90*dtr;
sigma_ex = Saturate(x(10),-sigmaMax,sigmaMax);
K = [0.5643, 1.2934, 2.4238]; %Bank Angle dynamics gains
M_align = 7;

% Current energy
E = 0.5*V^2 + planetModel.mu/planetModel.radiusEquatorial - planetModel.mu./r;

[g,L,D,~,M,~,rho,rhodot,~,C_D,C_D_dot] = EntryForces(x,planetModel,vehicleModel);

% Model uncertainty, control uncertainty, and "noise" acting in the input channel are applied here:
delta.CD = 0; % Constant offset
% delta.rho = .1; % Constant offset
delta.rho = 0*.1*rho; % Multiplicative factor
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
    phase = 4;
else
    % Nominal model drag dynamics
    % Could also used measured drag directly if we think its better than our observer
    D_dot = x(8);
    D_hat = x(7);
    [a,b] = DragFBL(g,L,D_hat,r,V,gamma,rho,rhodot,D_dot,C_D,C_D_dot); %g,L,D,r,V,gamma,rho,rho_dot,D_dot
    
    % Reference trajectory
    E_predict = 0.99*E ;
    D_des = interp1(ref.energy,ref.D,E_predict,'spline');
    D_des_dot = interp1(ref.energy,ref.D_dot,E_predict,'spline');
    D_des_ddot = interp1(ref.energy,ref.D_ddot,E_predict,'spline');
    
    % Error states
    e = [x(7)-D_des; x(8)-D_des_dot]; %Should we (always) use the observer states? Maybe wait for convergence like above?
    
    %Sliding Surface Computation
    s = S(e,t);
    
    dhat = x(9);
    dmax = abs(dhat);
    
    % Control
    if D_noisy/9.81 < 0.1  % Pre-entry control until .1 Earth g's are sensed
        u = ref.control(1);
        phase = 1;
        
    elseif D_noisy/9.81 >= 0.1 && M > M_align                  % Drag-tracking feedback control
        u0 = D_des_ddot - a - dmax*Saturate(s/.001,-1,1);
        u = controller(b,u0,e,t);
        U = 1;
        u = Saturate(u,0,U);
        phase = 2;
        
    else                % Heading Alignment
        phase = 3;
        sigma = -HeadingAlignment(x',L,ref.target,0.045);
    end
    
    % Compute the bank angle sign
    if phase < 3
        banksign = LateralGuidance(CR,D_noisy,banksign);
        %     banksign = LateralGuidance(CR,V,banksign);
        sigma = banksign*acos(u);
    end
    
    % Compute the dynamics
    if ~BANK_CONSTRAINTS
        sigma_ex = sigma;
    end
    dxhat = SMO(x(7:9),D_noisy,a,b,cos(sigma_ex),observerGains.k,observerGains.alpha); % Observer states
    dx = EntryDynamics(x(1:6),sigma_ex,g,L,D_perturbed);
    
    % Bank angle dynamic terms
    rateMax = 20*dtr;
    accMax = 5*dtr;
    lim.rate = rateMax;
    lim.acceleration = accMax;
    lim.angleMax = sigmaMax;
    lim.angleMin = sigmaMin;
    du = BankAngleDynamics(t,x(10:11),sigma,lim,K);
    dX = [dx;dxhat;du];
    
    
end

% Display transition info:
if currentPhase < phase
    currentPhase = phase;
    disp([phaseText{phase}, '(t = ',num2str(t), ' s)']);
end

end