function [dX,D_perturbed,D_noisy] = EntrySMO(t,x,sigma,planetModel,vehicleModel,observerGains)
r_eq = planetModel.radiusEquatorial;
dtr = pi/180;
sigmaMin = 18.19*dtr;
sigmaMax = 87.13*dtr;
sigma_ex = Saturate(x(10),-sigmaMax,sigmaMax);

[g,L,D,~,~,~,rho,rhodot,~,CD,C_D_dot] = EntryForces(x,planetModel,vehicleModel);

% Model uncertainty, control uncertainty, and "noise" acting in the input channel are applied here:
delta.CD = 0*(0.35 + .05*sin(.5*t)); % Constant offset
% delta.rho = .1; % Constant offset
delta.rho = 0*.3*rho; % Multiplicative factor
D_perturbed = D + delta.CD*D/CD + D*delta.rho/rho + D/CD*delta.rho/rho*delta.CD;

% Measurement noise
D_noisy = D_perturbed;%D_perturbed*(1 + .0001*sin(t)); % probably very small in truth


% Check parachute constraints - stop if slow enough or too low.
hmin = 6; %km
vmax = 480; %m/s
if x(4) < vmax || (x(1)-r_eq)/1000 < hmin
    dX = zeros(size(x));
    
else
    % Nominal model of forces and drag dynamics
    [r,~,~,V,gamma,~] = ParseState(x');
    [a,b] = DragFBL(g,L,x(7),r,V,gamma,rho,rhodot,x(8),CD,C_D_dot); %g,L,D,r,V,gamma,rho,rho_dot,D_dot,C_D,C_D_dot
%     [a,b] = DragFBL(g,L,D,r,V,gamma,rho,rhodot,[],CD,C_D_dot); %Purely model based

    % Compute the dynamics
    dxhat = SMO(x(7:9),D_noisy,a,b,cos(sigma_ex),observerGains.k,observerGains.alpha); % Observer states
    dx = EntryDynamics(x(1:6),sigma_ex,g,L,D_perturbed);
    
    % Bank Angle dynamics:
    rateMax = 20*dtr;
    accMax = 5*dtr;
    lim.rate = rateMax;
    lim.acceleration = accMax;
    lim.angleMax = sigmaMax;
    lim.angleMin = sigmaMin;
    K = [0.5643, 1.2934, 2.4238];
    du = BankAngleDynamics(t,x(10:11),sigma,lim,K);
    dX = [dx;dxhat;du];
end

end