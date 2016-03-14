function [dX,D_perturbed,D_noisy] = EntrySMO(t,x,sigma,planetModel,vehicleModel,observerGains)
r_eq = planetModel.radiusEquatorial;

[g,L,D,~,~,~,rho,rhodot,~,CD,C_D_dot] = EntryForces(x,planetModel,vehicleModel);

% Model uncertainty, control uncertainty, and "noise" acting in the input channel are applied here:
delta.CD = 0*.4; % Constant offset
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

    % Compute the dynamics
    dxhat = SMO(x(7:9),D_noisy,a,b,cos(sigma),observerGains.k,observerGains.alpha); % Observer states
    dx = EntryDynamics(x(1:6),sigma,g,L,D_perturbed);

    dX = [dx;dxhat];
end

end