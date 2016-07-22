function [t,x] = OpenLoopSim(x0, bank, ref)


%Ref needs target, deltas, and initial states

[t,x] = ode45(@(T,X) OLDyn(T,X,bank(T),ref.planet,ref.vehicle,ref.state,ref.target,ref.delta),[0,350],x0,[]);



end


function dx = OLDyn(~,x,sigma,planetModel,vehicleModel,state,target,delta)



% dtr = pi/180;
% sigmaMin = 0*dtr;
% sigmaMax = 90*dtr;

[r,theta,phi,V,gamma,~] = ParseState(x');

[g,L,D,~,~,~,rho,rhodot,~,CD,C_D_dot] = EntryForces(x,planetModel,vehicleModel);

% Model uncertainty, control uncertainty, and "noise" acting in the input channel are applied here:
% delta.CD = (0.35 + .05*sin(.5*t));
% delta.rho_offset = .1; % Constant offset
% delta.rho_mult = .3; % Multiplicative factor
D_perturbed = D + delta.CD*D/CD + D*delta.rho + D/CD*delta.rho*delta.CD;



% Check parachute constraints - stop if slow enough or too low.
hmin = 6.2; %km
vmax = 480; %m/s
vmin = 325; %m/s
[DR,CR] = Range(state(1,2),state(1,3),state(1,6),theta,phi);
if (V < vmax && (DR >= target.DR)) || ((r-planetModel.radiusEquatorial)/1000 < hmin) || (V < vmin)
    dx = zeros(size(x));
    
else
    % Compute the dynamics
    dx = EntryDynamics(x(1:6),sigma,g,L,D_perturbed);
    
   
    
end
end