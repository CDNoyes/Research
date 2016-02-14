function dX = EntrySMO(t,x,sigma,planetModel,vehicleModel)
r_eq = planetModel.radiusEquatorial;

% Check parachute constraints - stop if slow enough or too low.
hmin = 6; %km
vmin = 480; %m/s
if x(4) < vmin || (x(1)-r_eq)/1000 < hmin
    dX = zeros(size(x));
else
    [r,~,~,V,gamma,~] = ParseState(x');
    [g,L,D,~,~,~,rho,rhodot] = EntryForces(x,planetModel,vehicleModel);
    
    [a,b] = DragFBL(g,L,x(7),r,V,gamma,rho,rhodot,x(8)); %g,L,D,r,V,gamma,rho,rho_dot,D_dot
    [alpha, k] = computeSMOGains(4);
    dxhat = SMO(x(7:9),D,a,b,cos(sigma),k,alpha); % Observer states
    
    dx = EntryDynamics(x(1:6),sigma,g,L,D);
    
    
    dX = [dx;dxhat];
end

end

function [alpha, k] = computeSMOGains(lambda)

alpha = (factorial(3)./(factorial(1:3).*factorial(3-(1:3)))).*lambda.^(1:3);
k = [2.5852,14.2055,19.5150];
end