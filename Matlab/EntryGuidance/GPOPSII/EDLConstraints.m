function output = EDLConstraints( input )

S       = input.auxdata.vehicle.area;
m       = input.auxdata.vehicle.mass; % Initial mass 
mu      = input.auxdata.planet.mu;
rp      = input.auxdata.planet.radiusEquatorial;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entry Phase 
x = input.phase(1).state;
dSigma = input.phase(1).control;

r = x(:,1);
% theta = x(:,2);
phi = x(:,3);
v = x(:,4);
gamma = x(:,5);
psi = x(:,6);
Sigma = x(:,7);


g = mu./(r.^2);
h = r-rp;
[rho,a] = MarsAtmosphericDensity(h);
M = v./a;
[C_D,C_L] = AerodynamicCoefficients(M);
q = 0.5*rho.*v.^2*S/m;
drag = q.*(C_D);
lift = q.*C_L;


dr = v.*sin(gamma);
dtheta = v.*cos(gamma).*cos(psi)./(r.*cos(phi));
dphi = v.*cos(gamma).*sin(psi)./r;
dv = -drag - g.*sin(gamma);
dgamma = lift.*cos(Sigma)./v - (g./v - v./r).*cos(gamma);
dpsi = lift.*sin(Sigma)./(v.*cos(gamma)) - (v./r).*cos(gamma).*cos(psi).*tan(phi);

output(1).dynamics = [dr dtheta dphi dv dgamma dpsi dSigma];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SRP phase 
x = input.phase(2).state;
u = input.phase(2).control;

r = x(:,1);
% theta = x(:,2);
phi = x(:,3);
v = x(:,4);
gamma = x(:,5);
psi = x(:,6);
m = x(:,7);

T = u(:,1);
pitch = u(:,2);

Sigma = zeros(size(T)); % Level flight 

g = mu./(r.^2);
h = r-rp;
[rho, a] = MarsAtmosphericDensity(h);
M = v./a;
[C_D, C_L] = AerodynamicCoefficients(M);
q = zeros(size(T));%0.5*rho.*v.^2*S./m;
drag = q.*C_D;
lift = q.*C_L;


dr = v.*sin(gamma);
dtheta = v.*cos(gamma).*cos(psi)./(r.*cos(phi));
dphi = v.*cos(gamma).*sin(psi)./r;
dv = -drag - g.*sin(gamma) + T./m.*cos(pitch-gamma);
dgamma = lift.*cos(Sigma)./v - (g./v - v./r).*cos(gamma) + T./(v.*m).*sin(pitch-gamma);
dpsi = lift.*sin(Sigma)./(v.*cos(gamma)) - (v./r).*cos(gamma).*cos(psi).*tan(phi);
dm = -T/input.auxdata.vehicle.v_exit;

output(2).dynamics = [dr dtheta dphi dv dgamma dpsi dm];

end