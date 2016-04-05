function output = EntryConstraints( input )

ref     = input.auxdata.ref;
S       = input.auxdata.vehicle.area;
m       = input.auxdata.vehicle.mass;
mu      = input.auxdata.planet.mu;
rp      = input.auxdata.planet.radiusEquatorial;

x = input.phase.state;
dSigma = input.phase.control;

r = x(:,1);
theta = x(:,2);
phi = x(:,3);
v = x(:,4);
gamma = x(:,5);
psi = x(:,6);
Sigma = x(:,7);

E = 0.5*v.^2 + mu/rp - mu./r;
drag_ref = interp1(ref.energy,ref.D,E,'spline');

g = mu./(r.^2);
h = r-rp;
[rho,a] = MarsAtmosphericDensity(h);
M = v./a;
[C_D,C_L] = AerodynamicCoefficients(M);
q = 0.5*rho.*v.^2*S/m;
drag = q.*C_D;
lift = q.*C_L;


dr = v.*sin(gamma);
dtheta = v.*cos(gamma).*cos(psi)./(r.*cos(phi));
dphi = v.*cos(gamma).*sin(psi)./r;
dv = -drag - g.*sin(gamma);
dgamma = lift.*cos(Sigma)./v - (g./v - v./r).*cos(gamma);
dpsi = lift.*sin(Sigma)./(v.*cos(gamma)) - (v./r).*cos(gamma).*cos(psi).*tan(phi);

output.dynamics = [dr dtheta dphi dv dgamma dpsi dSigma];
output.integrand = (drag-drag_ref).^2;

end




