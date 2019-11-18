function output = EntryConstraints( input )

S       = input.auxdata.vehicle.area;
m       = input.auxdata.vehicle.mass;
mu      = input.auxdata.planet.mu;
rp      = input.auxdata.planet.radiusEquatorial;
d       = input.auxdata.delta;

x = input.phase.state;
nx = size(x,2);
if nx == 7
    dSigma = input.phase.control;
    Sigma = x(:,7);
else
    dSigma = [];
    Sigma = input.phase.control;
end

r = x(:,1);
theta = x(:,2);
phi = x(:,3);
v = x(:,4);
gamma = x(:,5);
psi = x(:,6);



g = mu./(r.^2);
h = r-rp;
[rho,a] = MarsAtmosphericDensity(h);
M = v./a;
[C_D,C_L] = AerodynamicCoefficients(M);
q = 0.5*rho.*v.^2*S/m;
drag = q.*(C_D+d.CD);
lift = q.*C_L;


dr = v.*sin(gamma);
dtheta = v.*cos(gamma).*cos(psi)./(r.*cos(phi));
dphi = v.*cos(gamma).*sin(psi)./r;
dv = -drag - g.*sin(gamma);
dgamma = lift.*cos(Sigma)./v - (g./v - v./r).*cos(gamma);
dpsi = -lift.*sin(Sigma)./(v.*cos(gamma)) - (v./r).*cos(gamma).*cos(psi).*tan(phi);

output.dynamics = [dr dtheta dphi dv dgamma dpsi dSigma];
output.integrand = zeros(size(r));
% output.integrand = (abs(Sigma)-pi/4).^2 ./(1+input.phase.time);


end




