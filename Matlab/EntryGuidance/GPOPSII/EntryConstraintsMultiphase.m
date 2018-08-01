function output = EntryConstraintsMultiphase( input )

S       = input.auxdata.vehicle.area;
m       = input.auxdata.vehicle.mass;
mu      = input.auxdata.planet.mu;
rp      = input.auxdata.planet.radiusEquatorial;
d       = input.auxdata.delta;

bankmin=0;
bankmax=90*pi/180;
bank_angles = [bankmin,-bankmax,bankmax,-bankmin];

for iphase = 1:input.auxdata.nPhases
    x = input.phase(iphase).state;
    Sigma = bank_angles(iphase)*ones(size(x(:,1)));
    
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
    
    output(iphase).dynamics = [dr dtheta dphi dv dgamma dpsi];
%     output(iphase).integrand = zeros(size(r));
end

end