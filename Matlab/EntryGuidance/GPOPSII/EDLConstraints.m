function output = EDLConstraints( input )

S       = input.auxdata.vehicle.area;
m       = input.auxdata.vehicle.mass; % Initial mass 
mu      = input.auxdata.planet.mu;
rp      = input.auxdata.planet.radiusEquatorial;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entry Phase 
x = input.phase(1).state;

if size(x, 2) == 7
dSigma = input.phase(1).control;
Sigma = x(:,7);
cosbank = cos(Sigma);
else
%    Sigma = input.phase(1).control;
cosbank = input.phase(1).control;
   dSigma = [];
   dpsi = zeros(size(cosbank));
end

r = x(:,1);
% theta = x(:,2);
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
drag = q.*(C_D);
lift = q.*C_L;


dr = v.*sin(gamma);
dtheta = v.*cos(gamma).*cos(psi)./(r.*cos(phi));
dphi = v.*cos(gamma).*sin(psi)./r;
dv = -drag - g.*sin(gamma);
dgamma = lift.*cosbank./v - (g./v - v./r).*cos(gamma);
if ~exist('dpsi', 'var')
    dpsi = lift.*sin(Sigma)./(v.*cos(gamma)) - (v./r).*cos(gamma).*cos(psi).*tan(phi);
end
output(1).dynamics = [dr dtheta dphi dv dgamma dpsi dSigma];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SRP phase 
% x = input.phase(2).state;
% u = input.phase(2).control;

% r = x(:,1);
% % theta = x(:,2);
% phi = x(:,3);
% v = x(:,4);
% gamma = x(:,5);
% psi = x(:,6);
% m = x(:,7);
% 
% T = u(:,1);
% pitch = u(:,2);
% 
% Sigma = zeros(size(T)); % Level flight 
% 
% g = mu./(r.^2);
% h = r-rp;
% [rho, a] = MarsAtmosphericDensity(h);
% M = v./a;
% [C_D, C_L] = AerodynamicCoefficients(M);
% q = zeros(size(T));%0.5*rho.*v.^2*S./m;
% drag = q.*C_D;
% lift = q.*C_L;
% 
% 
% dr = v.*sin(gamma);
% dtheta = v.*cos(gamma).*cos(psi)./(r.*cos(phi));
% dphi = v.*cos(gamma).*sin(psi)./r;
% dv = -drag - g.*sin(gamma) + T./m.*cos(pitch-gamma);
% dgamma = lift.*cos(Sigma)./v - (g./v - v./r).*cos(gamma) + T./(v.*m).*sin(pitch-gamma);
% dpsi = lift.*sin(Sigma)./(v.*cos(gamma)) - (v./r).*cos(gamma).*cos(psi).*tan(phi);
% dm = -T/input.auxdata.vehicle.v_exit;
% 
% output(2).dynamics = [dr dtheta dphi dv dgamma dpsi dm];

s1 = input.phase(2).state;
control1 = input.phase(2).control;

% Variables
% x1     = s1(:,1);
% y1     = s1(:,2);
% z1     = s1(:,3);
u1     = s1(:,4);
v1     = s1(:,5);
w1     = s1(:,6);
m1     = s1(:,7);


T     = control1(:,1);
ux    = control1(:,2);
uy    = control1(:,3);
uz    = control1(:,4);
aT    = T./m1;

gravity1 = 3.71;
Isp = 290;
g0 = 9.81;
ve = Isp*g0;

% EOM
xdot1 = u1;
ydot1 = v1;
zdot1 = w1;
udot1 = aT.*ux;
vdot1 = aT.*uy;
wdot1 = aT.*uz - gravity1;
mdot1 = -T/ve;

% Output
% output(1).integrand = T;
output(2).dynamics = [xdot1 ydot1 zdot1 udot1 vdot1 wdot1 mdot1];
output(2).path = (sum(control1(:,2:4).^2,2)); % Norm of the control directions

end