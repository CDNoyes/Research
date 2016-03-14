function [g,L,D,hs,M,a,rho,rho_dot,temp,C_D,C_D_dot] = EntryForces(x,planetModel,vehicleModel,ScaleFactor)

if nargin < 4 || isempty(ScaleFactor)
    ScaleFactor.radius = 1;
    ScaleFactor.velocity = 1;
end

%Constants:
r_eq = planetModel.radiusEquatorial;        % m
mu = planetModel.mu;                        % gravitational parameter, m^3/s^2
%Vehicle Parameters
S = vehicleModel.area;                      % reference wing area, m^2
m = vehicleModel.mass;                      % mass, kg

r = x(1);       %vehicle radius, m
V = x(4);       %velocity, m/s
gamma = x(5);
h_dot = V*sin(gamma);

g = mu/r^2/ScaleFactor.velocity^2/ScaleFactor.radius;      % gravity, m/s^2

h = r*ScaleFactor.radius-r_eq;                             % Altitude
[rho,a,hs,~,a_dot,rho_dot] = MarsAtmosphericDensity(h,h_dot);                       % Density and Speed of Sound
M = V*ScaleFactor.velocity/a;                              % Mach Number
[C_D,C_L,dC_DdM] = AerodynamicCoefficients(M); 
temp = .5*rho*V^2*S/m*ScaleFactor.radius;
L = temp*C_L;
D = temp*C_D;

V_dot = -D-g*sin(gamma);
dMdt = V_dot/a - V*a_dot/a^2;
C_D_dot = dC_DdM*dMdt;

end