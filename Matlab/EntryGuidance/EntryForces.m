function [g,L,D,hs,M,a,rho] = EntryForces(x,planetModel,vehicleModel,ScaleFactor)

if nargin < 4 || isempty(ScaleFactor)
    ScaleFactor.radius = 1;
    ScaleFactor.velocity = 1;
end

%Constants:
r_eq = planetModel.radiusEquatorial;      % m

%Vehicle Parameters
S = vehicleModel.area; %reference wing area, m^2
m = vehicleModel.mass; %mass, kg

r = x(1);       %vehicle radius, m
V = x(4);       %velocity, m/s

mu = 4.2830e13;     % gravitational parameter, m^3/s^2
g = mu/r^2/ScaleFactor.velocity^2/ScaleFactor.radius;      % gravity, m/s^2

h = r*ScaleFactor.radius-r_eq;                             % Altitude
[rho,a,hs] = MarsAtmosphericDensity(h);                       % Density and Speed of Sound
M = V*ScaleFactor.velocity/a;                              % Mach Number
[C_D,C_L] = AerodynamicCoefficients(M); 
temp = .5*rho*V^2*S/m*ScaleFactor.radius;
L = temp*C_L;
D = temp*C_D;

end