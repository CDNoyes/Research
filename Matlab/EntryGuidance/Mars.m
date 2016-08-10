%MARS Creates a basic model of Mars.
%   MARS returns a structure with information about the planets size and a
%   model of the density of its atmosphere. This function is meant to be an
%   easily extendable interface for all Mars-based analysis.

function mars = Mars()

mars.radiusEquatorial = 3396.2e3;           % meters
mars.radiusPolar = 3376.2e3;                % meters
mars.atmosphere = @MarsAtmosphericDensity;  % kg/m^3
% mars.mu = 0.04283e6;                      % km^3/s^2
mars.mu = 0.04283e15;                       % m^3/s^2
mars.omega = 7.095e-5;                      % rad/s, angular rate of planet rotation
mars.gm =  mars.mu/mars.radiusEquatorial^2; % m/s^2
mars.g0 = 9.81;                             % earth g, used for normalizing and in isp equations
end