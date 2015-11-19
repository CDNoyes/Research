%MARS Creates a basic model of Mars.
%   MARS returns a structure with information about the planets size and a
%   model of the density of its atmosphere. This function is meant to be an
%   easily extendable interface for all Mars-based analysis.

function mars = Mars()

mars.radiusEquatorial = 3396.2e3;
mars.radiusPolar = 3376.2e3;
mars.atmosphere = @MarsAtmosphericDensity;

end