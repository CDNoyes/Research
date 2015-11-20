%BANKANGLEDYNAMICS Computes the bank angle profile that tracks a given
%profile subject to constraints on maximum angle, angular rate, and angular
%acceleration.
%   LIMITS is a structure with the max rate and acceleration.
%   H is a prediction horizon parameter that can reduce the error resulting
%   from the constraints.
function dx = BankAngleDynamics(t,x,profile,limits,h)

Kp = 0.5;
Kd = 1.4;
% h = 5; %prediction horizon -  start maneuvers early
sigmaC = profile(t+h);
s = sign(x(1));
x(1) = s*Saturate(abs(x(1)),limits.angleMin, limits.angleMax);
dx(2,1) = Saturate(Kp*(sigmaC-x(1))-Kd*x(2),-limits.acceleration, limits.acceleration); % Acceleration
dx(1) = Saturate(x(2), -limits.rate, limits.rate); % Rate

end