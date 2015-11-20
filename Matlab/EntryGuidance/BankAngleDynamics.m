%BANKANGLEDYNAMICS
%LIMITS is a structure with the max rate and acceleration.

function dx = BankAngleDynamics(t,x,profile,limits)

Kp = 0.5;
Kd = 0.9;

sigmaC = profile(t);
s = sign(x(1));
x(1) = s*Saturate(abs(x(1)),limits.angleMin, limits.angleMax);
dx(2,1) = Saturate(Kp*(sigmaC-x(1))-Kd*x(2),-limits.acceleration, limits.acceleration); % Acceleration
dx(1) = Saturate(x(2), -limits.rate, limits.rate); % Rate

end