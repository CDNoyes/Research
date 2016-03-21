%BANKANGLEDYNAMICS Computes the bank angle profile that tracks a given
%profile subject to constraints on maximum angle, angular rate, and angular
%acceleration.
%   LIMITS is a structure with the max rate and acceleration.
%   K is an array of parameters that can reduce the error resulting
%   from the constraints; [Kp, Kd, h]

function dx = BankAngleDynamics(t,x,profile,limits,K)
% x - [sigma, sigma_dot] (executed)

Kp = K(1);
Kd = K(2);
h = K(3);
if isa(profile,'function_handle')
    sigmaC = profile(t+h); %Commanded bank angle
else
    sigmaC = profile;
end
x(1) = sign(x(1))*Saturate(abs(x(1)), limits.angleMin, limits.angleMax);

dx(2,1) = Saturate(Kp*(sigmaC-x(1))+Kd*(-x(2)), -limits.acceleration, limits.acceleration); % Acceleration
dx(1) = Saturate(x(2), -limits.rate, limits.rate); % Rate

end