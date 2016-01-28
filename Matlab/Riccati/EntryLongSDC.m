% ENTRYLONGSDC Computes the linear-like factorization of the longitudinal entry dynamics.
%   ENTRYLONGSDC() Gives the dynamic coefficient matrix A, the input coupling
%   matrix B, and the measurement sensitivity matrix C.
%   This particular implementation is for the entry equations with range as
%   the independent variable.

function [A,B,C] = EntryLongSDC(planetModel,vehicleModel)

% x = [r,V,gamma,sigma]

A = @(X) Amatrix(X,planetModel,vehicleModel);
B = [0;0;0;1];
C = @(X) Cmatrix(X,planetModel,vehicleModel);
end
function A = Amatrix(x,planetModel,vehicleModel)
r = x(1);
V = x(2);
gamma = x(3);
sigma = x(4);
[g,L,D] = EntryForces(x,planetModel,vehicleModel);
A = zeros(4,4);

%Radius:
A(1,3) = tan(gamma)/gamma;

%Velocity:
w1 = 1/3;
w2 = 1/3;
w3 = 1/3;
A(2,1) = -w1*(D/V/cos(gamma)/r + g*tan(gamma)/V/r);
A(2,2) = -w2*(D/cos(gamma)/V^2 + g*tan(gamma)/V^2);
A(2,3) = -w3*(D/cos(gamma)/V/gamma + g*tan(gamma)/V/gamma);

%FPA:
w1 = 0.25;
w2 = 0.75;
%This is discontinuous
% if abs(gamma) < .01;
% w3 = 0; %How would dynamic weights affect the outcome?
% else
%     w3 = .1;
% end
% if abs(sigma) < .01
%     w4 = 0;
% else
%     w4 = 1;
% end
w3 = .1+gamma^2;
w4 = 1+(sigma/pi)^6;
wt = w1+w2+w3+w4; %Total weight, for normalization

fpa_prime = (L*cos(sigma)/cos(gamma)-g)/V^2 + 1/r;
A(3,1) = w1/wt*fpa_prime/r;
A(3,2) = w2/wt*fpa_prime/V;
A(3,3) = w3/wt*fpa_prime/gamma; %This is not well-behaved for gamma = 0
A(3,4) = w4/wt*fpa_prime/sigma; %This is not well-behaved for sigma = 0

%Bank Angle:
%All zero since the bank angle is affine in the pseudo-control u
end

%% Drag-tracking output matrix:
function C = Cmatrix(x,planetModel,vehicleModel)

r = x(1);
V = x(2);
% gamma = x(3);
% sigma = x(4);
[~,~,D] = EntryForces(x,planetModel,vehicleModel);
C = zeros(1,4);
w1 = 0.5;
w2 = 0.5;
C(1) = w1*D/r;
C(2) = w2*D/V;

end