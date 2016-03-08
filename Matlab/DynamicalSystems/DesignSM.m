% DESIGNSM Designs a time-varying sliding manifold for a second order
% system.
% DESIGNSM(Error0, ControlMargin) Computes the parameters for a
% time-varying sliding line (or manifold) of three different types based on
% the initial error state (where error is x(t)-x_d(t)) (note that e_dot(0)
% is assumed to be 0) and the control margin, i.e. the difference between
% the control limit U and the feedforward control required in the presence
% of worst-case disturbances.
%
% The three types of sliding are:
% {constant acceleration, constant velocity, terminal}
%
% Reference: "Time Varying Sliding Modes for Second Order Systems"

function [T,A,B,C,beta] = DesignSM(e0, Umax)

se0 = sign(e0);


%Constant Velocity:
T = sqrt(2*e0/Umax);
A = Umax*se0;
B = -se0*sqrt(2*Umax*abs(e0));
C = [];
beta = sqrt(2*Umax/abs(e0));




end