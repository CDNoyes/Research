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

function [sig,U,T] = DesignSM(e0, Umax, type)

if e0 == 0
    sig = @(e,t) 0;
    U = @(b,u0,e,t) 1/b*u0;
    T = 0;
    return
end

se0 = sign(e0);


if strcmpi(type,'Velocity')     % Constant Velocity:
    T = sqrt(2*e0/Umax);
    A = Umax*se0;
    B = -se0*sqrt(2*Umax*abs(e0));
    beta = sqrt(2*Umax/abs(e0));
    sig = @(e,t) [beta,1]*e + (A*t+B)*(t<=T); 
    U = @(b,u0,e,t) 1/b*u0 + 1/b*(-beta*e(2)-A*(t<=T));
    
else                            % Terminal
    T = sqrt(6*abs(e0)/Umax);
    A = -0.5*se0*Umax;
    B = e0*sqrt(6*Umax/abs(e0));
    C = -3*e0;
    alpha = abs(fixPower(e0))*sqrt(6*Umax/abs(e0))/(4^(1/3));
    
    sig = @(e,t) SMT(T,A,B,C,alpha,e,t);
    U = @(b,u0,e,t) controllerT(T,A,B,C,alpha,b,u0,e,t);
end

end

function s = SMT(T,A,B,C,alpha,e,t)

% s = e(2) + real(2*A*t+B+alpha*Saturate((e(1)+A*t^2+B*t+C)/.01,-1,1)*(e(1)+A*t^2+B*t+C)^(2/3))*(t<=T) + real(alpha*Saturate(e(1)/.01,-1,1)*e(1)^(2/3))*(t>T);
s = e(2) + (2*A*t+B+alpha*Saturate((e(1)+A*t^2+B*t+C)/.01,-1,1)*fixPower(e(1)+A*t^2+B*t+C)^2)*(t<=T) + (alpha*Saturate(e(1)/.01,-1,1)*fixPower(e(1))^2)*(t>T);

end

function u = controllerT(T,A,B,C,alpha,b,u0,e,t)

% u = 1/b*u0 + 1/b*( real(-2*A+2/3*alpha^2*(e(1)+A*t^2+B*t+C)^(1/3))*(t<=T) + real(2/3*(t>T)*alpha^2*e(1)^(1/3)) );
u = u0/b + 1/b*( (-2*A+2/3*alpha^2*fixPower(e(1)+A*t^2+B*t+C))*(t<=T) + (2/3*(t>T)*alpha^2*fixPower(e(1))) );

end

function n = fixPower(m)

n = sign(m)*abs(m)^(1/3);

end