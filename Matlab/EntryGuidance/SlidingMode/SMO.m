% SMO Sliding mode observer dynamics for a second order system.
%	SMO(x,D,a,b,u,k,alpha) Computes the change in estimates of drag, drag
%	rate, and disturbance. The inputs are:
%     x - the vector of observer states
%     D - the measured drag acceleration (i.e. the reference variable)
%     a,b - scalar dynamic terms, D_ddot = a(x) + b(x)u
%     u - the control parameter (u=cos(gamma) for entry)
%     k - vector of linear gains
%     alpha - vector of nonlinear gains
%
%     The gains k and alpha can be computed using computeSMOGains.
%     Although the code itself and the help reference using this for drag
%     tracking entry guidance, this SMO is applicable to any second order
%     system. See TimeVaryingSMC for an example where it is applied to a
%     vanderpol oscillator.
function dx = SMO(x,D,a,b,u,k,alpha)

e = D-x(1);
% signe = sign(e);
signe = tanh(10*e); %smooth approximation to sign function
dx = [ x(2) + alpha(1)*e + k(1)*signe
       x(3) + alpha(2)*e + k(2)*signe + a + b*u
              alpha(3)*e + k(3)*signe ];

end


