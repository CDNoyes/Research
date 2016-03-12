% SMO Sliding mode observer dynamics for drag tracking.
%	SMO(x,D,a,b,u,k,alpha) Computes the change in estimates of drag, drag
%	rate, and disturbance. The inputs are:
%     x - the vector of observer states
%     D - the measured drag acceleration
%     a,b - Drag dynamics, D_ddot = a(x) + b(x)u
%     u - cos(gamma), the control parameter
%     k - vector of linear gains
%     alpha - vector of nonlinear gains

function dx = SMO(x,D,a,b,u,k,alpha)

e = D-x(1);
% signe = sign(e);
signe = tanh(10*e); %smooth approximation to sign function
dx = [ x(2) + alpha(1)*e + k(1)*signe
       x(3) + alpha(2)*e + k(2)*signe + a + b*u
              alpha(3)*e + k(3)*signe ];

end


