%Sliding mode observer for drag tracking
% x is the vector of observer states
% D is the measured drag acceleration
% a and b are the FBL terms for drag tracking
% k is a vector of linear gains
% alpha is a vector of nonlinear gains

function dx = SMO(x,D,a,b,u,k,alpha)

e = D-x(1);
% signe = sign(e);
signe = tanh(10*e); %smooth approximation to sign function
dx = [ x(2) + alpha(1)*e + k(1)*signe
       x(3) + alpha(2)*e + k(2)*signe + a + b*u
              alpha(3)*e + k(3)*signe ];

end


