function [dx,u,d_hat] = TerminalSMCDynamics(t,x)

%Extensive number of parameters to tune
k = 2500;
beta = 4;
epsilon = 0.5;
p0 = 5;
q0 = 9;
alpha1 = 50;
beta1 = 0.5;
p1 = 5;
p2 = 5;
q1 = 7;
q2 = 7;
delta = 60;
mu = 0.8;

bl = 0.8;

yd = 2*sin(0.2*t);
yd_dot = 0.2*2*cos(0.2*t);
yd_ddot = -0.2^2*2*sin(0.2*t);
d = 2*sin(0.1*pi*t)+3*sin(0.2*sqrt(t+1));
% d = 0;

s = x(3)-x(2);
s1 = x(1)-yd;
s2 = x(2)-yd_dot + alpha1*s1 + beta1*s1^(p1/q1)+s;
d_hat = -k*s-beta*Saturate(s/bl,-1,1)-epsilon*s^(p0/q0)-abs(-2*x(1)+3*(1-x(1)^2)*x(2))*Saturate(s/bl,-1,1)+2*x(1)-3*(1-x(1)^2)*x(2);
d_hat = real(d_hat);
vr = 2*x(1)-3*(1-x(1)^2)*x(2) + yd_ddot-alpha1*(x(2)-yd_dot)-beta1*(x(2)-yd_dot)^(p1/q1)-d_hat-delta*s2-mu*s2^(p2/q2);
u = Saturate(real(vr),-10,10);
du = u-vr;

dx = [  x(2)
        -2*x(1)+3*(1-x(1)^2)*x(2) + u + d
        -k*s-beta*Saturate(s/bl,-1,1)-epsilon*s^(p0/q0)-abs(-2*x(1)+3*(1-x(1)^2)*x(2))*Saturate(s/bl,-1,1)+vr];
dx = real(dx);
end