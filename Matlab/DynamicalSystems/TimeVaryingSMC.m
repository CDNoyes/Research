function TimeVaryingSMC()

x0 = [0.5;0.4];
e0 = x0(1);
limit = 8;
margin = limit-6;
[T,A,B] = DesignSM(e0,margin);

[t,x] = ode45(@dynamics,[0,100],x0,[],A,B,T);
for i = 1:length(t)
    [~,u(i),s(i),d(i)] = dynamics(t(i),x(i,:)',A,B,T);
end

figure
plot(t,[2*sin(0.2*t),0.2*2*cos(0.2*t)])
hold all
plot(t,d)
plot(t,x)
legend('Ref','Ref rate','Disturbance')

figure
plot(t,u)
ylabel('u')

figure
plot(t,s)
title('Sliding Surface')
end

function [dx,u,s,d] = dynamics(t,x,A,B,T)

%System format
f = -2*x(1)+3*(1-x(1)^2)*x(2);
b = 1;

% Sliding Surface Design Parameters
beta = 2;
c = [beta,1];

% Reference Trajectory
yd = 2*sin(0.2*t);
yd_dot = 0.2*2*cos(0.2*t);
yd_ddot = -0.2^2*2*sin(0.2*t);

% Disturbance
d = 2*sin(0.1*pi*t)+3*sin(0.2*sqrt(t+1));
% d = 0;
dmax = 5; %Estimated bound of the disturbance

% Error states
e = x-[yd;yd_dot];

%Sliding Surface Computation
s = c*e + (A*t+B)*(t<=T); 

% Control
u0 = yd_ddot-f - dmax*Saturate(s/.01,-1,1);
u = 1/b*u0 + 1/b*(-beta*e(2)-A*(t<=T));

dx = [  x(2)
        f + b*u + d];
    
end