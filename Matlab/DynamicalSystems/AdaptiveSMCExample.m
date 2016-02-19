clc; clear;
x0 = [0;0;-.35];
X0 = [x0;x0;5];

[t,X] = ode45(@AdaptiveSMCDynamics,[0,10],X0);

figure
plot(t,X(:,1:3))
hold all
plot(t,[sin(t),cos(t),-sin(t)])

figure
plot(t,X(:,7))
title('Gain convergence')

figure
plot(t,X(:,1:3) - [sin(t),cos(t),-sin(t)])
title('Tracking Error')

e = (X(:,1:3) - [sin(t),cos(t),-sin(t)])';
c = [12,7,1];
for i = 1:length(t)
    s(i) = c*e(:,i);
    
end

figure
plot(t,s)
title('Sliding Surface')