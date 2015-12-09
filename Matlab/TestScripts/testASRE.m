%% F8 Crusader
clear

A = @(x) [-0.877 + 0.47*x(1) + 3.846*x(1)^2-x(1)*x(3), -0.019*x(2), 1-0.088*x(1)
            0, 0, 1
            -4.208-0.47*x(1)-3.56*x(1)^2, 0, -0.396];
B = @(x,u) [-0.215+0.28*x(1)^2+0.47*x(1)*u + 0.63*u^2
                0
            -20.967 + 6.265*x(1)^2 + 46*x(1)*u + 61.4*u^2];

F = .1*eye(3);
Q = .01*eye(3);
R = 1;
x0 = [0.50;0;0];
tf = 20;
sol = ASRE(x0,tf,1,A,B,[],Q,R,F);

figure
plot(sol.time,sol.state{end})
figure
plot(sol.time,sol.control{end})

return
%% Inverted Pend

m = 1; %kg
l = 0.5; %m
b = .1;
I = m*l^2;
g = 9.81;
f = InvertedPendulum();
A = @(x) [0, 1; m*g*l*ReplaceNAN(sin(x(1))/x(1),0)/I, -b/I];
B = [0; 1/I];