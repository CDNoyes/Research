%% F8 Crusader Regulation Problem
clear; clc;

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
sol = ASRE(x0,tf,A,B,[],Q,R,F);

figure
plot(sol.time,sol.state)
figure
plot(sol.time,sol.control)

sol2 = SDRE(x0,tf,A,B,[],Q,R,F);

figure
plot(sol2.time,sol2.state)
figure
plot(sol2.time,sol2.control)

return

%% Inverted Pendulum Tracking
clc; clear;
m = 1; %kg
l = 0.5; %m
b = .1;
I = m*l^2;
g = 9.81;
[f,j] = InvertedPendulum();
A = @(x) [0, 1; m*g*l*ReplaceNAN(sin(x(1))/x(1),0)/I, -b/I];
B = [0; 1/I];
r = @(t) sin(t);
R = .07;
C = [1,0]; %Track the position
Q = @(x) (1000);
F = 1;
tf = 10;
x0 = [0;0];
sol = ASRE(x0,tf,A,B,C,Q,R,F,r);
sol2 = SDRE(x0,tf,A,B,C,Q,R,F,r);

figure
plot(sol.time,sol.state(:,1))
hold all 
plot(sol2.time,sol2.state(:,1))
plot(linspace(0,tf,50),r(linspace(0,tf,50)),'ko')
legend('ASRE Trajectory','SDRE Trajectory','Reference')
figure
plot(sol.time,sol.control)
hold all
plot(sol2.time,sol2.control)



%% Inverted Pendulum Swing up via Regulation

clc; clear;
m = 1; %kg
l = 0.5; %m
b = .1;
I = m*l^2;
g = 9.81;
[f,j] = InvertedPendulum();
A = @(x) [0, 1; m*g*l*ReplaceNAN(sin(x(1))/x(1),0)/I, -b/I];
B = [0; 1/I];
R = 300;
Q = @(x) 0*eye(2);
F = 100*eye(2);
tf = 3;
x0 = [pi;0];
sol = ASRE(x0,tf,A,B,[],Q,R,F,[]);
sol2 = SDRE(x0,tf,A,B,[],Q,R,F,[]);


figure
plot(sol.time,sol.state)
hold all
plot(sol2.time,sol2.state)
legend('Angle (ASRE)','Angular Velocity (ASRE)','Angle (SDRE)','Angular Velocity (SDRE)')
figure
plot(sol.time,sol.control)
hold all
plot(sol2.time,sol2.control)
legend('ASRE','SDRE')

%% Inverted Pendulum Swing up via final constraints

clc; clear;
m = 1; %kg
l = 0.5; %m
b = .1;
I = m*l^2;
g = 9.81;
[f,j] = InvertedPendulum();
A = @(x) [0, 1; m*g*l*ReplaceNAN(sin(x(1))/x(1),0)/I, -b/I];
B = [0; 1/I];
R = 1;
Q = @(x) 0*eye(2);
F = 0*eye(2);
tf = 3;
x0 = [pi;0];
xf = [0];
C = [1,0];
sol = ASRE(x0,tf,A,B,C,Q,R,F,xf);

figure
plot(sol.time,sol.state)

figure
plot(sol.time,sol.control)