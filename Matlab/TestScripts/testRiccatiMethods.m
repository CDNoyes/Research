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
A = @(x) [0, 1; m*g*l*ReplaceNAN(sin(x(1))/x(1),1)/I, -b/I];
B = [0; 1/I];
r = @(t) sin(t);
R = 1;
C = [1,0]; %Track the position
Q = @(x) (100000);
F = 10;
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
A = @(x) [0, 1; m*g*l*ReplaceNAN(sin(x(1))/x(1),1)/I, -b/I];
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
A = @(x) [0, 1; m*g*l*ReplaceNAN(sin(x(1))/x(1),1)/I, -b/I];
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


%% Inverted Pendulum Tracking with constrained input
clc; clear;
m = 1; %kg
l = 0.5; %m
b = .1;
I = m*l^2;
g = 9.81;
[f,j] = InvertedPendulum();
% A = @(x) [0, 1, 0; m*g*l*ReplaceNAN(sin(x(1))/x(1),1)/I, -b/I, ReplaceNAN(Saturate(x(3),-5,5)/x(3),0)/I;0,0,0];
% A = @(x) [0, 1, 0; m*g*l*ReplaceNAN(sin(x(1))/x(1),1)/I, -b/I, 1/I;0,0,0];
%Try the "integral servomechanism" method
A = @(x) [0, 1, 0, 0; m*g*l*ReplaceNAN(sin(x(1))/x(1),1)/I, -b/I, ReplaceNAN(Saturate(x(3),-5,5)/x(3),0)/I, 0; 0,0,0,0;1, 0, 0, 0];

B = [0; 0; 1e2;0];
r = @(t) [sin(t);cos(t);-cos(t)];
R = 1;
C = [1,0,0,0;0,1,0,0;0,0,0,1]; %Track the position, its derivative, and its integral
Q = @(x) diag([10000,1000,0]);
l = size(C,1);
F = diag([10,1,1]);
tf = 10;
x0 = [0;0;5;0];
sol = ASRE(x0,tf,A,B,C,Q,R,F,r);

figure
subplot 211
semilogy(sol.time,abs(sol.state(:,1:l)-r(sol.time)'))
ylabel('Error')
subplot 212
plot(sol.time,sol.state(:,1))
hold all
plot(sol.time,sol.state(:,2))
plot(sol.time,r(sol.time)')
legend('Position','Velocity','Desired')
figure
plot(sol.time,Saturate(sol.state(:,3),-5,5))
figure
plot(sol.time,sol.control)