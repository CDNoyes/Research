%% UT Demos
clear; clc;


% Simple sine function
theta = -15; % mean, degrees
std = 5; 
UnscentedTransformDemoFunction(theta, std);

% The conclusion is that linearization works quite well over a variety of
% mean values of x but performs increasingly poorly for large variance.

%% Covariance elements matter, right? (Answer is yes)
clear
x = [0,0]';
P = [1, 0.9; 0.9, 1];

[z,w] = UnscentedTransform(x, P, 3);
 figure
 plot(z(1,:),z(2,:), 'x')
 
 x = [0,0]';
P = [1, 0; 0, 1];

[z,w] = UnscentedTransform(x, P, 3);
hold all
plot(z(1,:),z(2,:), 'ko')
title('Sigma Points')