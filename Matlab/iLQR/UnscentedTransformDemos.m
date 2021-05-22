%% UT Demos
clear; clc;


% Simple sine function, linear variance produces poor variance results
theta = -15;
std = 1;

x = theta + std*randn(10000, 1);

y = sind(x);

[z,w] = UnscentedTransform(theta, std^2, 3);

mc_mean = mean(y);
mc_var = var(y);
linear_mean = sind(theta);
linear_var = cosd(theta)^2 * std^2;

yz = sind(z)';
ut_mean = sum(yz.*w);
ut_var = sum(w.*(yz-ut_mean).^2);
disp('      MC      Linear    UT')  
disp([mc_mean, linear_mean, ut_mean])
disp([mc_var, linear_var, ut_var].^0.5)

%% Covariance elements matter, right? (Answer is yes)
clear
x = [0,0]';
P = [1, 0.9; 0.9, 1];

[z,w] = UnscentedTransform(x, P, 3);
 figure
 plot(z(1,:),z(2,:), 'o')
 
 x = [0,0]';
P = [1, 0; 0, 1];

[z,w] = UnscentedTransform(x, P, 3);
hold all
plot(z(1,:),z(2,:), 'ko')