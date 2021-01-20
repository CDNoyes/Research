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

disp([mc_mean, linear_mean, ut_mean])
disp([mc_var, linear_var, ut_var].^0.5)