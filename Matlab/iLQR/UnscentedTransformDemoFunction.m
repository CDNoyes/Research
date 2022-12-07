function output = UnscentedTransformDemoFunction(theta, std)

N = 1e4; % number of monte carlo samples
x = theta + std*randn(N, 1); % Sample the inputs 
y = sind(x);  % MC output 

[z,w] = UnscentedTransform(theta, std^2, 3);

mc_mean = mean(y);
mc_var = var(y);
linear_mean = sind(theta);
linear_var = cosd(theta)^2 * std^2 * (pi/180)^2;

yz = sind(z)';
ut_mean = sum(yz.*w);
ut_var = sum(w.*(yz-ut_mean).^2);

output.mean = [mc_mean, linear_mean, ut_mean];
output.std = [mc_var, linear_var, ut_var].^0.5;

output.error_mean = abs([mc_mean-linear_mean, mc_mean-ut_mean]);
output.error_std = abs([mc_var.^0.5-linear_var.^0.5, mc_var.^0.5-ut_var.^0.5]);

output.error_mean_per = output.error_mean./output.mean(1) * 100;
output.error_std_per = output.error_std./output.std(1) * 100;



disp('Test Function:')
disp(['y = sind(x) with x ~ N(', num2str(theta),', ',num2str(std),')'])
disp(['N = ', num2str(N),' for monte carlo results'])
disp(" ")
disp('Results:')
disp(['     MC       Linear     UT'])  
disp(output.mean)
disp(output.std)

disp('Absolute Errors (relative to MC):')
disp(['     Linear     UT'])  
disp(output.error_mean)
disp(output.error_std)

disp('Percent Errors (relative to MC):')
disp(['     Linear     UT'])  
disp(output.error_mean_per)
disp(output.error_std_per)

