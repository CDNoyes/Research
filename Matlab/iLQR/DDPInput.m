function input = DDPInput(weights)

input.ddp = true;
% Initial State and Covariance 
% input.v0 = 5525;
% input.x0 = [0 , 54.5e3, -11.5*pi/180, 0]';
input.v0 = 5461.4;
input.x0 = [0, 39.4497e3, (-10.604*pi/180), 0]';   %

% input.P0 = diag([2500, 0.25*pi/180, 5e3]).^2; % SRL-like
% input.P0 = diag([3500, 0.35*pi/180, 10e3]).^2; % The added altitude variance is not suitable with the lower initial vel
input.P0 = diag([2500, 0.25*pi/180, 10e3]).^2;
input.P0_param = 1*diag(20/100/3*ones(3,1)).^2;
input.vf = 460; 

input.weights = weights; % objective function weights
input.ut_scale = 17; % Unscented Transform scale factor for computing sigma
input.optimize_gains = true;
input.closed_loop = true;
input.gains = [0.0725, -0.025, -4]; % D S fpa for the heavy vehicle
% input.gains = [0.1,-0.035, -5 ]; % for MSL-like vehicle 

% Vehicle Info
input.bounds = [0, 1]; % feed forward control limits 

% Algorithmic Parameters
input.terminal_plots = true;
input.running_plots = 1; % set to -1 to include deriv plots
input.max_iterations = 100;
input.horizon = 1000;
input.qn = false; % default for now 
input.parallel = true;
input.guess = [];
end