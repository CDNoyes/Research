function input = DDPInput(weights)

input.ddp = true;
input.full_hessian = false;
% Initial State and Covariance 
% input.v0 = 5525;
% input.x0 = [0 , 54.5e3, -11.5*pi/180, 0]';
if 0% heavy vehicle
    input.v0 = 5461.4;
    input.x0 = [0, 39.4497e3, (-10.604*pi/180), 0]';   
else % MSL-class
    input.v0 = 5500;
    input.x0 = [0, 53.5e3, -13*pi/180, 0]';   
end

% input.P0 = diag([1000, 0.25*pi/180, 10e3]).^2;
input.P0 = diag([2500, 0.25*pi/180, 10e3]).^2; % Values used for much analysis
input.P0_param = 1*diag([0.05, 0.05, 0.07]).^2;
input.vf = 460; 

input.weights = weights; % objective function weights
input.ut_scale = 15; % Unscented Transform scale factor for computing sigma
input.optimize_gains = true;
input.closed_loop = true;
% input.gains = [0.0725, -0.025, -4]; % D S fpa for the heavy vehicle
% input.gains = [0.1305   -0.0389   -2.5205 ]; % for MSL-like vehicle 
input.gains = [0.0636   -0.0307   -0.0744]; % for MSL from higher state
input.n_controls = 1;

% Vehicle Info
input.bounds = [0, 1]; % feed forward control limits 
input.lod = 0.24; % Nominal LoD at mach=24, scales the MSL profile 

% Algorithmic Parameters
input.terminal_plots = true;
input.running_plots = 1; % set to -1 to include deriv plots
input.max_iterations = 100;
input.horizon = 250;
input.qn = false; % default for now 
input.parallel = true;
input.guess = [];
end