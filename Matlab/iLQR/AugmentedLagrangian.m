%% First example solving a problem with augmented lagragian technique

inp = DDPInput([0,0,0]);
inp.max_iterations = 30; % want a fairly low number since we have to solve repeatedly
inp.terminal_plots = false;
inp.horizon = 1000;

lagrange = [0;0];
penalty = 0.001;
penalty_factor = 10; % multiply penalty by this value each time 

for i = 1:10
    disp(['Iteration ', num2str(i)])
    sol = entry_3d(inp, lagrange, penalty);

    inp.guess = sol.u; % so that next iteration uses prev iteration sol
    lagrange = lagrange + penalty*sol.constraint
    penalty = penalty * penalty_factor

    if max(abs(sol.constraint)) < 1e-4  % if the crossrange is under _ km, good enough 
        break
    end
    
end

%% GPOPS Solution:

% Trajectory Duration: 160.3306 s
% Final altitude:      5.6143 km
% Final velocity:      460 m/s
% Entry FPA:           -11.5 deg
% Max   FPA:           2.1167 deg
% Final FPA:           -9.5325 deg
% Final Downrange:     349.9176 km
% Final Crossrange:    0 km