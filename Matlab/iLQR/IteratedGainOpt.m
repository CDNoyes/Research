clear; clc;

% choose weights
%% Separate Reference Control and Gain Optimizations
wh = 1;
ws = 1;
wu = 0.3;

% choose alpha 1 and 2
for alpha_ref = 5:2:17
    % opt reference control
    input = DDPInput([wh,ws,wu]);
    input.ut_scale = alpha_ref;
    input.horizon = 250;
    input.n_controls = 1;
    input.max_iterations = 100;
    input.terminal_plots = 0;
    sol_ref = entry_stochastic_gains_params(input);
    
    for alpha_gain = 5:2:17
        % opt gains
        input.n_controls = 3;
        input.ut_scale = alpha_gain;
        input.guess = sol_ref.u;
        input.max_iterations = 100;
        sol = entry_stochastic_gains_params(input);
        sol.reference_solution = sol_ref;
        save(['./alpha_sweep/sol_',num2str(alpha_ref),'_',num2str(alpha_gain)],'sol')
    end
end


%% single optimization of gains and overcontrol
clear; 

wh = 1;
ws = 1;
wu = 0.3;
sols = {};
for alpha = 5:17
    input = DDPInput([wh,ws,wu]);
    input.horizon = 250;
    input.max_iterations = 150;
    input.terminal_plots = 0;
    input.n_controls = 2;
    input.ut_scale = alpha;
%     input.bounds = [0.1, 0.9]; %leave margins for feedback...
    sol = entry_stochastic_gains_params(input);
    sols{end+1} = sol;
end

save('./alpha_sweep/cov_overcontrol_sols','sols')
%% collect data in single file, also append gains to reference control
sols = {};
list = dir('./alpha_sweep/*.mat');
for i = 1:length(list)
    load([list(i).folder,'\', list(i).name])
    sol.u = [sol.reference_solution.u; sol.u];
    sol.ut_scale_pair = [sol.reference_solution.input.ut_scale, sol.input.ut_scale];
    sols{end+1} = sol;
end

save('./alpha_sweep/alpha_sweep.mat','sols')

figure(1)
hold all
plotted = [];
leg = {};
for i = 1:length(sols)
    sol = sols{i};
    if ~any(plotted == sol.ut_scale_pair(1))
        plot(sol.v(1:end-1), smoothdata(sol.u(1,:),'movmean', 5), 'linewidth', 2)
        
        plotted(end+1) = sol.ut_scale_pair(1);
        leg{end+1} = num2str(plotted(end));
    end
end
legend(leg{:})
% monte carlo in python for results