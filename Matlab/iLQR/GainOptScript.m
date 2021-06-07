%% Generate a single gain function opt to show bad results!
inp = DDPInput([1,1,0.2]);
inp.ddp = true;
inp.ut_scale = 15;
inp.n_controls = 4;
% inp.guess = sol.u;
sol = entry_stochastic_gains_params(inp);

save('MSLGainOptGoneWrongLowGain.mat', 'sol')

%% What happens if I use large weights? Or increased control cost near end
inp = DDPInput([3,3,3]);
inp.ddp = true;
inp.ut_scale = 9;
sol = entry_stochastic_gains_params(inp);
save('HighControlWeight.mat','sol')
sols{4} = sol;
inp.guess = sol.u;
for i = 0:2
    inp.weights(1) = i;
    sols{i+1} = entry_stochastic_gains_params(inp);

end
save('HighControlWeights.mat','sols')
%%
inp = DDPInput([1,1,0.2]);
inp.horizon = 250;
% inp.guess = ones(1,inp.horizon)*0.8;
% inp.guess(linspace(5505,460,inp.horizon) > 3000) = 0.2;

% Optimize gains for the new initial state/cov
inp.ut_scale = [0, 5, 10];
inp.terminal_plots = true;
sol_gains = FixedGainOpt(inp);


inp.ut_scale = 15;
inp.gains = sol_gains.input.gains;
sol = entry_stochastic_gains_params(inp);

inp.ut_scale = [0, 5, 10];
inp.terminal_plots = false;
inp.guess = sol.u;
sol_est = FixedGainOpt(inp, false);
sol_opt = FixedGainOpt(inp);
sols = {sol, sol_est, sol_opt};
save('MSLDetailedExample.mat', 'sols')

%% Use this block to re-optimize a solution in 'sol'
inp = DDPInput([1,1,0.2]);
inp.ut_scale = [0, 5, 10];
inp.horizon = 250;
inp.terminal_plots = false;
inp.guess = sol.u;
sol_opt = FixedGainOpt(inp);


%% Estimate the open loop dispersions for the heavy vehicle, guess traj
inp = DDPInput([0,0,0]);
inp.ut_scale = [0, 5, 10];
inp.ut_scale = 15;
inp.horizon = 250;
inp.gains = [0,0,0];
inp.terminal_plots = true;
sol = FixedGainOpt(inp, 0);

% SolutionPlots(sol)

% SolutionPlots(sol,'HeavyOpenLoop')

%% Call the Fixed GainOpt function with K=0 and K=input.gains
% sweep over alpha maybe? need to save one of each for use in python

% technically, for the closed loop especially, the linearization should use
% the nominal trajectory as the reference, set the weights on all other
% sigma points to zero
clear; clc
inp = DDPInput([0,0,0]);
inp.ut_scale = 12;
inp.horizon = 250;
inp.lod = 0.24;
sol = FixedGainOpt(inp);

% save('UTExample_ClosedLoop', 'sol')

% inp.gains = [0,0,0]';
% sol = FixedGainOpt(inp);

% save('UTExample_OpenLoop', 'sol')
%% Gather convergence data for three solutions, param unc only
inp = DDPInput([1,1,0.2]);
inp.ddp = true;
inp.full_hessian = true;
inp.ut_scale = 15;

sol_ddp = entry_stochastic_gains_params(inp);

inp.full_hessian = false;
sol_partial = entry_stochastic_gains_params(inp);

inp.ddp = false;
inp.max_iterations = 250;
sol_lqr = entry_stochastic_gains_params(inp);

sols = {sol_ddp, sol_partial, sol_lqr};
% save('ConvergenceExample', 'sols');

%% Convergence from a different guess?
inp = DDPInput([1,1,0.2]);
inp.ddp = true;
inp.ut_scale = 15;
inp.horizon = 1000;
inp.full_hessian = false;

inp.guess = 0.5 + zeros(1,inp.horizon);
sol_half = entry_stochastic_gains_params(inp);

inp.guess = ones(1,inp.horizon);
sol_full = entry_stochastic_gains_params(inp);

inp.guess = linspace(1,0,inp.horizon);
sol_rampdown = entry_stochastic_gains_params(inp);
sols = {sol_half, sol_full, sol_rampdown};
% save('ConvergenceGuessExample', 'sols');

% Make table of final cost and total iterations
%%
inp = DDPInput([1,1,0.2]);

load solutions_cl_ddp_max_guess

inp.terminal_plots = false;
inp.running_plots = true;
inp.horizon = 250;
inp.ut_scale = 15;
inp.guess = sols{6}.u;
clear sols
sols{1} = entry_stochastic_gains_params(inp);

%% Add the reoptimized gains

inp.ut_scale = 4;
inp.guess = sols{1}.u;
sols{2} = FixedGainOpt(inp);

save('DetailedExample', 'sols')

%% Load a file and it creates a new one with _optimized
% fname = 'margin_comparison_0p24';
fname = 'msl_weight_sweep_high_control';
fname_new = [fname,'_optimized'];
load(fname);
inp = DDPInput([0,0,0.2]);


for i = 1:length(sols)
    inp.terminal_plots = false;
    inp.running_plots = false;
    inp.horizon = 250;
    inp.ut_scale = [3,9,15];
    inp.weights = sols{i}.weights;
    inp.guess = sols{i}.u;
    sols{i} = FixedGainOpt(inp);
    gains(i,:) = sols{i}.input.gains;
end

save(fname_new, 'sols')
%% Load a file and it creates a new one with _reestimated
fname = 'msl_weight_sweep_open_loop';
fname_new = [fname,'_reestimated'];
load(fname);
inp = DDPInput([0,0,0.2]);


for i = 1:length(sols)
    inp.terminal_plots = false;
    inp.running_plots = false;
    inp.horizon = 250;
    inp.ut_scale = [0,5,10];
    inp.weights = sols{i}.weights;
    inp.guess = sols{i}.u;
    sols{i} = FixedGainOpt(inp,false);
    stats(i,:) = sols{i};
end

save(fname_new, 'sols')

%% Look at sensitivity to alpha in a single case
clear
load UTExample_ClosedLoop

inp = DDPInput([0,0,0.2]);
inp.gains = sol.input.gains;
alpha = 0:20;

for i = 1:length(alpha)
    inp.terminal_plots = false;
    inp.running_plots = false;
    inp.horizon = 250;
    inp.ut_scale = alpha(i);
    sol = FixedGainOpt(inp,false);
    stats(i,:) = [sol.mean(1,end)/1000, sol.mean(3,end), sol.var(1,end)^0.5/1000, sol.var(3,end)^0.5];
end
stats(:,2) = stats(:,2) - 300;
figure
plot(alpha, stats, 'linewidth', 2)
axis([0, alpha(end), 0, 13])
xlabel('UT Scaling Parameter, \alpha', 'fontsize', 14)
ylabel('Statistics (km)', 'fontsize', 14)
legend('Mean Altitude','Mean Range - 300 km','STD Altitude','STD Range','location','best')
grid on
saveas(gcf, 'C:\Users\Aero\Documents\EDL\Documents\Dissertation\Images\AlphaSensitivity.png')
