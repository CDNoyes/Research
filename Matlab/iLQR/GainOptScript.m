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
fname = 'msl_weight_sweep';
fname_new = [fname,'_optimized'];
load(fname);
inp = DDPInput([0,0,0.2]);


for i = 1:length(sols)
    inp.terminal_plots = false;
    inp.running_plots = false;
    inp.horizon = 250;
    inp.ut_scale = 4;
    inp.weights = sols{i}.weights;
    inp.guess = sols{i}.u;
    sols{i} = FixedGainOpt(inp);
    gains(i,:) = sols{i}.input.gains;
end

save(fname_new, 'sols')