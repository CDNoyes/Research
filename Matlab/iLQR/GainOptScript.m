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