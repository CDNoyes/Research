%%
% Conclusions:
% When reoptimizing gains, both low altitude and range errors were monotonic 
% wrt to the UT scaling parameter. The further out the sigma points are,
% the more variance is predicted. 

% Using Solution 6 for the guess and ut = 4 produced a MC with less than
% 2 km 3sigma range error!

%%
clear; clc;
inp = DDPInput([0,0,0]);
ut = 3:19;
inp.horizon = 250;
load solutions_cl_ddp_max_guess
% 6 = [1,1]
% 4 = [3,0]
sol_num = 6;
inp.guess = sols{sol_num}.u;
inp.terminal_plots = false;
inp.weights = sols{sol_num}.weights;
clear sols
for i = 1:length(ut)
    inp.ut_scale = ut(i);
    sols{i} = FixedGainOpt(inp);
    Kopt(i,:) = sols{i}.input.gains;
end
    close all
save('UTScale_max_alt.mat', 'sols')
%%
plot(ut, Kopt)

%% TODO:
% run MC on all of the sols. This represents a sweep over the UT parameter
% but with optimized gains for each run. If the results are always sufficiently
% incorrect, maybe we need to go to optimizing gains for one UT scale,
% and optimizing the trajectory for another one. 
% Should the gains be optimized with a low UT and the trajectory with a
% larger value, or vice versa. 

% From one test so far, optimizing the trajectory with a high UT (14) and
% optimizing the gain with low UT (4) worked really well. 


%%

load solutions_cl_ddp_max_guess
inp = DDPInput([0,0,0]);
inp.ut_scale = 4;
inp.horizon = 250;
inp.terminal_plots = false;

for sol_num = 1:length(sols)
    
    inp.guess = sols{sol_num}.u;
    inp.weights = sols{sol_num}.weights;
    inp.weights(2) = max(inp.weights(2), 0.01);
    sols_new{sol_num} = FixedGainOpt(inp);
end
sols = sols_new;
save('solutions_cl_ddp_optimized.mat', 'sols')

%% 
% Purpose: reopt the gains for case 6 to make the range errors really
% small, then maybe use those gains for 6 cases to keep them fixed for the
% comparison

% Note that we can also reoptimize the gains for each trajectory for a
% seprate comparison. Even the nominal optimization cases, since the gain
% opt uses the UT transform. Next section does this.

load margin_comparison
inp = DDPInput([0,0,0]);
inp.ut_scale = 4;
inp.horizon = 250;
inp.terminal_plots = false;

for sol_num = 6
    
    inp.guess = sols{sol_num}.u;
    inp.weights = sols{sol_num}.weights;
    inp.weights(2) = max(inp.weights(2), 0.01);
    sols{sol_num} = FixedGainOpt(inp);
end

for i = 1:5
   sols{i}.input.gains = sols{sol_num}.input.gains; 
end

save('reference_comparison_fixed_gains.mat', 'sols')

%%
load margin_comparison
inp = DDPInput([0,0,0]);
inp.ut_scale = 4;
inp.horizon = 250;
inp.terminal_plots = false;
for sol_num = 1:6
    
    v = sols{sol_num}.v;
    u = interp1(v(2:end), sols{sol_num}.u, linspace(v(1), v(end), inp.horizon),'spline');
    inp.guess = u;
    inp.weights = sols{sol_num}.weights;
    inp.weights(2) = max(inp.weights(2), 0.01);
    sols{sol_num} = FixedGainOpt(inp);
end

save('reference_comparison_opt_gains.mat', 'sols')