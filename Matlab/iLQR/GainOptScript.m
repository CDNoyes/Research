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
save('ConvergenceExample', 'sols');
%%
load ConvergenceExample
fields = fieldnames(sols{1}.trace);
for j = 2:length(fields)
    for i = 1:3
        trace = sols{i}.trace;
        iter = 1:trace(end).iter;
        y = [];
        for k = iter
           y(k) =  getfield(sols{i}.trace(k), fields{j});
        end
        figure(j)
        hold all
        plot(iter, y, 'linewidth', 2)
        xlabel('Iteration')
        ylabel(fields{j})
    end
    grid on
    legend('DDP','DDP-ControlHessian','iLQR','location','best')
end
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
fname = 'msl_weight_sweep';
fname_new = [fname,'_optimized'];
load(fname);
inp = DDPInput([0,0,0.2]);


for i = 1:length(sols)
    inp.terminal_plots = false;
    inp.running_plots = false;
    inp.horizon = 250;
    inp.ut_scale = [0,5,10];
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
%% Some solutions
% alpha = 5
% Kopt =
%
%     0.0582   -0.0215   -0.0226*200
%
% hf = 7.2347 +/- 3*0.80172 km (3s low = 4.8296)
% sf = 327.7621 +/- 3*0.51265 km
%
%
% alpha = 10
% Kopt =
%
%     0.0634   -0.0227   -0.0149*200
%
% hf = 7.2718 +/- 3*0.83534 km (3s low = 4.7657)
% sf = 327.2111 +/- 3*0.58646 km
%
% alpha = 15
% Kopt =
%
%     0.0703   -0.0244   -0.0091*200
%
% hf = 7.2849 +/- 3*0.88727 km (3s low = 4.6231)
% sf = 326.8169 +/- 3*0.7663 km
%
% alpha = [5, 15]
% Kopt =
%
%     0.0708   -0.0251   -0.0186*200
%
% hf = 7.2448 +/- 3*0.84436 km (3s low = 4.7117)
% sf = 327.776 +/- 3*0.86039 km
%
% alpha = [5, 10, 15]
% Kopt =
%
%     0.0684   -0.0244   -0.0176*200
%
% hf = 7.253 +/- 3*0.84166 km (3s low = 4.728)
% sf = 327.5989 +/- 3*0.79188 km
%
% alpha = [5,10,15] with initial guess K = 0
% Kopt =
%
%     0.0684   -0.0244    -3.5220
%
% hf = 7.253 +/- 3*0.84166 km (3s low = 4.728)
% sf = 327.5989 +/- 3*0.79188 km