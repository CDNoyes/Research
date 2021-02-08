function generate_solutions
sdir = 'E:\Documents\EDL\Documents\PropellantOptimalJournal\ddp\matlab\';

figs = {'Range','Altitude','Fpa','Control','Drag'};

%% Single Data Point UT Validation, find a decent solution with fixed gains, determine appropriate beta
wh = 1;
ws = 1;
inp = DDPInput([wh, ws, 0.3]);
inp.terminal_plots = false;
for i = 11:20
    i
    inp.ut_scale = i;
    sol = entry_stochastic(inp);
    [hm, fm, sm, hv, fv, sv]=stats(sol);
    sol.stats = [hm, fm, sm, hv, fv, sv]; % terminal stats
    if size(sol.u,1) == 4
        sol = UTSolve(sol, sol.u(2:4,:));
    else
        sol = UTSolve(sol, sol.input.gains);
    end
    save(['E:\Documents\GitHub\Research\Matlab\\iLQR\solutions\beta_sweep\MatlabBeta',num2str(inp.ut_scale)], 'sol');
    save(['E:\Documents\GitHub\Research\Matlab\\iLQR\solutions\beta_sweep\Beta',num2str(inp.ut_scale)], '-struct','sol');
end
disp('')
%% Unscented Transform "Validation"

if 1
    sols = {};
    ws = linspace(0, 0.4, 5);
    ws = [0.025, 0.05, 0.075];
    for i = 1:length(ws)
        i
        inp = DDPInput([2, ws(i), 0.3]);
        inp.running_plots = 0;
        inp.terminal_plots = 0;
        inp.horizon = 2000;
        sol = entry_stochastic_gains(inp);
        [hm, fm, sm, hv, fv, sv]=stats(sol);
        sol.stats = [hm, fm, sm, hv, fv, sv]; % terminal stats
        sols{end+1} = UTSolve(sol, sol.u(2:4,:));
    end
    
    save('UTValidationSet2000.mat','sols');
    %     save('UTValidationSetQuu.mat','sols'); % Retains only Quu hessian terms
else
    load UTValidationSet2000
end
disp('');

for i = 1:length(sols)
    S(i,:) = [sols{i}.weights(2), sols{i}.var(3,end)^0.5, sols{i}.ut.stats(3)/3];
    %     UTCompare(sols{i});
end
figure
semilogy(S(:,1), S(:,2),'o')

for i = 1:length(sols)
    sol = sols{i};
    save(['E:\Documents\GitHub\Research\Matlab\\iLQR\solutions\unscented\sol',num2str(i)], '-struct','sol');
end
%%
% wh = 1;
% ws = 2;
% wu = 0.2;
% inp = DDPInput([wh, ws, wu]);
% inp.ut_scale = 10;

% OL = entry_stochastic(0, [wh, ws, wu], [0, 1]);
% CL = entry_stochastic(1, [wh, ws, wu], [0, 1]);
% CLg = entry_stochastic_gains(inp);
% [hm, fm, sm, hv, fv, sv]=stats(CLg);
% CLg.stats = [hm, fm, sm, hv, fv, sv]; % terminal stats
% CLg = UTSolve(CLg, CLg.u(2:4,:));
% save(['Robust',num2str(wh),num2str(ws),'.mat'],'-struct','CLg')
% save(['Beta',num2str(inp.ut_scale),'.mat'],'-struct','CLg')
load ParametricUncertainty

% sols = {OL, CL, CLg};
for i = 1:4
    inp = DDPInput([i-1, i-1, 0.3]);
    %     sols{i} = entry_stochastic_gains_params(inp);
    [hm, fm, sm, hv, fv, sv]=stats(sols{i});
    sols{i}.stats = [hm, fm, sm, hv, fv, sv]; % terminal stats
    %     sols{i} = UTSolve(sols{i}, sols{i}.u(2:4,:));
end
% save('ParametricUncertainty.mat','sols');
% save('EqualWeightComparison.mat','sols')

% load EqualWeightComparison
% sols{3} = CLg;
% i=3;
% [hm, fm, sm, hv, fv, sv]=stats(sols{i});
%     sols{i}.stats = [hm, fm, sm, hv, fv, sv]; % terminal stats

load ParametricUncertainty
plot_sol(sols{1}, 'k', 0.6,0);
% plot_sol(sols{2}, 'c', 0.8,0);
% plot_sol(sols{3}, 'm', 0.4,0);
plot_sol(sols{4}, 'b', 0.6,0);

for i = 1:4
    %     figure(i)
    %     legend('Mean Opt', 'w=1', 'w=2', 'w=3')
    %     legend('Mean Opt', 'w=3')
    data(i,:) = [sols{i}.stats(1), sols{i}.stats(1)-3*sols{i}.stats(4)^0.5, sols{i}.stats(6)^0.5/1000];
    %     legend('Open Loop', 'Closed Loop', 'Closed Loop + Gain Opt')
end
figure(5)
plot(0:3, data)
legend('Mean Altitude', '3s Low Altitude', '3s Range Error')
close all
clc
kd = 0.1; %more lift up when too much drag
ks = -0.1;  % less lift up when too close
kf = -50.0;% left lift up when to shallow
K = [kd,ks,kf];

disp('Open loop optimization, flown open loop')
sols{1} = UTSolve(sols{1}); % Open loop flown open loop

disp('Open loop optimization, flown with fixed gains')
sol = UTSolve(sols{1}, K); % Open loop flown with fixed gains
% UTPlot(sol)
disp('Open loop optimization, flown with optimized gains')
sol = UTSolve(sols{1}, sols{3}.u(2:4,:)); % Open loop flown with optimized gains
% UTPlot(sol)
disp('Closed loop optimization with fixed gains')
sols{2} = UTSolve(sols{2}, K); % Closed loop flown with fixed gains
% UTPlot(sol)
disp('Closed loop optimized with fixed gains, flow with optimal gains')
sol = UTSolve(sols{2}, sols{3}.u(2:4,:)); % Closed loop flown with fixed gains

disp('Closed loop optimized jointly with gains, flown with fixed gains instead')
sol = UTSolve(sols{3}, K); % Optimal with fixed gains
% UTPlot(sol)
disp('Closed loop optimized jointly with gains')
sols{3} = UTSolve(sols{3}, sols{3}.u(2:4,:)); % Optimal with optimized gains, 1000 pts

% UTPlot(sol)
%
% UTPlot(sols{1})
% UTPlot(sols{1})

UTCompare(sols{1})
UTCompare(sols{2})
UTCompare(sols{3})
% UTCompare(sols{4})

%
% OL_alt = entry_stochastic(0, 0, 0); % no cost on variances terms
% CL_alt = entry_stochastic(1, 0, 0); % no cost on variances terms


%% Generate solutions for comparison
% three (open loop) altitude optimal solutions with different bank angle limits
% choose one compromised solution, may wh = 2, ws = 1, or just do mean opt
% again but closed loop to account for uncertainties

% Then, fly closed loop UT for all 4 profiles
if 1
    bounds = [[0,1];[0, 1]; [cosd(90-15), cosd(15)]; [cosd(90-30), cosd(30)]];
    W = zeros(4,3);
    
    W(1,:) = [1, 1, 0.2];
    sols = {};
    for i = 1:4
        inp = DDPInput(W(i,:));
        inp.terminal_plots = false;
        inp.running_plots=false;
        inp.bounds = bounds(i,:);
        if i == 1
            inp.ut_scale = 17;
            inp.horizon = 250;
            sols{i} = entry_stochastic_gains_params(inp);
            [hm, fm, sm, hv, fv, sv]=stats(sols{i});
            sols{i}.stats = [hm, fm, sm, hv, fv, sv]; % terminal stats
        else
            inp.horizon = 2000;
            sols{i} = entry_vel(inp); % deterministic optimization
            sols{i}.stats = zeros(1,6);
        end
        
    end
    sols{1}=UTSolve(sols{1}, sols{1}.input.gains);
    % margin_sols = sols;
    % load solutions_cl_ddp_all
    % margin_sols{1} = sols{6};
    % sols = margin_sols;
    save('margin_comparison.mat','sols') % just in case, allows for reloading quickly
    % load margin_comparison
    close all
    for i = 1:4
        if i >= 2
            sols{i}.Lm = sols{i}.L;
            sols{i}.Dm = sols{i}.D;
        end
        
        %     sols{i}=UTSolve(sols{i}, sols{1}.input.gains);
        
        %     stats(sols{i})
        %     sols{i} = UTSolve(sols{i}, sols{i}.input.gains);
        %     disp(sols{i}.ut.stats)
        %     UTPlot(sols{i})
    end
    UTCompare(sols{1}) % just a method to determine visually if we're discretizing accurately enough
    
    
    load solutions_cl_ddp_max
    sols = {sols{1}, sols{13}, sols{4}};
    
    close all
    plot_sol(sols{1}, 'k', 0.2, 0);
    plot_sol(sols{2}, 'c', 0.8, 0);
    plot_sol(sols{3}, 'b', 0.8, 0);
    plot_sol(sols{4}, 'r', 0.6, 0);
    
    j = 1;
    for i = 4:5
        figure(i)
        if 1
            legend('[0, 90]','[15, 75]', '[30, 60]', 'location','southwest')
            saveas(gcf, [sdir,'Nominal',figs{i},'.png'])
        else
            legend('w_h=0, w_s=0', 'w_h=3, w_s=0', 'w_h=0, w_s=3', 'location','best')
            saveas(gcf, [sdir,'Robust',figs{i},'.png'])
        end
        j = j+1;
        
    end
end
%% Sweep over the weights

Ws = 0:1:3;
Wh = 0:1:3;
% Ws = 0.5:1:2.5;
% Wh = 0.5:1:2.5;
% Ws = [0.5];
% Wh = [0,0.25,0.75:0.5:3, 3];
sols = {};
for j = 1:length(Ws)
    for i = 1:length(Wh)

        disp(length(sols))
        wh = Wh(i);
        ws = Ws(j);
        inp = DDPInput([wh, ws, 0.2]);
        inp.terminal_plots = false;
        inp.running_plots = true;
        inp.horizon = 250;
        inp.ut_scale = 15;
%         if ~isempty(sols)
%             inp.guess = sols{end}.u; % Use the robust mean solution as guess
%         end
        sols{end+1} = entry_stochastic_gains_params(inp);
        [hm, fm, sm, hv, fv, sv]=stats(sols{end});
        sols{end}.stats = [hm, fm, sm, hv, fv, sv]; % terminal stats
        if size(sols{end}.u,1) == 4
            sols{end} = UTSolve(sols{end}, sols{end}.u(2:4,:));
        else
            sols{end} = UTSolve(sols{end}, sols{end}.input.gains);
        end
        
    end
end


% save('solutions_cl_ddp.mat','sols') % just in case, allows for reloading quickly
save('solutions_cl_ddp_max_guess.mat','sols') % just in case, allows for reloading quickly
%%
for i = 1:length(sols)
    sol = sols{i};
    %     save(['E:\Documents\GitHub\Research\Matlab\iLQR\solutions\gainopt\','Sol',num2str(sol.weights(1)),num2str(sol.weights(2)),'.mat'],'-struct','sol')
    %     save(['E:\Documents\GitHub\Research\Matlab\iLQR\solutions\fixed_gain\','Sol',num2str(sol.weights(1)),num2str(sol.weights(2)),'.mat'],'-struct','sol')
    %     save(['E:\Documents\GitHub\Research\Matlab\iLQR\solutions\all_unc\','Sol',num2str(sol.weights(1)),num2str(sol.weights(2)),'.mat'],'-struct','sol')
    save(['E:\Documents\GitHub\Research\Matlab\iLQR\solutions\margin\','Sol',num2str(i),'.mat'],'-struct','sol')
    %     save(['E:\Documents\GitHub\Research\Matlab\iLQR\solutions\six_uncertainty\','Sol',num2str(i),'.mat'],'-struct','sol')
    
    disp(sol.mean(1,end))
end

disp(' ');
%%
% load solutions_cl_gains_ddp.mat
[n,m] = size(sols);
for i = 1:n*m
    W(i,:) = sols{i}.weights;
    S(i,:) = sols{i}.stats;
    U(i,:) = sols{i}.ut.stats;
end
keep = W(:,2) <= 3;

nstd = 3;
hlow = S(keep,1)-nstd*S(keep,4).^0.5;
hlow_ut = U(keep,2);

%%
% figure
% scatter(hlow, hlow_ut)
% figure
% scatter(3*S(:,6).^0.5, U(:,3))

N = 50 ;
xi = linspace(min(W(keep,1)), max(W(keep,1)), N) ;
yi = linspace(min(W(keep,2)), max(W(keep,2)), N) ;
[Xi, Yi] = meshgrid(xi,yi) ;
Hi = griddata(W(keep,1), W(keep,2), hlow, Xi, Yi) ;
Si = griddata(W(keep,1), W(keep,2), S(keep,6).^0.5, Xi, Yi) ;
DR = griddata(W(keep,1), W(keep,2), S(keep,3), Xi, Yi) ;

figure
contourf(Xi,Yi,Hi)
xlabel('w_h')
ylabel('w_s')
hold all
% title('3-sigma low altitude, km')
c = colorbar;
ylabel(c ,[num2str(nstd),'-sigma low altitude, km'],'fontsize',16)
% scatter(W(keep,1),W(keep,2),[],'r')
saveas(gcf, [sdir,'ClosedLoopAltitude.png'])

% figure
% surf(Xi,Yi,Hi)
% xlabel('w_h');
% ylabel('w_s');

figure
contourf(Xi, Yi, Si, 10);
hold all
xlabel('w_h');
ylabel('w_s');
% title('range error std, km')
c = colorbar;
ylabel(c ,'range error std, km','fontsize',16)
% scatter(W(keep,1),W(keep,2),[],'r')
saveas(gcf, [sdir,'ClosedLoopRangeError.png'])


figure
contourf(Xi, Yi, DR, 10);
hold all
xlabel('w_h');
ylabel('w_s');
% title('range error std, km')
c = colorbar;
ylabel(c ,'Mean Range, km','fontsize',16)
% scatter(W(keep,1),W(keep,2),[],'r')
saveas(gcf, [sdir,'ClosedLoopRangeMean.png'])
% It appears that more robust trajectories are longer

% figure
% surf(Xi, Yi, Si);
% colorbar


[hmax,imax] = max(hlow)
disp(W(imax,:))
[smin,imin] = min(S(keep,6).^0.5)
disp(W(imin,:))


% plot_sol(sols{1}, 'k', 0.2,0);
% plot_sol(sols{imax}, 'c', 0.8,0);
% plot_sol(sols{10}, 'b', 0.8,0);
%
% for i = 1:4
%     figure(i)
%     legend('Mean Alt Opt', 'w_h=2, w_s=0', 'w_h=1, w_s=1')
% end

%%
% Get the plots with all 4 for the control profile
% For init state unc or combined unc
plot_sol(sols{1}, 'k', 0.2,0);
plot_sol(sols{6}, 'c', 0.8,0);
plot_sol(sols{11}, 'b', 0.8,0);
plot_sol(sols{16}, 'r', 0.8,0);

% For param unc
% plot_sol(sols{1}, 'k', 0.2,0);
% plot_sol(sols{5}, 'c', 0.8,0);
% plot_sol(sols{9}, 'b', 0.8,0);
% plot_sol(sols{16}, 'r', 0.8,0);


for i = 1:4
    figure(i)
    legend('W=0', 'W=1', 'W=2', 'W=3')
end
%%
plot_sol(sols{1}, 'k', 0.8,0);
plot_sol(sols{16}, 'c', 0.6,0);
% fnames = {'Range_InitialStateUncertainty.png', 'Alt_InitialStateUncertainty.png', 'FPA_InitialStateUncertainty.png'};
fnames = {'Range_', 'Alt_', 'FPA_'};
suffix = 'AllUnc.png';
foldr = 'E:\Documents\EDL\Documents\Seminar\';

for i = 1:3
    figure(i)
    legend('W=0', 'W=3')
    saveas(gcf, [foldr,fnames{i},suffix])
end
%%
load solutions.mat

stats(OL_alt)
stats(OL)
stats(OL_dr)

stats(CL_alt)
stats(CL)
stats(CL_dr)


figure
plot([4.5, 3.78, 1.88],[5.459-3*0.58, 5.354-3*0.445, 5.021-3*0.576])
hold on
plot([21.7, 20.8, 17.7], [5.8762 - 3*0.079643, 5.7054 - 3*0.12262, 1.2051 - 3*0.36245])
xlabel('Range Error Standard Deviation, km')
ylabel('3-sigma low Altitude, km')
%% Rather than focus on how much closing the loop helps, examine the ability to reduce variance

close all
plot_sol(OL_alt, 'k', 0.2, 0);
plot_sol(OL    , 'c', 0.8, 0);
plot_sol(OL_dr , 'b', 0.8, 0);

j = 1;
for i = 1:4
    figure(i)
    legend('w_h=0.0, w_s=0.0', 'w_h=0.001, w_s=0.01','w_h=0.001, w_s=0.1')
    
    %     legend('Mean Alt Opt', 'w_h=0.001, w_s=0.01')
    savefig([sdir,'OpenLoop_',figs{j}])
    j = j+1;
    
end

plot_sol(CL_alt, 'k', 0.2, 4);
plot_sol(CL,     'c', 0.8, 4);
plot_sol(CL_dr,  'b', 0.8, 4);

j = 1;
for i = 5:8
    figure(i)
    legend('w_h=0.0, w_s=0.0', 'w_h=0.001, w_s=0.01','w_h=0.001, w_s=0.1')
    savefig([sdir,'ClosedLoop_',figs{j}])
    
    j = j+1;
    
end

function contour_plot()
figure
plot()

function [hm, fm, sm, hv, fv, sv]=stats(sol)
hm = sol.mean(1,end)/1000;
fm = sol.mean(2,end)*180/pi;
sm = sol.mean(3,end);
if sm(end) > 100e3
    k = 1000;
else
    k = 1;
end
sm = sm/k;

hv = sol.var(1,end)/1000^2;
fv = sol.var(2,end)*(180/pi)^2;
sv = sol.var(3,end)/k^2;

disp(['hf  = ',num2str(hm),'  +/- 3*',num2str(hv.^0.5),' km'])
disp(['fpa = ',num2str(fm),'  +/- 3*',num2str(fv.^0.5),' deg'])
disp(['sf  = ',num2str(sm),'  +/- 3*', num2str(sv.^0.5),' km'])
disp(' ')

function plot_sol(sol, color, alpha, figure_offset, nsmooth)
if nargin == 4
    nsmooth = 40;
end
v = sol.v;

hm = sol.mean(1,:);
fm = sol.mean(2,:)*180/pi;
sm = sol.mean(3,:);

if isfield(sol, 'var')
    hv = sol.var(1,:);
    fv = sol.var(2,:)*(180/pi)^2;
    sv = sol.var(3,:);
    if sm(end) > 2000
        sm = sm/1000;
        sv = sv/1000^2;
    end
    
    figure(1+figure_offset)
    hold all
    x2 = [v', fliplr(v')];
    inBetween = [(sm-3*sv.^0.5), fliplr((sm+3*sv.^0.5))];
    fill(x2, inBetween, color, 'Facealpha',alpha);
    
    xlabel('Velocity, m/s')
    ylabel('Downrange, km')
    grid on
    
    figure(2+figure_offset)
    hold all
    inBetween = [(hm-3*hv.^0.5)/1000, fliplr((hm+3*hv.^0.5)/1000)];
    fill(x2, inBetween, color, 'Facealpha',alpha);
    
    
    xlabel('Velocity, m/s')
    ylabel('Altitude, km')
    grid on
    
    figure(3+figure_offset)
    hold all
    inBetween = [(fm-3*fv.^0.5), fliplr((fm+3*fv.^0.5))];
    fill(x2, inBetween, color, 'Facealpha',alpha);
    
    xlabel('Velocity, m/s')
    ylabel('FPA, deg')
    grid on
end
figure(4+figure_offset)
hold all
if nsmooth > 0
    plot(v(1:end-1), smooth(sol.u(1,:), nsmooth), 'linewidth', 3)
else
    plot(v(1:end-1), sol.u(1,:), 'linewidth', 3)

end
ylim([0,1])
xlabel('Velocity, m/s')
ylabel('cos(bank), -')
grid on

figure(5+figure_offset)
hold all
plot(v, sol.Dm, 'linewidth', 3)
ylim([0, max(120, max(sol.Dm))])
xlabel('Velocity, m/s')
ylabel('Drag, m/s/s')
grid on

function sol = UTSolve(sol, K)
%%% Integrates equations of motion for each sigma point and computes stats
if nargin == 1
    K = [0,0,0];
    k = @(v) K;
else
    if length(K(:)) == 3 % fixed gain scenario
        k = @(v) K;
    else
        %        for i=1:3
        %            K(i,:) = smooth(K(i,:),5);
        %        end
        k = @(v) interp1(sol.v, [K, K(:,end)]', v);
    end
end
sol.u(1,end) = sol.u(1,end-1);
xr = @(v) interp1(sol.v, [sol.Dm; sol.mean(2:3,:)]', v);
% u = smooth([sol.u(1,:), sol.u(1,end)], 5);
u = [sol.u(1,:), sol.u(1,end)];

[V,X,U] = UnscentedEntry(sol.v, sol.X0, @(v)interp1(sol.v', u, v), k, sol.sigma_weights, xr);
sol.ut.v = V;
sol.ut.x = X;
sol.ut.u = U;

nx = length(sol.sigma_weights);
hmean = X(:,1:nx)*sol.sigma_weights;
fmean = X(:,nx+1:2*nx)*sol.sigma_weights;
smean = X(:,2*nx+1:3*nx)*sol.sigma_weights;

hvar = (X(:,1:nx)-hmean).^2*sol.sigma_weights;
fvar = (X(:,nx+1:2*nx)-fmean).^2*sol.sigma_weights;
svar = (X(:,2*nx+1:3*nx)-smean).^2*sol.sigma_weights;

sol.ut.mean = [hmean, fmean*180/pi, smean];
sol.ut.var = [hvar, fvar*(180/pi)^2, svar];
sol.ut.std = sol.ut.var.^0.5;
sol.ut.stats = [hmean(end)/1000, (hmean(end)-3*hvar(end)^0.5)/1000, 3*svar(end).^0.5]; % 3-sigma low altitude, range std
print_header()
disp(sol.ut.stats)

function print_header()
disp('Mean Alt, km  3sig low Alt, km  3sig DR (km)')

function UTPlot(sol)
nx = length(sol.sigma_weights);
v = sol.ut.v;

hm = sol.ut.mean(:,1)';
hv = sol.ut.var(:,1)';
sm = sol.ut.mean(:,3)';
sv = sol.ut.var(:,3)';

alpha = 0.75;

figure(1)
x2 = [v', fliplr(v')];
inBetween = [(hm-3*hv.^0.5)/1000,fliplr((hm+3*hv.^0.5)/1000)];
fill(x2, inBetween, 'c', 'Facealpha',alpha);
hold all
plot(sol.ut.v, sol.ut.x(:,1:nx)/1000, 'k--')
xlabel('Velocity, m/s')
ylabel('Altitude, km')
legend('3-sigma bounds','UT Sigma Point Trajectories', 'location','best')
grid on

% figure(2)
% plot(sol.ut.v, 180/pi*sol.ut.x(:,nx+1:nx+nx), 'k--')
% xlabel('Velocity, m/s')
% ylabel('FPA, deg')

figure(3)
inBetween = [(sm-3*sv.^0.5),fliplr((sm+3*sv.^0.5))];
fill(x2, inBetween, 'c', 'Facealpha',alpha);
hold all
plot(sol.ut.v, sol.ut.x(:,2*nx+1:3*nx), 'k--')
xlabel('Velocity')
ylabel('Downrange Distance, km')
legend('3-sigma bounds','UT Sigma Point Trajectories')
grid on


function UTCompare(sol)
%%% Compares unscented integration with output from DDP, assumes UTSolve
%%% has been called with appropriate gains for comparison
nx = size(sol.h,1);
% sol = UTSolve(sol);
disp('DDP')
disp([sol.stats(1), sol.stats(1)-3*sol.stats(4)^0.5, 3*sol.stats(end)^0.5])
disp("integration")
disp(sol.ut.stats)

figure
plot(sol.v, sol.h/1000)
hold all
plot(sol.ut.v, sol.ut.x(:,1:nx)/1000, '--')
xlabel('Velocity')
ylabel('Altitude, km')
title('Dashed Lines = Integration, Solid = DDP')

figure
plot(sol.v, sol.fpa)
hold all
plot(sol.ut.v, sol.ut.x(:,nx+1:nx+nx), '--')
title('Dashed Lines = Integration, Solid = DDP')
xlabel('Velocity')

figure
plot(sol.v, sol.s)
hold all
plot(sol.ut.v, sol.ut.x(:,2*nx+1:3*nx), '--')
xlabel('Velocity')
ylabel('DR')
title('Dashed Lines = Integration, Solid = DDP')


figure
plot(sol.v(1:end-1), sol.u(1,:), 'k','linewidth',2)
hold all
plot(sol.ut.v, sol.ut.u, '--')
xlabel('Velocity')
ylabel('Control')
title('Dashed Lines = Integration, Solid = DDP')