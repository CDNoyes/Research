function generate_solutions
sdir = 'E:\Documents\EDL\Documents\PropellantOptimalJournal\ddp\matlab\';

figs = {'Range','Altitude','Fpa','Control'};

wh = 3;
ws = 1;
wu = 0.2;

% OL = entry_stochastic(0, [wh, ws, wu], [0, 1]);
% CL = entry_stochastic(1, [wh, ws, wu], [0, 1]);
CLg = entry_stochastic_gains([wh, ws, wu], [0, 1]);
[hm, fm, sm, hv, fv, sv]=stats(CLg);
CLg.stats = [hm, fm, sm, hv, fv, sv]; % terminal stats 
CLg = UTSolve(CLg, CLg.u(2:4,:));
save(['Robust',num2str(wh),num2str(ws),'.mat'],'-struct','CLg')

% sols = {OL, CL, CLg};
% for i = 1:3
%     [hm, fm, sm, hv, fv, sv]=stats(sols{i});
%     sols{i}.stats = [hm, fm, sm, hv, fv, sv]; % terminal stats 
% end
% save('EqualWeightComparison.mat','sols')

load EqualWeightComparison 
% sols{3} = CLg;
% i=3;
% [hm, fm, sm, hv, fv, sv]=stats(sols{i});
%     sols{i}.stats = [hm, fm, sm, hv, fv, sv]; % terminal stats 

% plot_sol(sols{1}, 'k', 0.2,0);
% plot_sol(sols{2}, 'c', 0.8,0);
% plot_sol(sols{3}, 'b', 0.8,0);
% 
% for i = 1:4
%     figure(i)
%     legend('Open Loop', 'Closed Loop', 'Closed Loop + Gain Opt')
% end
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
if 0
% bounds = [[0,1];[0, 1]; [cosd(90-15), cosd(15)]; [cosd(90-30), cosd(30)]];
% W = zeros(4,3);
% 
% W(1,:) = [1, 0.25, 0.5];
% % W(1,3) = 0.5;
% cl = [1, 0, 0, 0];
% i=1;
% for i = 1:4
%     
%         sols{i} = entry_stochastic(cl(i), W(i,:), bounds(i,:));
%         [hm, fm, sm, hv, fv, sv]=stats(sols{i});
%         sols{i}.stats = [hm, fm, sm, hv, fv, sv]; % terminal stats 
% end

% save('margin_comparison.mat','sols') % just in case, allows for reloading quickly 
load margin_comparison
close all
for i = 1:4
    sols{i} = UTSolve(sols{i});
    disp(sols{i}.ut.stats)
%     UTPlot(sols{i})
end
UTCompare(sols{1}) % just a method to determine visually if we're discretizing accurately enough

close all
plot_sol(sols{1}, 'k', 0.2, 0);
plot_sol(sols{2}, 'c', 0.8, 0);
plot_sol(sols{3}, 'b', 0.8, 0);
plot_sol(sols{4}, 'r', 0.6, 0);

j = 1;
for i = 1:5
    figure(i)
    legend('Stochastic CL, w_h=1, w_s=0.25', '[cos(90), cos(0)]','[cos(75), cos(15)]', '[cos(60), cos(30)]', 'location','best')
%         saveas(gcf, [sdir,'Comparison',figs{j},'.png'])
    j = j+1;

end
end
%% Sweep over the weights

Ws = 0:1:3;
Wh = 0:1:3;
sols = {};
for i = 1:length(Wh)
    for j = 1:length(Ws)
        wh = Wh(i);
        ws = Ws(j);
%         sols{end+1} = entry_stochastic(1, [wh, ws, 0.2], [0,1]);
        sols{end+1} = entry_stochastic_gains([wh, ws, 0.2], [0,1]);

        [hm, fm, sm, hv, fv, sv]=stats(sols{end});
        sols{end}.stats = [hm, fm, sm, hv, fv, sv]; % terminal stats 
        close all
    end
end


% save('solutions_cl_ddp.mat','sols') % just in case, allows for reloading quickly 
save('solutions_cl_gains_ddp.mat','sols') % just in case, allows for reloading quickly 

disp(' ');
%%
% load solutions_cl.mat
[n,m] = size(sols);
for i = 1:n*m
    W(i,:) = sols{i}.weights;
    S(i,:) = sols{i}.stats;
    
end

nstd = 3;
hlow = S(:,1)-nstd*S(:,4).^0.5;

N = 50 ;
xi = linspace(min(W(:,1)), max(W(:,1)), N) ;
yi = linspace(min(W(:,2)), max(W(:,2)), N) ;
[Xi, Yi] = meshgrid(xi,yi) ;
Hi = griddata(W(:,1), W(:,2), hlow, Xi, Yi) ;
Si = griddata(W(:,1),W (:,2), S(:,6).^0.5, Xi, Yi) ;

figure
contourf(Xi,Yi,Hi)
xlabel('W_h');
ylabel('W_s');
% title('3-sigma low altitude, km')
c = colorbar;
ylabel(c ,[num2str(nstd),'-sigma low altitude, km'],'fontsize',16)
% saveas(gcf, [sdir,'ClosedLoopAltitude.png'])


levels = [linspace(0,2, 20),linspace(2,5,10)];
figure
contourf(Xi, Yi, Si,levels);
hold all
xlabel('W_h');
ylabel('W_s');
% title('range error std, km')
c = colorbar;
ylabel(c ,'range error std, km','fontsize',16)
scatter(W(:,1),W(:,2),[],'r')
% saveas(gcf, [sdir,'ClosedLoopRangeError.png'])


[hmax,imax] = max(hlow)
disp(W(imax,:))

plot_sol(sols{1}, 'k', 0.2,0);
plot_sol(sols{imax}, 'c', 0.8,0);
plot_sol(sols{10}, 'b', 0.8,0);

for i = 1:4
    figure(i)
    legend('Mean Alt Opt', 'w_h=2, w_s=0', 'w_h=1, w_s=1')
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


%%
% close all
plot_sol(OL, 'k', 0.2);
plot_sol(CL, 'c', 0.8);

for i = 11:14
    figure(i)
    legend('Open Loop', 'Closed Loop')
end

function stats_table(sols)

    for i = 1:length(sols)
        [hm, fm, sm, hv, fv, sv] = stats(sols{i});

    end

function [hm, fm, sm, hv, fv, sv]=stats(sol)
hm = sol.mean(1,end)/1000;
fm = sol.mean(2,end)*180/pi;
sm = sol.mean(3,end);

hv = sol.var(1,end)/1000^2;
fv = sol.var(2,end)*(180/pi)^2;
sv = sol.var(3,end);

disp(['hf  = ',num2str(hm),'  +/- 3*',num2str(hv.^0.5),' km'])
disp(['fpa = ',num2str(fm),'  +/- 3*',num2str(fv),' deg'])
disp(['sf  = ',num2str(sm),'  +/- 3*', num2str(sv.^0.5),' km'])
disp(' ')

function plot_sol(sol, color, alpha, figure_offset)

v = sol.v;

hm = sol.mean(1,:);
fm = sol.mean(2,:)*180/pi;
sm = sol.mean(3,:);

hv = sol.var(1,:);
fv = sol.var(2,:)*(180/pi)^2;
sv = sol.var(3,:);

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

figure(4+figure_offset)
hold all
plot(v(1:end-1), sol.u(1,:), 'linewidth', 3)
xlabel('Velocity, m/s')
ylabel('cos(bank), -')
grid on

figure(5+figure_offset)
hold all
plot(v, sol.Dm, 'linewidth', 3)
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
       k = @(v) interp1(sol.v, [K, K(:,end)]', v);
   end
end
sol.u(1,end) = sol.u(1,end-1);
xr = @(v) interp1(sol.v, [sol.Dm; sol.mean(2:3,:)]', v);
[V,X,U] = UnscentedEntry(sol.v, sol.X0, @(v)interp1(sol.v', [sol.u(1,:), sol.u(1,end)], v), k, sol.sigma_weights, xr);
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