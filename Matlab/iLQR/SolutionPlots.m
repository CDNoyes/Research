function SolutionPlots(sol, save_prefix)
if nargin > 1
    do_save = true;
else
    do_save = false;
end

savedir = 'E:\Documents\EDL\Documents\Dissertation\Images\Trajectory\';

lw = 2;
fs = 14;
ticksize = 10;
figsize = [500,300];
position = [100,100];

h = sol.h;
v = sol.v;
u = sol.u;
hm = sol.mean(1,:);
hv = sol.var(1,:);
D = sol.D;
Dm = sol.Dm;
Dv = sol.Dv;
fpa = sol.fpa;
fm = sol.mean(2,:);
fv = sol.var(2,:);
s = sol.s;
sm = sol.mean(3,:);
sv = sol.var(3,:);


f=figure;
set(f, 'Position',[position,figsize]);

x2 = [v', fliplr(v')];
inBetween = [(sm-3*sv.^0.5), fliplr((sm+3*sv.^0.5))];
fill(x2, inBetween, 'c');
hold all
plot(v, sm, 'm', 'linewidth', lw)
plot(v, s', 'k--','linewidth', 1)
set ( gca, 'xdir', 'reverse')
ax = gca;
ax.XAxis.FontSize = ticksize;
ax.YAxis.FontSize = ticksize;
xlabel('Velocity (m/s)', 'fontsize', fs)
ylabel('Downrange (km)', 'fontsize', fs)
legend('3\sigma region','Mean Trajectory', 'Sigma Trajectory','location','best')
grid on
if do_save
    saveas(gcf, [savedir,save_prefix,'Range.png'])
end

f=figure;
set(f, 'Position',[position,figsize]);
x2 = [v', fliplr(v')];
inBetween = [(hm-3*hv.^0.5)/1000, fliplr((hm+3*hv.^0.5)/1000)];
fill(x2, inBetween, 'c');
hold all
plot(v, hm/1000, 'm', 'linewidth', lw)
plot(v, h'/1000, 'k--','linewidth', 1)
set ( gca, 'xdir', 'reverse' )
ax = gca;
ax.XAxis.FontSize = ticksize;
ax.YAxis.FontSize = ticksize;
legend('3\sigma region','Mean Trajectory', 'Sigma Trajectory','location','best')
xlabel('Velocity (m/s)', 'fontsize', fs)
ylabel('Altitude (km)', 'fontsize', fs)
grid on
if do_save
    saveas(gcf, [savedir,save_prefix,'Altitude.png'])
end

f=figure;
set(f, 'Position',[position,figsize]);
inBetween = [(Dm-3*Dv.^0.5), fliplr((Dm+3*Dv.^0.5))];
fill(x2, inBetween, 'c');
hold all
plot(v, Dm, 'm', 'linewidth', lw)
plot(v, D', 'k--','linewidth', 1)
set ( gca, 'xdir', 'reverse' )
ax = gca;
ax.XAxis.FontSize = ticksize;
ax.YAxis.FontSize = ticksize;
xlabel('Velocity (m/s)', 'fontsize', fs)
ylabel('Drag (m/s^2)', 'fontsize', fs)
grid on
legend('3\sigma region','Mean Trajectory', 'Sigma Trajectory','location','best')
if do_save
    saveas(gcf, [savedir,save_prefix,'Drag.png'])
end

f=figure;
set(f, 'Position',[position,figsize]);
inBetween = [(fm-3*fv.^0.5), fliplr((fm+3*fv.^0.5))];
fill(x2, inBetween*180/pi, 'c');
hold all
plot(v, fm*180/pi, 'm', 'linewidth', lw)
plot(v, 180/pi*fpa', 'k--','linewidth', 1)

set ( gca, 'xdir', 'reverse' )
ax = gca;
ax.XAxis.FontSize = ticksize;
ax.YAxis.FontSize = ticksize;
xlabel('Velocity (m/s)', 'fontsize', fs)
ylabel('FPA (deg)', 'fontsize', fs)
grid on
legend('3\sigma region','Mean Trajectory', 'Sigma Trajectory','location','best')
if do_save
    saveas(gcf, [savedir,save_prefix,'FPA.png'])
end

f=figure;
set(f, 'Position', [position,figsize]);
plot(v(2:end-1), u(1,1:end-1), 'k', 'linewidth', lw)
hold all
set ( gca, 'xdir', 'reverse' )
ax = gca;
ax.XAxis.FontSize = ticksize;
ax.YAxis.FontSize = ticksize;
xlabel('Velocity (m/s)', 'fontsize', fs)
ylabel('Control', 'fontsize', fs)
grid on
if do_save
    saveas(gcf, [savedir,save_prefix,'Control.png'])
end

%% Four Grid
gridsize = figsize*2.5;
f=figure;
set(f, 'Position',[position,gridsize]);

subplot(2,2,1)
x2 = [v', fliplr(v')];
inBetween = [(hm-3*hv.^0.5)/1000, fliplr((hm+3*hv.^0.5)/1000)];
fill(x2, inBetween, 'c');
hold all
plot(v, hm/1000, 'm', 'linewidth', lw)
plot(v, h'/1000, 'k--','linewidth', 1)
% set ( gca, 'xdir', 'reverse' )
ax = gca;
ax.XAxis.FontSize = ticksize;
ax.YAxis.FontSize = ticksize;
legend('3\sigma region','Mean Trajectory', 'Sigma Trajectory','location','best')
xlabel('Velocity (m/s)', 'fontsize', fs)
ylabel('Altitude (km)', 'fontsize', fs)
grid on

subplot(2,2,2)
inBetween = [(fm-3*fv.^0.5), fliplr((fm+3*fv.^0.5))];
fill(x2, inBetween*180/pi, 'c');
hold all
plot(v, fm*180/pi, 'm', 'linewidth', lw)
plot(v, 180/pi*fpa', 'k--','linewidth', 1)
% set ( gca, 'xdir', 'reverse' )
ax = gca;
ax.XAxis.FontSize = ticksize;
ax.YAxis.FontSize = ticksize;
xlabel('Velocity (m/s)', 'fontsize', fs)
ylabel('FPA (deg)', 'fontsize', fs)
grid on
legend('3\sigma region','Mean Trajectory', 'Sigma Trajectory','location','best')

subplot(2,2,3)
inBetween = [(sm-3*sv.^0.5), fliplr((sm+3*sv.^0.5))];
fill(x2, inBetween, 'c');
hold all
plot(v, sm, 'm', 'linewidth', lw)
plot(v, s', 'k--','linewidth', 1)
% set ( gca, 'xdir', 'reverse')
ax = gca;
ax.XAxis.FontSize = ticksize;
ax.YAxis.FontSize = ticksize;
xlabel('Velocity (m/s)', 'fontsize', fs)
ylabel('Downrange (km)', 'fontsize', fs)
legend('3\sigma region','Mean Trajectory', 'Sigma Trajectory','location','best')
grid on

subplot(2,2,4)
plot(v(2:end-1), u(1,1:end-1), 'k', 'linewidth', lw)
hold all
% set ( gca, 'xdir', 'reverse' )
ax = gca;
ax.XAxis.FontSize = ticksize;
ax.YAxis.FontSize = ticksize;
xlabel('Velocity (m/s)', 'fontsize', fs)
ylabel('Control', 'fontsize', fs)
grid on
if do_save
    saveas(gcf, [savedir,save_prefix,'Grid.png'])
end
