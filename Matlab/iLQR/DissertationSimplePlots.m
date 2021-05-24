%% Convergence Example
close all; clear; clc;

load ConvergenceExampleSmoother
for i = 1:3
    [hm, fm, sm, hs, fs, ss] = get_stats(sols{i});
    stats(i,:) = [hm, fm, sm, hs, fs, ss];
end
fields = fieldnames(sols{1}.trace);
labels = {'Iter','Regularization, \lambda', 'dlamba','Total Cost','Linesearch Stepsize, \epsilon','max(||Q_u||)','Improvement (|\DeltaJ|)','Reduction Ratio', 'Time to Compute Derivatives (s)','Time Per Forward Pass','Time Per Backward Pass'};
for j = 2:length(fields)
    for i = 1:3
        cost(i) = sols{i}.cost;
        trace = sols{i}.trace;
        iter = 1:trace(end).iter;
        y = [];
        for k = iter
            y(k) =  getfield(sols{i}.trace(k), fields{j});
        end
        
        if j == 5
            y(y>1) = 1; %linesearch
        end
        if j == 9
            mean(y)
        end
        
        f=figure(j-1);
        
        if i > 1
            hold on
        else
            set(f, 'Position',[100,100,500,300])
        end
        if any(j == [2,3,5,6,7])
            semilogy(iter, y, 'linewidth', 2)
        else
            plot(iter, y, 'linewidth', 2)
        end
        xlabel('Iteration','fontsize',16)
        ylabel(labels{j},'fontsize',16)
    end
    grid on
    legend('DDP','DDP-ControlHessian','iLQR','location','northeast')
    %     saveas(gcf, ['C:\Users\Aero\Documents\EDL\Documents\Dissertation\Images\Convergence\',fields{j},'.png']) %laptop
    saveas(gcf, ['E:\Documents\EDL\Documents\Dissertation\Images\Convergence\',fields{j},'.png']) %pc
    
end
%% Control Profiles
f=figure;
set(f, 'Position',[100,100,500,300])
for i=1:3
    hold all
    plot(sols{i}.v(2:end), smooth(sols{i}.u,20), 'linewidth', 2)
end
grid on
legend('DDP','DDP-ControlHessian','iLQR','location','northwest')
xlabel('Velocity (m/s)','fontsize',14)
ylabel('Reference Control','fontsize',14)
set(gca, 'XDir','reverse')

%     saveas(gcf, ['C:\Users\Aero\Documents\EDL\Documents\Dissertation\Images\Convergence\',fields{j},'.png']) %laptop
saveas(gcf, 'E:\Documents\EDL\Documents\Dissertation\Images\Convergence\ControlProfiles.png') %pc

%%
fs = 14;

V = linspace(5500, 460, 500);
u = linspace(0.35, 1, 500);

figure
% subplot(1,2,1)
plot(V, u, 'linewidth',2)
xlabel('Velocity (m/s)', 'fontsize',fs)
ylabel('Reference Control', 'fontsize',fs)
grid on
set(gca, 'XDir','reverse')
saveas(gcf, 'E:\Documents\EDL\Documents\Dissertation\Images\GuessControl.png')

%% Aero/Atmo
fs = 14;

M = linspace(1,24,500);
[CD,CL] = AerodynamicCoefficients(M);

figure
% subplot(1,2,1)
plot(M, [CD,CL, CL./CD], 'linewidth',2)
legend('C_D', 'C_L', 'L/D', 'location','best')
xlabel('Mach', 'fontsize',fs)
grid on
set(gca, 'XDir','reverse')
saveas(gcf, 'C:\Users\Aero\Documents\EDL\Documents\Dissertation\Images\CoeffNominal.png')

h = linspace(0, 55e3, 10000);
rho = MarsAtmosphericDensity(h);

% subplot(1,2,2)
figure
plot(h/1000, rho, 'linewidth', 2)
ylabel('Atmospheric Density (kg/m^3)', 'fontsize',fs)
xlabel('Altitude (km)', 'fontsize', fs)
grid on
saveas(gcf, 'C:\Users\Aero\Documents\EDL\Documents\Dissertation\Images\DensityNominal.png')

%% Saturation Approx
fs = 14;
x = linspace(-0.5, 1.5, 5000);
y = Saturate(x, 0, 1);
figure(1)
% newcolors = {'#F00','#F80','#0BB'};
% colororder(newcolors)
hold all
grid on
for K = [1, 3, 20]
    ys = smooth_sat(x, K);
    plot(x, ys,'linewidth', 2)
    
end
plot(x, y, 'k--', 'linewidth', 2)
legend('M=1','M=3','M=20','Sat', 'location','northwest')
xlabel('x', 'fontsize', fs)
ylabel('Sat_{[0,1]}(x)', 'fontsize', fs)
saveas(gcf, 'C:\Users\Aero\Documents\EDL\Documents\Dissertation\Images\SmoothSat.png')

%% The derivative of smooth sat
% K = linspace(0.1, 20, 1000);

fs = 14;
x = linspace(-0.5, 1.5, 5000);
y = Saturate(x, 0, 1);
figure(1)
% newcolors = {'#F00','#F80','#0BB'};
% colororder(newcolors)
hold all
grid on
for K = [1, 3, 20]
    ys = smooth_sat(x, K);
    plot(x, ys,'linewidth', 2)
    
end
plot(x, y, 'k--', 'linewidth', 2)
legend('M=1','M=3','M=20','Sat', 'location','northwest')
xlabel('x', 'fontsize', fs)
ylabel('Sat_{[0,1]}(x)', 'fontsize', fs)
% saveas(gcf, 'C:\Users\Aero\Documents\EDL\Documents\Dissertation\Images\SmoothSatDeriv.png')