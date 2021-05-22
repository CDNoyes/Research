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