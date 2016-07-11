%ENTRYPLOTS Plot some standard graphs useful in examining entry
%trajectories.

function EntryPlots(trajectorySummary)
% close all
dtr = pi/180;
planet = Mars();
ts = trajectorySummary;
t = ts.time;
x = ts.state;
sig = ts.sigma;
u = ts.control;
r_eq = planet.radiusEquatorial;
hkm = (x(:,1)-r_eq)/1000;
e = ts.energy_norm;

% Find max Qdyn and L/D indices
% [Lmax,i_Lmax] = max(ts.L);
% [Lgmax,i_Lgmax] = max(ts.L.*cos(x(:,5))');
% [Dmax,i_Dmax] = max(ts.D);
[LoDmax,i_LoDmax] = max(ts.L./ts.D);
[Qmax,i_Qmax] = max(ts.D./ts.C_D);

[lineSpecs,textSpecs,figSpecs] = PlotSpecs();

% Altitude vs Velocity
n = 1;
figure(n)
grid on
box on
set(gcf,'name','Altitude vs Velocity', figSpecs{:})
ParachuteDeploymentConstraints(true);
hold on
plot(x(:,4),hkm, lineSpecs{:})
xlabel('Velocity (m/s)',textSpecs{:})
ylabel('Altitude (km)',textSpecs{:})
title(['Final altitude: ',num2str(hkm(end)), ' km'],textSpecs{:})

% Range vs Velocity
n=n+1;
figure(n)
grid on
box on
set(gcf,'name','Range To Go vs Velocity', figSpecs{:})
hold on
plot(x(:,4),ts.target.DR-ts.DR, lineSpecs{:})
xlabel('Velocity (m/s)',textSpecs{:})
ylabel('Range To Go (km)',textSpecs{:})
title(['Final altitude: ',num2str(hkm(end)), ' km'],textSpecs{:})

% Altitude vs Range:
n = n+1;
figure(n)
grid on
box on
set(gcf,'name','Altitude vs Range', figSpecs{:})
% ParachuteDeploymentConstraints(true);
hold on
plot(ts.DR,hkm, lineSpecs{:})
plot(ts.target.DR*ones(1,100),linspace(0,max(hkm)),'r--',lineSpecs{:})
xlabel('DownRange (km)',textSpecs{:})
ylabel('Altitude (km)',textSpecs{:})
title(['Final altitude: ',num2str(hkm(end)), ' km'],textSpecs{:})

%DR v CR
n = n+1;
figure(n)
plot(ts.CR,ts.DR, lineSpecs{:})
hold on
plot(ts.CR(i_Qmax),ts.DR(i_Qmax),'ko',ts.CR(i_LoDmax),ts.DR(i_LoDmax),'r*' )
ylabel('Downrange (km)',textSpecs{:})
xlabel('Crossrange (km)',textSpecs{:})
title(['Downrange: ',num2str(ts.DR(end)),' km, Crossrange: ',num2str(ts.CR(end)), ' km'],textSpecs{:})
axis([min(min(ts.CR),-5),max(max(ts.CR),5),0,100*ceil(max(ts.DR)/100)]) % will be ugly for DR > 1000
grid on
box on
legend('Trajectory','Peak Dynamic Pressure','Max L/D','Location','best')
set(gcf,'name','Downrange vs Crossrange', figSpecs{:})

%Lat/Lon
n = n+1;
figure(n)
plot(x(1,3)/dtr,x(1,2)/dtr,'ko' )
hold on
plot(ts.target.lat/dtr,ts.target.lon/dtr,'kx' )
plot(x(:,3)/dtr,x(:,2)/dtr, lineSpecs{:})
xlabel('Latitude (deg)',textSpecs{:})
ylabel('Longitude (deg)',textSpecs{:})
title('Ground Track', textSpecs{:})
legend('Initial point','Target','Location','Best')
grid on
box on
set(gcf,'name','Groundtrack', figSpecs{:})

%Flight path angle
n = n+1;
figure(n)
subplot 211
plot(t,x(:,5)/dtr, lineSpecs{:})
grid on
axis([0,max(t),min(x(:,5)/dtr)*(1-sign(min(x(:,5)))*.1),max(x(:,5)/dtr)*(1+sign(max(x(:,5)))*.1)])
xlabel('Time (s)',textSpecs{:})
ylabel('FPA (deg)',textSpecs{:})
subplot 212
plot(ts.energy_norm,x(:,5)/dtr, lineSpecs{:})
axis([0,1,min(x(:,5)/dtr)*(1-sign(min(x(:,5)))*.1),max(x(:,5)/dtr)*(1+sign(max(x(:,5)))*.1)])
xlabel('Normalized Energy',textSpecs{:})
ylabel('FPA (deg)',textSpecs{:})
grid on
box on
set(gcf,'name','FPA', figSpecs{:})

%Heading angle
n = n+1;
figure(n)
subplot 211
plot(t,x(:,6)/dtr, lineSpecs{:})
grid on
axis([0,max(t),min(x(:,6)/dtr)*(1-sign(min(x(:,6)))*.1),max(x(:,6)/dtr)*(1+sign(max(x(:,6)))*.1)])
xlabel('Time (s)',textSpecs{:})
ylabel('Heading (deg)',textSpecs{:})
subplot 212
plot(ts.energy_norm,x(:,6)/dtr, lineSpecs{:})
axis([0,1,min(x(:,6)/dtr)*(1-sign(min(x(:,6)))*.1),max(x(:,6)/dtr)*(1+sign(max(x(:,6)))*.1)])
xlabel('Normalized Energy',textSpecs{:})
ylabel('Heading (deg)',textSpecs{:})
grid on
box on
set(gcf,'name','Heading', figSpecs{:})

%CONTROL
n = n+1;
figure(n)
subplot 211
plot(t,sig/dtr, lineSpecs{:})
grid on
axis([0,max(t),-180,180])
xlabel('Time (s)',textSpecs{:})
ylabel('Bank Angle (deg)',textSpecs{:})
subplot 212
plot(ts.energy_norm,sig/dtr, lineSpecs{:})
axis([0,1,-180,180])
xlabel('Normalized Energy',textSpecs{:})
ylabel('Bank Angle (deg)',textSpecs{:})
grid on
box on
set(gcf,'name','Bank Profile', figSpecs{:})

% MACH
n = n+1;
figure(n)
subplot 211
plot(t,ts.M, lineSpecs{:})
grid on
axis([0,max(t),0,35])
xlabel('Time (s)',textSpecs{:})
ylabel('Mach Number (-)',textSpecs{:})
subplot 212
plot(ts.energy_norm,ts.M, lineSpecs{:})
axis([0,1,0,35])
xlabel('Normalized Energy',textSpecs{:})
ylabel('Mach Number (-)',textSpecs{:})
grid on
box on
set(gcf,'name','Mach Profile', figSpecs{:})

n = n+1;
figure(n)
subplot 211
plot(t,x(:,4), lineSpecs{:})
grid on
axis([0,max(t),0,6000])
xlabel('Time (s)',textSpecs{:})
ylabel('Velocity (m/s)',textSpecs{:})
subplot 212
plot(e,x(:,4), lineSpecs{:})
xlabel('Normalized Energy',textSpecs{:})
ylabel('Velocity (m/s)',textSpecs{:})
grid on
box on
set(gcf,'name','Velocity Profile', figSpecs{:})

gload = sqrt(ts.L.^2+ts.D.^2)/9.81;
n = n+1;
figure(n)
subplot 211
plot(t,gload, lineSpecs{:})
grid on
% axis([0,max(t),0,6000])
xlabel('Time (s)',textSpecs{:})
ylabel('Load (g)',textSpecs{:})
subplot 212
plot(e,gload, lineSpecs{:})
xlabel('Normalized Energy',textSpecs{:})
ylabel('Load (g)',textSpecs{:})
grid on
box on
set(gcf,'name','g-load profile', figSpecs{:})


if isfield(ts,'observer') && ts.observer
    n = n+1;
    figure(n)
    plot(e,ts.state(:,7),lineSpecs{:})
    hold all
    plot(e,ts.D,lineSpecs{:})
    plot(e,ts.D'./cos(x(:,5)), lineSpecs{:})
    
    legend('Observed Drag','Modeled Drag', 'Modeled Drag/cos(\gamma)')
    xlabel('Normalized Energy (-)',textSpecs{:})
    ylabel('Drag (m/s^2)',textSpecs{:})
    set(gcf,'name','Drag Profile', figSpecs{:})
    grid on
    
    n = n+1;
    figure(n)
    plot(ts.energy_norm,ts.state(:,9),lineSpecs{:})
    xlabel('Normalized Energy (-)',textSpecs{:})
    ylabel('Disturbance Estimate (m/s^4)',textSpecs{:})
    grid on
    set(gcf,'name','Disturbance Estimate', figSpecs{:})
end


% Various experimental crossrange plans
if true
%     n = n+1;
%     figure(n)
%     plot(ts.L./ts.D,abs(ts.CR),lineSpecs{:})
%     hold on
%     plot(LoDmax, abs(ts.CR(i_LoDmax)),'ko',lineSpecs{:} )
%     set(gcf,'name','Crossrange vs L/D', figSpecs{:})
%     legend(' ','Max L/D')
    
%     n = n+1;
%     figure(n)
%     plot(ts.D, abs(ts.CR),lineSpecs{:})
%     set(gcf,'name','Crossrange vs Drag', figSpecs{:})
%     grid on
    
    n = n+1;
    figure(n)
    plot(ts.L'.*cos(x(:,5)), abs(ts.CR),lineSpecs{:})
    hold on
    plot(ts.L(i_Qmax)*cos(x(i_Qmax,5)), abs(ts.CR(i_Qmax)),'ko',ts.L(i_LoDmax),abs(ts.CR(i_LoDmax)),'r*', lineSpecs{:} )
%     plot(Lmax, abs(ts.CR(i_Lmax)),'r*',lineSpecs{:} )
    set(gcf,'name','Crossrange vs Lcos\gamma', figSpecs{:})
    legend(' ','Max Dyn Press', 'Max L/D')
    grid on
    
end

% Drag dynamics, for analysis
if false
    n = n+1;
    figure(n)
    plot(ts.energy_norm,ts.D, lineSpecs{:})
    axis([0,1,0,1.1*max(ts.D)])
    xlabel('Normalized Energy',textSpecs{:})
    ylabel('Drag (m/s^2)',textSpecs{:})
    
    n = n+1;
    figure(n)
    subplot 211
    title('Drag Dynamics')
    plot(ts.time,[ts.a; ts.b ],'--','LineWidth',2)
    hold all
    plot(ts.time,ts.a + ts.b .* ts.control,'LineWidth',2)
    legend('a','b','a+bu')
    axis([0,max(t),-2,1])
    xlabel('Time (s)',textSpecs{:})
    
    subplot 212
    plot(ts.energy_norm,[ts.a; ts.b],'LineWidth',2)
    hold all
    plot(ts.energy_norm,ts.a + ts.b .* ts.control,'LineWidth',2)
    axis([0,1,-2,1])
    xlabel('Normalized Energy',textSpecs{:})
    
end