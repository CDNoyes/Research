%Compare variance of multiple similar solutions

clear


dir = 'EntryGuidance/HighElevationPlanner/';
sets = {'ControlMargin','ControlNoMargin','ControlNegLift1'};

dtr = pi/180;

T = linspace(0,300,500);

ref.state = [3540e3; -90.07*dtr; -43.90*dtr; 5505; -14.15*dtr; 4.99*dtr]';
ref.vehicle = VehicleModel();
ref.planet = Mars();
rp = ref.planet.radiusEquatorial;
delta.rho = linspace(0*-.20,.20, 29);
delta.CD = 0;%linspace(-.1,.1,30);
for i = 1:length(sets)
    
    load([dir,sets{i},'.mat'])
    control = @(t) BankAngleProfile(t,p(1),p(2),p(3),sigmaMin,sigmaMax);
    [ref.target.theta,ref.target.phi] = FinalLatLon(ref.state(2),ref.state(3),ref.state(6),DR,CR);
    u(:,i) = control(T)/dtr;
    xf{i} = [];
    % Conduct MonteCarlo Integrations Here
    for j = 1:length(delta.rho)
        for k = 1:length(delta.CD)
            ref.target.DR = DR;
            ref.target.CR = CR;
            ref.delta.CD = delta.CD(k);
            ref.delta.rho = delta.rho(j);
            [t,x] = OpenLoopSim(ref.state',control,ref);
            xf{i}(end+1,1:6) = x(end,:);
        end
    end 
    xf{i}(:,7) = (xf{i}(:,1)-rp)/1000; %Altitude

end


%% Figures and Analysis
% figure(1)
% plot(T,u)
% xlabel('Time (s)')
% ylabel('Bank Profile (deg)')
% set(gcf,'WindowStyle','docked','name','Bank Profiles','numbertitle','off')
clc
for i = 1:length(sets)
    Mean(i,:) = mean(xf{i});
    Std3(i,:) = 3*std(xf{i});
    disp(sets{i})
    disp('       lon    lat    alt    vel    fpa')
    disp(sprintf('Mean: %0.2f %0.2f %2.2f  %3.1f %2.1f',Mean(i,2)-ref.target.theta,Mean(i,3)-ref.target.phi,(Mean(i,1)-rp)/1000,Mean(i,4),Mean(i,5)))
    disp(sprintf('3-sig: %0.2f %0.2f %2.2f  %3.1f %2.1f',Std3(i,2),Std3(i,3),(Std3(i,1))/1000,Std3(i,4),Std3(i,5)))
    disp(' ');
    
end


figure(2)
subplot 131
plot(xf{1}(:,2)/dtr, xf{1}(:,3)/dtr,'x')
hold on
plot(ref.target.theta/dtr,ref.target.phi/dtr,'r^','MarkerSize',10)
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
subplot 132
plot(xf{2}(:,2)/dtr, xf{2}(:,3)/dtr,'o')
hold on
plot(ref.target.theta/dtr,ref.target.phi/dtr,'r^','MarkerSize',10)
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
subplot 133
plot(xf{3}(:,2)/dtr, xf{3}(:,3)/dtr,'*')
hold on
plot(ref.target.theta/dtr,ref.target.phi/dtr,'r^','MarkerSize',10)
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
%want lat/lon spread, altitude spread, FPA spread, V spread

figure(3)
subplot 131
hold all
ParachuteDeploymentConstraints(true)
plot(xf{1}(:,4), xf{1}(:,7),'x')
title(sets{1})
axis([290, 600, 4,18])

subplot 132
hold all
ParachuteDeploymentConstraints(true)
plot(xf{2}(:,4), xf{2}(:,7),'o')
title(sets{2})
axis([290, 600, 4,18])

subplot 133
hold all
ParachuteDeploymentConstraints(true)
plot(xf{3}(:,4), xf{3}(:,7),'*')
title(sets{3})
axis([290, 600, 4,18])
nBins = 35;
figure(4)
subplot 311
hist(xf{1}(:,5)/dtr,nBins)
subplot 312
hist(xf{1}(:,5)/dtr,nBins)
subplot 313
hist(xf{1}(:,5)/dtr,nBins)
