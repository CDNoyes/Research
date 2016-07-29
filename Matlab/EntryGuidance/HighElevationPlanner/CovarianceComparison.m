%Compare variance of multiple similar solutions

clear; close all; clc;


dir = 'EntryGuidance/HighElevationPlanner/Trajectories/';
% sets = {'ControlHEP2','ControlHEP','ControlHEP-PF','ControlHEP-UT'};
% labels = {'Standard Planner','Standard Planner - No Constraint','Robust Planner, PF','Robust Planner, UT'};
sets = {'ControlHEP2','ControlHEP-UT'};
labels = {'Standard Planner','Robust Planner, UT'};
dtr = pi/180;

T = linspace(0,300,500);
sigmaMin = 18.19*dtr;
sigmaMax = 87.13*dtr;

delta.rho = 0; %linspace(0*-.20,.20, 29);
delta.CD = (0.04)*randn(1,5000);%linspace(-.1,.1,30);
for i = 1:length(sets)
    
    load([dir,sets{i},'.mat'])
%     [CD_err,W] = SigmaPoints(0,0.04^2); % +-12% error in drag coefficient
    CD_err = linspace(-.20,.20,49);
    W = 0;
	ref.sp.delta.rho = zeros(size(delta.CD));
    ref.sp.delta.CD = CD_err;
    ref.sp.delta.p0 = normpdf(CD_err,0,0.04);
    x0 = repmat(ref.state(1,:)',1,length(CD_err));
    ref.P0 = zeros(9);
    ref.P0(7,7) = 0.04^2;
    ref.sp.state = x0;
%     ref.sp.delta = delta;
    ref.sp.W = W;
    ref.planet = Mars();
    ref.vehicle = VehicleModel();
    control = @(t) BankAngleProfile(t,p(1),p(2),p(3),sigmaMin,sigmaMax);
    u(:,i) = control(T)/dtr;
    xf{i} = [];
    ref.method = 'PF'; %force all estimates to use same method
    method{i} = ref.method;
    [costUT(i),xMeanUT{i},PfUT{i}] = HighElevationCostFunctionStochastic(p,ref.planet,ref.vehicle,ref.target.DR,ref.target.CR,ref); % Estimated cost
    % Conduct MonteCarlo Integrations Here
    for j = 1:length(delta.rho)
        for k = 1:length(delta.CD)
            ref.delta.CD = delta.CD(k);
            ref.delta.rho = delta.rho(j);
            [t,x] = OpenLoopSim(ref.state(1,:)',control,ref);
            xf{i}(end+1,1:6) = x(end,:);
        end
    end
    rp = ref.planet.radiusEquatorial;
%     [dr,cr] = Range(ref.state(1,2),ref.state(1,3),ref.state(1,6),mean(xf{i}(:,2)),mean(xf{i}(:,3)));
%     d = norm([dr-ref.target.DR,cr-ref.target.CR]);
%     costTrue(i) = d; % The "true" cost based on Monte Carlo results
    xf{i}(:,7) = (xf{i}(:,1)-rp)/1000; %Altitude
    
end


%% Figures and Analysis
figure(1)
plot(T,u)
xlabel('Time (s)')
ylabel('Bank Profile (deg)')
set(gcf,'WindowStyle','docked','name','Bank Profiles','numbertitle','off')
legend(labels{:},'Location','best')
clc
for i = 1:length(sets)
    Mean(i,:) = mean(xf{i});
    Std3(i,:) = 3*std(xf{i});
    Pf{i} = cov(xf{i}(:,1:6),1);
    errCov{i} = (Pf{i}-PfUT{i}(1:6,1:6))./Pf{i};
    Std3UT{i} = 3*(PfUT{i}.^0.5); % Should account for covariance
    costTrue(i) = (mean(xf{i}(:,2))-ref.target.lon)^2 + (mean(xf{i}(:,3))-ref.target.lat)^2 + trace(Pf{i}(2:3,2:3));
    disp(labels{i})
    disp('       lon    lat     alt    vel    fpa')
    disp(sprintf('Mean:  %0.2f  %0.2f  %2.2f  %3.1f   %2.1f',(Mean(i,2)-ref.target.lon)/dtr,(Mean(i,3)-ref.target.lat)/dtr,(Mean(i,1)-rp)/1000,Mean(i,4),Mean(i,5)))
    disp(sprintf('3-sig: %0.4f  %0.4f  %2.2f  %3.1f   %2.1f', Std3(i,2)/dtr, Std3(i,3)/dtr,(Std3(i,1))/1000,Std3(i,4),Std3(i,5)))
    disp(sprintf([method{i},' Mean:  %0.2f  %0.2f  %2.2f  %3.1f   %2.1f'],(xMeanUT{i}(2)-ref.target.lon)/dtr,(xMeanUT{i}(3)-ref.target.lat)/dtr,(xMeanUT{i}(1)-rp)/1000,xMeanUT{i}(4),xMeanUT{i}(5)))
    disp(sprintf([method{i},' 3-sig: %0.4f  %0.4f  %2.2f  %3.1f   %2.1f'], Std3UT{i}(2,2)/dtr, Std3UT{i}(3,3)/dtr,(Std3UT{i}(1,1))/1000,Std3UT{i}(4,4),Std3UT{i}(4,5)))
    disp(' ');
    
end


figure(2)
set(gcf,'WindowStyle','docked','name','Lat-Lon Footprint','numbertitle','off')
for i = 1:length(sets)
    subplot(1,length(sets),i)
    plot(xf{i}(:,2)/dtr, xf{i}(:,3)/dtr,'x')
    hold on
    plot(ref.target.lon/dtr,ref.target.lat/dtr,'r^','MarkerSize',10)
    plot(mean(xf{i}(:,2))/dtr, mean(xf{i}(:,3))/dtr,'rx','MarkerSize',10)
%     DrawNormalEllipse(Mean(i,2:3)/dtr,Pf{i}(2:3,2:3)/dtr,3,'r--');
    plot(xMeanUT{i}(2)/dtr,xMeanUT{i}(3)/dtr,'k*','MarkerSize',10)
    DrawNormalEllipse(xMeanUT{i}(2:3)/dtr,PfUT{i}(2:3,2:3)/dtr,3,'k--');
    axis([-72.7 -71.9 -41.4 -41.28])
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    title(labels{i})
    legend('MC','Target','True Mean','Est. Mean','Est. 3-\sigma')
end


figure(3)
set(gcf,'WindowStyle','docked','name','h-v','numbertitle','off')
for i = 1:length(sets)
    subplot(1,length(sets),i)
    hold all
    ParachuteDeploymentConstraints(true);
    plot(xf{i}(:,4), xf{i}(:,7),'x')
    title(labels{i})
    axis([290, 600, 4,18])
end

