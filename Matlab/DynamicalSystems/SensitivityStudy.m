function SensitivityStudy()
% Compare STM with finite difference approx
clear; clc; close all;
[lineSpecs,textSpecs,figSpecs] = PlotSpecs();

% Linear, scalar system
% F = @(t,x,u) -x+u(t);
% J = @(x,u) -1;
% x0 = 0;

% [F,J] = InvertedPendulum();
% x0 = [0;0];

[Fmu,Jmu,Hmu] = VDP();
mu = 5;
F = @(T,X,U) Fmu(T,X,U,mu);
J = @(X,U) Jmu(X,U,mu);
H = @(X,U) Hmu(X,U,mu);
x0 = [0;0];

% Design a nominal trajectory (arbitrarily, in this case)
v = @(t,u0) u0 + .5*t.*cos(t);
u0 = -1;
u = @(t) v(t,u0);
[~,~,x,t] = Integrate(x0,F,u);


figure
if length(x0) < 2
    plot(t,x)
else
    plot(x(:,1),x(:,2))
end
set(gcf,'name','Nominal Trajectory', 'numbertitle','off','WindowStyle','docked')

%% Sensitivity to initial state and control

% via FD:
tic
dxf_dx0_fd = ComplexDiff(@(X)Integrate(X,F,u),x0);
dxf_du0_fd = ComplexDiff(@(U)Integrate(x0,F,@(t) v(t,U)),u0);
tFD = toc;

% via transition matrices:
tic
dxf_dx0_stm = ComputeTransition(J,t,x,u(t),0);
tSTM = toc;
dxf_du0_stm = dxf_dx0_stm(1:2,3);

disp(['Complex Differencing:  ',num2str(tFD),' seconds.'])
disp(['STM-Based Sensitivity: ',num2str(tSTM),' seconds.'])

%% Now let's examine the "size" of the linear region
% I.e. let's find the set of all initial states such that the error in
% final state predicted by the linear sensitivity is accurate to some
% tolerance.

% Define the grid of initial conditions
dx1 = .02;
dx2 = .02;
x10 = linspace(x0(1)-dx1,x0(1)+dx1,19);
x20 = linspace(x0(2)-dx2,x0(2)+dx2,21);
% figure
% hold on
% plot(x0(1),x0(2),'ro',lineSpecs{:})
e = [];
for i = 1:length(x10)
    for j = 1:length(x20)
        xf = Integrate([x10(i);x20(j)],F,u);
%         plot(x10(i),x20(j),'c*',lineSpecs{:})
        e(:,end+1) = (xf-x(end,:)');
        e_pred = dxf_dx0_stm(1:2,1:2)*([x10(i);x20(j)]-x0);
        e_pred_cd = dxf_dx0_fd(1:2,1:2)*([x10(i);x20(j)]-x0);
%         plot(e(1,end),e(2,end),'kx',lineSpecs{:})
%         plot(e_pred(1),e_pred(2),'b^')
%         plot(e_pred_cd(1),e_pred_cd(2),'g^')

        e_fd(i,j) =  norm( e(:,end) - e_pred_cd );
        e_stm(i,j) = norm( e(:,end) - e_pred );
        
    end
end

% lin = polyfit(e(1,:),e(2,:),1);
% quad = polyfit(e(1,:),e(2,:),2);
% plot(e(1,:),polyval(lin,e(1,:)),'r--',e(1,:),polyval(quad,e(1,:)),'m*')
% legend('Nominal IC','\delta x_0','\delta x_f','predicted \delta x_f')
% box on
% axis equal
% set(gcf,'name','Initial and Final States', 'numbertitle','off','WindowStyle','docked')

e_methods = (e_stm-e_fd)';
e_binary = sign(e_methods);

lim = [1e-3,1e-2,5e-2,1e-1]*norm([dx1,dx2])/sqrt(2);


% Another nice visualization to add would be the inverse of the prediction
% error plot. I.e. a contour of dxf1,dxf1 vs norm(dx0) with viscircles
% showing the boundary defining our linear region (whatever it may be).


figure
subplot_tight(1,2,1,[0.1,0.05])
contourf(x10-x0(1),x20-x0(2),(e_stm'))
title('Error between STM sensitivity prediction and actual integrated trajectory',textSpecs{:})
xlabel('\delta x1_{0}',textSpecs{:})
ylabel('\delta x2_{0}',textSpecs{:})
colorbar
box on
set(gcf,'name','Prediction Errors', 'numbertitle','off','WindowStyle','docked')

subplot_tight(1,2,2,[0.1,0.12])
contourf(x10-x0(1),x20-x0(2),e_fd')
title('Error between CD sensitivity prediction and actual integrated trajectory',textSpecs{:})
xlabel('\delta x1_{0}',textSpecs{:})
ylabel('\delta x2_{0}',textSpecs{:})
box on
colorbar

% figure
% surf(x10-x0(1),x20-x0(2),(e_stm-e_fd)')
% title('Negative-> STM is better, Positive -> Complex/Finite Difference is better',textSpecs{:})
% xlabel('\delta x1_{0}',textSpecs{:})
% ylabel('\delta x2_{0}',textSpecs{:})
% box on
% set(gcf,'name','STM vs Differencing', 'numbertitle','off','WindowStyle','docked')
% colorbar

figure
contourf(x10-x0(1),x20-x0(2),e_binary,[-1,.9,1]) %,[-1,1]
title('Blue -> STM is better, Red -> Complex Difference is better',textSpecs{:})
xlabel('\delta x1_{0}',textSpecs{:})
ylabel('\delta x2_{0}',textSpecs{:})
box on
set(gcf,'name','STM vs Differencing Binary Map', 'numbertitle','off','WindowStyle','docked')


figure
for i = 1:4
    linear_lim = lim(i);
    e_linear = e_stm'>=linear_lim;
    per = 100*sum(sum(~e_linear))/(length(x10)*length(x20));
    
    subplot_tight(2,2,i,[.1,.05])
    contourf(x10-x0(1),x20-x0(2),e_linear,[0,1])
    title(['In Blue: \{ \delta x_0 such that ||\delta x_f-\delta x_{pred}(tf)|| <= ',num2str(linear_lim),' \}, ',num2str(per),'%'],textSpecs{:})
    if i>2
        xlabel('\delta x1_{0}',textSpecs{:})
    end
    if mod(i,2)
        ylabel('\delta x2_{0}',textSpecs{:})
    end
end
box on
set(gcf,'name','Linear Region Binary Maps - STM', 'numbertitle','off','WindowStyle','docked')

figure
for i = 1:4
    linear_lim = lim(i);
    e_linear = e_fd'>=linear_lim;
    per = 100*sum(sum(~e_linear))/(length(x10)*length(x20));
    subplot_tight(2,2,i,[.1,.05])
    contourf(x10-x0(1),x20-x0(2),e_linear,[0,1])
    title(['In Blue: \{ \delta x_0 such that ||\delta x_f|| <= ',num2str(linear_lim),' \}, ',num2str(per),'%'],textSpecs{:})
    if i>2
        xlabel('\delta x1_{0}',textSpecs{:})
    end
    if mod(i,2)
        ylabel('\delta x2_{0}',textSpecs{:})
    end
end
box on
set(gcf,'name','Linear Region Binary Maps - CD', 'numbertitle','off','WindowStyle','docked')

%% Control Sensitivity Plots
clear xf
du0 = .05;
du = linspace(-du0,du0,200);

for i = 1:length(du)
    xf(:,i) = Integrate(x0,F,@(t)v(t,u0+du(i)));
    ef(:,i) = xf(:,i) - x(end,:)';
    e_u_stm(:,i) = (ef(:,i) - dxf_du0_stm*du(i));
    e_u_cd(:,i) = (ef(:,i) - dxf_du0_fd*du(i));
    n_stm(i) = norm(e_u_stm(:,i));
    n_cd(i) = norm(e_u_cd(:,i));
    
end


figure
plot(e_u_stm(1,:),e_u_stm(2,:),'k--')
hold on
plot(e_u_cd(1,:),e_u_cd(2,:),'g--')
scatter(e_u_stm(1,:),e_u_stm(2,:),[],du)
scatter(e_u_cd(1,:),e_u_cd(2,:),[],du)
title('',textSpecs{:})
xlabel('\delta x1_{f}',textSpecs{:})
ylabel('\delta x2_{f}',textSpecs{:})
box on
set(gcf,'name','Control Sensitivity', 'numbertitle','off','WindowStyle','docked')
title('Final State Prediction Error Colored by Input Deviation')
legend('STM','CD')
colorbar

figure
plot(du,n_stm,'k',du,n_cd,'r',lineSpecs{:})
set(gcf,'name','Control Sensitivity 2', 'numbertitle','off','WindowStyle','docked')
legend('STM','CD')
title('Norm of Prediction Error versus Control Deviation')
end

function [xf,tf,x,t] = Integrate(x0,f,u)

[t,x] = ode45(f,linspace(0,20,1000),x0,[],u);
xf = x(end,:).';
tf = t(end);



end