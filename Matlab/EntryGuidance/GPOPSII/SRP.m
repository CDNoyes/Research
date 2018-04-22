function output = getOptTrajAccel()

output = SolveOCP();

sol = output.result.solution.phase(1);
mass = getMass(sol);
tSwitch = ThrustPlaneAnalysis(sol);
H = hamiltonian(sol);
[pitch,azimuth] = getAnglesFromUnitVector(sol.control(:,2:4));
pitch = pitch*180/pi;
azimuth = azimuth*180/pi;

disp(['Prop used: ',num2str(mass(1)-mass(end)), ' kg.'])

% Plots
figure(1)
plot(sol.time,sol.state(:,1:3))
legend('x','y','z')
title('Positions')
set(gcf,'name','Position vs Time','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(2)
plot(sol.time,sol.state(:,4:6))
legend('u','v','w')
title('Velocities')
set(gcf,'name','Velocity vs Time','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(3)
plot(sol.time,sol.control(:,2:4))
legend('u_x','u_v','u_z')
title('Control Directions')
set(gcf,'name','Control dirs vs Time','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(4)
legText = {'a_T'};
plot(sol.time,sol.control(:,1))
hold all
for i = 1:2
    if isreal(tSwitch(i))
        legText{end+1} = ['Theoretical Switch ',num2str(length(legText))];
        plot(tSwitch(i)*ones(size(sol.time)),linspace(min(sol.control(:,1)),max(sol.control(:,1)),length(sol.time)),'--')
    end
end
legend(legText{:})
title('Thrust Acceleration')
set(gcf,'name','Thrust Arcs','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(5)
plot(sol.time,mass)
title('Mass vs Time')
set(gcf,'name','Mass','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(6)
plot(sol.time,mass.*sol.control(:,1))
title('Thrust Magnitude')
set(gcf,'name','Thrust Force','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(7)
plot(sol.time,sol.costate)
legend('x','y','z','u','v','w','Location','best')
title('Costates')
set(gcf,'name','Costates','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(8)
plot(sol.time,H)
title('Hamiltonian')
set(gcf,'name','Hamiltonian','numbertitle','off')
set(gcf,'WindowStyle','docked')

figure(9)
plot(sol.time,pitch)
hold on
plot(sol.time,azimuth)
legend('Pitch','Yaw','Location','best')
title('Control Angles')
set(gcf,'name','Control Angles','numbertitle','off')
set(gcf,'WindowStyle','docked')


x = -sum(sol.state(:,1).^2 + sol.state(:,2).^2, 2).^0.5;
z = sol.state(:,3);
theta = atan2(z, x);
% figure(666)
% plot(sol.time, theta*180/pi)
% ylabel('Angle from target')
% figure(667)
E = sum(sol.state(:,4:6).^2,2)*0.5 - 3.71*sol.state(:,3)/1e5;
% Edot = sum(sol.state(:,4:6).*repmat(sol.control(:,1),[1,3]).*sol.control(:,2:4),2);
% plot(sol.time,Edot)
% 
% figure(10)
% % plot(sol.state(:,3),sol.state(:,1:2))
% plot(E,sol.state(:,1:3))
% legend('x','y','z')
% % axis('equal')
% title('Positions vs Energy')
% set(gcf,'name','Position vs E','numbertitle','off')
% set(gcf,'WindowStyle','docked')
% 
% figure(11)
% plot(E,sol.state(:,4:6))
% legend('x','y','z')
% title('Velocity vs Energy')
% set(gcf,'name','Velocity vs Energy','numbertitle','off')
% set(gcf,'WindowStyle','docked')
% 
% figure(13)
% plot(E,mass)
% title('Mass vs Energy')
% set(gcf,'name','Mass vs E','numbertitle','off')
% set(gcf,'WindowStyle','docked')

figure(12)
plot3(sol.state(:,1),sol.state(:,2),sol.state(:,3),'o-')
hold on 
plot3(sol.state(:,1),sol.state(:,2),0*sol.state(:,3),'k--') % Ground track 
title('Trajectory')
axis('equal')
set(gcf,'name','3D Trajectory','numbertitle','off')
set(gcf,'WindowStyle','docked')


% 
% 
% d = flipud(sum(sol.state(:,1:3).^2,2))/1e6;
% figure(15)
% plot(d,sol.state(:,1:3))
% legend('x','y','z')
% title('Positions vs d^2')
% set(gcf,'name','Position vs d^2','numbertitle','off')
% set(gcf,'WindowStyle','docked')
% 
% figure(16)
% plot(d,sol.state(:,4:6))
% legend('x','y','z')
% title('Velocity vs d^2')
% set(gcf,'name','Velocity vs d^2','numbertitle','off')
% set(gcf,'WindowStyle','docked')
% 
% figure(17)
% plot(d,mass)
% title('Mass vs d^2')
% set(gcf,'name','Mass vs d^2','numbertitle','off')
% set(gcf,'WindowStyle','docked')

output = output.result.solution;

end

function output = SolveOCP()
t0 = 0;
tmin = 1;
tmax = 40;

x0 = -3200;
y0 = 400;
z0 = 3200;
u0 = 625;
v0 = 0;
w0 = -270;

xmin = -15e3;
xmax = 15e3;
ymin = xmin;
ymax = xmax;
zmin = 0;
zmax = 100e3;

umin = -2e3;
umax = 2e3;
vmin = umin;
vmax = umax;
wmin = -3e3;
wmax = 100;

xf = 0;
yf = 0;
zf = 0;
uf = 0;
vf = 0;
wfmin = -10;
wfmax = -10;

% Thrust acceleration limits
amin = 40;
amax = 70;

% For slightly less conservatism, we could split the problem into three
% phases, and adjust the thrust bounds based on the mass at the
% transitions. Let's gauge how suboptimal these solutions are relative to
% constant thrust solutions.

%% Bounds Definition
iphase = 1;

bounds.phase(iphase).initialtime.lower    = t0;
bounds.phase(iphase).initialtime.upper    = t0;
bounds.phase(iphase).finaltime.lower      = tmin;
bounds.phase(iphase).finaltime.upper      = tmax;

bounds.phase(iphase).initialstate.lower          = [x0 y0 z0 u0 v0 w0];
bounds.phase(iphase).initialstate.upper          = [x0 y0 z0 u0 v0 w0];

bounds.phase(iphase).state.lower          = [xmin ymin zmin umin vmin wmin];
bounds.phase(iphase).state.upper          = [xmax ymax zmax umax vmax wmax];

bounds.phase(iphase).finalstate.lower          = [xf yf zf uf vf wfmin];
bounds.phase(iphase).finalstate.upper          = [xf yf zf uf vf wfmax];

bounds.phase(iphase).control.lower = [amin, -1,-1,-1];
bounds.phase(iphase).control.upper = [amax,  1, 1, 1];

bounds.phase(iphase).path.lower = 1;
bounds.phase(iphase).path.upper = 1;

bounds.phase(iphase).integral.lower = 0;
bounds.phase(iphase).integral.upper = tmax*amax^2;

%% Guess
guess.phase(iphase).time           = [0; 13];
guess.phase(iphase).state     = [bounds.phase(iphase).initialstate.lower; bounds.phase(iphase).finalstate.lower];
guess.phase(iphase).control   = [amax 1 0 0; amin 0 0 -1];
guess.phase(iphase).integral = [ 0.5*(amin+amax)*15];

%% NLP Parameters and Optimization Call
setup.name = 'SRP Optimization';
setup.functions.continuous = @Dynamics;
setup.functions.endpoint = @Cost;
setup.nlp.solver = 'snopt';
if strcmp(setup.nlp.solver,'snopt')
    setup.nlp.snoptoptions.tolerance = 1e-6; % Default 1e-6
    setup.nlp.snoptoptions.maxiterations = 2000; % Default 2000
elseif strcmp(setup.nlp.solver,'ipopt')
    setup.nlp.ipoptoptions.linear_solver = 'mumps'; % Default 'mumps' ('ma57')
    setup.nlp.ipoptoptions.tolerance = 1e-6; % Default 1e-7
    setup.nlp.ipoptoptions.maxiterations = 2000; % Default 2000
end
setup.bounds = bounds;
setup.guess = guess;
setup.derivatives.supplier = 'sparseFD'; % 'sparseFD' or 'sparseBD' or 'sparseCD' or 'adigator'
setup.derivatives.derivativelevel = 'first'; % 'first' or 'second'
setup.derivatives.dependencies = 'sparseNaN'; % 'full', 'sparse' or 'sparseNaN'
setup.scales.method = 'automatic-guessUpdate'; % 'none' or 'automatic-bounds' or 'automatic-guess' or 'automatic-guessUpdate' or 'automatic-hybrid' or 'automatic-hybridUpdate' or 'defined'
setup.mesh.method = 'hp-PattersonRao'; % 'hp-PattersonRao' or 'hp-DarbyRao' or 'hp-LiuRao'
setup.mesh.tolerance = 1e-3; % Default 1e-3
setup.mesh.maxiterations = 50; % Default 10
setup.method = 'RPM-Differentiation'; % 'RPM-Differentiation' or 'RPM-Integration'
setup.displaylevel = 0; % 0 = no output. 1 = only mesh refinement. 2 = NLP solver iteration output and mesh refinement

output = gpops2(setup);

end

function output = Cost(input)

output.objective = input.phase(1).integral;

end

function output = Dynamics(input)


%---------------------%
% Dynamics in Phase 1 %
%        SRP          %
%---------------------%
t1 = input.phase(1).time;
s1 = input.phase(1).state;
control1 = input.phase(1).control;

% Variables
x1     = s1(:,1);
y1     = s1(:,2);
z1     = s1(:,3);
u1     = s1(:,4);
v1     = s1(:,5);
w1     = s1(:,6);

aT    = control1(:,1);
ux    = control1(:,2);
uy    = control1(:,3);
uz    = control1(:,4);

gravity1 = 3.71;

% EOM
xdot1 = u1;
ydot1 = v1;
zdot1 = w1;
udot1 = aT.*ux;
vdot1 = aT.*uy;
wdot1 = aT.*uz - gravity1;

% Output
output(1).integrand = aT;
output(1).dynamics = [xdot1 ydot1 zdot1 udot1 vdot1 wdot1];
output(1).path = sqrt(sum(control1(:,2:4).^2,2)); %Norm of the control directions

end


function m = getMass(sol)

m0 = 8500;
Isp = 290;
g0 = 9.81;
ve = Isp*g0;

m = m0*exp(-cumtrapz(sol.time,sol.control(:,1))/ve);


end

function tSwitch = ThrustPlaneAnalysis(sol)

% Compute the true quantities involved in the thrust plane solution

state = sol.state;
time = sol.time;
tau0 = time(end);
amin = min(sol.control(:,1));
amax = max(sol.control(:,1));

CR = sol.costate(1,1:3);
CV = sol.costate(end,4:6);

g = 3.71;
gvec = g*[0,0,-1];
Vbar = state(end,4:6)-state(1,4:6) - gvec*time(end);
Rbar = state(end,1:3)-state(1,1:3) - state(1,4:6)*time(end)-0.5*gvec*time(end)^2;

CRu = CR/norm(CR);
CVu = CV/norm(CV);
Vbaru = Vbar/norm(Vbar);
Rbaru = Rbar/norm(Rbar);
vecs = [CRu; CVu; Vbaru; Rbaru];
% figure
% hold all
% for i = 1:4
%     plot3([0,vecs(i,1)],[0,vecs(i,2)],[0,vecs(i,2)])
% end
% legend('CR','CV','Vtilde','Rtilde')

D = dot(CV,CV);
E = dot(CR,CV);
F = dot(CR,CR);

tSwitch = time(end) - roots([F,2*E,D-1]);

xTP = Vbar./norm(Vbar);
zTP = cross(Vbar,Rbar)./norm(cross(Vbar,Rbar));
yTP = cross(zTP,xTP);
RTP = [xTP',yTP',zTP']; % The transpose rotates from I to TP

CVTP = RTP'*CV';
CRTP = RTP'*CR';
gTP = RTP'*gvec';
VbarTP = RTP'*Vbar';
RbarTP = RTP'*Rbar';

% uTP = (RTP'*sol.control(:,2:4)')';

% get angles in thrust plane representation
[phiV,thetaV] = getAnglesFromUnitVector(VbarTP');
[phiR,thetaR] = getAnglesFromUnitVector(RbarTP');
[phiCV,thetaCV] = getAnglesFromUnitVector(CVTP');
[phiCR,thetaCR] = getAnglesFromUnitVector(CRTP');
[phig,~] = getAnglesFromUnitVector(gTP');

tanThetaCV = tan(thetaCV);
cosThetaCV = cos(thetaCV);

% phiCR = asin(CRTP(3)); % Pitch angle of CR, should be 0
% phiCV = asin(CVTP(3)); % Pitch angle of CV, should be 0

PHI = cos(thetaCR)*cos(thetaCV)+sin(thetaCR)*sin(thetaCV); % = dot(CR,CV)/norm(CR)/norm(CV)

E_check = norm(CR)*norm(CV)*PHI;

% PGK computed the general way for comparison
% These turn out to be much better than those below
P = dot(CR,Vbar);
G = sqrt(F)*(F*dot(CV,Vbar)-E*dot(CR,Vbar))/(D*F-E^2);
K = -dot(CR,Rbar)-0.5*dot(CV,Vbar);

[Pi,Gi,Ki] = getPGK(5, tau0, amin, amax, norm(CR), norm(CV), PHI);
% This function returns something similar enough to believe it

% Compute Vbar and Rbar from D,E,F,P,G,K and CR and CV (eq 4.1)
Vbar_ = (P/F - E*G/sqrt(F^3))*CR' + G/sqrt(F)*CV';
Rbar_ = 0.5*( (3*E^2/sqrt(F^5) - D/sqrt(F^3))*G - 2*K/F - 3*E*P/F^2)*CR + ...
    (P/F-E*G/sqrt(F^3))*CV;

disp('Velocity Error (m/s):')
norm(Vbar-Vbar_')
disp('Position Error (m):')
norm(Rbar-Rbar_)

end

function H = hamiltonian(sol)
aT = sol.control(:,1);
CR = sol.costate(1,1:3);
gvec = [0;0;-3.7].';

for i = 1:length(aT)
    H(i) = aT(i) + dot(CR,sol.state(i,4:6)) + dot(sol.costate(i,4:6),aT(i)*sol.control(i,2:4)+gvec);
    
end
end

function [pitch,azimuth] = getAnglesFromUnitVector(u)

for i = 1:size(u,1)
%     pitch(i) = 180-asin(u(i,3))*180/pi;
    pitch(i) = atan2(u(i,3),norm(u(i,1:2)));
    azimuth(i) = atan(u(i,2)/u(i,1));
end
end

% function [CTCR2,CTCR1,CTCR0] = Ci(a0,af,R,V) 
% % V is magnitude of Vtilde
% % R is magnitude of Rtilde

% CA1 = V^2*sin(thetaV)^2-(a0*tf)^2;
% CA2 = 2*V^2*sin(thetaV)*cos(thetaV);
% CA3 = V^2*cos(thetaV)^2-(a0*tf)^2;
% 
% CB1 = -2*a0^2*tf;
% CB2 = -2*af*V*sin(thetaV);
% CB3 = -2*af*V*cos(thetaV);
% 
% CC1 = 2*(af-a0)*V*sin(thetaV);
% CC2 = 2*(af-a0)*V*cos(thetaV);
% 
% CD1 = (af^2-a0^2);
% CD2 = -2*af*(af-a0);
% CD3 = (af-a0)^2;
% 
% 
% CT1 = V*Vfmag*cos(phiVf)*sin(thetaVf)*sin(thetaV) + R*g*cos(thetag)*cos(phig)*cos(thetaR);
% CT2 = -V*Vfmag*cos(phiVf)*sin(thetaVf)*cos(thetaV) + R*g*cos(phig)*sin(thetag)*cos(thetaR);
% CT3 = -af*R*cos(thetaR);
% CT4 = V*Vfmag*cos(phiVf)*cos(thetaVf)*sin(thetaV)-R*g*cos(phig)*cos(thetag)*sin(thetaR);
% CT5 = -V*Vfmag*cos(phiVf)*cos(thetaVf)*cos(thetaV) - R*g*cos(phig)*sin(thetag)*sin(thetaR);
% CT6 = af*R*sin(thetaR);
% CT7 = -af*V*sin(thetaV);                 % sin(thetaV)
% CT8 = af*V*cos(thetaV);
% CT9 = -af*R*cos(thetaR);
% CT10 = af*R*sin(thetaR); % Same as CT6? Possible typo in thesis
% 
% CTV1 = CT1*cos(thetaCV) + CT2*sin(thetaCV) + CT3;
% CTV2 = CT4*cos(thetaCV) + CT5*sin(thetaCV) + CT6;
% CTV3 = CT7*cos(thetaCV) + CT8*sin(thetaCV);
% CTV4 = CB1*sin(thetaCV) + CB2;
% CTV5 = CB1*cos(thetaCV) + CB3;
% 
% 
% CTCR2 = CTV3^2*CA1 + CTV1^2*CD3 + ...
%         CTV1*CTV3*CC1 + CTV3*CTV4*CT9 + ...
%         CTV1*CD2*CT9 + CD1*CT9^2;
%     
% CTCR1 = CTV3*CA2 + CTV3*(CTV4*CT10+CTV5*CT9) + ...
%         CTV1*(CD2*CT10 + CTV2*CD3 + CTV3*CC2) + ...
%         CTV2*(CD2*CT9 + CTV1*CD3 + CTV3*CC1) + ...
%         2*CD1*CT9*CT10;
% 
% CTCR0 = CTV3^2*CA3 + CTV2^2*CD3 + ...
%         CTV2*CTV3*CC2 + CTV3*CTV5*CT10 + ...
%         CTV2*CD2*CT10 + CD1*CT10^2;
% 
% 
% end

function [P,G,K] = getPGK(scenario, tau0, amin, amax, CR, CV, phi)
switch scenario
    case 1
        rootTerm = sqrt(CV^2+2*CR*CV*phi*tau0 + CR^2*tau0^2);
        P = amin*(CV-rootTerm);
        G = -amin*log((rootTerm + CR*tau0 + CV*phi)/(CV+CV*phi));
        K = 0.5*amin*tau0*rootTerm;
    case 2
        
    case 3
        
    case 4
        
    case 5
        P = amax*(CV-sqrt(CV^2+2*CR*CV*phi*tau0 + CR^2*tau0^2));
        G = -log( ((sqrt(CV^2+2*CR*CV*phi*tau0 + CR^2*tau0^2) + CR*tau0 + CV*phi)/(CV+CV*phi))^amax) ...
            -log( ( (1-sqrt(CV^2*phi^2-(CV^2-1)))/(1+sqrt(CV^2*phi^2-(CV^2-1))) )^(amax-amin));
        K = 0.5*amax*tau0*sqrt(CV^2+2*CR*CV*phi*tau0 + CR^2*tau0^2) - ...
            (amax-amin)*sqrt(CV^2*phi^2-(CV^2-1))/CR;
end
end