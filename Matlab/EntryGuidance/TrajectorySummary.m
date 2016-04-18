%TRAJECTORY Creates a structure with all the data possibly
%required for trajectory tracking computations.

function output = TrajectorySummary(t,x,sigma,DR_t,CR_t)

tf = EntryAnalysis(t,x,DR_t,CR_t);
t = t(t<=tf);
x = x(1:length(t),:);
sigma = sigma(1:length(t));
% dtr = pi/180;
planet = Mars();
vm = VehicleModel();
% r_eq = planet.radiusEquatorial;
% hkm = (x(:,1)-r_eq)/1000;
E = 0.5*x(:,4).^2 + planet.mu/planet.radiusEquatorial - planet.mu./x(:,1); %Actual energy
En = (E-E(1))/(E(end)-E(1)); %Normalized Energy, 0 at entry and 1 at deployment



[theta_t,phi_t] = FinalLatLon(x(1,2),x(1,3),x(1,6),DR_t,CR_t);
output.target.lat = phi_t;
output.target.lon = theta_t;
output.target.DR = DR_t;
output.target.CR = CR_t;
output.energy = E;
output.energy_norm = En;
output.time = t;
output.state = x;
output.sigma = sigma;
output.control = cos(sigma);
[DR,CR] = Range(x(1,2),x(1,3),x(1,6),x(:,2),x(:,3));
output.DR = DR;
output.CR = CR;

for i = 1:length(t)
psi_d = DesiredHeading(x(i,2),x(i,3),theta_t,phi_t);
output.headingError(i) = psi_d-x(i,6);
end

observer = size(x,2) >= 9;
for i = 1:length(t)
    [output.g(i),output.L(i),output.D(i),~,output.M(i),output.a(i),output.rho(i),output.rho_dot(i),~,output.C_D(i),output.C_D_dot(i)] = EntryForces(x(i,:),planet,vm);
    if  observer
        [output.a(i),output.b(i),output.D_dot(i)] = DragFBL(output.g(i),output.L(i),output.D(i),x(i,1),x(i,4),x(i,5),output.rho(i),output.rho_dot(i),x(i,8),output.C_D(i),output.C_D_dot(i));
    else %Purely model
         [output.a(i),output.b(i),output.D_dot(i)] = DragFBL(output.g(i),output.L(i),output.D(i),x(i,1),x(i,4),x(i,5),output.rho(i),output.rho_dot(i),[],output.C_D(i),output.C_D_dot(i));
    end
    output.D_ddot(i) = output.a(i)+output.b(i)*output.control(i);
end

end