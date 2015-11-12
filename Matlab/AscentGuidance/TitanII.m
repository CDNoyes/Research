function vehicleModel = TitanII()

vehicleModel.name = 'Titan II';
vehicleModel.nStages = 2;
vehicleModel.isp = [258,316];
vehicleModel.burnTime = [156,180];
vehicleModel.thrust = [1913e3,445e3]; % N
vehicleModel.mass.payload = 3580; % kg
vehicleModel.mass.fuel = [117910,25839]; % kg
vehicleModel.mass.structural = [4122.6,2748.4];
vehicleModel.mass.total = 154200;
vehicleModel.aero = @aerodynamicsSimple;
end

function [CL,CD] = aerodynamicsSimple(alpha)

CL = -0.01 + 0.018*alpha;
CD = 0.01 + 3.4*CL^2;

end

function [CL,CD] = aerodynamics(alpha, M)
CL0 = -0.01;
CD0 = 0.01;

if M >= 0 && M <= 0.8
    CDalpha = 0.29;
elseif M >= 0.8 && M <= 1.068
    CDalpha = M-0.51;
else
    CDalpha = 0.091 + 0.5/M;
end

if M >= 0 && M <= 0.25
    CLalpha = 2.8;
elseif M > 0.25 && M <= 1.1
    CLalpha = 2.8 + 0.447*(M-0.25);
elseif M > 1.1 && M <= 1.6
    CLalpha = 3.18 - 0.66*(M-1.1);
elseif M > 1.6 && M <= 3.6
    CLalpha = 2.85 + 0.350*(M-1.6);
else
    CLalpha = 3.55;
end

CL = CL0 + CLalpha*alpha;
CD = CD0 + CDalpha*alpha;
end