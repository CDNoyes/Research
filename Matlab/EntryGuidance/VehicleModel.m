%VEHICLEMODEL Creates a model of an entry vehicle.
%   VEHICLEMODEL returns a structure with the parameters of an entry
%   vehicle. This function allows one to change a parameter such as the
%   vehicle mass and have it update everywhere it is used.

function VM = VehicleModel(mass)

%MSL-like vehicle:
VM.area = 15.8; % reference area, m^2
% VM.mass = 2804; % mass, kg MSL
if nargin == 0
    VM.mass = 5000; %  vehicle for SRL

else
    VM.mass = mass;
end
% VM.mass = 8500; % heavy vehicle for SRP 
VM.aerodynamics = @AerodynamicCoefficients;
VM.thrust = VM.mass*3.71*10; % T/W of 10 at ignition on Mars
VM.v_exit = 295*9.81; % Exit velocity
[Cd,Cl] = VM.aerodynamics(10);
VM.bc = VM.mass/VM.area/Cd;
end