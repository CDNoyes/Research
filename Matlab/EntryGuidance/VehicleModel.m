%VEHICLEMODEL Creates a model of an entry vehicle.
%   VEHICLEMODEL returns a structure with the parameters of an entry
%   vehicle. This function allows one to change a parameter such as the
%   vehicle mass and have it update everywhere it is used.

function VM = VehicleModel()

%MSL-like vehicle:
VM.area = 15.8; % reference wing area, m^2
VM.mass = 2804; % mass, kg
VM.aerodynamics = @AerodynamicCoefficients;
VM.thrust = VM.mass*3.71*10; % T/W of 10 at ignition on Mars
VM.v_exit = 295*9.81; % Exit velocity
end