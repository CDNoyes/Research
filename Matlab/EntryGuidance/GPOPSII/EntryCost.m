function output = EntryCost( input )
xf = input.phase.finalstate;
% tf = input.phase.finaltime;

% hf = ((xf(1)-input.auxdata.planet.radiusEquatorial)/1000 - 8)^2;
% output.objective = hf;

output.objective = -(xf(1)-input.auxdata.planet.radiusEquatorial)/1000;%+ 50*xf(5).^2;
% output.objective = -xf(5);
% output.objective = input.phase.integral;
% output.objective = xf(4); % Minimize final velocity 
end