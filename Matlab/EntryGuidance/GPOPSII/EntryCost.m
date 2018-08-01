function output = EntryCost( input )
xf = input.phase.finalstate;
% tf = input.phase.finaltime;


output.objective = -(xf(1)-input.auxdata.planet.radiusEquatorial)/1000;%+ 50*xf(5).^2;
% output.objective = xf(5).^2;
% output.objective = input.phase.integral;
end