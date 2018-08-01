function output = EntryCostMultiphase( input )
xf = input.phase(end).finalstate;
p = input.parameter;

output.objective = -(xf(1)-input.auxdata.planet.radiusEquatorial)/1000 + 0*xf(5).^2;
% output.objective = -xf(3); % Max cross range when starting on equator and heading east or west
% output.objective = -xf(2); % Max down range when starting on equator and heading east

% Linkage constraints
for i = 1:(input.auxdata.nPhases-1)
    % Variables at Terminus of Phase i
    tf1 = input.phase(i).finaltime;
    xf1 = input.phase(i).finalstate;
    % Variables at Start of Phase i+1
    t02 = input.phase(i+1).initialtime;
    x02 = input.phase(i+1).initialstate;
    
    output.eventgroup(i).event = [x02-xf1, t02-tf1, p(i)-tf1];
    
end

end