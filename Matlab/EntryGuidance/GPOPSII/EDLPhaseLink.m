function output = EDLPhaseLink(input)

tf1 = input.phase(1).finaltime;
xf1 = input.phase(1).finalstate;

t02 = input.phase(2).initialtime;
x02 = input.phase(2).initialstate;
xf2 = input.phase(2).finalstate;

x0_srp = targetless_transform(xf1, input.auxdata.planet.radiusEquatorial+input.auxdata.target.altitude);

%TODO: Transform from entry to srp 
output.eventgroup(1).event = [x02(1:6)-x0_srp];
output.objective = -xf2(7);
end


function x_srp = targetless_transform(x_entry, Rtarget)

R = x_entry(:,1);
V = x_entry(:,4);
fpa = x_entry(:,5);
vx = -V.*cos(fpa);
vz = V.*sin(fpa);

x_srp = [0, 0, R-Rtarget, vx, 0, vz];

end