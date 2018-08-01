function output = EDLPhaseLink(input)

tf1 = input.phase(1).finaltime;
xf1 = input.phase(1).finalstate;

t02 = input.phase(2).initialtime;
x02 = input.phase(2).initialstate;
xf2 = input.phase(2).finalstate;

output.eventgroup(1).event = [x02(1:6)-xf1(1:6), t02-tf1];
output.objective = -xf2(7);