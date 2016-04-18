%ENTRYANALYSIS Performs rudimentary analysis of an entry trajectory.
%   ENTRYANALYSIS(TIME,STATE,DOWNRANGE,CROSSRANGE) displays the total
%   trajectory length, final altitude, flight path angle, ranges etc. If
%   the target downrange and crossrange are provided the range error is
%   computed.

function [tf,d] = EntryAnalysis(t,x,dr,cr)

if nargin < 2
    error('Trajectory Analysis requires at least the time and state histories.')
end

dtr = pi/180;                           % deg to rad
planet = Mars();
r_eq = planet.radiusEquatorial;
hkm = (x(:,1)-r_eq)/1000;
[DR,CR] = Range(x(1,2),x(1,3),x(1,6),x(:,2),x(:,3));
tf = findTrajLength(t,x);

disp(['Trajectory Duration: ',num2str(tf), ' s'])
disp(['Final altitude:      ',num2str(hkm(end)), ' km'])
disp(['Final FPA:           ',num2str(x(end,5)/dtr), ' deg'])
disp(['Final Downrange:     ',num2str(DR(end)),' km'])
disp(['Final Crossrange:    ',num2str(CR(end)), ' km'])
if nargin >= 4 && ~isempty([dr,cr])
    d = norm([dr-DR(end),cr-CR(end)]);
    if d > 1
        disp(['Total range error:   ',num2str(d), ' km'])
    else
        disp(['Total range error:   ',num2str(d*1000), ' m'])
    end
else
    d = [];
end

end

function tf = findTrajLength(t,x)

d = diff(x(:,1));
idx = find(d==0,1);
if isempty(idx)
    tf = t(end);
else
    tf = t(idx);
end
end