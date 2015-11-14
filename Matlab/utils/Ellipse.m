%ELLIPSE Plots an ellipse.
%   ELLIPSE(CENTER,A,B,THETA,LINESPECS) draws an ellipse centered at CENTER
%   with semimajor axis A, semiminor axis B, rotated counter-clockwise by
%   an angle THETA (in radians) with optional argument LINESPECS consisting
%   of typical Matlab linespec arguments.

function Ellipse(x_center,a,b,theta,LineSpecs)
if nargin < 5 || isempty(LineSpecs)
    LineSpecs = 'b';
end

x_center = x_center(:);
t = -pi:.01:pi;
x = a*cos(t);
y = b*sin(t);
R = [cos(theta) -sin(theta)
    sin(theta) cos(theta)];
XY = (R*[x;y] + repmat(x_center,1,length(t)))';

plot(XY(:,1),XY(:,2),LineSpecs)
end