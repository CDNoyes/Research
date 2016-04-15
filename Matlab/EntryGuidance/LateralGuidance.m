%Reference: "Drag-Based Predictive Tracking Guidance for Mars Precision
%Landing"

function [s,headingErrorLimit] = LateralGuidance(CR,IV,currentsign)


% Velocity Based Profile:
% V = IV;
% lim0 = 3;
% limf = 0.25;
% v0 = 5505;
% vf = 2000;
% M = (lim0-limf)/(v0-vf);
% CR0 = lim0-M*v0;
% 
% headingErrorLimit = limf*(V<vf) + (CR0 + V*M).*(V>=vf);

% Drag Based Profile:
D = IV;
[cr_coeff,s] = polyfit([0,80.62, 7.416, 120],[0,5.05-1, 0.8, 6-1],3);
headingErrorLimit = polyval(cr_coeff,D);


if abs(CR)>=headingErrorLimit && currentsign == sign(CR)
    s = -currentsign;
else
    s = currentsign;
end

end