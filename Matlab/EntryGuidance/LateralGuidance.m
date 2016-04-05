%Reference: "Drag-Based Predictive Tracking Guidance for Mars Precision
%Landing"

function [s,headingErrorLimit] = LateralGuidance(CR,V,currentsign)

lim0 = 3;
limf = 0.25;
v0 = 5505;
vf = 2000;
M = (lim0-limf)/(v0-vf);
CR0 = lim0-M*v0;

headingErrorLimit = limf*(V<vf) + (CR0 + V*M).*(V>=vf);

if abs(CR)>=headingErrorLimit && currentsign == sign(CR)
    s = -currentsign;
else
    s = currentsign;
end

end