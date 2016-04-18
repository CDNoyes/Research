% From Joel's dissertation, see eqns 2.40-2.44

function psi_d = DesiredHeading(theta,phi,theta_t,phi_t)

d = acos(cos(phi)*cos(phi_t)*cos(theta-theta_t) + sin(phi)*sin(phi_t));
if abs(theta-theta_t) < 1e-8
    if phi_t-phi > 0
        PHI = 0;
    else
        PHI = pi;
    end
else
    PHI = sign(theta_t-theta)*acos( (sin(phi_t) - sin(phi)*cos(d))/(cos(phi)*sin(d)) );
end
psi_d = pi/2 - PHI;
end

