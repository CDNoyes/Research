function sigma = HeadingAlignment(state,L,target,gain)

[r,theta,phi,V,gamma,psi] = ParseState(state);

delta = 2*asin( sin(phi-target.lat)^2 + cos(phi)*cos(target.lat)*sin(0.5*(theta-target.lon))^2 );

psi_d = pi/2 - acos( (sin(target.lat)-sin(phi)*cos(delta))/(sin(delta)*cos(phi)) );
psi_dot = 0;
sinSigma = Saturate(-V*cos(gamma)/L*(-gain*(psi-psi_d)+psi_dot), -1,1);
sigma = asin(sinSigma);

end