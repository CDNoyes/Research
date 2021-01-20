function [m,S,cl,cd] = aero_const()

% L/D ~ 0.26
% cl = 0.357;
% cd = 1.408;

% L/D ~ 0.27
% cl = 0.3763;
% cd = 1.392;

% L/D ~ 0.28
cl = 0.3879;
cd = 1.3824;

m = 5000;
BC = 225;
S = m/BC/cd; % this keeps both mass and BC constant as cd changes 

% [cd,cl] = srl_constant_aero(M, alpha)