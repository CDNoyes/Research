function [m, S, cl, cd] = aero_const()

% L/D ~ 0.24
% cl = 0.351;
% cd = 1.46;
% S = 15.8;
% m = 2804;

% L/D ~ 0.26
% cl = 0.357;
% cd = 1.408;

% L/D ~ 0.27
% cl = 0.3763;
% cd = 1.392;

% L/D ~ 0.28
% cl = 0.3879;
% cd = 1.3824;

if 0
    % Old heavy vehicle used these:
    % L/D ~ 0.29
    cl = 0.3986;
    cd = 1.3736;
    
    m = 5000;
    % BC = 225;
    BC = 185;
    S = m/BC/cd; % this keeps both mass and BC constant as cd changes
else
    % MSL-class examples use these
    cl = 0.351;
    cd = 1.46;
    S = 15.8;
    m = 2804;
    BC = m/(S*cd);
end

% [cd,cl] = srl_constant_aero(M, alpha)