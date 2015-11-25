%Everything related to state transition matrices

%% 1st, 2nd order STM propagations vs numerical differences
%Should be exact for linear systems, like VDP with mu = 0. 

clear; clc;

mu = 0.25;
[fVDP,jVDP,hVDP] = VDP(mu);

x0 = [3;5];
[dim,ind] = Dimension(2,1,0);

[f,x0] = 