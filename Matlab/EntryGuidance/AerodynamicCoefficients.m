%Calculates the coefficient of drag and coefficient of lift for a given
%Mach number for the MSL-like vehicle described in Benito's dissertation.

function [CD,CL,dCDdM] = AerodynamicCoefficients(Mvector)

pD = [2.598e4; -1022; -2904; 678.6; -44.33; 1.373];
qD = [1.505e4; 1687; -2651; 544.1; -34.11; 1];

pL = [1.172e4;-3654; 485.6; -14.61;0.4192];
qL = [2.53e4; -7846; 1086; -28.35; 1];

for i = 1:length(Mvector)
    M = Mvector(i);
    mD = M.^(0:5);
    CD(i,1) = (mD*pD)/(mD*qD);

    dCDdM(i,1) = ( ((1:5).*mD(1:5))*pD(2:6) )/(mD*qD) - CD(i)/(mD*qD)*( ((1:5).*mD(1:5))*qD(2:6) );


    mL = M.^(0:4);
    CL(i,1) = (mL*pL)/(mL*qL);

end


%% Old, ugly code that I didn't write...

% Drag Coefficient parameters
% RD = 5;
% SD = 5;
% pD0 =  2.598E4;
% pD1 = -1022;
% pD2 = -2904;
% pD3 =  678.6;
% pD4 = -44.33;
% pD5 =  1.373;
% pD = [pD0 pD1 pD2 pD3 pD4 pD5].';
% qD0 =  1.505E4;
% qD1 =  1687;
% qD2 = -2651;
% qD3 =  544.1;
% qD4 = -34.11;
% qD = [qD0 qD1 qD2 qD3 qD4].';

% Lift Coefficient parameters
% RL = 4;
% SL = 4;
% pL0 =  1.172E4;
% pL1 = -3654;
% pL2 =  485.6;
% pL3 = -14.61;
% pL4 =  0.4192;
% pL = [pL0 pL1 pL2 pL3 pL4].';
% qL0 =  2.53E4;
% qL1 = -7846;
% qL2 =  1086;
% qL3 = -28.35;
% qL = [qL0 qL1 qL2 qL3].';


% % CD
% numD = 0;
% denD = 0;
% for i=1:RD+1
%     numD_new = pD(i,1)*M.^(i-1);
%     numD = numD + numD_new;
% end
% for j=1:SD
%     denD_new = qD(j,1)*M.^(j-1);
%     denD = denD + denD_new;
% end
% denD = M.^SD + denD;
% CD = numD./denD;

% CL
% numL = 0;
% denL = 0;
% for i=1:RL+1
%     numL_new = pL(i)*M.^(i-1);
%     numL = numL + numL_new;
% end
% for j=1:SL
%     denL_new = qL(j)*M.^(j-1);
%     denL = denL + denL_new;
% end
% denL = M.^SL + denL;
% CL = numL./denL;