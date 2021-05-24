%Calculates the coefficient of drag and coefficient of lift for a given
%Mach number for the MSL-like vehicle described in Benito's dissertation.

function [CD,CL,dCDdM,dCLdM] = AerodynamicCoefficients(Mvector)

partials = (nargout > 2);
% 
% pD = [2.598e4; -1022; -2904; 678.6; -44.33; 1.373];
% qD = [1.505e4; 1687; -2651; 544.1; -34.11; 1];
% 
% pL = [1.172e4;-3654; 485.6; -14.61;0.4192];
% qL = [2.53e4; -7846; 1086; -28.35; 1];
% 
% for i = 1:length(Mvector)
%     M = Mvector(i);
%     mD = M.^(0:5);
%     CD(i,1) = (mD*pD)/(mD*qD);
% 
% 
% 
%     mL = M.^(0:4);
%     CL(i,1) = (mL*pL)/(mL*qL);
%     
%     if partials
%     dCDdM(i,1) = ( ((1:5).*mD(1:5))*pD(2:6) )/(mD*qD) - CD(i)/(mD*qD)*( ((1:5).*mD(1:5))*qD(2:6) );
%     dCLdM(i,1) = ( ((1:4).*mL(1:4))*pL(2:5) )/(mL*qL) - CL(i)/(mL*qL)*( ((1:4).*mL(1:4))*qL(2:5) );
%     end
% end

CD = 1.408;
CL = 0.357;
    if partials
    dCDdM = 0;
    dCLdM = 0;
    end