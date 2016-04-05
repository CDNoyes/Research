function output = EntryCost( input )
ref = input.auxdata.ref;
xf = input.phase.finalstate;
% tf = input.phase.finaltime;

[DR,CR] = Range(ref.state(1,2),ref.state(1,3),ref.state(1,6),xf(2),xf(3));


output.objective =  CR.^2 + 0*input.phase.integral; %(DR-780).^2 ; % Final CR + Drag Tracking

end