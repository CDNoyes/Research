function output = EDLCost( input )
xf = input.phase(2).finalstate;
output.objective = -xf(7); % Maximized final mass 
end