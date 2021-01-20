function sol = DDPSolve(inp)


if inp.optimize_gains
    sol = entry_stochastic_gains(inp);
else
    sol = entry_stochastic(inp);    
end