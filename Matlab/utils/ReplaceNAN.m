function value = ReplaceNAN(inputValues, replace)

value = inputValues;
value(isnan(value)) = replace;

end