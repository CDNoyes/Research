function Submatrix(M,ind)

%Should we deduce whether the matrix is a jacobian or hessian and whether
%it includes state/control/parameters?

d = size(M);
isHessian = (length(d) == 3 || isa(M,'cell'));

if isHessian % Like the second output of Hessian.m

else
    
end


end