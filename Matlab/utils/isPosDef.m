%ISPOSDEF Determines if a matrix is positive definite.
%   ISPOSDEF(MATRIX,SEMI) returns true if MATRIX is positive definite based
%   on its eigenvalues. If SEMI is true, the function checks for positive
%   semidefiniteness of MATRIX instead.
%   
%   This function can also be used to check for negative
%   (semi)definiteness very easily. E.g.,
%   isNegDef     = @(M) isPosDef(-M,true);
%   isNegSemiDef = @(M) isPosDef(-M);
function value = isPosDef(M,semi)
if nargin == 1 || isempty(semi)
    semi = false;
end
v = eig(M);
if semi
    tol = 0;
else
    tol = 1e-10; %An eigenvalue must be greater than or equal to this value to be considered nonzero
end
if all(v>=tol)
    value = true;
else
    value = false;
end



end