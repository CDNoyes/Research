%MATRIXINTERP Interpolates a series of matrices.
%   MATRIXINTERP(X,M,XI) Interpolates the matrices M given at the values X
%   at the new point(s) XI. If X has length n, then M should be a cell
%   array with n cells, each containing an m-by-p matrix. The result MI is
%   a cell array of m-by-p matrices when the same number of cells as the
%   length of XI.

function Mi = MatrixInterp(x,M,xi)

if isa(M,'cell')
    M = cat(3,M{:});   
end
s = size(M);

%Method one: griddedInterpolant

%Method two: reshaping
V = reshape(M,[],2).';
VI = interp1(x,V,xi)';
if ~isscalar(xi)
    Mi =  squeeze(mat2cell(reshape(VI,[s(1:2),length(xi)]),s(1),s(2),ones(1,length(xi))));
else
    Mi = reshape(VI,[s(1:2),length(xi)]);
end

end