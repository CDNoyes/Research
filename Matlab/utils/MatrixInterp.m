%MATRIXINTERP Interpolates a series of matrices.
%   MATRIXINTERP(X,M,XI) Interpolates the matrices M given at the values X
%   at the new point(s) XI. If X has length n, then M should be a cell
%   array with n cells, each containing an m-by-p matrix. The result MI is
%   a cell array of m-by-p matrices when the same number of cells as the
%   length of XI.

function Mi = MatrixInterp(x,M,xi,method)

if nargin < 4 || isempty(method)
    method = 'linear';
end

if isa(M,'cell')
    M = cat(3,M{:});
end
s = size(M);

%The two methods are actually very competitive.
%Method two: reshaping
% M = squeeze(M);
% if length(s) ~= length(size(M))
%     M = M';
% end

% V = reshape(M,[],s(1)*s(2));
% VI = interp1(x,V,xi,method);

%This method works but is somewhat slow because of the loop.
VI = zeros([s(1:2),length(xi)]);
for i = 1:s(1)
    VI(i,1:s(2),:) = interp1(x,squeeze(M(i,:,:)).',xi,method).';
end

if ~isscalar(xi)
%     VI = reshape(VI',[s(1:2),length(xi)]);
    Mi =  squeeze(mat2cell(VI,s(1),s(2),ones(1,length(xi))));
else
    Mi = reshape(VI,[s(1:2),length(xi)]);
end

end