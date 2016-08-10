function y = ndspace(d1,d2,n,f)

% NDSPACE N-dimensionally spaced points
%   NDSPACE(X1, X2) generates a column matrix of 10^n linearly
%   equally spaced points in the hypercube defined by the two diametrically
%   opposite cornerpoints X1 and X2 both n-by-1 vectors. If either is a
%   scalar it will be expanded to repmat(X,n,1).
%
%   NDSPACE(X1, X2, P) generates P^n points if P is a scalar or prod(P)
%   points if P is a n-by-1 vector.
%
%   NDSPACE(X1, X2, P, F) generates points using a function defined by the
%   function handle or string F which needs to have the calling syntax of
%   LINSPACE(X1, X2, N). If F is a n-by-1 cell of function handles or 
%   strings different generating functions can be used per dimension.
%
%   Example
%
%     x = ndspace([0,0],[10,100],5,@linspace)
%     y = ndspace([0,1],[1,2],6,@logspace)
%     z = ndspace([10,1],[20,2],[5,8],{@linspace,@logspace})
%
%     plot(x(:,1),x(:,2),'bo'), hold on
%     plot(y(:,1),y(:,2),'gx')
%     plot(z(:,1),z(:,2),'r*'), hold off
%
% Author: Christophe Lauwerys
% Release date: 12/08/2011

error(nargchk(2,4,nargin))
assert(isvector(d1),'D1 must be a vector')
assert(isvector(d2),'D2 must be a vector')

if nargin<3 || isempty(n)
    n = 10;
else
    assert(isvector(n),'N must be a vector')
end

if nargin<4 || isempty(f)
    f = {@linspace};
elseif iscell(f)
    assert(all(cellfun(@(x)ischar(x) || isa(x,'function_handle'),f)),'F must be a cell of strings and function handles')
else
    assert(ischar(f) || isa(f,'function_handle'),'F must be a string or a function handle')
    f = {f};
end

nd = max([numel(d1),numel(d2),numel(n),numel(f)]); % Number of input dimensions

assert(numel(d1)==nd || isscalar(d1),'D1 must be a scalar or a vector of compatible length with D2, N and F')
assert(numel(d2)==nd || isscalar(d2),'D2 must be a scalar or a vector of compatible length with D1, N and F')
assert(numel(n)==nd || isscalar(n),'N must be a scalar or a vector of compatible length with D1, D2 and F')
assert(numel(f)==nd || isscalar(f),'F must be a scalar or a vector of compatible length with D1, D2 and N')

d1 = expand(num2cell(d1),nd);
d2 = expand(num2cell(d2),nd);
n = expand(num2cell(n),nd);
f = expand(f,nd);

y = cellfun(@feval,f,d1,d2,n,'Uni',0);
if nd>1
    [y{:}]=ndgrid(y{:});
end
y = cellfun(@(x)x(:),y,'Uni',0);
y = [y{:}];

function x = expand(x,nd)

if isscalar(x)
    x=repmat(x,nd,1);
else
    x=x(:);
end