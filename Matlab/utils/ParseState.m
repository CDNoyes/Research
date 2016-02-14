%PARSESTATE splits a generic state vector x into its component states.
%   PARSESTATE(STATE) returns min{outputs, states in STATE} states taken
%   columnwise from STATE. Input state may be a vector or matrix.
function varargout = ParseState(x)

if iscolumn(x)
    x = x';
end
n = size(x,2);

for i = 1:n
    varargout{i} = x(:,i);
end
end
