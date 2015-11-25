function [dim,ind] = Dimension(nState,nControl,nParameter)

dim.state = nState;
dim.control = nControl;
if nargin == 3 && ~isempty(nParameter)
    dim.parameter = nParameter;
else
    dim.parameter = 0;
end

dim.total = dim.state+dim.control+dim.parameter;

ind = Indices(dim);

end

function ind = Indices(dim)

ind.state = 1:dim.state;
ind.control = dim.state + 1:dim.control;
ind.parameter = dim.state + dim.control + 1:dim.parameter;
ind.total = 1:dim.total;
end