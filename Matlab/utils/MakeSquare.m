% MAKESQUARE Pads a non-square matrix with zeros.
%     MAKESQUARE(M) Returns a square version of M such that the smaller of
%     its two dimensions is padded with zeros. For square M, MAKESQUARE
%     returns M.

function Ms = MakeSquare(M)

s = size(M);
d = abs(diff(s));
[~,dim] = min(s);
s(dim) = d;                 % Change size to necessary size of zeros to append
Ms = cat(dim,M,zeros(s));

end