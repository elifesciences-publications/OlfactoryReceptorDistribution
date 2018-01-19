function x = ipr(v)
% IPR Calculate inverse participation ratio.
%   x = IPR(v) calculates the inverse participation ratio for the vector v.
%   This is defined as sum(v)^2 / sum(v.^2).
%
%   When applied to a matrix, the IPR is calculated for every row.

x = sum(v, 2).^2 ./ sum(v .^ 2, 2);

end