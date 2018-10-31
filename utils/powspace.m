function v = powspace(x1, x2, n, p)
% POWSPACE Generate a powerlaw interpolation between two points.
%   v = POWSPACE(x1, x2, n, p) generates an `n`-step powerlaw interpolation
%   between `x1` and `x2`, using exponent `p`:
%       v(i) = x1 + (x2-x1) * ((i-1)/(n-1))^p
%
%   This reduces to Matlab's linspace when `p = 1`.
%
%   XXX But it would be nice if
%           powspace(a, b, n, p) == fliplr(powspace(b, a, n, p))
%       which isn't the case currently...
%
%   See also: LINSPACE.

r = 0:(1/(n-1)):1;
v = x1 + (x2 - x1)*r.^p;

end