function c = divergent(c1, c2, n, cm, p)
% DIVERGENT Create a divergent color map, going from one color, to white in
% the middle, to another color.
%   c = DIVERGENT(c1, c2, n) creates an n-step color map going from color
%   `c1` at one end, to white in the middle, to `c2` at the other end. The
%   colors `c1` and `c2` must be three-component vectors.
%
%   c = DIVERGENT(c1, c2, n, cm) uses `cm` instead of white for the middle
%   color.
%
%   c = DIVERGENT(c1, c2, n, [], p)
%   c = DIVERGENT(c1, c2, n, cm, p)
%       use a nonlinear function with exponent `p` for interpolating
%       between the colors.
%
%   See also: REDBLUE.

% handle default
if nargin < 4 || isempty(cm)
    cm = [1 1 1];
end
if nargin < 5 || isempty(p)
    p = 1;
end

% normalize shapes of inputs
c1 = c1(:)';
c2 = c2(:)';

% prepare the output
c = zeros(n, 3);

% draw first half
n2 = ceil(n/2);
for i = 1:3
    c(1:n2, i) = powspace(c1(i), cm(i), n2, p);
end

% draw second half
for i = 1:3
    c(end-n2+1:end, i) = powspace(cm(i), c2(i), n2, p);
end

end