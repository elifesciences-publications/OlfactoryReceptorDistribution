function lum = getlum(c, mode)
% GETLUM Convert RGB colors to luminance.
%   lum = GETLUM(c) returns a luminance estimate for every row of the
%   3-column vector `c`, containing RGB colors. This uses the 'lumsrgb'
%   mode below.
%
%   lum = GETLUM(c, mode) uses particular ways for calculating the
%   luminance. Mode can be
%     * 'lumlin':    lum = 0.2126*R + 0.7152*G + 0.0722*B
%     * 'lumsrgb':   lum = 0.2126*r + 0.7152*g + 0.0722*b,
%                      where r, g, b are gamma corrected values,
%                    x = X/12.92 if x <= 0.03928 else ((X+0.055)/1.055)^2.4
%                      with x \in {r, g, b}, X \in {R, G, B}
%     * 'lumapprox': lum = sqrt(0.299*R^2 + 0.587*G^2 + 0.114*B^2).

if size(c, 2) ~= 3
    error([mfilename ':badsize'], 'The color argument should have 3 columns.');
end

if nargin < 2
    mode = 'lumsrgb';
end

R = c(:, 1);
G = c(:, 2);
B = c(:, 3);

switch mode
    case 'lumlin'
        lum = 0.2126*R + 0.7152*G + 0.0722*B;
    case 'lumsrgb'
        cr = c;
        mask = (cr <= 0.03928);
        cr(mask) = cr(mask)/12.92;
        cr(~mask) = ((cr(~mask) + 0.055)/1.055).^2.4;
        lum = 0.2126*cr(:, 1) + 0.7152*cr(:, 2) + 0.0722*cr(:, 3);
    case 'lumapprox'
        lum = sqrt(0.299*R.^2 + 0.587*G.^2 + 0.114*B.^2);
end

end