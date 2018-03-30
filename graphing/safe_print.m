function safe_print(fname, type, varargin)
% SAFE_PRINT Print image to file, avoiding overwrite.
%   SAFE_PRINT(fname) saves the current figure to a PDF.
%
%   SAFE_PRINT(fname, 'png') saves the current figure to a PNG file. If no
%   additional arguments are passed (see below), this includes the '-r300'
%   option to use 300 DPI.
%
%   Additional arguments are passed to Matlab's PRINT command.
%
%   The extension is added automatically if not provided (which is
%   recommended).

if nargin < 2
    type = 'pdf';
end

if strcmp(type, 'png') && isempty(varargin)
    varargin = {'-r300'};
end

[~, ~, ext] = fileparts(fname);
if isempty(ext)
    fname = [fname '.' type];
end

if exist(fname, 'file')
    warning([mfilename ':fexist'], 'File already exists.');
    return;
end

switch type
    case 'pdf'
        print('-dpdf', varargin{:}, fname);
    case 'png'
        print('-dpng', varargin{:}, fname);
    otherwise
        error([mfilename ':badtype'], 'Unrecognized file type.');
end

end