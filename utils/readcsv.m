function data = readcsv(file, varargin)
% READCSV Read CSV file with arbitrary data types.
%   data = READCSV(file) reads a CSV file, treating all data as strings.
%
%   Options:
%    'delimiter'
%       Change the delimiter.
%       (default: ',')

% parse options
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('delimiter', ',', @(c) ischar(c) && (numel(c) == 1 || (numel(c) == 2 && c(1) == '\')));

parser.parse(varargin{:});
params = parser.Results;

% load raw data
f = fopen(file);
% split into lines
data_raw = textscan(f, '%s', 'Delimiter', '');
fclose(f);
data_raw = data_raw{1};

% split into columns
data0 = {};
for i = 1:length(data_raw)
    crt_line = strtrim(data_raw{i});
    if ~isempty(crt_line)
        crt_tokens = strsplit(crt_line, params.delimiter, 'CollapseDelimiters', false);
        if isempty(data0)
            data0 = cell(length(data_raw), length(crt_tokens));
        end
        data0(i, :) = crt_tokens; %#ok<AGROW>
    end
end
% skip lines that are completely empty
mask = arrayfun(@(i) ~all(cellfun(@isempty, data0(i, :))), 1:size(data0, 1));
data = data0(mask, :);

end