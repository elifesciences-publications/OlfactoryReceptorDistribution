function h = plotOlfactoryRF(responses)
% plotOlfactoryRF Plot a receptive field for an olfactory receptor, as in
% Hallem&Carlson 2006.
%   plotOlfactoryRF(responses) makes a symmetrized bar plot of the entries
%   in the `responses` vector, to depict an olfactory receptive field in a
%   manner similar to visual receptive fields.

sorted = sort(responses, 'descend');

locations = zeros(size(sorted));
locations(2:2:end) = 1:floor(length(locations) / 2);
locations(3:2:end) = -(1:floor((length(locations)-1) / 2));

h = bar(locations, sorted, 1, 'k');

end