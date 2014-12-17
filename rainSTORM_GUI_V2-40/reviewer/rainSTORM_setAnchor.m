function params = rainSTORM_setAnchor(params)
% function flagMarkAnchor = rainSTORM_setAnchor(flagMarkAnchor)
% Identifies the mean position of a fiducial marker, 
% So that other marker positions (even in another stack of im data) can be 
% subtracted to give required drift corrections.

% markerPosits = params.reviewer.results.markerPosits;

params.reviewer.results.markerAnchor = mean(params.reviewer.results.markerPosits,1);

params.flags.MarkAnchor = 1; % Indicates a fiducial marker "anchor position" is set
end