function params = rainSTORM_undrift(params)
% function flagUndrifted = rainSTORM_undrift(flagUndrifted)
% Copyright 2012. Refer to 00_license.txt for details.

% 0. Load required information from base workspace
flagSavedSupResData = params.flags.SavedSupResData;
markerAnchor = params.reviewer.results.markerAnchor;
markerFrames = params.reviewer.results.markerFrames; % List of frame #s with mark
% markerParams = params.markerParams;   % Presently not needed
markerPosits   = params.reviewer.results.markerPosits; 
reviewedParams = params.reviewer.results.reviewedSupResParams;
% reviewedPosits = params.reviewedPosits;


% 1. Save the uncorrected localisation data, if necessary
% To reload, run the next two lines (without commenting them out):
%    SupResPosits = SavedSupResPosits;
%    SupResParams = SavedSupResParams;
%    flagSavedSupResData is set=0 on running the rainSTORM.m Search
%    This Saved Data flag MAY also be set by the optical offset panel
if (flagSavedSupResData == 0) % If this is the first try at drift correction, 
    % Then save the uncorrected data, for user reference / reload
%     SavedSupResPosits = params.SupResPosits; 
%     SavedSupResParams = params.SupResParams;  
%     
%     params.SavedSupResPosits = SavedSupResPosits;
%     params.SavedSupResParams = SavedSupResParams;
    params.filtered_undrift = params.localization;
    params.flags.SavedSupResData = 1;%Flag unedited results exist
end

% 2. PROCESS: Evaluate required drift correction:
driftCorrection = ones(size(markerPosits,1),1)*markerAnchor - markerPosits;
% Note that a vector outer product generates a N-row by 2-Col matrix 
% 'mean(markerPosits,1)' gives the [Row,Col] mean marker postion, by taking
% the mean value down each row of data 


% 3. PROCESS: CORRECT DRIFT
% Use the identified fiducial marker data to correct drift
% This method uses a simple translation - affine TFORM needs real thought
% This method is elementwise; and ideally should be vectorised!
for lpLc = 1:size(reviewedParams,1) % For each quality-controlled position

    % Logical test: "Does a marker position come from this frame?"
   thisFrameMark = (markerFrames == reviewedParams(lpLc).frame_idx ); 
   
   if(sum(thisFrameMark) == 0 )     % If this frame has no marker position
       reviewedParams(lpLc).I = -1; %...And clear error-coded data at end
   else
       reviewedParams(lpLc).x = reviewedParams(lpLc).x + ...
           driftCorrection(thisFrameMark, 1); 
       reviewedParams(lpLc).y = reviewedParams(lpLc).y + ...
           driftCorrection(thisFrameMark, 2); 
   end
end

reviewedParams([reviewedParams.I] == -1) = [];



% 4. OUTPUT
% Write results to the base workspace,
% Copy reviewed data (now corrected) over SupResPosits and SupResParams, 
%... so that subsequent uses of the Reviewer refer to corrected data.

params.localization.results.SupResParams = reviewedParams;
% params.SupResPosits = reviewedPosits;
params.reviewer.results.reviewedSupResParams = reviewedParams;
% params.reviewedPosits = reviewedPosits;

params.flags.Undrifted = 1; % Return 1 to indicate SupResPosits have been adjusted
end