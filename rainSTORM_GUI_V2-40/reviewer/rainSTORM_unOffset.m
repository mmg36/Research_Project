function params = rainSTORM_unOffset( params )
% function flagOpoffCorrected = rainSTORM_unOffset( ~ )
% rainSTORM_unOffset
% Copyright 2012. Refer to 00_license.txt for details.
%   Corrects optical offset of localised CH2 data
% 

flagSavedSupResData = params.flags.SavedSupResData;
% reviewedDeltaX = params.reviewedDeltaX;
reviewedParams = params.reviewer.results.reviewedSupResParams;
reviewedPosits = [[reviewedParams.x]' [reviewedParams.y]'];
optOffTFORM    = params.reviewer.results.optOffTFORM;

% Copied from _UNDRIFT, save data for reversion if necessary
if (flagSavedSupResData == 0) % 
%     % Then save the uncorrected data, for user reference / reload

%     SavedSupResPosits = params.SupResPosits; 
%     SavedSupResParams = params.SupResParams;  
%     
%     params.SavedSupResPosits = SavedSupResPosits;
%     params.SavedSupResParams = SavedSupResParams;

    params.filtered_offset_correction = params.localization;
    params.flags.SavedSupResData = 1;%Flag unedited results exist
end

% SHOULD THIS BE FORWARDS, OR BACKWARDS? Am trying Inverse.
reviewedPosits = tforminv(optOffTFORM, reviewedPosits);

X = num2cell(reviewedPosits(:,1));
Y = num2cell(reviewedPosits(:,2));
[reviewedParams.x] = deal(X{:});
[reviewedParams.y] = deal(Y{:});
% Write the offset-corrected positions to the base workspace
% For both immediate analysis (rev...) and image reconstruction (Sup...)

params.localization.results.SupResParams = reviewedParams;
% params.SupResPosits = reviewedPosits;
params.reviewer.results.reviewedSupResParams = reviewedParams;
% params.reviewedPosits = reviewedPosits;

params.flags.OpOffCorrected = 1;

end

