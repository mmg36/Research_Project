function params = rainSTORM_DeleteBoxed(params)
% function flagSomeCuts = rainSTORM_DeleteBoxed(flagSomeCuts)
% Copyright 2012. Refer to 00_license.txt for details.
% 0. 
% Load required inputs from base workspace
boxCols = params.reviewer.settings.boxtrack_params.boxCols;
boxRows = params.reviewer.settings.boxtrack_params.boxRows;
flagSavedSupResData = params.flags.SavedSupResData;
reviewedSupResParams = params.reviewer.results.reviewedSupResParams;
% reviewedPosits = params.reviewedPosits;


% 1. 
% Save unedited raw localisation results if they aren't yet saved
if (flagSavedSupResData == 0) % If unedited results aren't yet saved
                              % Then save them:    
                              
    % TODO: filtered_<algo_name>-be elmenteni 
%     SavedSupResPosits = params.SupResPosits;
%     SavedSupResParams = params.SupResParams;  
%     params.SavedSupResPosits = SavedSupResPosits;
%     params.SavedSupResParams = SavedSupResParams;

    params.filtered_deleteboxed = params.localization;
    params.flags.SavedSupResData = 1;%Flag unedited results exist
end


% 2. 
% Remove boxed positions from the reviewed data 
% This can prevent spuriously accurate blobs appearing in reconstruction
% reviewedParams(reviewedPosits(:,1) > boxRows(1) & ...
%                reviewedPosits(:,1) < boxRows(2) & ...
%                reviewedPosits(:,2) > boxCols(1) & ...
%                reviewedPosits(:,2) < boxCols(2) ...
%                , :) = []; % CLEAR PARAMS FIRST - NEED POSIT DATA TO DO IT!
% 
% reviewedPosits(reviewedPosits(:,1) > boxRows(1) & ...
%                reviewedPosits(:,1) < boxRows(2) & ...
%                reviewedPosits(:,2) > boxCols(1) & ...
%                reviewedPosits(:,2) < boxCols(2) ...
%                , :) = []; % CLEAR POSITIONS LAST
reviewedSupResParams(...
               [reviewedSupResParams.x] > boxRows(1) & ...
               [reviewedSupResParams.x] < boxRows(2) & ...
               [reviewedSupResParams.y] > boxCols(1) & ...
               [reviewedSupResParams.y] < boxColss(2) ...
               )=[];
% 3. 
% Write results to the base workspace,
% Copy reviewed data (with empty box) over SupResPosits and SupResParams, 
%... so that subsequent uses of the Reviewer refer to corrected data.
params.localization.results.SupResParams = reviewedSupResParams;
% params.SupResPosits = reviewedPosits;
params.reviewer.results.reviewedSupResParams = reviewedSupResParams;
% params.reviewedPosits = reviewedPosits;
           
% SomeCuts flag visszat?r?si ?rt?ke?
end