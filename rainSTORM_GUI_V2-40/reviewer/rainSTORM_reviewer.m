% rainSTORM_reviewer 
%   Copyright 2012. Refer to 00_license.txt for details.
%   Eric Rees, and GUI developed initially by Mark Deimund
%
% FUNCTION
%   A function to apply Quality Control to Localisation Microscopy data
%   And then create the required Visualisation of the data
%   This function can be run via the "Reviewer" GUI
%
% USING THIS SCRIPT:
%   Use this function after running rainSTORM_main on some data,
%   Which may produce ~10000+ SupResPosits and corresponding SupResParams.
%
% NOTES
%   This script reviews Localisation data, then calls an image
%     reconstruction
%   Reviewing parameters permits tightening of acceptance parameters
%     (such as threshold, tolerance, allowSig).
%   This obviates re-running the (5 min) image analysis algorithm
%   But this script only reviews fitted data, 
%     so only higher Threshold       (new Thresh)
%     or smaller tolerance           (newTol)
%     or a narrower range of widths  (newSig)
%     will have any effect on the SupResPosits that are accepted.
%
%   A new sub-pixel resolution is defined for the visualisation  (linMag)
%   A scale bar length (in CCD pixels) is defined (scaleBarLn)
%
% FURTHER READING
%   Rees et al. Optical Nanoscopy 1:12
%
% REQUIRED INPUTS - 
% these must be available:
% 
% SupResPosits  % [Rows, Cols] of localisations, with sub-pixel resolution
% SupResParams  % [Threshold, Tolerance, Total counts, SigX] 
% scaleBarLn    % e.g. scaleBarLn = 8.33; % For 2 micro-m, at 240 nm/Px
% linMag        % e.g. linMag = 10; % But we could redefine this here
% size(myFrame) % e.g. [64 64]

function params =rainSTORM_reviewer(params)

newThresh = params.reviewer.settings.filter_settings.newThresh;
newTol = params.reviewer.settings.filter_settings.newTol;
newSig = params.reviewer.settings.filter_settings.newSigma;
newPrecision = params.reviewer.settings.filter_settings.newPrecision;
% SupResPosits = params.SupResPosits;
frame_size = params.rawdata_mgr.myImInfo.frame_size;
flagSB = params.flags.SB;
newFrames = params.reviewer.settings.filter_settings.newFrames;
reconstructionScaleFactor = params.reviewer.settings.linMag;
algVisual = params.reviewer.settings.algVisual;

% newThresh = New threshold 
% newTol = New tolerance (*almost* identical to _fitLM method)
% newSig = New acceptable range of widths for Gaussian fits
% flagSB = If true, plot scalebar


% Compute precision in fitted positions, in [row-direction, col-direction]
[params,nPhotons] = rainSTORM_precision(params);
SupResParams = params.localization.results.SupResParams;
% SupResDeltaX = params.SupResDeltaX;
% Use the mean Thompson precision for 2D quality control
% For 3D, let us change this and apply the precision limit to both axes
% deltaX = mean(SupResDeltaX,2); 

% Review localisation parameters to choose good localisations
% 
% SupResParams(:,1): Is a 3x3 spot brighter than a threshold? AND
% SupResParams(:,2): Is the residual of the Gaussian fit small enough? AND
% SupResParams(:,4): Is the Gaussian Std Dev in the row direction in range,
% SupResParams(:,5): Is the Gaussian Std Dev in the col direction in range,
% SupResDeltaX(:,1): Is the Thompson Precision estimate precise enough- row
% SupResDeltaX(:,2): Is the Thompson Precision estimate precise enough- col
% SupResParams(:,7): Is the localisation from a wanted CCD frame number?
% 

qualityApprovedRows = ( ([SupResParams.I] > newThresh) &...
                        ([SupResParams.res] < newTol) &...
                        ([SupResParams.sig_x] > newSig(1) ) &...
                        ([SupResParams.sig_x] < newSig(2) ) &...
                        ([SupResParams.sig_y] > newSig(1) ) &...
                        ([SupResParams.sig_y] < newSig(2) ) &...
                        ([SupResParams.x_std] < newPrecision) &... 
                        ([SupResParams.y_std] < newPrecision) &...
                        ([SupResParams.frame_idx] >= newFrames(1) ) &...
                        ([SupResParams.frame_idx] <= newFrames(2) )  ... 
                       );

% reviewedPosits = SupResPosits( qualityApprovedRows,: ); % All Columns
reviewedParams = SupResParams( qualityApprovedRows );
% reviewedDeltaX = SupResDeltaX( qualityApprovedRows,: );

% reviewedPhotonNums = nPhotons( qualityApprovedRows,: );


% params.reviewedPosits = reviewedPosits;
params.reviewer.results.reviewedSupResParams = reviewedParams;
% params.reviewedDeltaX = reviewedDeltaX;
% params.reviewedPhotonNums = reviewedPhotonNums;
 reviewedDeltaX = [[reviewedParams.x_std]' [reviewedParams.y_std]'];
reviewedPosits = [[reviewedParams.x]' [reviewedParams.y]'];

meanRevDeltaX = mean(reviewedDeltaX); % Mean precision in reconstruction
stdRevDeltaX  = std(reviewedDeltaX);  % Std Dev of precisions- poor metric?
SparrowThompsonLimit = 2*sqrt(mean(reviewedDeltaX.^2)) % Apprx 'resolution'

% Reconstruct an image using the reviewed localisations
% The following line creates a "Simple Histogram Image"
[SupResIm] = rainSTORM_recon(reviewedPosits,reviewedParams, ...
               reconstructionScaleFactor,frame_size); 

params.reviewer.results.SupResIm = SupResIm;
% Now either plot the "Simple Histogram Image"
% Or create and plot the "Jittered Histogram Image" 
if(algVisual == 1)         
    params = rainSTORM_display(params);
elseif(algVisual == 2)
  params = rainSTORM_recon_JH(params);
%   reconstructionScaleFactor = params.reviewer.results.jhLinMag; % Keep recent value updated in base
end
   
numberOfFits = size(SupResParams,1);
densestPoint = max(SupResIm(:));

params.reviewer.results.SparrowThompsonLimit = SparrowThompsonLimit;
% params.SupResDeltaX = SupResDeltaX; % Check this works
params.reviewer.results.densestPoint = densestPoint;
% params.figNewReconHandle = figNewReconHandle;
params.reviewer.results.meanRevDeltaX = meanRevDeltaX;
% params.newThresh = newThresh;
% params.newTol = newTol;
% params.newSigma = newSig;
% params.newPrecision = newPrecision;
% params.newFrames = newFrames;
params.localization.numberOfLocs = numberOfFits;
params.reviewer.results.qualityApprovedLocs = qualityApprovedRows;
% params.reconstructionScaleFactor = reconstructionScaleFactor;
params.reviewer.results.stdRevDeltaX = stdRevDeltaX;

end