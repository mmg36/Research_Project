function params = rainSTORM_LocsPerFrame(params)
% Copyright 2012. Refer to 00_license.txt for details.

numberOfFrames  = params.localization.results.numberOfFrames;
reviewedSupResParams = params.reviewer.results.reviewedSupResParams;
SupResParams   = params.localization.results.SupResParams;

numberAcceptedPerFrame = zeros(1,numberOfFrames);
numberCandidatesPerFrame = zeros(1,numberOfFrames);

for lpLoc = 1:size(reviewedSupResParams,1)
  numberAcceptedPerFrame(reviewedSupResParams(lpLoc).frame_idx) = ...
      numberAcceptedPerFrame(reviewedSupResParams(lpLoc).frame_idx)+1;
end

for lpLoc = 1:size(SupResParams,1)
  numberCandidatesPerFrame(SupResParams(lpLoc).frame_idx) = ...
      numberCandidatesPerFrame(SupResParams(lpLoc).frame_idx)+1;
end

numberRejectedPerFrame = numberCandidatesPerFrame - numberAcceptedPerFrame;

numberAcceptedPerFrameSm = smooth(numberAcceptedPerFrame,100);
numberRejectedPerFrameSm = smooth(numberRejectedPerFrame,100);

plot(1:numberOfFrames, numberRejectedPerFrameSm, 'r')
xlabel('CCD frame number', 'fontSize', 12, 'fontWeight', 'bold');
ylabel('Number per Frame', 'fontSize', 12, 'fontWeight', 'bold');
title('Localisations per Frame (smoothed)', ... 
    'fontSize', 12,'fontWeight','bold');
set(gca, 'fontSize', 12, 'fontWeight', 'bold')
  hold on
  plot(1:numberOfFrames, numberAcceptedPerFrameSm,'b')
  legend('Rejected','Accepted');
  hold off

params.numberAcceptedPerFrame = numberAcceptedPerFrame;
params.numberRejectedPerFrame = numberRejectedPerFrame;

params.flags.LocsPerFrame = 1;
end