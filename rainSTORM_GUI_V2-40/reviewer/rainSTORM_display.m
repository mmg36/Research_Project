function params = rainSTORM_display(params, SupResIm, linMag)
% rainSTORM_display
% Copyright 2012. Refer to 00_license.txt for details.
% This script displays a localisation density histogram in a MATLAB figure
% Optionally, a scalebar is plotted.
% To manually select a scalebar (or not), edit flagSB in the base workspace

% SupResIm is the Simple Histogram visualsation of localisation density
% linMag is the factor for sub-dividing pixel widths
%     linMag is the "reconstructionScaleFactor" for reviewed images
%     linMag is set by "prevSF" for the preview (5 is hardcoded - often OK)
if nargin == 1
    linMag = params.localization.settings.prevSF;
    SupResIm = params.reviewer.results.SupResIm;
end

flagSB = params.flags.SB;
pixelWidth = params.rawdata_mgr.myImInfo.pixelWidth;

scaleBarLn = 1000/pixelWidth; % Length for a 1 micron scalebar

% BRACKET THIS WITH IF OR CASE STATEMENTS, FOR ALTERNATIVE VISUALISATION
% For a Simple Histogram visualisation (as the on-screen figure)
figNewReconHandle = figure;
imshow(SupResIm, 'border', 'tight')
% For a Jittered Histogram, call rainSTORM_reconJH() and plot
% Need to thread linMag through as an argument - variable origins
%
     
hold on
caxis([min(SupResIm(:)) max(SupResIm(:))]);
colormap(hot); % Default colormap - can be changed in Reviewer
    
if(flagSB) 
  plot([max(xlim) (max(xlim)-scaleBarLn*linMag)]-10, ...
        [max(ylim)*0.9 max(ylim)*0.9],'w-','LineWidth',3);
  text((max(xlim)-1.0*scaleBarLn*linMag - 14), ...
      (max(ylim)*0.90 -18),'1 ?m','FontSize',12,'Color','w');
end

hold off

params.reviewer.settings.linMag=linMag;
params.reviewer.results.scaleBarLn=scaleBarLn;
params.reviewer.results.figNewReconHandle = figNewReconHandle;
% Function returns figNewReconHandle; this the visualisation figure number 
end