% rainSTORM_extras_simulate_3D_BW
%
% 2011 November 29  pseudo-3D_STORM_SIMULATION from Black and White object
%  Simulate some frames for astigmatic 3D storm...  
% See    C:\Documents and Settings\ejr36\My
% Documents\Projects\2011_3D_Anisotropic_STORM\2011_11_29_3Dstorm

% Eric Rees
% 2D work on ERF()-matrix method by Mark Deimund
%
% Simulate a stack of anisotropic (3D) localisation microscopy data.
% This simulation uses our usual dSTORM assumptions: temporally independent
% frames, and randomly-activated fluorophores that persist for a duration
% of exactly one frame.
% It also assumes elliptical Gaussian PSFs, after Holtzer-2007:
%
% Holtzer-2007 Gaussian beams [APPLIED PHYSICS LETTERS 90, 053902 2007]
%
% The units in this script (used implicitely) are nanometres
% 16 nm per true-object pixel, 160 nm per CCD pixel
% Beware if other scripts use microns (mine often do - EJR)
%
% Outline:
% 1. Import "crossed lines" as the "true object"
%    This script assumes a black (0) and white (>0) true object
%    Greyscale true objects can be simulated using small edits, if needed
% 2. Determine the "true" fluorophore positions in x-y
%    Assign some z-positions (first try a pseudo-2D inclined plane)
% 3. Simulate some noise-free images using a 2D Gaussian allocation
% 4. Add (optionally) Gaussian and/or Poisson noise
% 5. Save each image sequentially
% 
% ImageJ can then be used to compile the stack into a multipage .tif
%
% NOTES FOR SUCCESSFUL USE OF THIS SCRIPT
%
% Memo: in 2D we could use simple conv2, but for a Z-dep Gaussian, we'll
% need something more sophisticated. One approach that could work would be
% to allocate individual PSFs.
% 
% Avoid simulating objects whose psfs will overfill their input image size!
% This script will fail if it tries to allocate signal beyond the border o
% the CCD (which equals the true object area). Will try to fix later.
%
% New users will want to edit file output directory (do find: imwrite )


% 0. Set flags to control the algorithm
rng('shuffle')               % Start a random number sequence (see 'rng')

flag2D               = 1;    % True: all z-posits = 0. False: 'inclined' Z

flagGaussianNoise    = 1;    % Add Gaussian "camera noise" to ccdSignal
flagPoissonNoise     = 0;    % Add Poisson "background noise" to ccdSignal
flagQuantumPoisson   = 1;    % Convert expected photon number to Poisson variate. Simulates quantum behaviour


numberOfImages     = 1000;      % How many simulated images do we want?
% dyeFraction        = 0.0004;   % Probability a given dye is "on". Indep!! 
dyeFraction   =  0.1;
dyePhotons         = 1000;
noiseGaussianMean  = 100;
noiseGaussianSigma = 10;
noisePoissonMean   = 0;

simSig2D   = 160;      % (2D) Define PSF shape for simulated data...


% USE THESE PARAMETERS IN THE MODEL EQUATION (Holtzer 2007):
% sigXsqr = ((simSigZero/simZR)^2) *(yData+simGamma).^2 + simSigZero^2;
% sigYsqr = ((simSigZero/simZR)^2) *(yData-simGamma).^2 + simSigZero^2;

simGridSize  = [64, 64];       % size of a "Simulation Grid" space 
  pxSimGrid    = 160;            % nm per unit of "Simulation Grid"

  rescale  = 1;                  % rescale simulation grid to camera pixels
  pxCCD    = pxSimGrid*rescale;  % nm width of CCD pixels
  ccdGridSize = simGridSize./rescale;
  
  ccdEdgesRow = 0 : pxCCD : pxSimGrid*simGridSize(1);
  ccdEdgesCol = 0 : pxCCD : pxSimGrid*simGridSize(2);  

% 2. Determine (or define) the dye locations in [Row,Col], and Z
% Assume dyes stick in the centre of the hi-res pixels

%XXXXX

% 3. Simulate the images, using subsets of dyes, then get ccd
% images. Also need some way to visualise this as an image stack

figure(1)
colormap(gray)

tic % Start timing how long it takes to simulate a stack of images
for lpIm = 1:numberOfImages    
  
  dyeIsActive = rand(size(SimData,1));     % "Activate" a subset of dyes
  dyeIsActive = (dyeIsActive < dyeFraction);
  dyeIndexes  = find(dyeIsActive);
  
  imageCCD     = zeros( ccdGridSize ); % Setup empty camera measurement

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for lpIm = 1:numberOfImages    
  
  dyeChance   = rand(size(dyeObPosRow));  % Random number for each dye
 % Determines which fluoresce
 
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For each dye glowing in the current frame:
  for loopDyes = 1:length(dyeIndexes); % Add photons for each dye in frame

      if(flag2D)
       mySigRow = simSig2D; % Circular Gaussian PSF, for 2D
       mySigCol = simSig2D;
      else
       myZ = dyeObPosZ(dyeIndexes(loopDyes));
       sigXsqr = ((simSigZero/simZR)^2) *(myZ+simGamma).^2 + simSigZero^2;
       sigYsqr = ((simSigZero/simZR)^2) *(myZ-simGamma).^2 + simSigZero^2;
       mySigRow = sqrt(sigYsqr);
       mySigCol = sqrt(sigXsqr); % This determines the PSF ellipticity
      end
      
      % Allocate signal from this dye onto the CCD grid
      % Use a reasonable analytical approximation...
      % First find relative CCD grid distances
      % "-0.0 used to be -0.5, for registration..."
      myPosRow = (dyeObPosRow(dyeIndexes(loopDyes)) -0.0)*pxObj;
      myPosCol = (dyeObPosCol(dyeIndexes(loopDyes)) -0.0)*pxObj; % Dye position in nm
      myCcdEdgeRow = find(ccdEdgesRow == max(ccdEdgesRow(ccdEdgesRow < myPosRow)) );
      myCcdEdgeCol = find(ccdEdgesCol == max(ccdEdgesCol(ccdEdgesCol < myPosCol)) );
      relCcdEdgesRow = (myCcdEdgeRow-4:myCcdEdgeRow+5)*pxCCD - myPosRow;
      relCcdEdgesCol = (myCcdEdgeCol-4:myCcdEdgeCol+5)*pxCCD - myPosCol;
      
      % Determine amount of psf in each pixel using vectorised format
      thisDyeSignal = (10000).*...  % Try handling numbers near 1
            ( erf(relCcdEdgesRow(2:10)/(mySigRow)) - erf(relCcdEdgesRow(1:9)/(mySigRow)) )'*...
            ( erf(relCcdEdgesCol(2:10)/(mySigCol)) - erf(relCcdEdgesCol(1:9)/(mySigCol)) );
       
      % Allocate signal photons to each pixel in proportion to psf fraction
      totalPsfSignal = sum(thisDyeSignal(:));
      thisDyeSignal = thisDyeSignal./totalPsfSignal; % Normalise
      thisDyeSignal = thisDyeSignal*dyePhotons;      % Photons in ccd pixels
 
      thisDyeSignal = floor(thisDyeSignal);
      
      if(flagQuantumPoisson)
       % Thompson-2002 converts expected photon number to a Poisson variate
       % This simulates Quantum Mechanical photon statistics 
       % It worsens localisation precision (physically realistically)
       thisDyeSignal = uint16(thisDyeSignal);
       thisDyeSignal = imnoise(thisDyeSignal, 'poisson');
       thisDyeSignal = double(thisDyeSignal);
      end
      
      imageCCD(myCcdEdgeRow-4:myCcdEdgeRow+4, ...
               myCcdEdgeCol-4:myCcdEdgeCol+4) = ...
      imageCCD(myCcdEdgeRow-4:myCcdEdgeRow+4, ...
               myCcdEdgeCol-4:myCcdEdgeCol+4) + thisDyeSignal;
      
  end % Have now simulated image of all "active" dyes in this frame
  
% Now add camera (Gaussian) and shot (Poisson) noise to simulated data  
  if(flagGaussianNoise)
    % Note the use of "floor" to quantise integer numbers of photons
    W = ones(size(imageCCD))*noiseGaussianMean + floor(noiseGaussianSigma*randn(size(imageCCD)) );
    imageCCD = imageCCD + W;
  end
    
  if(flagPoissonNoise) 
    % Note the use of "uint16" to set integer scale for noise
    B = imnoise(uint16( noisePoissonMean*ones(size(imageCCD)) ),'poisson');
    imageCCD = uint16(imageCCD) + B;
  end
  
    
 imageCCD(imageCCD<0) = 0;        % Let negative noisy intensities = 0
  
  % Save image as a tif 
  if(flagSaveTifs)
  imageCCD = uint16(imageCCD); % Vista Preview dislikes 16-bit.. ImageJ OK
  tifName = ['C:\Documents and Settings\ejr36\My Documents\Work_Papers\2013JoOptics\Matlab\image',int2str(lpIm),'.tif'];
  % Write as tif. Note that "append" mode is N^2 slow; avoid if possible.
  imwrite(imageCCD,tifName,'tif'); %,'writemode', 'append'); 
  % Don't use 'writemode', 'append' - stacking in ImageJ is faster
  end
  
  if(flagShowCCDImages)
  figure(1)            % NOTE THAT DISPLAYING IMAGES ON SCREEN IS SLOW!
  imagesc(imageCCD)
  end
  
  waitbar(lpIm/numberOfImages)
end     % Have now simulated one CCD image
toc     % Have now simulated a stack of dSTORM images. Stack them in ImageJ
