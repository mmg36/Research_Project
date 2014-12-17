function params=rainSTORM_main(params)

    global parCompTool

    filename = params.rawdata_mgr.filename;
    ext = params.rawdata_mgr.ext;
    initX0 = params.localization.settings.initX0;
    initSig = params.localization.settings.initSig;
    rad = params.localization.settings.rad;
    tol = params.localization.settings.tol;
    Thresh = params.localization.settings.Thresh;
    maxIts = params.localization.settings.maxIts;
    alg = params.localization.algo_id;
    allowSig = params.localization.settings.allowSig; 
    allowX = params.localization.settings.allowX; 
    prevSF = params.localization.settings.prevSF;

    flagSum = params.flags.Sum;
    
    tic % See % http://www.matlabtips.com/waiting-for-the-waitbar/

    image = myImage(params.rawdata_mgr, parCompTool);
    frameSize = image.getFrameSize();
    numberOfImages = image.getNumberOfImages();
    
    SupResParams = [];
%     SupResParams = params.localization.results.SupResParams;
   
if (parCompTool)
%     SupResPosits = [];
%      SupResParams = [];
%     SupResNumber = 0;
    
    for lpIm=1:numberOfImages;
        SRP = rainSTORM_main_loop(lpIm, image.getFrame(lpIm),initX0, initSig, rad, tol, Thresh, maxIts, alg, allowSig, allowX);
        %[myFits, myParams] = rainSTORM_main_loop(lpIm, image.getFrame(lpIm),initX0, initSig, rad, tol, Thresh, maxIts, alg, allowSig, allowX);
        % Add frame number to myParams.
%         myParams(:,7) = lpIm;
        % Store super-resolution positions, and corresponding accuracy indicators
%         SupResPosits = [SupResPosits; myFits];
        SupResParams = [SupResParams, SRP];
%         SupResNumber = SupResNumber + size(myFits,1); % Track write row
        
    end
else
%     SupResPosits = zeros(numberOfImages*params.estNum,2); % To hold [row,col] of fits
%     SupResParams = zeros(numberOfImages*params.estNum,7); % And [Str,Res,Int,Sigma]
%     SupResNumber = 1;
    
    for lpIm=1:numberOfImages;
        SRP = rainSTORM_main_loop(lpIm, image.getFrame(lpIm),initX0, initSig, rad, tol, Thresh, maxIts, alg, allowSig, allowX);
%         [myFits, myParams] = rainSTORM_main_loop(lpIm, image.getFrame(lpIm),initX0, initSig, rad, tol, Thresh, maxIts, alg, allowSig, allowX);
        % Add frame number to myParams.
%         myParams(:,7) = lpIm;
%         % Store super-resolution positions, and corresponding accuracy indicators
%         SupResPosits(SupResNumber:SupResNumber + size(myFits,1)-1,:) = myFits;
%         SupResParams(SupResNumber:SupResNumber + size(myFits,1)-1,:) = myParams;
%         SupResNumber = SupResNumber + size(myFits,1); % Track write row
        SupResParams = [SupResParams, SRP];
        
    end
end

% close(progress)
% After image analysis, close files.

toc

% Remove any empty placeholders from the list of localised positions
% Remove empty first row from SupResPosits and SupResParams (by starting at row 2).
% SupResPosits = SupResPosits(1:SupResNumber-1,:);
% SupResParams = SupResParams(1:SupResNumber-1,:);


% 4. Reconstruct localised positions into a super-resolution image
%    Inputs are [positions, pixel-densification factor, [nRows nCols]]
%    SupResIm is the directly binned image.
%    Im2 is SupResIm with an image closure applied.
% To display a reconstuction immediately, without running the reviewer, 
% uncomment the following:
    
       
% Find sum of fluorescence images if desired.
if(flagSum) 
  sumFrame = uint32( zeros(frameSize) );
  sumProgress = waitbar(0,'Generating sum image');
    
  for lpIm = 1:numberOfImages;
      myFrame = uint32(image.getFrame(lpIm));
      sumFrame = sumFrame + myFrame;
      waitbar(lpIm/numberOfImages)
  end
 
  close(sumProgress)
  params.localization.results.sumFrame = sumFrame;
    
  figSum = figure;
  % Scale Sum image display to match preview reconstruction - uses % size
  imshow(sumFrame,'border','tight','InitialMagnification',prevSF*100)
  hold on
  caxis([min(sumFrame(:)) max(sumFrame(:))]);
  colormap(gray)    
  hold off
end

% params.localization.results.numberOfLocs = SupResNumber;
params.localization.results.SupResParams = SupResParams;
params.localization.results.numberOfFrames = numberOfImages; % (Means number of frames)
params.rawdata_mgr.myImInfo.frame_size = frameSize; % Size of CCD image

end