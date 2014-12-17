% Fitting algorithm which determines the "centre of mass" of the ROI.
function [SupResParams] = rainSTORM_fitCoM(frameIdx, myFrame,myPixels,rad,bkgdSig)
% function [myFits,myParams] = rainSTORM_fitCoM(frameIdx, myFrame,myPixels,rad,bkgdSig)
% Calculate centroid for the ROI and use these coordinates as the
% localisation.

%     myFits = zeros(size(myPixels,1),2);   % [row col] matrix for each candidate
%     myParams = zeros(size(myPixels,1),7); % Parameters for accepted fits
    SupResParams = [];

    for lpPx = 1:size(myPixels,1); % For local maxima in descending order
        myRow = myPixels(lpPx,1);
        myCol = myPixels(lpPx,2);
        myROI = myFrame(myRow-rad:myRow+rad,myCol-rad:myCol+rad); 
        myROI = myROI - min(myROI(:));  % Square region to fit. Subtract minimum.
        % Note that myROI = ones(7) gives sigX = NaN: real data peaks obviate this 
        weight = (-rad:rad); % Basis. Let the ROI centre be (0,0).                

        % Calculate total intensity of the ROI.
        totalIntensity = sum(myROI(:));

        % Calculate x-value for centroid.
        xSums = sum(myROI,2);              % Determine the sum of each row.
        xValues = xSums' * weight';        % Weighted sum of rows
        xCoord = xValues / totalIntensity; % Calculate x-coordinate.
        sigX = sqrt( (xSums' * (weight.^2)')/totalIntensity - xCoord^2); % CHECK!

        % Calculate y-value for centroid.
        ySums = sum(myROI,1);      % Determine the sum of each column.
        yValues = ySums * weight'; % Weighted sum of columns
        yCoord = yValues / totalIntensity; % Centre of Mass y-coordinate.
        sigY = sqrt( (ySums * (weight.^2)')/totalIntensity - yCoord^2); % CHECK!

%         % Add localisation to myFits. Note -(0.5,0.5) for pixel registration.
%         myFits(lpPx,:) = [(double(myRow)+xCoord- 0.5),(double(myCol)+yCoord- 0.5)];
% 
%         myParams(lpPx,1)=myPixels(lpPx,3); % Copy signal counts to parameters 
%         myParams(lpPx,2)=0;                % Bodge tolerance as perfect: no fit!
%         myParams(lpPx,3)=totalIntensity;   % = Counts - (estimated) background
%         myParams(lpPx,4)=sigX;             % Sample Std dev in row-direction
%         myParams(lpPx,5)=sigY;             % Sample Std dev in col-direction
%         myParams(lpPx,6)=bkgdSig;          % Background estimate for this frame
%         % The 7th column, frame number, is applied in the MAIN script
        
        SupResParams(lpPx).frame_idx = frameIdx;
        SupResParams(lpPx).x = (double(myRow)+xCoord- 0.5);
        SupResParams(lpPx).y = (double(myCol)+yCoord- 0.5);
        SupResParams(lpPx).z = 0;
        SupResParams(lpPx).I = myPixels(lpPx,3); % Averaged magnitude of this signal
        SupResParams(lpPx).sig_x = sigX;  % X-width (sigma, rows, fitted) of this Gaussian
        SupResParams(lpPx).sig_y = sigY;  % Y-width (sigma, cols, fitted) of this Gaussian
        SupResParams(lpPx).avg_brigthness = bkgdSig; % Background for each ROI
        SupResParams(lpPx).res = 0; % Mean critical tol for fit
        SupResParams(lpPx).res_Row = 0;
        SupResParams(lpPx).res_Col = 0;
        SupResParams(lpPx).Sum_signal = totalIntensity; % Sum of signal (counts) for this fit

    end
end