%rainSTORM_2BPolarizationAnalysis
%This script uses the main RainSTORM Software to analyse 2 stacks of x-y
%polarized images

% Made by Daan van Kleef & Mehdi Goudarzi

% Only change to false if analysing multiple angles for plotting
flagInitialRun = false;

if flagInitialRun 
    clear
    flagInitialRun = true;
end

%Real angle of polarization
Phireal = 45;

% Flagsetzero sets missing observations to count 10^-10
% Flagsetbackground sets missing observations to background level
% Flagbackgroundsubtraction subtracts background level from observed photon
% counts and puts missing observations to 0
% NOTE, do not use backgroundsubtraction and setbackground together!

Flagsetzero = 0;
Flagsetbackground = 0;
Flagbackgroundsubtraction = 0;
FlagrefitMin = 0;
FlagrefitAverage= 1;

if sum(Flagsetzero+Flagsetbackground+Flagbackgroundsubtraction+FlagrefitMin+FlagrefitAverage) > 1
    error('Make sure there is only one active flag!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Use rainSTORM to analyse raw image data
%    Then read the data from the Matlab struct into a useful array...
startup  % Run rainSTORM once (manually), on the X-polarisation data
uiwait 

% Read in as an array!
Nxtable = struct2table((params.localization.results.SupResParams));

for count = 1:size(Nxtable,2)
    Nxdummy = double(table2array(Nxtable(:,count)));
    if count == 1
        Nxarray = Nxdummy;
    else
        Nxarray = cat(2,Nxarray,Nxdummy);
    end
end

% Read the path of the stored stack. 
XPath = params.rawdata_mgr.filename;


startup  % Run rainSTORM once (manually), on the Y-polarisation data
uiwait

% Read in as an array!
Nytable = struct2table((params.localization.results.SupResParams));
for count = 1:size(Nytable,2)
    Nydummy = double(table2array(Nytable(:,count)));
    if count == 1
        Nyarray = Nydummy;
    else
        Nyarray = cat(2,Nyarray,Nydummy);
    end
end

% Read the path of the stored stack. 
YPath = params.rawdata_mgr.filename;


% 2. Identify matching images
%    Method 1: assume fluorophore position is known (from simulation
%    params)

% Method 1:
% Remove spurious localizations from Nxarray and place in Nxarraycorrect
% Place rejected localizations into Nxreject
NxarrayCorrect = Nxarray;
NxDistance     = sqrt( (Nxarray(:,2)-32).^2 + (Nxarray(:,3)-30).^2 );
Nxreject = NxarrayCorrect(NxDistance > 1, :);
NxarrayCorrect(NxDistance > 1, :) = [];

%Idem for Nyarray
NyarrayCorrect = Nyarray;
NyDistance     = sqrt( (Nyarray(:,2)-32).^2 + (Nyarray(:,3)-30).^2 );
Nyreject = NyarrayCorrect(NyDistance > 1, :);
NyarrayCorrect(NyDistance > 1, :) = [];

numberOfFrames = params.localization.results.numberOfFrames;

%Initiate matrix to hold combined count list
% ComCount will hold the combined count list.
% Column Structure:
    %  Frame ID  | Nx | Ny | Nx present? | Ny present? | Estimated angle
    % Last 2 column will hold flags to indicate succesful counts
        % 0 = absent, 1 = present
        % Can be used to relocalize missing observations later
ComCount = zeros(numberOfFrames,6);

%Initiate frame ID list
for count = 1:numberOfFrames
    ComCount(count,1) = count;
end

%Copy data from Nxarraycorrect into ComCount
for count = 1:size(NxarrayCorrect)
    FrameID = NxarrayCorrect(count,1);
    ComCount(FrameID,2) = NxarrayCorrect(count, 12); % Sum-Signal used, rather than I.
    ComCount(FrameID,4) = 1; 
end

%Copy data from Nyarraycorrect into Comcount
for count = 1:size(NyarrayCorrect)
    FrameID = NyarrayCorrect(count,1);
    ComCount(FrameID,3) = NyarrayCorrect(count, 12); % Sum-Signal used, rather than I.
    ComCount(FrameID,5) = 1; 
end

%2b. (optional) Determine average residual background level from rejected
%localisations

ResBackground = 0.5*(mean(Nxreject(:,5)) + mean(Nyreject(:,5)));
% What is there is no rejected values or very few
% background = 49*std(pixelBrightnesses < mean (pixelBrightnesses) ????


% 3. Complete missing paired observations

% Method 1: Assume missing observations are count 10^-10
if(Flagsetzero == 1);

for count = 1:numberOfFrames
    if ComCount(count, 4) == 0
        ComCount (count , 2) = 10^(-10);
    end
    if ComCount(count,5) == 0
        ComCount(count,3) = 10^(-10);
    end
end
end

% Method 2: Assume missing observations are background count level
if(Flagsetbackground == 1);
for count = 1:numberOfFrames
    if ComCount(count, 4) == 0
        ComCount (count , 2) = ResBackground;
    end
    if ComCount(count,5) == 0
        ComCount(count,3) = ResBackground;
    end
end
end

% Method 3. Subtract average residual background level from measured
%accepted intensities and put missing observations to 0
if Flagbackgroundsubtraction == 1
    ComCount(:,2:3) = ComCount(:,2:3)-ResBackground;
    for count = 1:numberOfFrames
    if ComCount(count, 4) == 0
        ComCount (count , 2) = 10^(-10);
    end
    if ComCount(count,5) == 0
        ComCount(count,3) = 10^(-10);
    end
end
end

% Method 4: Re-inspect camera data for missing observations
% (a) Find out which frames of which data file need to be inspected
% (b) Process these and find region of interest
% (c) Then use minimum of ROI as background noise value
% (d) Use sum of ROI as missing observation
if FlagrefitMin == 1
    for count = 1:size(ComCount, 1)
        if ComCount(count, 4) == 0 
               myFrame = double(imread([XPath,'.tif'], count));
               myROI = myFrame(29:35,27:33); % For a single, simulated spot.
               myROI = myROI - min(min(myROI)); % Existing, crude background subtraction. 
               sum_signal = sum(myROI(:));
               % Then pass sum_signal this into the relevant place (ComCount(number, X=2, Y=3) )
               ComCount(count, 2) = sum_signal;
               ComCount (count, 4) = 2; % To check later if the values were moderated 
               
        end
        if ComCount(count, 5) == 0 
             myFrame = double(imread([YPath,'.tif'], count));
             myROI = myFrame(29:35,27:33); % For a single, simulated spot.
             myROI = myROI - min(min(myROI)); % Existing, crude background subtraction. 
             sum_signal = sum(myROI(:)); 
             % Then pass sum_signal this into the relevant place (ComCount(number, X=2, Y=3) ) 
             ComCount(count, 3) = sum_signal;
             ComCount (count, 5) = 2; % To check later if the values were moderated 
        end
    end 
end

% Method 5: Re-inspect camera data for missing observations
% (a) Find out which frames of which data file need to be inspected
% (b) Process these and find region of interest
% (c) Then use average value outside of ROI as background noise value
% (d) Use sum of ROI as missing observation
if FlagrefitAverage == 1
    for count = 1:size(ComCount, 1)
        if ComCount(count, 4) == 0 
               myFrame = double(imread([XPath,'.tif'], count));
               myROI = myFrame(29:35,27:33); % For a single, simulated spot.
               ResBackground =(sum(sum(myFrame))-sum(sum(myROI)))/((size(myFrame,1)^2)-(size(myROI,1)^2)); % Background is average count outside ROI.
               myROI = myROI - ResBackground;
               sum_signal = sum(myROI(:));
               % Then pass sum_signal this into the relevant place (ComCount(number, X=2, Y=3) )
               ComCount(count, 2) = sum_signal;
               ComCount (count, 4) = 2; % To check later if the values were moderated 
               
        end
        if ComCount(count, 5) == 0 
             myFrame = double(imread([YPath,'.tif'], count));
             myROI = myFrame(29:35,27:33); % For a single, simulated spot.
             ResBackground =(sum(sum(myFrame))-sum(sum(myROI)))/((size(myFrame,1)^2)-(size(myROI,1)^2)); % Background is average count outside ROI.
             myROI = myROI - ResBackground;
             sum_signal = sum(myROI(:)); 
             % Then pass sum_signal this into the relevant place (ComCount(number, X=2, Y=3) ) 
             ComCount(count, 3) = sum_signal;
             ComCount (count, 5) = 2; % To check later if the values were moderated 
        end
    end 
end


%Method 6: Re-inspect camera data for missing observations
% (a) Find out which frames of which data file need to be inspected
% (b) Process these using rainSTORM_fit_LocGF3_version2
%    ... (This is likely to be tricky)
% (c) Then use this fitted data, and subtract the background above,
%      to estimate orientation...


% 4. Caclulate polarization angle
% Take negative angle if corrected counts are negative

ComCount(:,6) = (180/pi)*acot(sqrt(abs(ComCount(:,2)./ComCount(:,3)))).*sign((ComCount(:,2)./ComCount(:,3)));



%5. Conclusions
% Results is a vector with the avergae estimated angle and the standard
% deviation of the phi column (not error of the mean). 
Phi = [Phireal,mean(ComCount(:,6)), std(ComCount(:,6))];


%6. Store data for plotting
%{
FileExists = exist('Polang.mat');
if FileExists == 2
    save ('Polang.mat',Phi45,'-append');
else
    save ('Polang.mat',['Phi',int2str(Phireal)]);
end
%}
if flagInitialRun
        PhiPlot = Phi;
else
        PhiPlot = cat(1, PhiPlot, Phi);
end

