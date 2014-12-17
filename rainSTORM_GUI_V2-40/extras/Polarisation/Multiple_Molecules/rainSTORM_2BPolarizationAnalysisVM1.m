%rainSTORM_2BPolarizationAnalysis
%This script uses the main RainSTORM Software to analyse 2 stacks of x-y
%polarized images

% Made by Daan van Kleef & Mehdi Goudarzi

clear 


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

% Setting the noise threshold. The number of succesful localisations per
% pixel should exceed this value in order to be marked as a fluorophore
% position.
noiseThresh = 1;

% Maximimum search radius to find local maximima in countmap to mark
% fluorophore positions.
LocRad = 1;

% Maximim distance any localisation can be from a temporary molecule
% location to be associated with that molecule.
LocRange = 1;

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

% 2.Quality control: This creates a density map for the frame and 
frameSize = params.rawdata_mgr.myImInfo.frame_size (1);
CountMapX = zeros (frameSize);

NxCorrect = [Nxarray(:,1) Nxarray(:,2) Nxarray(:,3) Nxarray(:,12) -1*ones(size(Nxarray,1),2)];

for countNarr = 1 : size (Nxarray, 1)
    xtemp = floor(Nxarray(countNarr, 3)); 
    ytemp = floor(Nxarray(countNarr, 2));
    % All of the localisation positions are rounded down. 
    CountMapX (xtemp, ytemp) = CountMapX (xtemp, ytemp) +1; 
    % Creating the number of counts per pixel in the image to create a
    % density map later. 
end
    
%To remove most of the background noise. 
CountMapX = CountMapX  - noiseThresh;
CountMapX(CountMapX<0) = 0;

% Determine preliminary locations of fluorophores, using avNMWS
PrelimLoc = rainSTORM_avNMS(CountMapX,LocRad);

%Now look for localisations within 2 pixels of the preliminary locations
%(5x5 grid) to find the real average of the localisations. Then assign
%localisations to the corresponding molecule

for LocCount = 1:size(PrelimLoc,1);
    for ArrayCount = 1:size(Nxarray,1);
        Locdistance = sqrt((Nxarray(ArrayCount,3)-PrelimLoc(LocCount,1))^2+(Nxarray(ArrayCount,2)-PrelimLoc(LocCount,2))^2);
        if Locdistance < LocRange
         NxCorrect(ArrayCount,5) = PrelimLoc(LocCount,2);
         NxCorrect(ArrayCount,6) = PrelimLoc(LocCount,1);
        end
    end
end

NxCorrect(NxCorrect(:, 5) == -1, :) = [];

% Idem for Y-stack
CountMapY = zeros (frameSize);
NyCorrect = [Nyarray(:,1) Nyarray(:,2) Nyarray(:,3) Nyarray(:,12) -1*ones(size(Nyarray,1),2)];

for countNarr = 1 : size (Nyarray, 1)
    xtemp = floor(Nyarray(countNarr, 3)); 
    ytemp = floor(Nyarray(countNarr, 2));
    % All of the localisation positions are rounded down. 
    CountMapY (xtemp, ytemp) = CountMapY (xtemp, ytemp) +1; 
    % Creating the number of counts per pixel in the image to create a
    % density map later. 
end
    
%To remove most of the background noise. 
CountMapY = CountMapY  - noiseThresh;
CountMapY(CountMapY<0) = 0;

% Determine preliminary locations of fluorophores, using avNMWS
PrelimLoc = rainSTORM_avNMS(CountMapY,LocRad);

%Now look for localisations within 2 pixels of the preliminary locations
%(5x5 grid) to find the real average of the localisations. Then assign
%localisations to the corresponding molecule

for LocCount = 1:size(PrelimLoc,1);
    for ArrayCount = 1:size(Nyarray,1);
        Locdistance = sqrt((Nyarray(ArrayCount,3)-PrelimLoc(LocCount,1))^2+(Nyarray(ArrayCount,2)-PrelimLoc(LocCount,2))^2);
        if Locdistance < LocRange
         NyCorrect(ArrayCount,5) = PrelimLoc(LocCount,2);
         NyCorrect(ArrayCount,6) = PrelimLoc(LocCount,1);
        end
    end
end

NyCorrect(NyCorrect(:, 5) == -1, :) = [];


% 2. Identify matching images


%Initiate matrix to hold combined count list
% ComCount will hold the combined count list.
% Column Structure:
    %  Frame ID  | Nx | Ny | Nx present? | Ny present? | Estimated angle
    % Last 2 column will hold flags to indicate succesful counts
        % 0 = absent, 1 = present
        % Can be used to relocalize missing observations later

numberOfFrames = params.localization.results.numberOfFrames;
ComCount = []; 

for count = 1 : numberOfFrames
        dumx = NxCorrect(find (NxCorrect(: , 1) == count),:);
        dumy = NyCorrect(find (NyCorrect(: , 1) == count),:);
        if (~isempty(dumx) && ~isempty(dumy)) 
            dumx(:,[5,6,7]) = zeros(size(dumx,1),3);
            dumy(:,[5,6]) = [];
            for countX = 1:size(dumx,1);
                for countY = 1:size(dumy,1)
                    DistanceXY = sqrt((dumx(countX,2)-dumy(countY,2))^2+(dumx(countX,3)-dumy(countY,3))^2);
                    if DistanceXY <= LocRange;
                        dumx(countX, [5,6,7] ) = dumy(countY, [2, 3, 4]); 
                        dumy(countY,1) = -1;
                    end;
                end;
            end;
            dumy(dumy(:, 1) == -1, :) = [];
            dumy = [dumy(:,1) zeros(size(dumy,1),3) dumy(:,[2,3,4])];
            dumx = cat(1, dumx, dumy);
            if isempty(ComCount) 
                ComCount = dumx;
            else 
                ComCount = cat (1, ComCount, dumx); 
            end;
        elseif isempty(dumy)
            dumx(:,[5,6,7]) = zeros(size(dumx,1),3);
            if isempty(ComCount) 
                ComCount = dumx;
            else 
                ComCount = cat (1, ComCount, dumx); 
            end;
        elseif isempty(dumy)
            dumy = [dumy(:,1) zeros(size(dumy,1),3) dumy(:,[2,3,4])];
            if isempty(ComCount) 
                ComCount = dumy;
            else 
                ComCount = cat (1, ComCount, dumy); 
            end;
        end;
end;

% 3. Complete missing paired observations

% Method 1: Assume missing observations are count 10^-10
if(Flagsetzero == 1);

for count = 1:size(ComCount,1);
    if ComCount(count, 4) == 0
        ComCount (count , 4) = 10^(-10);
         ComCount (count , [2,3]) = ComCount (count , [5,6]);
    end
    if ComCount(count,7) == 0
        ComCount(count,7) = 10^(-10);
        ComCount (count , [5,6]) = ComCount (count , [2,3]);
    end
end
end

% Method 2: Assume missing observations are background count level
%Legacy 

% Method 3. Subtract average residual background level from measured
%accepted intensities and put missing observations to 0
% Legacy, currently not in use


% Method 4: Re-inspect camera data for missing observations
% (a) Find out which frames of which data file need to be inspected
% (b) Process these and find region of interest
% (c) Then use minimum of ROI as background noise value
% (d) Use sum of ROI as missing observation
% Legacy, currently not in use
    

% Method 5: Re-inspect camera data for missing observations
% (a) Find out which frames of which data file need to be inspected
% (b) Process these and find region of interest
% (c) Then use average value outside of ROI as background noise value
% (d) Use sum of ROI as missing observation
if FlagrefitAverage == 1
    for count = 1:size(ComCount, 1)
        BGN = mean([params.localization.results.SupResParams.noiseBGN]);
        if ComCount(count, 4) == 0 
               myFrame = double(imread([XPath,'.tif'], ComCount(count, 1)));
               ComCount(count, [2,3]) = ComCount(count, [5,6]);
               myROI= myFrame(ceil(ComCount(count,2))-3:ceil(ComCount(count,2))+3,ceil(ComCount(count,3))-3:ceil(ComCount(count,3))+3) - BGN;
               sum_signal = sum(myROI(:));
             %  if sum_signal < 0 
             %      sum_signal = 0;
             %  end;
               % Then pass sum_signal this into the relevant place (ComCount(number, X=2, Y=3) )
               ComCount(count, 4) = sum_signal; 
               
        end
        if ComCount(count, 7) == 0 
             myFrame = double(imread([YPath,'.tif'], ComCount(count, 1)));
             ComCount(count, [5,6]) = ComCount(count, [2, 3]);
             myROI= myFrame(ceil(ComCount(count,2))-3:ceil(ComCount(count,2))+3,ceil(ComCount(count,3))-3:ceil(ComCount(count,3))+3) - BGN;
             sum_signal = sum(myROI(:));
             % Then pass sum_signal this into the relevant place (ComCount(number, X=2, Y=3) )
            % if sum_signal < 0 
            %       sum_signal = 0;
           %  end;
             ComCount(count, 7) = sum_signal; 
                
        end
    end 
end

% 4. Caclulate polarization angle
% Take negative angle if corrected counts are negative

ComCount(:,8) = (180/pi)*acot(sqrt(abs(ComCount(:,4)./ComCount(:,7)))).*sign((ComCount(:,4)./ComCount(:,7)));



%5. Conclusions
ComCountAv = ComCount;
ComCountAv (:,2) = 0.5*(ComCountAv(:,2)+ComCountAv(:,5));
ComCountAv (:,3) = 0.5*(ComCountAv(:,3)+ComCountAv(:,6));
ComCountAv(:, [5,6]) = [];
CeilComCount  = [ ComCountAv(:,1)  ceil(ComCountAv(:,[2,3])) ComCountAv(:,[4,5,6])] ;
CountMapPhi = zeros (frameSize);
PhiMap = nan (frameSize);
StdevMap = nan (frameSize);
for countRow = 1 : frameSize
    for countCol = 1 : frameSize
        for MapCount = 1 : size (CeilComCount, 1) 
            dumPhi =[];
            if (CeilComCount(MapCount, 2) == countRow && CeilComCount(MapCount, 3) == countCol)
               dumPhi= cat(1,CeilComCount(MapCount, 6),dumPhi); 
               CountMapPhi(countRow, countCol) = CountMapPhi(countRow, countCol)+1;
            end
            
            PhiMap(countRow, countCol) = mean(dumPhi);
            StdevMap(countRow, countCol) = std(dumPhi);
        end
    end
end


% Results is a vector with the avergae estimated angle and the standard
% deviation of the phi column (not error of the mean). 

