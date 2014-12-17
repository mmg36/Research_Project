
%rainSTORM_noiseCancelling
%This script removes the background noise from the polarised images capture
%straight from the camera before any analysis. 

% Made by Daan van Kleef & Mehdi Goudarzi.

function [CountMapX, CountMapY, NxCorrect, NyCorrect]  = rainSTORM_noiseCancelling(Nxarray, Nyarray, NxCorrect, NyCorrect, CountMapX, CountMapY, noiseThresh, LocRad, LocRange)
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
end
