% Trajectory Identifying Program
%
% DRAFT!
%   Likely bugs: Merging trajectories
%                Consequences of Naive jump-radius method of analysis
%
% Feb 2013 EJR
%
% NOTES
%   Start with boxedPosits, boxedParams
%   Plot some extraneous 
%
% FURTHER WORK
%   Integrate "flagPrintFriendly" methods from Hyungsun Kim, for graphics
%


% 0. FLOW CONTROL
flagShowDurationsHist = 1; 
flagPlotTrajects = 1;

% 1. INPUT
% Assuming data exists: boxedPosits, boxedParams, pixelWidth
% To use all reviewed localisations, uncomment the following:
% boxedPosits = reviewedPosits;
% boxedParams = reviewedParams;

jumpRad = 2.5; % Units: CCD pixel widths
               % Temporally-consequtive points 
               % which are spatially this close,
               % are identified as part of the same trajectory.
               % IF multiple points match, only the closest is matched
timeCycle = 0.04; % Camera readout cycle time (e.g. 20 ms)

trajectN = -ones(size(boxedPosits,1),1); % Trajectory identifier
               % Stores a trajectory for each localisation (boxedPosits)
               % Initialise as "-1" which is an error code if not replaced

trajectCount = 1; %Trajectory identifier for the next trajectory
numberOfPoints = size(boxedPosits,1);

boxedSequence = (1:1:numberOfPoints)'; % Column vector needed below

% 2. PROCESS
% Identify trajectories by "Nearest Successor Threading"
% The method is a loop. For each boxed position:
for lpPoints = 1:numberOfPoints

  myPosit= boxedPosits(lpPoints,:); % Identify this boxed position
  nowFrm = boxedParams(lpPoints,7); % And its frame number
    
  if(trajectN(lpPoints) == -1)
    % If this point does not yet belong to a trajectory
    trajectN(lpPoints) = trajectCount; % Assign it as the next trajectory
    trajectCount       = trajectCount+1;
  end
  myTrajectN = trajectN(lpPoints); % The new or existing trajectory number

  % Look ahead by 1 frame to see if this trajectory continues...
  % And, if so, label this continuation.
  nextFramePosits = boxedPosits(boxedParams(:,7)==nowFrm+1,:);
  nextFrameSequen = boxedSequence(boxedParams(:,7)==nowFrm+1);
  nextFrameDists  = 100*ones( size(nextFramePosits,1) ,1 );
  for lpNFPs = 1:size(nextFramePosits,1)
    nfDspmt = nextFramePosits(lpNFPs,:)-myPosit; % Displacement
    
    nextFrameDists(lpNFPs) = sqrt( nfDspmt(1)^2 + nfDspmt(2)^2 );
  end
  
  if(size(nextFramePosits,1) >0)
    % If some data exist, making sortrows possible...
    % Are any of these points close enough to be the same fluorophore?
    % Catenate columns into a matrix for simple sorting by rows...
    nextData = [nextFramePosits,nextFrameSequen,nextFrameDists];
    nextData = sortrows(nextData,4); % Sort by distance
    if(nextData(1,4) < jumpRad )     % If the closest point is close enough
      matchSequen = nextData(1,3);   % Identify the point (boxedPosits Row)
      % Annote the current trajectory number to the future position...
      trajectN(matchSequen) = myTrajectN;
    end
  end 
  % And then continue to process the next localisation
  
  mybar = waitbar(lpPoints/numberOfPoints);
end
% Have now processed through all the boxed localisations
close(mybar)


% 3. OUTPUT
% All the boxedPosits should now be threaded by trajectory number 
% So extract trajectories and plot a suitable image

numberOfTrajectories = max(trajectN(:)); % Or ( trajectCount - 1 )

figure(3)
set(gca,'Color',[0,0,0])
hold on
diffCmap = colormap(hsv(48));
diffCmapTop = 0.3E-12; % iE-12 is for six microns^2 / second. Try 4.

for lpTr = 1:numberOfTrajectories
   
    aTrajPosits = boxedPosits(trajectN == lpTr,:); % Get the trajectory
    
    % Colour code for plotting trajectories by MSD (of a one-frame step)
    if (size(aTrajPosits,1)>1) % If the trajecory is not a single point..
      deltaCoords = aTrajPosits(2:end,1:2) - aTrajPosits(1:end-1,1:2);
      squaredDisplacement = sum(deltaCoords.^2,2);
      msdOneStep = mean(squaredDisplacement);
      
      aTrajD = (pixelWidth^2)*msdOneStep/(4*timeCycle); % Diffusivity
      aTrajD = aTrajD*1E-18; % Convert to m^2/s
      
      aTrajDcolor = ceil(48*aTrajD/diffCmapTop); % Place on Colour scale
        if(aTrajDcolor>48)
          aTrajDcolor=48;
        end
        
      aTrajColor = diffCmap(aTrajDcolor,:); % Assign colour
    else
      aTrajColor = [0,0,0]; % If the "trajectory" is a single point
    end
    
    if(flagPlotTrajects)
     if (size(aTrajPosits,1) > 3) % Plot trajectories at least this long
     plot(aTrajPosits(:,2)*pixelWidth,-aTrajPosits(:,1)*pixelWidth, ...
         'lineWidth',1.5, 'color', aTrajColor );
     end
    end
%     if (rem(lpTr,20) == 0) % Save frames for a movie
%     myframe = getframe(gcf);
%     myframe = myframe.cdata;
%     imwrite(myframe,['C:\sims\sims\', int2str(lpTr),'.png'], 'png');
%     end
     myTbar = waitbar(lpTr/numberOfTrajectories);
end
close(myTbar)

% figure(3) % Absolute references may mess up V 2-34 rainSTORM
figure
axis square

xlabel('x, nm', 'fontsize',14)
ylabel('y, nm', 'fontsize',14)
title('Fluorophore Trajectories', 'fontsize',14)
set(gca,'FontSize',14,'fontweight','bold');

hold off

% A histogram to show Spot lifetimes:
% Process trajectN...
trajectoryDurations = zeros(numberOfTrajectories,1);
trajectoryPhotonNums= zeros(numberOfTrajectories,1);

for lpTr = 1:numberOfTrajectories
    
    trajectoryDurations(lpTr) = sum(trajectN==lpTr);
    trajectoryPhotonNums(lpTr) = sum( reviewedPhotonNums(trajectN==lpTr) );
    
end

if (flagShowDurationsHist)
  % figure(4) % Absolute references may mess up V 2-34 rainSTORM
  figure 
   hist(trajectoryDurations,1:1:50)
   xlabel('Number of Frames Duration')
   ylabel('Number of Trajectories')
  figure
   hist(trajectoryPhotonNums,100);
   xlabel('Photons detected in trajectory (estimate)')
   ylabel('Number of trajectories')
  % Multiply the number of frames by the camera cycle time (e.g. 40 ms)
  % to determine the "on duration" or "persistence time" 
  % of the identified trajectories
end

