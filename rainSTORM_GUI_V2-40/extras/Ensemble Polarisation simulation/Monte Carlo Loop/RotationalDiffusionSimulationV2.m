% Forward simulation of expected single molecule polarisations
%
% Version 1: Monte Carlo wobbling axis simulation
%            Consider rotation in 2D only
%            Model detection as Iy = I0_cos^2(phi) and Ix = I0_sin^2(phi)
%
% References: Harms, BioPhysJ, 1999, vol 77, p. 2864
%             And Debye-Stokes-Einstein relation in Cichos
%
% IIB project aims: 
%      1. See if the physics seems reasonable. 
%          It seems to match Harms 1999, so hopefully it is.
%      2. Generate histograms for key cases.
%           Solid, our biogel, our melted biogel, glycerol, water
%           Show these as curves, showing progression over time
%           Possibly animate this for powerpoint.
%      3. Predict the expected distribution for our biogel (solid and melt)
%           Which has viscoity ~ 240 Pa s, GFP radius 3 nm, 20 ms exposure
%           time
%      4. Combine results into the figure for this
%
%      5. Generate a table of simulated results, and look for a correlation
%           between std(listPolarisations) and camera exposure time / T_C


numberOfDyes = 100000;    % This many separate fluorophore simulations
% KILLED: numberOfTimePoints = 100; % Number of time points during camera exposure
                         % One point means an instantaneous measurement
                         % Number increments = number of points - 1

photonsPerFrame = 1000; % Expected number of photons collected per frame
                        % Assuming infinite repeats, there is no point 
                        %  making these into random numbers and averaging
                        %  

D_rot = 6E-60; % k_B T / (8 pi eta a^3)... 6000 per second in 1 Pa s, 3 nm
            % Debye-Stokes-Einstein relation
            % Use 100 Pa s (gel) and 3 nm radius for GFP.

T_C = 1/(6*D_rot); % rotational correlation time
tIncrement = 5E-4; % 1 ms % Use a time increment smaller than T_C

sigma = sqrt( 2 * D_rot * tIncrement ); % See Harms 1999;

generalpath          = 'C:\Users\user\Desktop\results\MonteCarlo\';
mkdir(generalpath);

listIy = zeros(numberOfDyes,1);
listIx = zeros(numberOfDyes,1);

for lpDye = 1:numberOfDyes
    
    phi0 = pi*rand(); % Random initial orientation in 2D
    phiT = phi0;
    
    for lpTime = 1:numberOfTimePoints
        
        Iy = photonsPerFrame * (cos(phiT))^2; % photons on X channel
        Ix = photonsPerFrame - Iy;            % photons on Y channel
        
        listIx(lpDye) = listIx(lpDye) + Ix;
        listIy(lpDye) = listIy(lpDye) + Iy;
        
        phiT = phiT + sigma*randn;  % Simulate axis wobble
    end
    
end

listPol = (listIx - listIy)./(listIx + listIy);

fg = figure(1);
hist(listPol,400)
  xlim([-1 1])
  xlabel('Polarisation', 'fontSize', 18)
  ylabel('Number of observations', 'fontSize', 18)
  title([num2str(0.5*numberOfTimePoints),' ms'], 'fontSize', 18)
  set(gca, 'fontSize', 18)
  saveas(fg, fullfile(generalpath, int2str(numberOfTimePoints)), 'png')
close(fg)