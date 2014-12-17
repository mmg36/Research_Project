% Script by MG and DvK to perform a Monte-Carlo type simulation to estimate
% the overall polarisation distribution of an ensemble of fixed dipole
% emitters in a 2D plane. [23/11/2014]
clear

% Structure of program: 
% 1. Define distribution of dipole angles in 2D-plane.
% 2. Calculate Probability vs. polarisation histograms for each molecule.
% 3. Sum histograms for each molecule to obtain overall histogram.

TauD = 5;
AllMap(:,1) = [-1:0.001:1];
AllMap(:,2) = 0;
% 1. Define distribution of dipole angles in 2D-plane.
% Assumed to be a uniform distribution from 0 to pi.



for Phi = pi.*(0:0.01:1);
    Phi
    % 2. Calculate Probability vs. polarisation histograms for each molecule
    % 2.0 Initiate lists.
    Timestep = 0.1; %(fraction of tau between each frame).
    Discretisation = 1000; % Number of datapoints.
    Alpha = (0:Timestep:1);
    Alpha(1,1) = 0.001; % Fudging problems with inifity when stdev = 0.
    PolMap =zeros((Discretisation+1),size(Alpha,2)); 

    % For Alpha loop thingy
    for AlphaCount = 1:size(Alpha,2)
        clear DeltaPhiMap PhiMap PMap   
        
        AlphaTime = Alpha(AlphaCount);
        TMMap(:,1) = [-1:0.001:1];
        TMMap(:,2) = 0;
        % 2.1 For each point in time (value of alpha), simulate Probability distribution of dPhi.    
        
        Stdev = sqrt(AlphaTime*(1/3)*TauD);
        %DeltaPhiMap = zeros((Discretisation),2);
        DeltaPhiMap(:,1) = [-4*Stdev:0.001:4*Stdev];
        DeltaPhiMap(:,2) = (1/(sqrt(2*pi)*Stdev))*exp(-(DeltaPhiMap(:,1).^2)/(2*Stdev^2));
       % 2.2 For each point in time, simulate Probability distribution of Phi.
        %PhiMap = zeros((Discretisation),2);
        PhiMap(:,1) = DeltaPhiMap(:,1)+ ones(size(DeltaPhiMap,1),1).*Phi; 
        PhiMap(:,2) = DeltaPhiMap(:,2);

       % 2.3 For each point in time, simulate Probability distribution of
       % Polarisation.
        %PMap = zeros((Discretisation),2);
        PMap(:,1) = cos(2*PhiMap(:,1)); 
        PMap(:,2) = PhiMap(:,2);
        
       % 2.4 Discretise Pmap and plot histogram.
       PMap(:,1) = (1/1000)*round(1000*PMap(:,1));
       % Add Polarisation values to master Polarisation list. To do so, determine the correct bucket for each 
       %polarisation value to be added to. In addition, keep track of the
       %number of counts in each bucket. Then determine the average function value per bucket. 
       BPMap(:,1) = [-1:0.001:1];
       BPMap(:,2) = 0;
       CountMap = zeros(size(BPMap,1),1);
       for countRow = 1:size(PMap,1)
            RowNumMPM = (PMap(countRow,1)+1)*1000+1;
            RowNumMPM =int64(RowNumMPM); % Determine RowNumber in BPMap and CountMap for each localisation
            BPMap(RowNumMPM,2) =  BPMap(RowNumMPM,2)+ PMap(countRow,2);
            CountMap(RowNumMPM,1) =  CountMap(RowNumMPM,1) + 1;
       end
       
       % Now divide BPMap by CountMap to calculate average in each bucket.
       for countRow = 1:size(BPMap,1)
        if~ (BPMap(countRow,2) == 0)
           BPMap(countRow,2) = BPMap(countRow,2) ./ CountMap(countRow,1);
        end
       end
       
       % Now normalise each BPmap before adding it.
       NormFactor = trapz( BPMap(:,1),BPMap(:,2));
       NormBPMap = BPMap;
        NormBPMap(:,2) = BPMap(:,2)/NormFactor;
      NormFactor;
       
       % 2.4 Sum histograms for each alpha.
       TMMap(:,2) =  TMMap(:,2) + NormBPMap(:,2);
       trapz(TMMap(:,1),TMMap(:,2))
    end
    
       AllMap(:,2) =  AllMap(:,2) + TMMap(:,2);
       scatter(AllMap(:,1),AllMap(:,2))
end

    AllMap(AllMap(:,2) == 0) = [];
    