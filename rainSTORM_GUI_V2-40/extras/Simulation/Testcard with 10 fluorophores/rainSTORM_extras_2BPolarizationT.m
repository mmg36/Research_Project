% rainSTORM_extras_2BPolarization
% Mehdi Goudarzi & Daan van Kleef
% 20/10/14
%
% Later Edits
%   None yet
%
% Function
%   Simulates molecular diffusion imaging data for polarisation-resolved Localisation Microscopy
%   Modified from rainSTORM_extras_simulate.m in V2.32
% This is the parent script which set the value of intrinsic variables.
%    Near the end of this script another script is called which deals with
%   the creating the images of photons. 
%
% Motivation
%   Demonstrate and explore forward-simulation of diffusion imaging
%   Generate data for validation and testing of diffusion analysis methods
%
% Further Reading
%   Einstein 1905 Brownian Motion
%   Thompson 2002 Biophys J 
%   Savin and Doyle 2005 Biophys J
%   Crocker and Grier 199(2)
%
% Outline of this program
%   0. Set the initial conditions of the simulation which is a single
%       fixed molecule. 
%   1. Generate Polarised photon count fractions on x and y image 
%   2. Call a second script to generate the images for eac plane of 
%       polarised light.
% 
% Notes
%   The folder of individual TIFs can be assembled into a stack by ...
%     dropping the folder onto the ImageJ bar.
%   Molecules can drift off-camera in a slightly non-physical way.
%     For now, avoid generating this kind of data by keeping dyes mid-CCD.
%   The dimensions used in this simualtion are nanometres and seconds
%   A later simulation will consider time-correlated fluorescence
%   The hard-coded boxcar window in the PSF simulation is a flaw
%     But it seems insignificant for the default widths !!

clear

% 0. Sets the initial conditions of the simulation
close all                    % Close figures
rng('shuffle')               % Start a random number sequence (see 'rng')

% Initial conditions
flagOneMol           = 1;    % Setup simulation with one dye molecule
flagTwoMol           = 0;    % Alternative. Simulation of 2 molecules
flagMultiMol         = 0;    % Alternative. Simulation of many molecules
flagMultiMolBrowDis  = 0;    % Alternative. Multi mols, distinct D's
flagTimeCorrBlinks   = 0;    % Alternative. See below for details. 
numberOfImages       = 1000;    % How many frames of data to simulate?

% Purely in-loop simulation conditions
flagBrownianConst    = 0;    % Molecule(s) move with Brownian diffusion

% Polarization options
phiAngle             = [0:5:45]';  % in degrees

% Positions of fluorophores [Number, PosX,PosY]
posVec               = [1 10 10; 2 30 10; 3 50 10; 4 50 20; 5 10 40; 6 30 35; 7 30 40; 8 30 45; 9 50 30; 10 50 50];

% 3D options
flag2D               = 1;    % True: all z-positions are zero

% Output options
flagSaveTifs         = 1;    % Save grayscale tif frames
flagShowCCDImages    = 1;    % Show simulated frames on screen (IS SLOW!)
generalpath          = 'C:\Users\user\Desktop\results\Timage\';

% General path stores the location of the results folder which all of the
% folders containing the images of this code are located there. 
% based on the angle of the molecule in question, a subfolder sharing 
% the same name is created.
% Polarisation is measured in X and Y planes and images corresponding to
% each direction is stored in the corresponding subfolder of the above subfolder. 
writeFileX = strcat(generalpath,'\X\'); % Save images here
writeFileY = strcat(generalpath,'\Y\');
mkdir(writeFileX);
mkdir(writeFileY);

% 1. Generate Polarised photon count fractions on x and y image 
% Not delayed evaluation, might need to update later
dyePhotonsX = (cos(phiAngle.*pi/180)).^2;
dyePhotonsY = (sin(phiAngle.*pi/180)).^2;

% 2. Combine data into matrix to pass on to generation script
% [Number, Xpos, Ypos, photonfracX, photonfracY]
SimDatax = [posVec dyePhotonsX];
SimDatay = [posVec dyePhotonsY];

%3. Generate random number matrix to determine which fluorophores are
%active in each frame.
dyeRand = rand(size(SimDatax,1),numberOfImages);

% 4. Call a second script to do the rest of image processing
% Analysing images with polarised light in the X-direction
SimData = SimDatax;
writeFile = writeFileX;
% A second script which generates the images is called here
rainSTORM_extras_2BSimV3T;

% Analysing images with polarised light in the Y-direction
SimData = SimDatay;
writeFile = writeFileY;
% A second script which generates the images is called here
rainSTORM_extras_2BSimV3T;
