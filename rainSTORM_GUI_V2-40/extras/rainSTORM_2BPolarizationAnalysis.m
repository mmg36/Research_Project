%rainSTORM_2BPolarizationAnalysis
%This script uses the main RainSTORM Software to analyse 2 stacks of x-y
%polarized images

% Made by Daan van Kleef & Mehdi Goudarzi


clear
startup  % Run rainSTORM once (manually), on the X-polarisation data
uiwait 

Nxtable = struct2table(params.localization.results.SupResParams);
Nx = double(table2array(Nxtable(:,5)));

'xdone'

startup  % Run rainSTORM once (manually), on the Y-polarisation data
uiwait


Nytable = struct2table(params.localization.results.SupResParams);
Ny =double(table2array(Nytable(:,5)));


for count = 1:size(Nx)
    phiTable(count)= acot(sqrt(Nx(count)./Ny(count)));
end

