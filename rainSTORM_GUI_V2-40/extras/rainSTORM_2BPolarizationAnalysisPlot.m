%rainSTORM_2BPolarizationAnalysisPlot
%This script uses the results from the PolarizationAnalysisScript to plot
%the observed and real angles for the processed values of phi
%Run the analysis script first!

% Made by Daan van Kleef & Mehdi Goudarzi

%{
S = load('Polang.mat');
trial_array = struct2array(S);
trial_size = size((trial_array));
final_size = trial_size(2)/3;
PhiPlot  = zeros(final_size,3);
for count=1:final_size
    PhiPlot(count,1) = trial_array(3*(count-1)+1);  
    PhiPlot(count,2) = trial_array(3*(count-1)+2);
    PhiPlot(count,3) = trial_array(3*(count-1)+3);
end
%}

hold on
axis([-5 50 -15 60])
xlabel('Real polarisation angle (degrees)');
ylabel('Calculated polarisation angle (degrees)');
errorbar(PhiPlot1(:,1),PhiPlot1(:,2),PhiPlot1(:,3),'--*b');
plot (PhiPlot1(:,1), PhiPlot1(:,1),'m');
errorbar(PhiPlot2(:,1),PhiPlot2(:,2),PhiPlot2(:,3),'--*r');
plot (PhiPlot2(:,1), PhiPlot2(:,1),'m');
hold off