
count = size(Dfactorx,1);
pos = [-3:1:3]';
bar (pos, yRows+7*Dfactorx(count,4),1), colormap([0.4,0.4,0.4]);
   hold on
    x = [-3.5:0.1:3.5]';
    xlabel('CCD pixel(cols)') % x-axis label
    ylabel('Counts in column') % y-axis label
    plot (x, Dfactorx(count, 1)*exp(-((x-Dfactorx(count, 2)).^2)/(2*Dfactorx(count,3).^2))+7*Dfactorx(count,4), 'r');
savefig('Xplot');
    %close;

        