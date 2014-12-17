count = size(Dfactory,1);
pos = [-3:1:3]';
bar (pos, yCols+7*Dfactory(count,4),1), colormap([0.4,0.4,0.4]);
   hold on
    x = [-3.5:0.1:3.5]';
    xlabel('CCD pixel(rows)') % x-axis label
    ylabel('Counts in rows') % y-axis label
    plot (x, Dfactory(count, 1)*exp(-((x-Dfactory(count, 2)).^2)/(2*Dfactory(count,3).^2))+7*Dfactory(count,4), 'r');
savefig('Yplot');
    %close;

        