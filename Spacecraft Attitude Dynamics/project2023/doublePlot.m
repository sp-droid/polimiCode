function doublePlot(x1, y1, x2, y2, xLabel, yLabel1, yLabel2, plotTitle)
    colors = {[0.9 0.7 0.1],'b','r'};
    figure;
    hold on;
    yyaxis left;
    plot(x1, y1,'LineWidth',1.5,'color',colors{1})
    ylabel(yLabel1);
    set(gca,'YColor', [0.9 0.7 0.1])
    yyaxis right;

    plot(x2, y2,'LineWidth',1.5,'color',colors{2})
    xlabel(xLabel);
    ylabel(yLabel2);
    set(gca,'YColor', 'b')
    title(plotTitle);
    grid on;
    set(gca,'fontsize', 15)
end