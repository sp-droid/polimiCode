function simplePlot(x, y, xLabel, yLabel, plotTitle, legends)
    colors = {[0.9 0.7 0.1],'b','r'};
    figure;
    hold on;
    for i=1:length(y)
        plot(x, y{i},'LineWidth',1.5,'color',colors{i},'DisplayName',legends{i})
    end
    xlabel(xLabel);
    ylabel(yLabel);
    title(plotTitle);
    grid on;
    legend;
    set(gca,'fontsize', 15)
end