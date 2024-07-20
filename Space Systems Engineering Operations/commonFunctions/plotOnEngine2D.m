function plotOnEngine2D(graphx,graphy,titleStr,yLabelLeft,yLabelRight,xLabelStr,y2logscale)
fnx = fieldnames(graphx);
fny = fieldnames(graphy);

figure;
if numel(fny)==1
    scatter(graphx.(fnx{1}),graphy.(fny{1}))
    hold on
else
    C1 = [0 0.4470 0.7410];  % blue
    C2 = [0.8500 0.3250 0.0980];  % orange
    C3 = [0.9290 0.6940 0.1250];  % yellow
    C4 = [0.4940 0.1840 0.5560];  % purple
    C5 = [0.4660 0.6740 0.1880];  % green
    C6 = [0.3010 0.7450 0.9330];  % light blue
    C7 = [0.6350 0.0780 0.1840];  % red
    C8 = [1.0000 0.4000 0.6000];  % pink
    C9 = [0.6350 0.0780 0.1840];  % brown
    C10 = [0.0000 0.0000 0.0000];  % black
    for i=1:numel(fny)
        if numel(fnx)==1
            fnxi = fnx{1};
        else
            fnxi = fnx{i};
        end
        if i==1
            yyaxis left
            plot(graphx.(fnxi),graphy.(fny{i}),'LineWidth',2,'DisplayName',fny{i},'Color',eval(['C' num2str(i)]))
            title(titleStr,'Interpreter','latex');
            xlabel(xLabelStr,'Interpreter','latex'); ylabel(yLabelLeft,'Interpreter','latex');
            % ylim([8,32])
            yyaxis right
            continue
        end
        %scatter(graphx.(fnxi),graphy.(fny{i}),2,'DisplayName',fny{i},'MarkerEdgeColor',eval(['C' num2str(i)]))
        plot(graphx.(fnxi),graphy.(fny{i}),'LineWidth',2,'LineStyle','-','DisplayName',fny{i},'Color',eval(['C' num2str(i)]))
        hold on
    end
end
if exist('y2logscale', 'var')
    if (y2logscale == true)
        set(gca, 'YScale', 'log');
    end
end
ylabel(yLabelRight,'Interpreter','latex');
grid on;
xlim([min(graphx.(fnx{1})), max(graphx.(fnx{1}))])
if numel(fny)~=1
    legend('Location','best');
end
set(gca,'fontsize', 16) 
hold off
end