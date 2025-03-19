clear
close all
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [50 100 1000 500]);
addpath(genpath('../shared'))

load("f107.mat");

%%

t = f107(:,1)./365.25+2000;

tsecs = linspace(0, 17.05, 300)*365.25*24*3600;
projected = arrayfun(@(x) f107estimation(x, 7287.5), tsecs);

plot(t,f107(:,2),'LineWidth',1.5)
hold on
plot(tsecs/365.25/86400+2019.95,projected,'LineWidth',1.5)
xlim([2004.8,2037])
ylim([0,400])
ylabel('SFU')
xlabel('Year')
legend('Real','Projected')
grid on;

print(gcf, 'output/f107.png', '-dpng', '-r300');