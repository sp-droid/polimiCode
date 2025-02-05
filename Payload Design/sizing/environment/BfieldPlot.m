clear
close all
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [100 100 700 500]);

load("Bfield.mat");

%%

t = linspace(0, 1.1714, length(Bfield));

plot(t,Bfield,'LineWidth',1.5)
%xlim([2004.8,2037])
%ylim([0,400])
xlabel('Orbit time [h]')
ylabel('B [Gauss]')
grid on;

print(gcf, 'output/Bfield.png', '-dpng', '-r300');