%% Test on attenuation functions
clear
close all
clc

addpath(genpath('../shared'))

%% Figure on sea-level specific attenuation
T = 15;
P = 101325;
rhoVaporStandard = 7.5;
rhoVaporDry = 0;

freqs = 0:1000;
specificAttDRY = arrayfun(@(freq) attenuationITUR676( freq, T, P, rhoVaporDry ), freqs);
specificAttSTD = arrayfun(@(freq) attenuationITUR676( freq, T, P, rhoVaporStandard ), freqs);

%%
figure;
semilogy(freqs, specificAttDRY, 'LineWidth', 2, 'DisplayName', 'Dry')
hold on;
semilogy(freqs, specificAttSTD, 'LineWidth', 2, 'DisplayName', 'Standard')
grid on;
xlabel('Frequency [GHz]')
ylabel('Specific attenuation [dB/km]')
legend;
% set(gca, 'fontsize', 14)
% print(gcf, 'output/attenuationSeaLevel.png', '-dpng', '-r200');

%% Figure on zenith attenuation
elevAngle = 90;

freqs = 0:350;
attDRY = arrayfun(@(freq) pathAttenuationITUR676( elevAngle, freq, T, P, rhoVaporDry ), freqs);
attSTD = arrayfun(@(freq) pathAttenuationITUR676( elevAngle, freq, T, P, rhoVaporStandard ), freqs);

figure;
semilogy(freqs, attDRY, 'LineWidth', 2, 'DisplayName', 'Dry')
hold on;
semilogy(freqs, attSTD, 'LineWidth', 2, 'DisplayName', 'Standard')
grid on;
title('Zenith attenuation due to atmospheric gases')
xlabel('Frequency [GHz]')
ylabel('Attenuation [dB]')
legend;
set(gca, 'fontsize', 14)
print(gcf, 'output/attenuationZenith.png', '-dpng', '-r200');

%% Figure on slanted attenuation
freqs = 0:350;
elevAngles = 10:20:90;

figure;
for i=1:length(elevAngles)
    elevAngle = elevAngles(i);
    att = arrayfun(@(freq) pathAttenuationITUR676( elevAngle, freq, T, P, rhoVaporStandard ), freqs);

    semilogy(freqs, att, 'LineWidth', 2, 'DisplayName', strcat('Slant=',string(elevAngle),'º'))
    hold on
end

grid on;
title('Slanted attenuation due to atmospheric gases')
xlabel('Frequency [GHz]')
ylabel('Attenuation [dB]')
legend('Location','best');
print(gcf, 'output/attenuationSlanted.png', '-dpng', '-r200');

%%
close all
clc

freqs = 10:100;
att = arrayfun(@(freq) pathAttenuationITUR676(90 - 4.2068, freq, T, P, rhoVaporStandard ), freqs);

figure()
semilogy(freqs, att, 'LineWidth', 2)
hold on

xlabel('Frequency [GHz]')
ylabel('Atmospheric Attenuation [dB]')
grid minor
xlim([freqs(1) freqs(end)])
% text(35, log10(att(35)), ['(' num2str(35) ', ' num2str(log10(att(35))) ')'], 'FontSize', 18, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')

%%
T = 30;
P = 101325;
rhoVaporStandard = 30;

freqs = linspace(20,50,200);
specificAttSTD = arrayfun(@(freq) attenuationITUR676( freq, T, P, rhoVaporStandard ), freqs);

figure;
plot(freqs, specificAttSTD, 'LineWidth', 2)
grid minor;
xlabel('Frequency [GHz]')
ylabel('Specific attenuation [dB/km]')
set(gca, 'fontsize', 20)
title('ITU-R P.676-12 attenuation model', 'FontSize', 23)
% print(gcf, 'output/attenuationZoom.png', '-dpng', '-r200');
exportgraphics(gcf, 'spec_att.pdf', 'ContentType', 'vector')