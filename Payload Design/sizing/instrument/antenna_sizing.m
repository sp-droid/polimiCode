clear
close all
clc

addpath(genpath('../shared'))

c = 3e8; % speed of light [m/s]
T_sys = 290; % [K]
k = 1.380649E-23; % Boltzmann constant

delta_r = 0.75; % range resolution for sub-cm precision [m]
f = 35e9; % center frequency
lambda = c/f;
B = c / (2*delta_r);

aperture = 10; % [m]

% array spacing multiple of wavelength to achieve peak performance 
% and maximize directivity while maintaining the impact of grating lobes
% limited
d = lambda / 2; 

% beam half-angle typical value for high resolution
theta_beam = 0.1*pi/180; % [rad]

% under the hypoteses of narrow beam, the following expr can be used for
% array length
L = lambda/theta_beam;

% the width of the array allows for more precise control over the beam 
% steering in the elevation direction (across the swath), ensuring accurate 
% topographic data --> higher angular resolution

% typical dimension for an antenna element in this frequency band
width = 0.3;

% number of elements for each antenna array
N = L/d + 1; % L = (N-1)*d

% antenna gains
G_t = 4*pi*(L*width)/(lambda^2); % Lesson 07 - slide 12
G_r = G_t;

%% Attenuation
% slant path distance
R = 901.83e3; % [m]
h_atm = 100e3; % height of the atmosphere (assumption)
slant_d = R/cosd(4.2062); % [m]
eta = 0.7;

RCS = 1;%4*pi*width^2*L^2/lambda^2; % from Lesson07 - slide 15
sigma0 = 10; % 10.^([-30, 5] .* 1e-1);

L_atm = pathAttenuationITUR676(90 - 4.2068, 35, 15, 101325, 40); % [dB]
L_fs_dB = 0; % 20*log10(lambda/(4*pi*slant_d)); % 20*log10(slant_d) + 20*log10(f) + 20*log10(4*pi/c) - 10*log10(G_t) - 10*log10(G_r);
L_tot_dB  = L_atm + L_fs_dB;

L_tot = 10^(L_tot_dB * 1e-1);

P_N = k*T_sys*B/1100; % thermal noise [W]
P_r = P_N; % 1; % [W] 
P_t = P_r * (4*pi)^3 * slant_d^4 / (G_t*G_r * lambda^2 * sigma0 * L_tot)

%% 
widths = linspace(0.1, 1, 1000);
PP_t   = zeros(length(widths), 2);

for i = 1:length(widths) 
    G_t = 4*pi*(L*widths(i))/(lambda^2);
    G_r = G_t;

    RCS = 4*pi*widths(i)^2*L^2/lambda^2; % from Lesson07 - slide 15

    L_atm = pathAttenuationITUR676(90 - 4.2068, 35, 15, 101325, 7.5); % [dB]
    % L_fs_dB = 20*log10(lambda/(4*pi*slant_d));
    L_fs_dB = 20*log10(slant_d) + 20*log10(f) + 20*log10(4*pi/c) - 10*log10(G_t) - 10*log10(G_r);
    L_tot_dB = L_atm + L_fs_dB;
    
    L_tot = 10^(L_tot_dB * 1e-1);
    
    PP_t(i, :) = P_r * (4*pi)^3 * R^4 ./ (G_t*G_r * lambda^2 * RCS.*sigma0 * L_tot);
end

close all

figure()
hold on

% Plot the data
plot(widths, PP_t(:, 1), 'LineWidth', 2)
% Define the y-values for the two horizontal lines
y1 = 1500;
y2 = 3000;

% Plot the horizontal lines
yline(y1, 'LineWidth', 1.3, 'Color', '#D95319') 
yline(y2, 'LineWidth', 1.3, 'Color', '#77AC30')

% Shade the area between the two lines
x_fill = [widths, flip(widths)];               % x values for the shaded area (same as plot)
y_fill = [ones(size(widths)) * y1, ones(size(widths)) * y2];  % y values for the shaded area

% Fill the area between the two horizontal lines
fill(x_fill, y_fill, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

grid minor
xlim([widths(1) widths(end)])
xlabel('Antenna Width [m]'), ylabel('Trasmitted Power $P_{\mathrm{TX}}$ [W]')
l = legend('', '$P_{\mathrm{min}} = 1500$ W', '$P_{\mathrm{max}} = 3000$ W');
t = get(l, 'Title');
set(t, 'String', 'Power Thresholds')

% exportgraphics(gcf, 'antenna_width.pdf', 'ContentType', 'vector')