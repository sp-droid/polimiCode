close all

mu_E = astroConstants(13);
a = 7000;
E0 = 0;
t0 = 0;

e = [0 0.2 0.4 0.6 0.8 0.95];

tOrb = 2*pi*sqrt( a^3/mu_E );
t = linspace(0,tOrb,100);

E = zeros(length(e),100);

for j=1:length(e)
    for i=1:100
        E(j,i) = keplerSolver(t(i), e(j), a, mu_E, t0, E0);
    end
end

Edeg = rad2deg(E);

%% Plots
figure()
for j=1:length(e)
    plot(t, Edeg(j,:), '-', LineWidth=2)
    hold on
end
xlabel('t [s]'); ylabel('T [deg]');
ylim([0,360])
title('2d plot');
legend('e = 0', 'e = 0.2', 'e = 0.4', 'e = 0.6', 'e = 0.8', 'e = 0.95', 'Location', 'southeast');
grid on;
hold off

figure()
surf(t, e, Edeg)
xlabel('t [s]'); ylabel('e [-]'); zlabel('E [deg]')
title('Surface plot');
grid on;
hold off

%% Functions
% Solver for Kepler's equation
function E = keplerSolver(t, e, a, mu, t0, E0)
%timescaling for relevant time information in plots
%
% PROTOTYPE
% E = keplerSolver(t, e, a, mu, t0, E0)
%
% INPUT:
% t[1] Time vector
% e[1] Eccentricity
% a[1] Semi-major axis
% mu[1] Gravitational parameter
% t0[1] Reference time
% E0[1] Reference eccentric anomaly
%
% OUTPUT:
% E[1] Eccentric anomaly
% tname String of the relevant time scale
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% VERSIONS
% 2022-09-21: v1
%
% -------------------------------------------------------------------------
% Mean motion and initial guess for E
n = sqrt(mu/a^3);
Eguess = n*t+e*sin(n*t)/(1-sin(n*t+e)+sin(n*t));

% Define function to solve for
func = @(E) E-e*sin(E)-n*t;

% Solve
E = fzero(func, Eguess);
end