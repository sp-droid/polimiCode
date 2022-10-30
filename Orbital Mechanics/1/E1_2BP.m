close all
% 2 BODY PROBLEM - 2BP

% Physical parameters
%astroConstants(13) outputs the Earth's gravitational parameter [km^3/s^2]
mu_E = astroConstants(13);

% Initial condition
r0 = [ 26578.137; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0; v0 ];

% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
Torb = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]

tspan = linspace( 0, Torb, 1000 );

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );

%% Compute parameters
[Ts, tname] = timescaling(T);
Norbit = T/Torb;

Rnorm = vecnorm(Y(:,1:3).');
Vnorm = vecnorm(Y(:,4:6).');

h = cross(Y(:,1:3),Y(:,4:6));
hnorm = vecnorm(h.');

e = cross(Y(:,4:6),h)/mu_E-Y(:,1:3)./Rnorm.';
enorm = vecnorm(e.');

ehCheck = dot(h',e');

% Specific energy
specE = Vnorm.^2/2-mu_E./Rnorm;

% Radial and transversal velocity
vr = dot(Y(:,4:6)',Y(:,1:3)')./Rnorm;
vt = sqrt(Vnorm.^2-vr.^2);%dot(Y(:,1:3)',cross(h',Y(:,1:3)')./hnorm./Rnorm);

%% Plots
% Plot the orbit
figure()
plot3( Y(:,1), Y(:,2), Y(:,3), 'red-', LineWidth=2)
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;

% Plot r, v
figure()
plot( Ts, Y(:,1), LineWidth=2)
hold on
plot( Ts, Y(:,2), LineWidth=2)
hold on
plot( Ts, Y(:,3), LineWidth=2)
hold on
plot( Ts, Rnorm, 'black', LineWidth=2)
hold on
xlabel(strcat('Time [',tname,']')); ylabel('R [km]');
title('Distance');
grid on;
legend('rx', 'ry', 'rz', 'r')
hold off

figure()
plot( Ts, Y(:,4), LineWidth=2)
hold on
plot( Ts, Y(:,5), LineWidth=2)
hold on
plot( Ts, Y(:,6), LineWidth=2)
hold on
plot( Ts, Vnorm, 'black', LineWidth=2)
hold on
xlabel(strcat('Time [',tname,']')); ylabel('V [km/s]');
title('Velocity');
grid on;
legend('vx', 'vy', 'vz', 'v')
hold off

% Plot specific energy
figure()
plot( Norbit, specE, 'blue', LineWidth=2)
xlabel('Orbit number'); ylabel('€ [km^2/s^2]');
title('Specific Energy');
grid on;

% Plot h, e
figure()
plot( Norbit, h(:,1), '--', LineWidth=2)
hold on
plot( Norbit, h(:,2), '--', LineWidth=2)
hold on
plot( Norbit, h(:,3), '--', LineWidth=2)
hold on
plot( Norbit, hnorm, 'black--', LineWidth=2)
hold on
xlabel('Orbit number'); ylabel('h [km^2/s]');
title('Angular momentum');
grid on;
legend('hx', 'hy', 'hz', 'h')
hold off

figure()
plot( Norbit, e(:,1), '--', LineWidth=2)
hold on
plot( Norbit, e(:,2), '--', LineWidth=2)
hold on
plot( Norbit, e(:,3), '--', LineWidth=2)
hold on
plot( Norbit, enorm, 'black--', LineWidth=2)
hold on
xlabel('Orbit number'); ylabel('e [-]');
title('Eccentricity');
grid on;
legend('ex', 'ey', 'ez', 'e')
hold off

% Plot e-h dot product
figure()
plot( Norbit, ehCheck, LineWidth=2)
xlabel('Orbit number'); ylabel('e*h [km^2/s]');
title('e-h dot product');
grid on;

% Plot vr & vt
figure()
plot( Norbit, vr, '--', LineWidth=2)
hold on
plot( Norbit, vt, '--', LineWidth=2)
hold on
xlabel('Orbit number'); ylabel('v [km/s]');
title('Radial and transversal velocity');
grid on;
legend('vr', 'vt')
hold off

%% Functions
% Function ode_2bp
function dy = ode_2bp( ~, y, mu )
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
% Juan Luis Gonzalo Gomez
%
% VERSIONS
% 2018-09-26: First version
%
% -------------------------------------------------------------------------
% Position and velocity
r = y(1:3);
v = y(4:6);
% Distance from the primary
rnorm = norm(r);
% Set the derivatives of the state
dy = [ v
(-mu/rnorm^3)*r ];
end

% Function timescaling
function [Y, tname] = timescaling( T )
%timescaling for relevant time information in plots
%
% PROTOTYPE
% T, tname = timescaling( T )
%
% INPUT:
% T[nx1] Time vector
%
% OUTPUT:
% Y[nx1] Scaled time vector
% tname String of the relevant time scale
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% VERSIONS
% 2022-09-21: v1
%
% -------------------------------------------------------------------------
maxT = max(T);
if maxT<181
    Y = T;
    tname = 's';
elseif maxT<10860
    Y = T/60;
    tname = 'min';
elseif maxT<259200
    Y = T/3600;
    tname = 'h';
elseif maxT<7776000
    Y = T/86400;
    tname = 'days';
elseif maxT<62208000
    Y = T/2592000;
    tname = 'months';
else
    Y = T/31557600;
    tname = 'years';
end
end