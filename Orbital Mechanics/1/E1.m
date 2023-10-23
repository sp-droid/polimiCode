close all

% Physical parameters
%astroConstants(13) outputs the Earth's gravitational parameter [km^3/s^2]
mu_E = astroConstants(13);

% Initial condition
r0 = [ 1599.4; 5859.1; 3019.2 ]; % [km]
v0 = [ -5.9909; -2.3882; 7.8083 ]; % [km/s]
y0 = [ r0; v0 ];

% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
Torb = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 2*Torb, 1000 );

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );

% Compute parameters
Rnorm = vecnorm(Y(:,1:3).');
Vnorm = vecnorm(Y(:,4:6).');

Norbit = T/2/Torb;

h = cross(Y(:,1:3),Y(:,4:6));
hnorm = vecnorm(h.');

% Specific energy
specE = Vnorm.^2/2-mu_E./Rnorm;

% Plot the orbit
figure()
plot3( Y(:,1), Y(:,2), Y(:,3), 'red-', LineWidth=2)
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;

% Plot specific energy
figure()
plot( T, specE, 'blue', LineWidth=2)
xlabel('Time [s]'); ylabel('€ [km^2/s^2]');
title('Specific Energy');
grid on;

% Plot r, v
figure()
plot( T, Y(:,1), LineWidth=2)
hold on
plot( T, Y(:,2), LineWidth=2)
hold on
plot( T, Y(:,3), LineWidth=2)
hold on
plot( T, Rnorm, 'black', LineWidth=2)
hold on
xlabel('Time [s]'); ylabel('R [km]');
title('Distance');
grid on;
legend('rx', 'ry', 'rz', 'r')
hold off

figure()
plot( T, Y(:,4), LineWidth=2)
hold on
plot( T, Y(:,5), LineWidth=2)
hold on
plot( T, Y(:,6), LineWidth=2)
hold on
plot( T, Vnorm, 'black', LineWidth=2)
hold on
xlabel('Time [s]'); ylabel('V [km/s]');
title('Velocity');
grid on;
legend('vx', 'vy', 'vz', 'v')
hold off

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