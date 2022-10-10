close all
% J2-perturbed 2BP

% Physical parameters
%astroConstants(13) outputs the Earth's gravitational parameter [km^3/s^2]
mu_E = astroConstants(13);
J2_E = astroConstants(9);
R_E = astroConstants(23);

% Initial condition 1 (circular orbit)
%r0 = [ 26578.137; 0; 0 ]; % [km]
%v0 = [ 0; 2.221; 3.173 ]; % [km/s]
%y0 = [ r0; v0 ];

% Initial condition 2
r0 = [ 6495; -970; -3622]; % [km]
v0 = [ 4.752; 2.130; 7.950]; % [km/s]
y0 = [ r0; v0 ];

% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
Torb = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
norbits = 800;
tspan = linspace( 0, norbits*Torb, 30000 );

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

%% Compute unperturbed
% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E, 0, R_E), tspan, y0, options );

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

%% Compute perturbed
% Perform the integration
[ ~, Y_p ] = ode113( @(t,y) ode_2bp(t,y,mu_E, J2_E, R_E), tspan, y0, options );

Rnorm_p = vecnorm(Y_p(:,1:3).');
Vnorm_p = vecnorm(Y_p(:,4:6).');

h_p = cross(Y_p(:,1:3),Y_p(:,4:6));
hnorm_p = vecnorm(h_p.');

e_p = cross(Y_p(:,4:6),h_p)/mu_E-Y_p(:,1:3)./Rnorm_p.';
enorm_p = vecnorm(e_p.');

ehCheck_p = dot(h_p',e_p');

% Specific energy
specE_p = Vnorm_p.^2/2-mu_E./Rnorm_p;

% Radial and transversal velocity
vr_p = dot(Y_p(:,4:6)',Y_p(:,1:3)')./Rnorm_p;
vt_p = sqrt(Vnorm_p.^2-vr_p.^2);

%% Plots
% Plot the orbit
figure()
plot3( Y(:,1), Y(:,2), Y(:,3), 'red-', LineWidth=2)
hold on
scatter3( Y_p(:,1), Y_p(:,2), Y_p(:,3), 1, Norbit)
hold on
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
cbar = colorbar;
cbar.Title.String = 'Orbit period';
axis equal;
grid on;
legend('2BP, J2-perturbed', 'Location','north');
hold off

% Plot specific energy
figure()
plot( Norbit, specE, 'blue', LineWidth=2)
hold on
plot( Norbit, specE_p, 'red--', LineWidth=2)
hold on
xlabel('Orbit number'); ylabel('€ [km^2/s^2]');
title('Specific Energy');
legend('2BP', 'J2-perturbed');
grid on;
hold off

% Plot h, e
figure()
plot( Norbit, h(:,1), 'blue', LineWidth=2)
hold on
plot( Norbit, h(:,2), 'red', LineWidth=2)
hold on
plot( Norbit, h(:,3), 'green', LineWidth=2)
hold on
plot( Norbit, hnorm, 'black', LineWidth=2)
hold on
plot( Norbit, h_p(:,1), 'blue--', LineWidth=2)
hold on
plot( Norbit, h_p(:,2), 'red--', LineWidth=2)
hold on
plot( Norbit, h_p(:,3), 'green--', LineWidth=2)
hold on
plot( Norbit, hnorm_p, 'black--', LineWidth=2)
hold on
xlabel('Orbit number'); ylabel('h [km^2/s]');
title('Angular momentum');
grid on;
legend('hx', 'hy', 'hz', 'h', 'hx_p', 'hy_p', 'hz_p', 'h_p')
hold off

figure()
plot( Norbit, e(:,1), 'blue', LineWidth=2)
hold on
plot( Norbit, e(:,2), 'red', LineWidth=2)
hold on
plot( Norbit, e(:,3), 'green', LineWidth=2)
hold on
plot( Norbit, enorm, 'black', LineWidth=2)
hold on
plot( Norbit, e_p(:,1), 'blue--', LineWidth=2)
hold on
plot( Norbit, e_p(:,2), 'red--', LineWidth=2)
hold on
plot( Norbit, e_p(:,3), 'green--', LineWidth=2)
hold on
plot( Norbit, enorm_p, 'black--', LineWidth=2)
hold on
xlabel('Orbit number'); ylabel('e [-]');
title('Eccentricity');
grid on;
legend('ex', 'ey', 'ez', 'e', 'ex_p', 'ey_p', 'ez_p', 'e_p')
hold off

% Plot e-h dot product
figure()
plot( Norbit, ehCheck, 'blue', LineWidth=2)
hold on
plot( Norbit, ehCheck_p, 'red--', LineWidth=2)
xlabel('Orbit number'); ylabel('e*h [km^2/s]');
title('e-h dot product');
grid on;
legend('2BP', 'J2-perturbed')
hold off

% Plot vr & vt
figure()
plot( Norbit, vr, 'blue--', LineWidth=2)
hold on
plot( Norbit, vt, 'red--', LineWidth=2)
hold on
plot( Norbit, vr_p, 'green:', LineWidth=2)
hold on
plot( Norbit, vt_p, 'black:', LineWidth=2)
hold on
xlabel('Orbit number'); ylabel('v [km/s]');
title('Radial and transversal velocity');
grid on;
legend('vr', 'vt', 'vr_p', 'vt_p')
hold off

%% Functions
% Function ode_2bp
function dy = ode_2bp( ~, y, mu, J2, R )
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

% aJ2 term
aJ2 = 5*r(3)^2/rnorm^2;
aJ2 = [r(1)/rnorm*(aJ2-1)
    r(2)/rnorm*(aJ2-1)
    r(3)/rnorm*(aJ2-3)];
aJ2 = aJ2*1.5*J2*mu*R^2/rnorm^4;

% Set the derivatives of the state
dy = [ v
(-mu/rnorm^3)*r+aJ2 ];
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