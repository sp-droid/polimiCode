clc
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.3 0.3 0.4 0.4]);

%% Given data
% Spacecraft (inertia, initial angle, torque)
I = [0.01; 0.05; 0.09];
w0 = [0.001; 0.001; 0.02];
M = [0;0;0];

% Wheel (inertia, direction, initial angle, torque)
Ir = 0.005;
k = [0;0;1];
wr0 = 0;
w0 = cat(1,w0,wr0);
Mr = 0;

% Make matrices operate vectorially
I = diag(I);
I_inv = pinv(I);

%% Initial DCM, quaternion & Euler Angles of INERTIAL FRAME N
A0=[1 0 0; 0 1 0; 0 0 1];

%Compute q0
q0 = zeros(4,1);
q0(4) = 0.5*sqrt(1+trace(A0));
q0(1) = 0.25/q0(4)*(A0(2,3)-A0(3,2));
q0(2) = 0.25/q0(4)*(A0(3,1)-A0(1,3));
q0(3) = 0.25/q0(4)*(A0(1,2)-A0(2,1));

%Compute Euler Angles for 312
euler312_0 = [-atan2(A0(2,1),A0(2,2)); asin(A0(2,3)); -atan2(A0(1,3),A0(3,3))];
%Compute Euler Angles for 313
euler313_0 = [-atan2(A0(3,1),A0(3,2)); asin(A0(3,3)); atan2(A0(1,3),A0(2,3))];

%% Pointing error
gamma0 = [0;0;1];

%% DCM of MOVING REFERENCE FRAME L
height = 200; %km
mu_E = 398600; %km^3/s^2
R_E = 6371; %km

n = sqrt(mu_E/(height+R_E)^3); %rd/s
%n = 0.02;
w_LN = [0;0;n]; %Rotation of L wrt N

%A_LN defined in simulink

%% Simulation
% Params
simfile = 'lab6task2';
sim_time = 20;
max_dt = .01; %[s]
abs_tol = 1e-7;
rel_tol = 1e-7;

if bdIsLoaded('vdp')==1 %Checks if model is open so we can edit these params
    set_param(simfile, 'StopTime', num2str(sim_time),...
    'MaxStep', num2str(max_dt), 'AbsTol', num2str(abs_tol), 'RelTol', num2str(rel_tol))
end

% Simulation
simu = sim(simfile);
w = simu.w;
w_d = simu.w_d;
time = simu.time;
nsteps = length(time);
A = simu.A;
q = simu.q;
Aeuler = simu.Aeuler;
gamma = simu.gamma;
pointingError = simu.pointingError;
w_BL = simu.w_BL;
attitudeError = simu.attitudeError(1,:);

%% Post-processing
% Calculate h and T
Idiag = diag(I);
T = 0.5*(Idiag(1)*w(1,:).^2+Idiag(2)*w(2,:).^2+Idiag(3)*w(3,:).^2);
h = sqrt((Idiag(1)*w(1,:)).^2+(Idiag(2)*w(2,:)).^2+(Idiag(3)*w(3,:)).^2);

% Calculate w in the original frame, from w (body frame) and A
wInert = zeros(3,1,nsteps);
for i=1:nsteps
    wInert(:,:,i) = A(:,:,i)*w(:,:,i);
end

% Calculate w in the original frame, from w (body frame) and q
%https://gamedev.stackexchange.com/questions/28395/rotating-vector3-by-a-quaternion
wInertq = zeros(3,1,nsteps);
for i=1:nsteps
    wInertq(:,:,i) = rotQuaternion(w(:,:,i),q(:,:,i));
end

% Calculate w in the original frame, from w (body frame) and Aeuler
%https://gamedev.stackexchange.com/questions/28395/rotating-vector3-by-a-quaternion
wInertE = zeros(3,1,nsteps);
for i=1:nsteps
    wInertE(:,:,i) = Aeuler(:,:,i)*w(:,:,i);
end

%% Plots
figure()
plot( time, w(1,:), 'blue', LineWidth=2)
hold on
plot( time, w(2,:), 'red', LineWidth=2)
hold on
plot( time, w(3,:), 'green', LineWidth=2)
hold on
xlabel('Time [s]'); ylabel('w [rd/s]');
title('w in body frame');
grid on;
legend('wx', 'wy', 'wz')
hold off

figure()
plot( time, wInert(1,:), 'blue', LineWidth=2)
hold on
plot( time, wInert(2,:), 'red', LineWidth=2)
hold on
plot( time, wInert(3,:), 'green', LineWidth=2)
hold on
xlabel('Time [s]'); ylabel('w [rd/s]');
title('w in inertial frame (DCM)');
grid on;
legend('wx', 'wy', 'wz')
hold off

figure()
plot( time, wInertq(1,:), 'blue', LineWidth=2)
hold on
plot( time, wInertq(2,:), 'red', LineWidth=2)
hold on
plot( time, wInertq(3,:), 'green', LineWidth=2)
hold on
xlabel('Time [s]'); ylabel('w [rd/s]');
title('w in inertial frame (quaternion)');
grid on;
legend('wx', 'wy', 'wz')
hold off

figure()
plot( time, wInertE(1,:), 'blue', LineWidth=2)
hold on
plot( time, wInertE(2,:), 'red', LineWidth=2)
hold on
plot( time, wInertE(3,:), 'green', LineWidth=2)
hold on
xlabel('Time [s]'); ylabel('w [rd/s]');
title('w in inertial frame (Euler angles)');
grid on;
legend('wx', 'wy', 'wz')
hold off

figure()
plot(time, T)
title('Energy')

figure()
plot(time, h)
title('Angular momentum')

figure()
%Energy ellipsoid
% wx^2*Ix/2/T + wy^2*Iy/2/T + wz^2*Iz/2/T = 1
[X,Y,Z] = ellipsoid(0,0,0,sqrt(2*T(1)/Idiag(1)),sqrt(2*T(1)/Idiag(2)),sqrt(2*T(1)/Idiag(3)));
energyEllipsoid = surf(X,Y,Z,'FaceAlpha', 1, 'EdgeColor', 'none');
hold on

%Ang. momentum ellipsoid
% wx^2/h^2*Ix^2 + wy^2/h^2*Iy^2 + wz^2/h^2*Iz^2 = 1
[X,Y,Z] = ellipsoid(0,0,0,h(1)/Idiag(1),h(1)/Idiag(2),h(1)/Idiag(3));
momentEllipsoid = surf(X,Y,Z,'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on

%Plot followed path
plot3(w(1,:),w(2,:),w(3,:),'k')
axis equal
hold off

%Pointing direction
figure()
plot( time, gamma(1,:), 'blue', LineWidth=2)
hold on
plot( time, gamma(2,:), 'red', LineWidth=2)
hold on
plot( time, gamma(3,:), 'green', LineWidth=2)
hold on
xlabel('Time [s]'); ylabel('[]');
title('Pointing direction of the s/c');
grid on;
legend('gammax', 'gammay', 'gammaz')
hold off

%Pointing error
figure()
plot(time, rad2deg(pointingError))
title('Pointing error [deg]')

%Attitude error
figure()
plot(time, attitudeError)
title('Attitude error abs(trace($A_{B/L} - \mathcal{I}$))')

figure()
plot( time, w_BL(1,:), 'blue', LineWidth=2)
hold on
plot( time, w_BL(2,:), 'red', LineWidth=2)
hold on
plot( time, w_BL(3,:), 'green', LineWidth=2)
hold on
xlabel('Time [s]'); ylabel('w [rd/s]');
title('Body vs moving reference frame');
grid on;
legend('w_BLx', 'w_BLy', 'w_BLz')
hold off

%% Functions
function rotated = rotQuaternion(v,q)
%https://gamedev.stackexchange.com/questions/28395/rotating-vector3-by-a-quaternion
u = q(1:3);
s = q(4);
rotated = 2*(dot(v,u))*u +(s^2-dot(u,u))*v +2*s*cross(u,v);
end