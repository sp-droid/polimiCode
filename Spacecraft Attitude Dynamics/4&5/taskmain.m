clc
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.3 0.3 0.4 0.4]);

stoptime = 20;

% Spacecraft (inertia, initial angle, torque)
I = [0.07; 0.0109; 0.0504];
w0 = [0.01; 0.01; 2];
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

% Initial DCM & quaternion
A0=[1 0 0; 0 1 0; 0 0 1];
q0=[0;0;0;1];

simu = sim("task.slx");
w = simu.w;
w_d = simu.w_d;
time = simu.time;

% Calculate h and T
Idiag = diag(I);
T = 0.5*(Idiag(1)*w(1,:).^2+Idiag(2)*w(2,:).^2+Idiag(3)*w(3,:).^2);
h = sqrt((Idiag(1)*w(1,:)).^2+(Idiag(2)*w(2,:)).^2+(Idiag(3)*w(3,:)).^2);

%% Plots
figure()
plot( time, w(1,:), 'blue', LineWidth=2)
hold on
plot( time, w(2,:), 'red', LineWidth=2)
hold on
plot( time, w(3,:), 'green', LineWidth=2)
hold on
xlabel('Time [s]'); ylabel('w [rd/s]');
title('Angular velocity');
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
% Energy ellipsoid
% wx^2*Ix/2/T + wy^2*Iy/2/T + wz^2*Iz/2/T = 1
[X,Y,Z] = ellipsoid(0,0,0,sqrt(2*T(1)/Idiag(1)),sqrt(2*T(1)/Idiag(2)),sqrt(2*T(1)/Idiag(3)));
energyEllipsoid = surf(X,Y,Z,'FaceAlpha', 1, 'EdgeColor', 'none');
hold on

% Ang. momentum ellipsoid
% wx^2/h^2*Ix^2 + wy^2/h^2*Iy^2 + wz^2/h^2*Iz^2 = 1
[X,Y,Z] = ellipsoid(0,0,0,h(1)/Idiag(1),h(1)/Idiag(2),h(1)/Idiag(3));
momentEllipsoid = surf(X,Y,Z,'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on

% Plot followed path
plot3(w(1,:),w(2,:),w(3,:),'k')
axis equal
hold off