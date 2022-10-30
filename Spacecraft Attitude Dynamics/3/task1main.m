% FIXES WILL ONLY BE APPLIED TO THE LAST SCRIPT OF THE LAB SESSION, UNLESS
% THEY ARE FUNDAMENTALLY DIFFERENT
clc
close all

stoptime = 10;

% Inertia and initial angles
I = [0.07, 0.0504, 0.0109];
w0 = [0.45, 0.52, 0.55];

simulation = sim("task1.slx");
w = simulation.w;
w_d = simulation.w_d;
time = simulation.time;

% Calculate h and T
T = 0.5*(I(1)*w(:,1).^2+I(2)*w(:,2).^2+I(3)*w(:,3).^2);
h = sqrt((I(1)*w(:,1)).^2+(I(2)*w(:,2)).^2+(I(3)*w(:,3)).^2);

%% Plots
figure()
plot( time, w(:,1), 'blue', LineWidth=2)
hold on
plot( time, w(:,2), 'red', LineWidth=2)
hold on
plot( time, w(:,3), 'green', LineWidth=2)
hold on
plot( time, w_d(:,1), 'blue--', LineWidth=2)
hold on
plot( time, w_d(:,2), 'red--', LineWidth=2)
hold on
plot( time, w_d(:,3), 'green--', LineWidth=2)
hold on
xlabel('Time [s]'); ylabel('w [rd], w_d [rd(s)]');
title('Angular velocity & acceleration');
grid on;
legend('wx', 'wy', 'wz', 'wx_d', 'wy_d', 'wz_d')
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
[X,Y,Z] = ellipsoid(0,0,0,sqrt(2*T(1)/I(1)),sqrt(2*T(1)/I(2)),sqrt(2*T(1)/I(3)));
energyEllipsoid = surf(X,Y,Z,'FaceAlpha', 1, 'EdgeColor', 'none');
hold on

% Ang. momentum ellipsoid
% wx^2/h^2*Ix^2 + wy^2/h^2*Iy^2 + wz^2/h^2*Iz^2 = 1
[X,Y,Z] = ellipsoid(0,0,0,h(1)/I(1),h(1)/I(2),h(1)/I(3));
momentEllipsoid = surf(X,Y,Z,'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on

% Plot followed path
plot3(w(:,1),w(:,2),w(:,3),'k')
axis equal
hold off