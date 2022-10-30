% FIXES WILL ONLY BE APPLIED TO THE LAST SCRIPT OF THE LAB SESSION, UNLESS
% THEY ARE FUNDAMENTALLY DIFFERENT
clc
close all

stoptime = 30;

% Inertia and initial angles
I = [0.01; 0.05; 0.07];
w0 = [0.45; 0.52; 0.55];

% We are going to be working vectorially this time
I = diag(I);
I_inv = pinv(I);

simulation = sim("task3.slx");
w = simulation.w;
w_d = simulation.w_d;
time = simulation.time;

figure()
plot( time, w(1,:), 'blue', LineWidth=2)
hold on
plot( time, w(2,:), 'red', LineWidth=2)
hold on
plot( time, w(3,:), 'green', LineWidth=2)
hold on
% plot( time, w_d(1,:), 'blue--', LineWidth=2)
% hold on
% plot( time, w_d(2,:), 'red--', LineWidth=2)
% hold on
% plot( time, w_d(3,:), 'green--', LineWidth=2)
% hold on
xlabel('Time [s]'); ylabel('w [rd], w_d [rd/s]');
title('Angle & Angular velocity');
grid on;
legend('wx', 'wy', 'wz')%, 'wx_d', 'wy_d', 'wz_d')
hold off