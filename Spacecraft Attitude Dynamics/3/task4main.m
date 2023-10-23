% FIXES WILL ONLY BE APPLIED TO THE LAST SCRIPT OF THE LAB SESSION, UNLESS
% THEY ARE FUNDAMENTALLY DIFFERENT
clc
close all

stoptime = 20;

% Principal inertia, initial angle, torque
I = [0.07; 0.0109; 0.0504];
w0 = [0.01; 0.01; 2];
M = [0;0;0];

% Wheel inertia, direction, initial angle, torque
Ir = 0.005;
k = [0;0;1];
wr0 = 2*pi;
w0 = cat(1,w0,wr0);
Mr = 0.001;%0.001;

% Make matrices to operate vectorially
I = diag(I);
I_inv = pinv(I);

simulation = sim("task4.slx");
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
xlabel('Time [s]'); ylabel('w [rd/s]');
title('Angle & Angular velocity');
grid on;
legend('wx', 'wy', 'wz')%, 'wx_d', 'wy_d', 'wz_d')
hold off