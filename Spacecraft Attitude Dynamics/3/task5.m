% FIXES WILL ONLY BE APPLIED TO THE LAST SCRIPT OF THE LAB SESSION, UNLESS
% THEY ARE FUNDAMENTALLY DIFFERENT
clc
close all

% Principal inertia, initial angle, torque
I = [0.07; 0.0109; 0.0504];
w0 = [0.01; 0.01; 2];

T = 0.5*(I(1)*w0(1)^2+I(2)*w0(2)^2+I(3)*w0(3)^2);
h = sqrt((I(1)*w0(1))^2+(I(2)*w0(2))^2+(I(3)*w0(3))^2);

% Energy ellipsoid
% wx^2*Ix/2/T + wy^2*Iy/2/T + wz^2*Iz/2/T = 1
[X,Y,Z] = ellipsoid(0,0,0,sqrt(2*T/I(1)),sqrt(2*T/I(2)),sqrt(2*T/I(3)));
energyEllipsoid = surf(X,Y,Z,'FaceAlpha', 1, 'EdgeColor', 'none');

hold on
% Ang. momentum ellipsoid
% wx^2/h^2*Ix^2 + wy^2/h^2*Iy^2 + wz^2/h^2*Iz^2 = 1
[X,Y,Z] = ellipsoid(0,0,0,h/I(1),h/I(2),h/I(3));
momentEllipsoid = surf(X,Y,Z,'FaceAlpha', 0.5, 'EdgeColor', 'none');

axis equal
hold off