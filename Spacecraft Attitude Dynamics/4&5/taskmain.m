clc
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.3 0.3 0.4 0.4]);

stoptime = 20;

%% Given data
% Spacecraft (inertia, initial angle, torque)
I = [0.07; 0.0504; 0.0109];
w0 = [1; 0.1; 0.1];
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

%% Initial DCM, quaternion & Euler Angles
%Other method:
%h2 = cross(H0,[0 1 0]');
%h3=cross(H0,h2);
%A_0(:,1)=H0
%A_0(:,2)=h2
%A_0(:,3)=h3

%A0=[1 0 0; 0 1 0; 0 0 1];
%q0=[0;0;0;1];

%The initial attitude parameters must be (in this example) such that
% h(0)inertial = [1;0;0]. For that, we are only permitted to change wx
h0inert = [1;0;0];

%Compute wx such that hbody has the same norm as hinertial
w0(1) = sqrt(1-(I(2,2)*w0(2))^2+(I(3,3)*w0(3))^2)/I(1,1);

%Compute A0
h0 = I*w0(1:3);
axisangle = vrrotvec(h0,h0inert);   %Axis-angle that rotates h0 to h0inert
A0 = vrrotvec2mat(axisangle);       %DCM that fullfills that rotation

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
gamma0 = [1;0;0];

%% Simulation
simu = sim("task.slx");
w = simu.w;
w_d = simu.w_d;
time = simu.time;
nsteps = length(time);
A = simu.A;
AnonNorm = simu.nonNormA;
q = simu.q;
qnonNorm = simu.nonNormq;
Aeuler = simu.Aeuler;
gamma = simu.gamma;
pointingError = simu.pointingError;

%% Post-processing
% Calculate h and T
Idiag = diag(I);
T = 0.5*(Idiag(1)*w(1,:).^2+Idiag(2)*w(2,:).^2+Idiag(3)*w(3,:).^2);
h = sqrt((Idiag(1)*w(1,:)).^2+(Idiag(2)*w(2,:)).^2+(Idiag(3)*w(3,:)).^2);

% Calculate w in the original frame, from w (body frame) and A
wInert = zeros(3,1,nsteps);
wInertnonNorm = zeros(3,1,nsteps);
for i=1:nsteps
    wInert(:,:,i) = A(:,:,i)*w(:,:,i);
    wInertnonNorm(:,:,i) = AnonNorm(:,:,i)*w(:,:,i);
end

%Calculate degree of non-orthonormality in AnonNorm
AnonNormvalue = zeros(1,nsteps);
temp = permute(AnonNorm,[2 1 3]);
Identity = eye(3);
for i=1:nsteps
    AnonNormvalue(i) = norm(Identity-temp(:,:,i)*AnonNorm(:,:,i));
end

% Calculate w in the original frame, from w (body frame) and q
%https://gamedev.stackexchange.com/questions/28395/rotating-vector3-by-a-quaternion
wInertq = zeros(3,1,nsteps);
wInertnonNormq = zeros(3,1,nsteps);
for i=1:nsteps
    wInertq(:,:,i) = rotQuaternion(w(:,:,i),q(:,:,i));
    wInertnonNormq(:,:,i) = rotQuaternion(w(:,:,i),qnonNorm(:,:,i));
end

%Calculate degree of non-normality in qnonNorm
qnonNormvalue = zeros(1,nsteps);
for i=1:nsteps
    qnonNormvalue(i) = norm(qnonNorm(:,:,i));
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
title('Angular velocity in body frame');
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
title('Angular velocity in inertial frame (DCM)');
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
title('Angular velocity in inertial frame (quaternion)');
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
title('Angular velocity in inertial frame (Euler angles)');
grid on;
legend('wx', 'wy', 'wz')
hold off

figure()
plot( time, abs(wInert(3,:)-wInertnonNorm(3,:)), 'blue', LineWidth=2)
hold on
plot( time, abs(wInert(2,:)-wInertnonNorm(2,:)), 'red', LineWidth=2)
hold on
plot( time, abs(wInert(1,:)-wInertnonNorm(1,:)), 'green', LineWidth=2)
hold on
xlabel('Time [s]'); ylabel('w [rd/s]');
title('DCM vs non normalized DCM effect on wInertial');
grid on;
legend('wzDiff','wyDiff','wxDiff')
hold off

figure()
plot( time, AnonNormvalue, LineWidth=2)
xlabel('Time [s]'); ylabel('Norm(I-A.T*A)');
title('Non normalized DCM divergence from orthonormality');
grid on;

figure()
plot( time, abs(wInertq(3,:)-wInertnonNormq(3,:)), 'blue', LineWidth=2)
hold on
plot( time, abs(wInertq(2,:)-wInertnonNormq(2,:)), 'red', LineWidth=2)
hold on
plot( time, abs(wInertq(1,:)-wInertnonNormq(1,:)), 'green', LineWidth=2)
hold on
xlabel('Time [s]'); ylabel('w [rd/s]');
title('q vs non normalized q effect on wInertial');
grid on;
legend('wzDiff','wyDiff','wxDiff')
hold off

figure()
plot( time, qnonNormvalue, LineWidth=2)
xlabel('Time [s]'); ylabel('Norm(q)');
title('Non normalized q divergence from orthonormality');
grid on;

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

%% Functions
function rotated = rotQuaternion(v,q)
%https://gamedev.stackexchange.com/questions/28395/rotating-vector3-by-a-quaternion
u = q(1:3);
s = q(4);
rotated = 2*(dot(v,u))*u +(s^2-dot(u,u))*v +2*s*cross(u,v);
end