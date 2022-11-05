clc
clear
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

%simu.A_LN defined in simulink

%% Simulation
% Params
simfile = 'lab7task1';
sim_time = 2*pi/n;
%max_dt = .01; %[s] 'MaxStep', num2str(max_dt)
abs_tol = 1e-7;
rel_tol = 1e-7;

if bdIsLoaded(simfile)==1 %Checks if model is open so we can edit these params
    set_param(simfile, 'StopTime', 'sim_time',...
    'AbsTol', num2str(abs_tol), 'RelTol', num2str(rel_tol))
end

% Simulation
simu = sim(simfile);
time = simu.time;
nsteps = length(time);

%% Post-processing
% Calculate h and T
Idiag = diag(I);
T = 0.5*(Idiag(1)*simu.w(1,:).^2+Idiag(2)*simu.w(2,:).^2+Idiag(3)*simu.w(3,:).^2);
h = sqrt((Idiag(1)*simu.w(1,:)).^2+(Idiag(2)*simu.w(2,:)).^2+(Idiag(3)*simu.w(3,:)).^2);

%% Plots
figure()
for i=1:3
    for j=1:3
        plot(time,reshape(simu.A_BN(i,j,:),1,[]))
        hold on
    end
end
xlabel('t [s]');
title('$A^{BN}$');
grid on; axis tight;
hold off

figure()
plot( time, simu.w(1,:), 'blue', LineWidth=2)
hold on
plot( time, simu.w(2,:), 'red', LineWidth=2)
hold on
plot( time, simu.w(3,:), 'green', LineWidth=2)
hold on
xlabel('t [s]'); ylabel('$\omega$ [rd/s]');
title('$\omega^B$');
grid on; axis tight;
legend('\omega^B_x', '\omega^B_y', '\omega^B_z')
hold off

kinematicMethod = 'DCM';
switch kinematicMethod
    case 'DCM'
        % Calculate simu.w in the original frame, from simu.w (body frame) and simu.A_BN
        wInert = zeros(3,1,nsteps);
        for i=1:nsteps
            wInert(:,:,i) = simu.A_BN(:,:,i)*simu.w(:,:,i);
        end

        figure()
        plot( time, wInert(1,:), 'blue', LineWidth=2)
        hold on
        plot( time, wInert(2,:), 'red', LineWidth=2)
        hold on
        plot( time, wInert(3,:), 'green', LineWidth=2)
        hold on
        xlabel('t [s]'); ylabel('$\omega$ [rd/s]');
        title('$\omega^{BN}$ (DCM)');
        grid on; axis tight;
        legend('\omega^{BN}_x', '\omega^{BN}_y', '\omega^{BN}_z')
        hold off
    case 'quaternion'
        % Calculate simu.w in the original frame, from simu.w (body frame) and simu.q
        %https://gamedev.stackexchange.com/questions/28395/rotating-vector3-by-a-quaternion
        wInertq = zeros(3,1,nsteps);
        for i=1:nsteps
            wInertq(:,:,i) = rotQuaternion(simu.w(:,:,i),simu.q(:,:,i));
        end

        figure()
        plot( time, wInertq(1,:), 'blue', LineWidth=2)
        hold on
        plot( time, wInertq(2,:), 'red', LineWidth=2)
        hold on
        plot( time, wInertq(3,:), 'green', LineWidth=2)
        hold on
        xlabel('t [s]'); ylabel('$\omega$ [rd/s]');
        title('$\omega^{BN}$ (quaternion)');
        grid on; axis tight;
        legend('\omega^{BN}_x', '\omega^{BN}_y', '\omega^{BN}_z')
        hold off
    case 'eulerAngles'
        % Calculate simu.w in the original frame, from simu.w (body frame) and simu.Aeuler
        %https://gamedev.stackexchange.com/questions/28395/rotating-vector3-by-a-quaternion
        wInertE = zeros(3,1,nsteps);
        for i=1:nsteps
            wInertE(:,:,i) = simu.Aeuler(:,:,i)*simu.w(:,:,i);
        end

        figure()
        plot( time, wInertE(1,:), 'blue', LineWidth=2)
        hold on
        plot( time, wInertE(2,:), 'red', LineWidth=2)
        hold on
        plot( time, wInertE(3,:), 'green', LineWidth=2)
        hold on
        xlabel('t [s]'); ylabel('$\omega$ [rd/s]');
        title('$\omega^{BN}$ (Euler angles)');
        grid on; axis tight;
        legend('\omega^{BN}_x', '\omega^{BN}_y', '\omega^{BN}_z')
        hold off
end

figure()
plot(time, T)
xlabel('t [s]');
grid on; axis tight;
title('Energy')

figure()
plot(time, h)
xlabel('t [s]');
grid on; axis tight;
title('Angular momentum')

plotEllipsoid = 0;
if plotEllipsoid==1
    % Energy ellipsoid, ang. momentum ellipsoid, followed path
    figure()
    % wx^2*Ix/2/T + wy^2*Iy/2/T + wz^2*Iz/2/T = 1
    [X,Y,Z] = ellipsoid(0,0,0,sqrt(2*T(1)/Idiag(1)),sqrt(2*T(1)/Idiag(2)),sqrt(2*T(1)/Idiag(3)));
    energyEllipsoid = surf(X,Y,Z,'FaceAlpha', 1, 'EdgeColor', 'none');
    hold on
    
    % wx^2/h^2*Ix^2 + wy^2/h^2*Iy^2 + wz^2/h^2*Iz^2 = 1
    [X,Y,Z] = ellipsoid(0,0,0,h(1)/Idiag(1),h(1)/Idiag(2),h(1)/Idiag(3));
    momentEllipsoid = surf(X,Y,Z,'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on
    
    plot3(simu.w(1,:),simu.w(2,:),simu.w(3,:),'k')
    axis equal
    hold off
end

plotPointing = 0;
if plotPointing==1
    %Pointing direction
    figure()
    plot( time, simu.gamma(1,:), 'blue', LineWidth=2)
    hold on
    plot( time, simu.gamma(2,:), 'red', LineWidth=2)
    hold on
    plot( time, simu.gamma(3,:), 'green', LineWidth=2)
    hold on
    xlabel('t [s]'); ylabel('$\Gamma$ [-]');
    title('Pointing direction: $\Gamma^B$');
    grid on; axis tight;
    legend('\Gamma^B_x', '\Gamma^B_y', '\Gamma^B_z')
    hold off
    
    %Pointing error
    figure()
    plot(time, rad2deg(simu.pointingError(:)))
    xlabel('t [s]'); ylabel('$\theta$ [deg]');
    grid on; axis tight;
    title('Pointing error $\theta = \arccos(\Gamma \cdot \Gamma_0)$')
end

%Attitude error
figure()
plot(time, simu.attitudeError(:))
xlabel('t [s]'); 
grid on; axis tight;
title('Attitude error: trace($A_{B/L} - \mathcal{I}$)')

figure()
for i=1:3
    for j=1:3
        plot(time,reshape(simu.A_BL(i,j,:),1,[]))
        hold on
    end
end
xlabel('t [s]');
title('$A^{BL}$');
grid on; axis tight;
hold off

figure()
plot(time, simu.w_BL(1,:), 'blue', LineWidth=2)
hold on
plot(time, simu.w_BL(2,:), 'red', LineWidth=2)
hold on
plot(time, simu.w_BL(3,:), 'green', LineWidth=2)
hold on
xlabel('t [s]'); ylabel('$\omega$ [rd/s]');
title('$\omega^{BL}$');
grid on; axis tight;
legend('\omega^{BL}_x', '\omega^{BL}_y', '\omega^{BL}_z')
hold off

figure()
for i=1:3
    for j=1:3
        plot(time,reshape(simu.A_LN(i,j,:),1,[]))
        hold on
    end
end
xlabel('t [s]');
title('$A^{LN}$');
grid on; axis tight;
hold off

%Perturbations
figure()
plot(time, simu.Mtotal(1,:), 'blue', LineWidth=2)
hold on
plot(time, simu.Mtotal(2,:), 'red', LineWidth=2)
hold on
plot(time, simu.Mtotal(3,:), 'green', LineWidth=2)
hold on
xlabel('t [s]'); ylabel('M');
title('External torque');
grid on; axis tight;
legend('M_x', 'M_y', 'M_z')
hold off

figure()
plot(time, simu.Mgg(1,:), 'blue', LineWidth=2)
hold on
plot(time, simu.Mgg(2,:), 'red', LineWidth=2)
hold on
plot(time, simu.Mgg(3,:), 'green', LineWidth=2)
hold on
xlabel('t [s]'); ylabel('M');
title('Gravity gradient perturbation');
grid on; axis tight;
legend('M^{GG}_x', 'M^{GG}_y', 'M^{GG}_z')
hold off

%% Functions
function rotated = rotQuaternion(v,q)
%https://gamedev.stackexchange.com/questions/28395/rotating-vector3-by-a-quaternion
u = q(1:3);
s = q(4);
rotated = 2*(dot(v,u))*u +(s^2-dot(u,u))*v +2*s*cross(u,v);
end