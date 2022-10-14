clc
clear all
%%
% inertia matrix (Principal inertial axis)
I_x=0.07; %[kg*m^2]
I_y=0.0504;
I_z=0.0109;
I_rotor=0;
J=[I_x 0 0 ; 0 I_y 0; 0 0 I_z];
% external torque
M_x=0;
M_y=0;
M_z=0;
M_rotor=0;
M=[M_x M_y M_z];
% initial conditions
C=0.2; %parameter 0.2<=C<=2*pi
omega_x=C;
omega_y=0.1;
omega_z=0.1;
omega_rotor=0;
omega=[omega_x, omega_y, omega_z];

k_rotor=[0, 0, 1]; %direction of the angular momentum of the rotor

h0=[1 0 0]';

H0=[I_x*omega_x I_y*omega_y I_z*omega_z]';

%{
A0=[1 0 0; 0 1 0; 0 0 1]; %Initial cosine matrix
q=[0 0 0 1]'; %Intial quaternion
%}

A0=[I_x*omega_x 1 0; I_y*omega_y 0.15 -263.2388451; I_z*omega_z 0.35 1];
q=[0 0 0 1]';

A = sim('lab4_5','SimulationMode','normal');
b = A.get('simout');
assignin('base','b',b);

A1 = sim('lab4_5','SimulationMode','normal');
b1 = A1.get('simout1');
assignin('base','b1',b1);

A2 = sim('lab4_5','SimulationMode','normal');
b2 = A2.get('simout2');
assignin('base','b2',b2);

A3 = sim('lab4_5','SimulationMode','normal');
b3 = A3.get('simout3');
assignin('base','b3',b3);

%%
% a=vector of omega in time
a=A.simout.Data;
a1=A1.simout1.Data;
a3=A3.simout3.Data;

omega_X=a(:,1);
omega_Y=a(:,2);
omega_Z=a(:,3);

[time, omega]=size(a);
T=1:time;
%ortho_controll=zeros(time,3);
%magn_controll=zeros(time,3);
%om_matrix=zeros(3,3);
%detA=zeros(time,1);
detA1=zeros(time,1);
detA_norm=zeros(time,1);

%{
for (i=1:time)
    ortho_controll(i,1)=A0(1,:)*A0(2,:)';
    ortho_controll(i,2)=A0(2,:)*A0(3,:)';
    ortho_controll(i,3)=A0(3,:)*A0(1,:)';
    magn_controll(i,1)=norm(A0(1,:))^2;
    magn_controll(i,2)=norm(A0(2,:))^2;
    magn_controll(i,3)=norm(A0(3,:))^2;
    om_matrix=[0 -omega_Z(i,1) omega_Y(i,1); omega_Z(i,1) 0 -omega_X(i,1); -omega_Y(i,1) omega_X(i,1) 0];
    detA(i,1)=det(A0);
    A0=-om_matrix*A0;
end
%}
for (i=1:time)
    detA1(i,1)=det(a1(:,:,i));
    detA_norm(i,1)=det(a3(:,:,i));
end 
%{
figure()
plot(T, ortho_controll(:,1))
hold on
plot(T, ortho_controll(:,2))
plot(T, ortho_controll(:,3))
xlabel('t'); ylabel('-');
line([min(T) max(T)],[0 0])
title('Orthogonality');
grid on;

figure()
plot(T, magn_controll(:,1))
hold on
plot(T, magn_controll(:,2))
plot(T, magn_controll(:,3))
xlabel('t'); ylabel('-');
line([min(T) max(T)],[1 1])
title('Magnitude');
grid on;
%}

figure()
plot(T, detA_norm(:,1),'Color','b')
hold on
plot(T, detA1(:,1),'Color','r')
xlabel('t'); ylabel('-');
%line([min(T) max(T)],[1 1],'Color','b')
grid on
xlabel('Iteration [-]')
ylabel('Determinant [-]')
title('Determinant Cosine Matrix');
legend('DCM normalized','DCM not normalized')
grid on;


%{
%A0=[1 0 0; 0 1 0; 0 0 1];
for (i=1:time)
    ortho_controll(i,1)=A0(1,:)*A0(2,:)';
    ortho_controll(i,2)=A0(2,:)*A0(3,:)';
    ortho_controll(i,3)=A0(3,:)*A0(1,:)';
    magn_controll(i,1)=norm(A0(1,:))^2;
    magn_controll(i,2)=norm(A0(2,:))^2;
    magn_controll(i,3)=norm(A0(3,:))^2;
    detA(i,1)=det(A0);
    A0=A0*1.5-A0*A0'*A0*0.5;
end

figure()
plot(T, ortho_controll(:,1))
hold on
plot(T, ortho_controll(:,2))
plot(T, ortho_controll(:,3))
xlabel('t'); ylabel('-');
line([min(T) max(T)],[0 0])
title('Orthogonality');
grid on;

figure()
plot(T, magn_controll(:,1))
hold on
plot(T, magn_controll(:,2))
plot(T, magn_controll(:,3))
xlabel('t'); ylabel('-');
line([min(T) max(T)],[1 1])
title('Magnitude');
grid on;

figure()
plot(T, detA(:,1),'Color','r')
hold on
xlabel('t'); ylabel('-');
line([min(T) max(T)],[1 1],'Color','b')
title('Determinante');
grid on;
%}
%% Quaternions
a2=A2.simout2.Data;

%{
magn_controll_q=zeros(time,1);
OMEGA=zeros(4,4);
for (i=1:time)
    magn_controll_q(i,1)=norm(q)^2;
    OMEGA(1,:)=[0 omega_Z(i,1) -omega_Y(i,1) omega_X(i,1)];
    OMEGA(2,:)=[-omega_Z(i,1) 0 omega_X(i,1) omega_Y(i,1)];
    OMEGA(3,:)=[omega_Y(i,1) -omega_X(i,1) 0 omega_Z(i,1)];
    OMEGA(4,:)=[-omega_X(i,1) -omega_Y(i,1) -omega_Z(i,1) 0];
    q=0.5.*OMEGA*q;
end

figure()
plot(T, magn_controll_q(:,1))
xlabel('t'); ylabel('-');
line([min(T) max(T)],[1 1], 'Color', 'r')
title('Magnitude');
grid on;
%}
qq=zeros(time,1);
q_norm=zeros(time,1);

for (i=1:time)
   qq(i,1)=norm([a2(i,1) a2(i,2) a2(i,3) a2(i,4)]);
   q_norm(i,1)=norm([a2(i,5) a2(i,6) a2(i,7) a2(i,8)]);
end

figure()
plot(T, qq(:,1),'Color','r')
hold on
plot(T, q_norm(:,1),'Color','b')
xlabel('t'); ylabel('-');
%line([min(T) max(T)],[1 1],'Color','b')
%ylim([-10,10]);
grid on
xlabel('Iteration [-]')
ylabel('Norm [-]')
title('Quaternion norm');
legend('Not normalized','Normalized')
grid on;

%% Task2

omega_in=zeros(time, 3); %omega vector in inertial frame
h_in=zeros(time,3);
angle=zeros(time,1);

for (i=1:time)
    omega_in(i,:)=inv(a3(:,:,i))*a(i,:)';
    h_in(i,:)=[I_x*omega_X(i,1) I_y*omega_Y(i,1) I_z*omega_Z(i,1)];
    angle(i,1)=acosd((dot(h0, h_in(i,:)'))/(norm(h0)*norm(h_in(i,:))));
end 

figure()
plot(T, angle(:,1),'Color','r')
xlabel('t'); ylabel('deg');
grid on
xlabel('Iteration [-]')
ylabel('Angle [deg]')
title('Angle');
