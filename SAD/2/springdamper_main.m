clc
close all

stoptime = 10;

% Mass, damping, stiffness
m = 1;
c = 1;
k = 1;
omega = 1;

%Initial conditions
x0 = 0.5;
x_p0 = 0.5;

simulation = sim("springdamper.slx");
x = simulation.x;
time = simulation.time;

figure()
plot(time, x(:,1), 'LineWidth', 2)
hold on
plot(time, x(:,2), 'LineWidth', 2)
hold on
plot(time, x(:,3), 'LineWidth', 2)
hold on
legend('Displacement', 'Velocity', 'Acceleration')
grid on
hold off