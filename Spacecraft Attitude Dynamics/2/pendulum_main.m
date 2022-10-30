clc
close all

stoptime = 10;
d_theta0 = 0;
theta0 = 1.57;

simulation = sim("pendulum.slx");
out = simulation.simout;

figure()
plot(out(:,1), out(:,2), 'LineWidth', 2)
hold on
plot(out(:,1), out(:,3), 'LineWidth', 2)
hold on
plot(out(:,1), out(:,4), 'LineWidth', 2)
hold on
legend('angle', 'ang. velocity', 'ang. acceleration')
grid on
hold off