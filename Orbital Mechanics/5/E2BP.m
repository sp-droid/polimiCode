clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);

%% Inputs
muEarth = astroConstants(13);
radEarth = astroConstants(23);
J2Earth = astroConstants(9);
r0 = [ 7495.3;0.0000;0.0000 ];
v0 = [ 0.0000;0.2686;-7.3239 ];
norbits = 100;
ngrid = 10000;

y0 = [ r0; v0 ];
[a,e,i,bOmega,sOmega,theta] = car2kep(r0,v0,muEarth,'rad');
kep0 = [a,e,i,bOmega,sOmega,theta]';

%% Orbit propagation
opts = struct;
opts.TinPeriods = true;
opts.RelTol = 1e-13;
opts.AbsTol = 1e-14;
opts.J2 = J2Earth;
opts.R = radEarth;

[Ycar, T] = timed2BP(y0, muEarth, opts, ngrid, [], norbits);

opts.keplerian = true;
[Y2kep, ~] = timed2BP(kep0, muEarth, opts, ngrid, [], norbits);

Ykep = zeros(ngrid,6);
Y2car = zeros(ngrid,6);
for j=1:ngrid
    [a,e,i,bOmega,sOmega,theta] = car2kep(Ycar(j,1:3),Ycar(j,4:6),muEarth,'rad');
    Ykep(j,1) = a; Ykep(j,2) = e; Ykep(j,3) = i;
    Ykep(j,4) = bOmega; Ykep(j,5) = sOmega; Ykep(j,6) = theta;

    [r,v] = kep2car(Y2kep(j,1),Y2kep(j,2),Y2kep(j,3),Y2kep(j,4),...
        Y2kep(j,5),Y2kep(j,6),muEarth,'rad');
    Y2car(j,1) = r(1); Y2car(j,2) = r(2); Y2car(j,3) = r(3);
    Y2car(j,4) = v(1); Y2car(j,5) = v(2); Y2car(j,6) = v(3);
end

%% Plots
% Plot the orbit
figure;
plot3( Ycar(:,1), Ycar(:,2), Ycar(:,3),'LineWidth',2)
hold on
plot3( Y2car(:,1), Y2car(:,2), Y2car(:,3),'LineWidth',2)
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
legend('Cartesian','Gauss equations')
axis equal;
grid on;

triPlot(T, Ykep(:,1), Y2kep(:,1), 'a')

triPlot(T, Ykep(:,2), Y2kep(:,2), 'e')

triPlot(T, Ykep(:,3), Y2kep(:,3), 'i')

triPlot(T, Ykep(:,4), Y2kep(:,4), 'bOmega')

triPlot(T, Ykep(:,5), Y2kep(:,5), 'sOmega')

triPlot(T, Ykep(:,6), Y2kep(:,6), 'theta')

% Plot percentage difference in 100T, value in 100T and value in first 10%
function triPlot(T, var1, var2, name)

if (name=='a')
    percentRes = abs(var1-var2)/var1(1);
elseif (name=='e')
    percentRes = abs(var1-var2);
else
    percentRes = abs(var1-var2)/2/pi;
end
tenN = ceil(length(var1)/10);
secular = movmean( var1, 100 );

figure;
semilogy(T, percentRes)
xlabel('time [T]'); ylabel('(varCar - varGauss)/var0 [-]');
title(strcat(name,{' '},'Difference'));
grid on;

figure;
plot(T, var1)
hold on
plot(T, var2)
plot(T, secular)
xlabel('time [T]'); ylabel('var value');
title(strcat(name,{' '},'Value'));
legend('Cartesian','Gauss equations','Secular')
grid on;
axis tight;
hold off

figure;
plot(T(1:tenN), var1(1:tenN))
hold on
plot(T(1:tenN), var2(1:tenN))
plot(T(1:tenN), secular(1:tenN))
xlabel('time [T]'); ylabel('var value');
title(strcat(name,{' '},'Value'));
legend('Cartesian','Gauss equations','Secular')
grid on;
axis tight;
hold off
end