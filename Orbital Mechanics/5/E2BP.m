clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);

%% Inputs
muEarth = astroConstants(13);
radEarth = astroConstants(23);
J2Earth = astroConstants(9);
muMoon = astroConstants(20);
r0 = [ 7425.5;900.4;1774.2 ];
v0 = [ -2.4135;4.0763;8.032 ];
tWindow = [date2mjd2000([2022;1;1;0;0;0]);0];
norbits = 1000;
ngrid = 100000;

y0 = [ r0; v0 ];
[a,e,i,bOmega,sOmega,theta] = car2kep(r0,v0,muEarth,'rad');
kep0 = [a,e,i,bOmega,sOmega,theta]';

%% Orbit propagation
opts = struct;
opts.TinPeriods = true;
opts.RelTol = 1e-13;
opts.AbsTol = 1e-14;
opts.moonPos = @(t) relativeMoon(t, tWindow(1));
opts.muMoon = muMoon;

[Ycar, T] = timed2BP(y0, muEarth, opts, ngrid, [], norbits);

%opts.keplerian = true;
%[Y2kep, ~] = timed2BP(kep0, muEarth, opts, ngrid, [], norbits);

Ykep = zeros(ngrid,6);
%Y2car = zeros(ngrid,6);
for j=1:ngrid
    [a,e,i,bOmega,sOmega,theta] = car2kep(Ycar(j,1:3),Ycar(j,4:6),muEarth,'rad');
    Ykep(j,1) = a; Ykep(j,2) = e; Ykep(j,3) = i;
    Ykep(j,4) = bOmega; Ykep(j,5) = sOmega; Ykep(j,6) = theta;

    %[r,v] = kep2car(Y2kep(j,1),Y2kep(j,2),Y2kep(j,3),Y2kep(j,4),...
    %    Y2kep(j,5),Y2kep(j,6),muEarth,'rad');
    %Y2car(j,1) = r(1); Y2car(j,2) = r(2); Y2car(j,3) = r(3);
    %Y2car(j,4) = v(1); Y2car(j,5) = v(2); Y2car(j,6) = v(3);
end

% Unwrap keplerian 4, 6 and 6
Ykep(:,4) = unwrap(Ykep(:,4));
Ykep(:,5) = unwrap(Ykep(:,5));
Ykep(:,6) = unwrap(Ykep(:,6));

%% Plots
% Plot the orbit
figure;
plot3( Ycar(:,1), Ycar(:,2), Ycar(:,3),'LineWidth',2)
%hold on
%plot3( Y2car(:,1), Y2car(:,2), Y2car(:,3),'LineWidth',2)
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
legend('Cartesian')
axis equal;
grid on;
hold off

triPlot(T, Ykep(:,1), 'a', norbits)

triPlot(T, Ykep(:,2), 'e', norbits)

triPlot(T, Ykep(:,3), 'i', norbits)

triPlot(T, Ykep(:,4), 'bOmega', norbits)

triPlot(T, Ykep(:,5), 'sOmega', norbits)

triPlot(T, Ykep(:,6), 'theta', norbits)

% Plot percentage difference in 100T, value in 100T and value in first 10%
function triPlot(T, var, name, norbs)

% if (name=='a')
%     percentRes = abs(var-var2)/var(1);
% elseif (name=='e')
%     percentRes = abs(var-var2);
% else
%     percentRes = abs(var-var2)/2/pi;
% end
tenN = ceil(length(var)/10);
secular = movmean( var, ceil(length(var)/norbs) ); % 1 orbit
longterm = movmean( var, ceil(0.25*length(var)/norbs) ); % 0.25 orbits
shortterm = movmean( var, ceil(0.05*length(var)/norbs) ); % 0.05 orbits

% figure;
% semilogy(T, percentRes)
% xlabel('time [T]'); ylabel('(varCar - varGauss)/var0 [-]');
% title(strcat(name,{' '},'Difference'));
% grid on;

figure;
plot(T, var)
hold on
plot(T, secular)
plot(T, longterm)
plot(T, shortterm)
xlabel('time [T]'); ylabel('var value');
title(strcat(name,{' '},'Value'));
legend('Unfiltered','Secular','Long term','Short term')
grid on;
axis tight;
hold off

figure;
plot(T(1:tenN), var(1:tenN))
hold on
plot(T(1:tenN), secular(1:tenN))
plot(T(1:tenN), longterm(1:tenN))
plot(T(1:tenN), shortterm(1:tenN))
xlabel('time [T]'); ylabel('var value');
title(strcat(name,{' '},'Value'));
legend('Unfiltered','Secular','Long term','Short term')
grid on;
axis tight;
hold off
end

function r = relativeMoon(tseconds, initialMJD)
T = initialMJD + tseconds/24/3600;
[r,~] = ephMoon(T);
r = r'; % 3x1
end