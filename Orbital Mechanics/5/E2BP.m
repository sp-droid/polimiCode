clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);

%% Inputs
muEarth = astroConstants(13);
r0 = [ 26578.137; 0; 0 ];
v0 = [ 0; 2.221; 3.173 ];
norbits = 0.8;
ngrid = 300;

y0 = [ r0; v0 ];
[a,e,i,bOmega,sOmega,theta] = car2kep(r0,v0,muEarth,'rad');
kep0 = [a,e,i,bOmega,sOmega,theta]';

%% Orbit propagation
opts = struct;
Ycar = timed2BP(y0, muEarth, opts, ngrid, [], norbits);

opts.keplerian = true;
Y2kep = timed2BP(kep0, muEarth, opts, ngrid, [], norbits);
Y2car = zeros(ngrid,6);
for i=1:ngrid
    [r,v] = kep2car(Y2kep(i,1),Y2kep(i,2),Y2kep(i,3),Y2kep(i,4),...
        Y2kep(i,5),Y2kep(i,6),muEarth,'rad');
    Y2car(i,1) = r(1); Y2car(i,2) = r(2); Y2car(i,3) = r(3);
    Y2car(i,4) = v(1); Y2car(i,5) = v(2); Y2car(i,6) = v(3);
end

%% Plots
% Plot the orbit
figure()
plot3( Ycar(:,1), Ycar(:,2), Ycar(:,3),'LineWidth',2)
hold on
plot3( Y2car(:,1), Y2car(:,2), Y2car(:,3),'LineWidth',2)
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
legend('Cartesian','Gauss equations')
axis equal;
grid on;