clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);

%% Input data
G = astroConstants(1);
radEarth = astroConstants(23);

muEarth = astroConstants(13);
muSun = astroConstants(4);
AU = astroConstants(2);
vinc = [15.1;0;0]; % Incoming velocity at Earth's SOI's boundary
%impact = [0;-9200;0]; % Points to the incoming asymptote from Earth
rEarth = [AU;0;0]; % Position of the Earth

massSun = muSun/G;
massEarth = muEarth/G;
vEarth = [0;sqrt(muSun/rEarth(1));0]; % Assuming circular orbit
vBefore = vEarth+vinc;
%%
ngrid = 5;

%% Fly-by grid plots
impact = linspace(9200, 13200, ngrid);

% Near Earth
figure;
planet3dOptions.Units = 'km';
planet3dOptions.Size = 1;
planet3D('Earth', planet3dOptions);
hold on

hx = gobjects(ngrid);
hxlgnd = strings(ngrid);
for i=1:ngrid
    impactVector = [0;impact(i);0];
    [vout, rper, vper, turnAngle] = hyperb2BP(vinc, impactVector, muEarth);

    vAfter = vEarth+vout;

    time = 2000;
    Yplot = timed2BP(rper, vper, muEarth, 100, [-time,time]);
    hx(i) = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'LineWidth',3);
    scatter3(rper(1),rper(2),rper(3),60,'filled');
    hxlgnd(i) = strcat('Impact =',{' '},num2str(impact(i)),{' '},'[km]');
end

h6 = quiver3(0,0,0,0,3*radEarth,0,'LineWidth',3,'Color','black');
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Fly-by near Earth')
lgnd = legend([h6,hx(1),hx(2),hx(3),hx(4),hx(5)],'Planet motion',...
    hxlgnd(1),hxlgnd(2),hxlgnd(3),hxlgnd(4),hxlgnd(5));
lgnd.FontSize = 12;
grid on;
hold off

% Flight path plot
figure;
planet3dOptions.Units = 'km';
planet3dOptions.Size = 20;
planet3D('Sun', planet3dOptions);
hold on
Yplot = timed2BP(rEarth, vBefore, muSun, 100, [], -0.5); % Before encounter
h1 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0, 0.4470, 0.7410],'LineWidth',3);

hx = gobjects(ngrid);
hxlgnd = strings(ngrid);
for i=1:ngrid
    impactVector = [0;impact(i);0];
    [vout, rper, vper, turnAngle] = hyperb2BP(vinc, impactVector, muEarth);

    vAfter = vEarth+vout;

    Yplot = timed2BP(rEarth, vAfter, muSun, 100, [], 0.5); % After encounter
    hx(i) = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'LineWidth',3);
    hxlgnd(i) = strcat('Impact =',{' '},num2str(impact(i)),{' '},'[km]');
end

title('Fly-by on Sun frame')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
lgnd = legend([h1,hx(1),hx(2),hx(3),hx(4),hx(5)],'Before encounter',...
    hxlgnd(1),hxlgnd(2),hxlgnd(3),hxlgnd(4),hxlgnd(5));
lgnd.FontSize = 12;
grid on;
hold off

%% Functions;
function [vout, rper, vper, turnAngle] = hyperb2BP(vinc, impact, mu)

% Solve hyperbola
vinf = vecnorm(vinc); % Excess velocity
a = -mu/vinf^2;
turnAngle = atan(-a/vecnorm(impact))*2;
%beta = (pi-turnAngle)/2;
e = 1/sin(turnAngle/2);
rper = a*(1-e); % Distance at perigee
%deltaV = 2*vinf*sin(turnAngle/2);
u = normalize(cross(impact,vinc),'norm');

vout = vinc*cos(turnAngle)+cross(u,vinc)*sin(turnAngle)+u*dot(u,vinc)*(1-cos(turnAngle));

nodeline = normalize(vinc-vout,'norm');
%centerPoint = (rper-a)*nodeline;
rper = rper*nodeline;

% Max velocity
h = sqrt(a*mu*(1-e^2));
vper = mu/h*(1+e);

% Fly-by flight path
bVector = normalize(vinc+vout,'norm');
vper = bVector*vper;

end