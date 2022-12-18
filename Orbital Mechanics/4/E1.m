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
impactVector = [0;0;-9200]; % Points to the incoming asymptote from Earth
rEarth = [AU;0;0]; % Position of the Earth

massSun = muSun/G;
massEarth = muEarth/G;
impact = vecnorm(impactVector); % Impact parameter or capital delta

%% Fly-by
% Solve hyperbola
vinf = vecnorm(vinc); % Excess velocity
a = -muEarth/vinf^2;
turnAngle = atan(-a/impact)*2;
beta = (pi-turnAngle)/2;
e = 1/sin(turnAngle/2);
rper = a*(1-e); % Distance at perigee
deltaV = 2*vinf*sin(turnAngle/2);
u = normalize(cross(impactVector,vinc),'norm');

vout = vinc*cos(turnAngle)+cross(u,vinc)*sin(turnAngle)+u*dot(u,vinc)*(1-cos(turnAngle));

nodeline = normalize(vinc-vout,'norm');
rperPoint = rper*nodeline;
centerPoint = (rper-a)*nodeline;

% Max velocity
h = sqrt(a*muEarth*(1-e^2));
vper = muEarth/h*(1+e);

% Fly-by flight path
bVector = normalize(vinc+vout,'norm');
vperVector = bVector*vper;

%% Flight path
% We aren't given the speed of the Earth so I'm gonna assume it's circular
vEarth = [0;sqrt(muSun/rEarth(1));0];
vBefore = vEarth+vinc;
vAfter = vEarth+vout;

%% Plots
% Flight path near Earth
figure;
% Flight paths
time = 2000;
Yplot = timed2BP(rperPoint, vperVector, muEarth, 100, -time); % Before encounter
h1 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0, 0.4470, 0.7410],'LineWidth',3);
hold on
Yplot = timed2BP(rperPoint, vperVector, muEarth, 100, time); % After encounter
plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0, 0.4470, 0.7410],'LineWidth',3);
% Perigee
h3 = scatter3(rperPoint(1),rperPoint(2),rperPoint(3),60,'filled');
% Center
h4 = scatter3(centerPoint(1),centerPoint(2),centerPoint(3),60,'filled');
% Incoming and outgoing asymptotes
Yplot = [(centerPoint-6*radEarth*normalize(vinc,'norm')) centerPoint]';
h2 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'--','Color',[0, 0.4470, 0.7410],'LineWidth',2);
Yplot = [(centerPoint+6*radEarth*normalize(vout,'norm')) centerPoint]';
plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'--','Color',[0, 0.4470, 0.7410],'LineWidth',2);
% Node line
Yplot = 3*radEarth*[-nodeline nodeline]';
h5 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'--','Color','black','LineWidth',3);
planet3dOptions.Units = 'km';
planet3dOptions.Size = 1;
planet3D('Earth', planet3dOptions);
% Planet motion
h6 = quiver3(0,0,0,0,3*radEarth,0,'LineWidth',3,'Color','black');
title('Fly-by on Earth frame')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
lgnd = legend([h1,h2,h3,h4,h5,h6],'Hyperbolic orbit','Asymptotes','Perigee',...
        'Center','Node line','Planet motion','Location','northoutside','NumColumns',3);
lgnd.FontSize = 12;
grid on;
hold off

% Flight path
figure;
planet3dOptions.Units = 'km';
planet3dOptions.Size = 20;
planet3D('Sun', planet3dOptions);
hold on
Yplot = timed2BP(rEarth, vBefore, muSun, 100, [], -0.5); % Before encounter
h1 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0, 0.4470, 0.7410],'LineWidth',3);
Yplot = timed2BP(rEarth, vAfter, muSun, 100, [], 0.5); % After encounter
h2 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0.9290, 0.6940, 0.1250],'LineWidth',3);
title('Fly-by on Sun frame')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
lgnd = legend([h1,h2],'Before encounter','After encounter','Location','northoutside','NumColumns',2);
lgnd.FontSize = 12;
grid on;
hold off