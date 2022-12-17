clear
close all

%% Input data
G = astroConstants(1);
radEarth = astroConstants(23);

muEarth = astroConstants(13);
muSun = astroConstants(4);
AU = astroConstants(2);
vinc = [15.1;0;0]; % Incoming velocity at Earth's SOI's boundary
impactVector = [0;9200;0]; % Points to the incoming asymptote from Earth
rEarth = [AU;0;0]; % Position of the Earth

massSun = muSun/G;
massEarth = muEarth/G;
impact = vecnorm(impactVector); % Impact parameter or capital delta

%% Fly-by
% Solve hyperbola
vinf = vecnorm(vinc); % Excess velocity
a = -muEarth/vinf^2;
turnAngle = atan(-a/impact)*2;
e = 1/sin(turnAngle/2);
rper = a*(1-e); % Distance at perigee
deltaV = 2*vinf*sin(turnAngle/2);
u = normalize(cross(impactVector,vinc),'norm');

vout = vinc*cos(turnAngle)+cross(u,vinc)*sin(turnAngle)+u*dot(u,vinc)*(1-cos(turnAngle));

% Earth SOI and orbit nodeline for the perigee point
% Revisar esto
planetoSOIearth = radEarth*(massSun/massEarth)^(2/5);
nodeline = normalize(vinc-vout,'norm');
rperVector = rper*nodeline;

% Max velocity
h = sqrt(a*muEarth*(1-e^2));
vper = muEarth/h*(1+e);

% Path from SOI's boundary
thetainc = -acos((h^2/muEarth/planetoSOIearth-1)/e);
rinc = planetoSOIearth*[cos(thetainc);sin(thetainc);0];
%rad2deg(atan2(e*sin(-thetainc),1+e*cos(-thetainc)))

%% Flight path
% We aren't given the speed of the Earth so I'm gonna assume it's circular
vEarth = [0;sqrt(muSun/rEarth(1));0];
vBefore = vEarth+vinc;
vAfter = vEarth+vout;

%% Plots
% Flight path near Earth
figure;
%Yplot = timed2BP(rEarth, vBefore, muSun, 100, [], -0.5); % Before encounter
%h1 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0, 0.4470, 0.7410],'LineWidth',3);
%hold on
h2 = scatter3(rperVector(1),rperVector(2),rperVector(3),'Color',[0.9290, 0.6940, 0.1250],'LineWidth',3);
title('Fly-by on Earth frame')
%lgnd = legend([h1,h2],'Before encounter','After encounter','Location','northoutside','NumColumns',2);
%lgnd.FontSize = 12;
grid on;
hold off

% Flight path
figure;
Yplot = timed2BP(rEarth, vBefore, muSun, 100, [], -0.5); % Before encounter
h1 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0, 0.4470, 0.7410],'LineWidth',3);
hold on
Yplot = timed2BP(rEarth, vAfter, muSun, 100, [], 0.5); % After encounter
h2 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0.9290, 0.6940, 0.1250],'LineWidth',3);
title('Fly-by on Sun frame')
lgnd = legend([h1,h2],'Before encounter','After encounter','Location','northoutside','NumColumns',2);
lgnd.FontSize = 12;
grid on;
hold off