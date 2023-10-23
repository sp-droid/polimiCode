clear
close all
set(groot, 'defaultFigureUnits', 'normalized', 'defaultFigurePosition', [0.2 0.2 0.6 0.6]);

%% Input data
G = astroConstants(1);
radEarth = astroConstants(23);

muEarth = astroConstants(13);
muSun = astroConstants(4);
AU = astroConstants(2);
rEarth = [0;-AU;0]; % Position of the Earth

massSun = muSun/G;
massEarth = muEarth/G;
vEarth = [sqrt(muSun/AU);0;0]; % Assuming circular orbit

vBefore = [ 31.5; 5.2; 0.0 ];
vAfter = [ 36.0; 0.0; 0.0 ];

vinc = vBefore-vEarth;
vout = vAfter-vEarth;
%% Fly-by

[deltaVper, rper, vpers, centers, turnAngle] = gravAssist(vinc, vout, muEarth, radEarth);

%% Plots
% Flight path near Earth
figure;
% Flight paths
time = 4000;
opts.show = true;
Yplot = timed2BP([rper;vpers(:,1)], muEarth, opts, 100, -time); % Before encounter
h1 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'LineWidth',3);
hold on
Yplot = timed2BP([rper;vpers(:,2)], muEarth, [], 100, time); % After encounter
h2 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'LineWidth',3);
% Perigee
h3 = scatter3(rper(1),rper(2),rper(3),60,'filled');
% Centers
scatter3(centers(1,1),centers(2,1),centers(3,1),60,'filled');
scatter3(centers(1,2),centers(2,2),centers(3,2),60,'filled');
% Incoming and outgoing asymptotes
Yplot = [(centers(:,1)-6*radEarth*normalize(vinc,'norm')) centers(:,1)]';
h4 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'--','LineWidth',2);
Yplot = [(centers(:,2)+6*radEarth*normalize(vout,'norm')) centers(:,2)]';
h5 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'--','LineWidth',2);
% Node line
Yplot = 4*radEarth*[-normalize(rper,'norm') normalize(rper,'norm')]';
plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'--','Color','black','LineWidth',3);
planet3dOptions.Units = 'km';
planet3dOptions.Size = 1;
planet3D('Earth', planet3dOptions);
% Planet motion
h6 = quiver3(0,0,0,0,3*radEarth,0,'LineWidth',3,'Color','black');
title('Fly-by on Earth frame')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
lgnd = legend([h1,h2,h3,h4,h5,h6],'Incoming hyperbola','Outgoing hyperbola','Perigee',...
        'Incoming asymptote','Outgoing asymptote','Planet motion',...
        'Location','northoutside','NumColumns',3);
lgnd.FontSize = 12;
grid on;
hold off

% Flight path
figure;
planet3dOptions.Units = 'km';
planet3dOptions.Size = 20;
planet3D('Sun', planet3dOptions);
hold on
Yplot = timed2BP([rEarth;vBefore], muSun, [], 100, [], -0.5); % Before encounter
h1 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0, 0.4470, 0.7410],'LineWidth',3);
Yplot = timed2BP([rEarth;vAfter], muSun, [], 100, [], 0.5); % After encounter
h2 = plot3(Yplot(:,1), Yplot(:,2), Yplot(:,3),'Color',[0.9290, 0.6940, 0.1250],'LineWidth',3);
title('Fly-by on Sun frame')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
lgnd = legend([h1,h2],'Before encounter','After encounter','Location','northoutside','NumColumns',2);
lgnd.FontSize = 12;
grid on;
hold off

%% Functions;
function [deltaVper, rper, vpers, centers, turnAngle] = gravAssist(vinc, vout, mu, iGuess)

% Total turn angle
turnAngle = atan2(norm(cross(vout,vinc)),dot(vout,vinc));

% Define function to solve for
vi = vecnorm(vinc); vo = vecnorm(vout);
func = @(r) turnAngle-asin(1/(1+r*vi^2/mu))-asin(1/(1+r*vo^2/mu));

% Solve
rper = fzero(func, iGuess);

as = [-mu/vi^2; -mu/vo^2];
es = [1+rper*vi^2/mu; 1+rper*vo^2/mu];

nodeline = normalize(vinc-vout,'norm');
centers = [(rper-as(1))*nodeline (rper-as(2))*nodeline];
rper = rper*nodeline;

% Max velocity
hs = [sqrt(as(1)*mu*(1-es(1)^2)); sqrt(as(2)*mu*(1-es(2)^2))];
vpers = [mu/hs(1)*(1+es(1)); mu/hs(2)*(1+es(2))];
deltaVper = vpers(2)-vpers(1);

% Fly-by flight path
bVector = normalize(vinc+vout,'norm');
vpers = [bVector*vpers(1) bVector*vpers(2)];
end