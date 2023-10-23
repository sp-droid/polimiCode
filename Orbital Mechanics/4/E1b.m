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

%%
ngrid = 100;

%% Fly-by grid plots
impact = linspace(0, 40*radEarth, ngrid);
Y = zeros(3,ngrid);

for i=1:ngrid
    impactVector = [0;impact(i);0];
    [vout, rper, vper, turnAngle] = hyperb2BP(vinc, impactVector, muEarth);
    Y(1,i) = vecnorm(rper);
    Y(2,i) = vecnorm(vper);
    Y(3,i) = turnAngle;
end

figure;
impact = impact/radEarth;
plot(impact,Y(1,:)/radEarth,'LineWidth',2)
title('Minimum altitude')
xlabel('Impact parameter [R]')
ylabel('Fly-by minimum altitude');
grid on;

figure;
plot(impact,Y(2,:),'LineWidth',2)
title('Maximum velocity')
xlabel('Impact parameter [R]')
ylabel('Fly-by maximum velocity');
grid on;

figure;
plot(impact,rad2deg(Y(3,:)),'LineWidth',2)
title('Turning angle')
xlabel('Impact parameter [R]')
ylabel('Fly-by turning angle');
grid on;

%% Flight path
% We aren't given the speed of the Earth so I'm gonna assume it's circular
vEarth = [0;sqrt(muSun/rEarth(1));0];
vBefore = vEarth+vinc;
vAfter = vEarth+vout;

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