clear
close all

load('horizons.mat')
JD = Yreal(1,1);
ngrid = length(Yreal);
for j=1:ngrid
    Yreal(j,1) = jd2mjd2000(Yreal(j,1));
end
% Yreal has vars in this order: (mjd2000,e,i,bOmega,sOmega,theta,a)

%% Inputs
% Chosen
nmax = 360;

% Rocket body parameters
cD = 2.2;
AoverM = 7.88241304748947 / 900; % m^2/kg
cR = 0.6; %Radiation pressure coefficient

% Celestial parameters
muEarth = astroConstants(13);
Rearth = astroConstants(23);
Tearth = 23*3600+56*60+4.09053;
wEarth = 2*pi/Tearth;
J2 = astroConstants(9);
CS = load('egm96/egm96_to360.ascii', '-ascii');
greenwich = 0;

muSun = astroConstants(4);
Rsun = astroConstants(3);
TTsun = 5778;

muMoon = astroConstants(20);

% Orbital parameters
a = Yreal(1,7);
e = Yreal(1,2);
i = Yreal(1,3);
bOmega = Yreal(1,4);
sOmega = Yreal(1,5);
theta = Yreal(1,6);

[r0, v0] = kep2car(a, e, i, bOmega, sOmega, theta, muEarth, 'deg'); y0 = [r0; v0];
tWindow = [date2mjd2000([2021;11;1;7;19;20.285]);];
tspan = Yreal(end,1)-Yreal(1,1); tWindow(2)=tWindow(1)+tspan; tspan = tspan*86400;

%% Orbit propagation J2+drag
densityModel = @(r) densitySimplified(norm(r)-Rearth);
opts.j2Pert = @(r) j2Pert(r,J2,Rearth,muEarth);
opts.drag = @(r,v) drag(r, v, densityModel(r), wEarth, cD, AoverM);

opts.RelTol = 1e-13;
opts.AbsTol = 1e-14;
[Y1, T] = timed2BP(y0, muEarth, opts, ngrid, tspan);

[scaledT, Tname] = timescaling(T);
Y1kep = zeros(ngrid,6);
for j=1:ngrid
    [a,e,i,bOmega,sOmega,theta] = car2kep(Y1(j,1:3),Y1(j,4:6),muEarth,'deg');
    Y1kep(j,1) = a; Y1kep(j,2) = e; Y1kep(j,3) = i;
    Y1kep(j,4) = bOmega; Y1kep(j,5) = sOmega; Y1kep(j,6) = theta;
    T(j) = T(j)/86400;
end

%% Orbit propagation, full model
% Initial greenwich sidereal time
jdT = (JD - 2451545.0 ) / 36525;
opts2.theta0 = wrapTo2Pi(deg2rad(280.46061837 + 360.98564736629*(JD-2451545.0) ...
    + 0.000387933*jdT^2 - jdT^3/38710000.0));

sunPos = @(t) relativeSun(t, tWindow(1), muSun);
moonPos = @(t) relativeMoon(t, tWindow(1));
opts2.wEarth = wEarth;
densityModel = @(r) densitySimplified(norm(r)-Rearth);
[A,B] = legendreAB(nmax);

opts2.RelTol = 1e-13;
opts2.AbsTol = 1e-14;
opts2.OutputFcn = @odeplot;
opts2.sunThirdBody = @(r,t) thirdBodyPert(r, sunPos(t)-r, muSun);
opts2.moonThirdBody = @(r,t) thirdBodyPert(r, moonPos(t)-r, muMoon);
opts2.egm96 = @(r,thetaG) egm96(r, thetaG, Rearth, muEarth, nmax, CS, A, B);
opts2.relativEffect = @(r,v) relativEffect(r, v, muEarth);
opts2.drag = @(r,v) drag(r, v, densityModel(r), wEarth, cD, AoverM);
opts2.srp = @(r,t) srp(r, sunPos(t)-r, TTsun, Rsun, Rearth, cR, AoverM);
opts2.perturbShow = true;

[ Y2, T ] = timed2BP(y0, muEarth, opts2, ngrid, tspan);

Y2kep = zeros(ngrid,6);
for j=1:ngrid
    [a,e,i,bOmega,sOmega,theta] = car2kep(Y2(j,1:3),Y2(j,4:6),muEarth,'deg');
    Y2kep(j,1) = a; Y2kep(j,2) = e; Y2kep(j,3) = i;
    Y2kep(j,4) = bOmega; Y2kep(j,5) = sOmega; Y2kep(j,6) = theta;
    T(j) = T(j)/86400;
end

%% Plotting fixes
Y1kep(:,6) = unwrap(Y1kep(:,6));
Y2kep(:,6) = unwrap(Y2kep(:,6));
Yreal(:,6) = unwrap(Yreal(:,6));
Yreal(:,1) = Yreal(:,1)-Yreal(1,1);
limits = tWindow-tWindow(1);

%% Plots
figure;
subplot(3,2,1)
plot(T,Y1kep(:,1),'LineWidth',1)
hold on
plot(T,Y2kep(:,1),'LineWidth',1)
plot(Yreal(:,1),Yreal(:,7),'LineWidth',1,'Color','k'); title('Semi-major axis'); axis tight
subplot(3,2,2)
plot(T,Y1kep(:,2),'LineWidth',1)
hold on
plot(T,Y2kep(:,2),'LineWidth',1)
plot(Yreal(:,1),Yreal(:,2),'LineWidth',2,'Color','k'); title('Eccentricity'); axis tight
subplot(3,2,3)
plot(T,Y1kep(:,3),'LineWidth',2)
hold on
plot(T,Y2kep(:,3),'LineWidth',2)
plot(Yreal(:,1),Yreal(:,3),'LineWidth',2,'Color','k'); title('Inclination'); axis tight
subplot(3,2,4)
plot(T,Y1kep(:,4),'LineWidth',2)
hold on
plot(T,Y2kep(:,4),'LineWidth',2)
plot(Yreal(:,1),Yreal(:,4),'LineWidth',2,'Color','k'); title('RAAN'); axis tight
subplot(3,2,5)
plot(T,Y1kep(:,5),'LineWidth',2)
hold on
plot(T,Y2kep(:,5),'LineWidth',2)
plot(Yreal(:,1),Yreal(:,5),'LineWidth',2,'Color','k'); title('Argument of periapsis'); axis tight
subplot(3,2,6)
plot(T,Y1kep(:,6),'LineWidth',2)
hold on
plot(T,Y2kep(:,6),'LineWidth',2)
plot(Yreal(:,1),Yreal(:,6),'LineWidth',2,'Color','k'); title('True anomaly'); axis tight
legend('J2+drag','Full model','NORAD 25850U','Location','northwest')

figure;
subplot(3,2,1)
plot(T,(Yreal(:,7)-Y1kep(:,1))/Yreal(1,7),'LineWidth',1)
hold on
plot(T,(Yreal(:,7)-Y2kep(:,1))/Yreal(1,7),'LineWidth',1)
title('Semi-major axis'); axis tight
subplot(3,2,2)
plot(T,(Yreal(:,2)-Y1kep(:,2)),'LineWidth',1)
hold on
plot(T,(Yreal(:,2)-Y2kep(:,2)),'LineWidth',1)
title('Eccentricity'); axis tight
subplot(3,2,3)
plot(T,(Yreal(:,3)-Y1kep(:,3))/180,'LineWidth',2)
hold on
plot(T,(Yreal(:,3)-Y2kep(:,3))/180,'LineWidth',2)
title('Inclination'); axis tight
subplot(3,2,4)
plot(T,(Yreal(:,4)-Y1kep(:,4))/360,'LineWidth',2)
hold on
plot(T,(Yreal(:,4)-Y2kep(:,4))/360,'LineWidth',2)
title('RAAN'); axis tight
subplot(3,2,5)
plot(T,(Yreal(:,5)-Y1kep(:,5))/360,'LineWidth',2)
hold on
plot(T,(Yreal(:,5)-Y2kep(:,5))/360,'LineWidth',2)
title('Argument of periapsis'); axis tight
subplot(3,2,6)
plot(T,(Yreal(:,6)-Y1kep(:,6))/360,'LineWidth',2)
hold on
plot(T,(Yreal(:,6)-Y2kep(:,6))/360,'LineWidth',2)
title('True anomaly'); axis tight
legend('vs.J2+drag','vs.Full model','Location','northwest')

figure;
subplot(3,2,1)
plot(T,(Yreal(:,7)-Y2kep(:,1))/Yreal(1,7),'LineWidth',1)
title('Semi-major axis'); axis tight
subplot(3,2,2)
plot(T,(Yreal(:,2)-Y2kep(:,2)),'LineWidth',1)
title('Eccentricity'); axis tight
subplot(3,2,3)
plot(T,(Yreal(:,3)-Y2kep(:,3))/180,'LineWidth',2)
title('Inclination'); axis tight
subplot(3,2,4)
plot(T,(Yreal(:,4)-Y2kep(:,4))/360,'LineWidth',2)
title('RAAN'); axis tight
subplot(3,2,5)
plot(T,(Yreal(:,5)-Y2kep(:,5))/360,'LineWidth',2)
title('Argument of periapsis'); axis tight
subplot(3,2,6)
plot(T,(Yreal(:,6)-Y2kep(:,6))/360,'LineWidth',2)
title('True anomaly'); axis tight
legend('vs.Full model','Location','northwest')

%% Secular variations
secularvals(T, Y2kep(:,1), 'Semi-major axis', ngrid)
secularvals(T, Y2kep(:,2), 'Eccentricity', ngrid)
secularvals(T, Y2kep(:,3), 'Inclination', ngrid)
secularvals(T, Y2kep(:,4), 'RAAN', ngrid)
secularvals(T, Y2kep(:,5), 'Argument of periapsis', ngrid)
secularvals(T, Y2kep(:,6), 'True anomaly', ngrid)

function secularvals(T, var, name, ngrid)

T = T*86400;
secular = movmean( var, ceil(ngrid/2) ); % half a year, full ngrid is 1 year

name
X = [ones(length(T),1) T];
b = X\secular
R2 = 1-sum((secular - X*b).^2)/sum((secular - mean(secular)).^2)
end

%% Functions
function r = relativeSun(tseconds, initialMJD, muSun)
T = initialMJD + tseconds/24/3600;
[kep,~] = uplanet(T, 3);
[r,~] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6), muSun, 'rad');
% This is a position vector pointing from the Sun to the Earth, but in
% ecliptic frame. We need to transform back to equatorial
r = -r;

tilt = deg2rad(23.4365472133); %Ecliptic obliquity in January 2021
r = rotRx(tilt)'*r;
end

function r = relativeMoon(tseconds, initialMJD)
T = initialMJD + tseconds/24/3600;
[r,~] = ephMoon(T);
r = r'; % 3x1
end