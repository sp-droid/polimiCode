% Satellite groundtrack
clear
close all
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [50 100 1000 500]);
addpath(genpath('../shared'))

% Physical parameters
muEarth = astroConstants(13);
Rearth = astroConstants(23);
Tearth = 23*3600+56*60+4.1;
wEarth = 2*pi/Tearth;
J2 = astroConstants(9);
greenwich = 0;

% Orbit parameters
a = 7271;                           % Semi-major axis
e = 0.0;                            % Eccentricity
i = 78;                             % Inclination
bOmega = 0;                         % Right ascension of ascending node
sOmega = 0;                         % Argument of pericentre
theta = 0;                          % Initial true anomaly

date = [2032;1;1;0;0;0];
tWindow = [date2mjd2000(date);0];

%koverm = Tearth/2/pi*sqrt(muEarth/a^3)
%%
% Repeating GT
koverm = 14;
N = 13;
P = 20;
Q = 21;
nOrbits = N*Q+P; nOrbits = 299;
nDays = Q; nDays = 21;

%a = ((Tearth/koverm/2/pi)^2*muEarth)^(1/3);
a = ((nDays*Tearth/nOrbits/2/pi)^2*muEarth)^(1/3);
%a = Rearth+885.26;

[r0, v0] = kep2car(a, e, i, bOmega, sOmega, theta, muEarth, 'deg'); y0 = [r0; v0];
Torb = 2*pi*sqrt( a^3/muEarth );          % Orbital period
lambda = Torb*wEarth;                     % Ground track drift
nOrb = 1;
nPoints = 100000;
ttime = nOrbits*Torb;

%% Computation
% Perform the integration
opts.RelTol = 1e-12;
opts.AbsTol = 1e-13;
%opts.j2Pert = @(r) j2Pert(r,J2,Rearth,muEarth);

[ Y, t ] = timed2BP(y0,muEarth,opts,nPoints,ttime);
t = t';

r = Y(:,1:3);
v = Y(:,4:6);
Rnorm = vecnorm(r');
Vnorm = vecnorm(v');
tOrbs = t/Torb;


delta = asin(r(:,3)'./Rnorm);                  % Declination
alpha = atan2(r(:,2)',r(:,1)');                 % Right ascension
long = wrapTo180(rad2deg(alpha-greenwich-wEarth*t));
lat = rad2deg(delta);

%%


distances = zeros(1,nPoints);
for i=1:nPoints
    distances(i) = haversine(long(1),lat(1),long(i),lat(i));
end

maxDistance = 2*2*pi*Rearth/nOrbits*sind(78);
revisits = isolateRevisits(distances, maxDistance);
distances(revisits)
length(distances(revisits))

figure;
plot((t-tOrbs(1))./86400, distances, 'LineWidth', 1)
hold on
scatter((t(revisits)-t(1))./86400, distances(revisits), 30, 'MarkerEdgeColor','red','MarkerFaceColor','red')
grid on
xlabel('Day')
ylabel('Distance')
xlim([0,(max(t)-t(1))/86400])
ylim([0,max(distances)])
fontsize(14,"points")
hold off

function indices = isolateRevisits(distances, maxDistance)

is_real = distances < maxDistance;
% Find the start and end of each island
island_boundaries = diff([0, is_real, 0]);
start_indices = find(island_boundaries == 1); % start of each island
end_indices = find(island_boundaries == -1) - 1; % end of each island
% Initialize an array to hold the minimum values
indices = NaN(1, length(start_indices));
% Loop over each island and find the minimum value
for i = 1:length(start_indices)
    % Extract the current island
    island = distances(start_indices(i):end_indices(i));
    
    % Find and store the minimum of this island
    [~, idx] = min(island);
    indices(i) = idx + start_indices(i) - 1;
end

end

function distance = haversine(lat1, lon1, lat2, lon2)

R = 6371;
a = sind((lat2-lat1)/2)^2 + cosd(lat2)*cosd(lat1)*sind((lon2-lon1)/2)^2;
c = 2*atan2(sqrt(a),sqrt(1-a));
distance = R*c;

end