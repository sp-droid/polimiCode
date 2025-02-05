%% Tidal aliasing
clear
close all
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [50 100 1000 500]);
addpath(genpath('../shared'))

%% Tidal components
tideNames = ["S1","S2","K1","K2","M2","N2","P1","O1","Q1","Mf","Mm","Ssa","Sa"];
tidePeriods = [1,0.5,0.99727,0.498635,0.517525,0.527431,1.002745,1.075806,1.119515,13.66,27.55,182.621095,365.2594];  % days
tideGroup1 = [1,7]; % For filtering retrograde orbits
tideGroup2 = [1,2,3,4,7]; % For filtering higher inclinations
tideGroup3 = [1,2,4,7]; % Only for filter
tideGroup4 = [5,6,8,9,10,11]; % For filtering over the 2 cpy threshold

for i=1:length(tidePeriods)
    fprintf(strcat(tideNames(i),': ',string(tidePeriods(i)),'\t'))
end
fprintf('\n')

%% Load candidates and calculate aliasing frequencies [cpy]
load("candidates.mat");

%% Calculate aliasing frequencies and component separability
aliasCPYs = zeros(length(orbits),length(tidePeriods)); % cycles per year
tidalSeparability = zeros(length(orbits),1); % years

for i=1:length(orbits)
    aliasCPYs(i,:) = arrayfun(@(x) 365.25/aliasingPeriods(orbits(i,5), x), tidePeriods);

    for m=1:length(tidePeriods)-1
        for n=m+1:length(tidePeriods)
            separability = 1/abs(aliasCPYs(i,m)-aliasCPYs(i,n)); % Rayleigh criterion of separability
            if separability>tidalSeparability(i)
                tidalSeparability(i) = separability;
            end
        end
    end
end

%% Plotting and filtering by CPYs

C1 = [0 0.4470 0.7410];  % blue
C2 = [0.8500 0.3250 0.0980];  % orange
C3 = [0.9290 0.6940 0.1250];  % yellow
C4 = [0.4940 0.1840 0.5560];  % purple
C5 = [0.4660 0.6740 0.1880];  % green
C6 = [0.3010 0.7450 0.9330];  % light blue
C7 = [0.6350 0.0780 0.1840];  % red
C8 = [1.0000 0.4000 0.6000];  % pink
C9 = [0.6350 0.0780 0.1840];  % brown
C10 = [0.0000 0.0000 0.0000];  % black
C11 = [0.0, 0.5, 0.5]; % dark teal
C12 = [0.333, 0.420, 0.184]; % olive green
C13 = [0.0, 0.5, 0.5]; % dark teal
C14 = [0.333, 0.420, 0.184]; % olive green

for option=1:5
    switch option
        case 1
            xplot = 90:115;
            tideGroup = tideGroup1;
            xlabelName = 'Inclination [deg]';
            ylims = [0,3.5];
            name = "aliasingRetrograde";
            filtered = ones(length(orbits),1);
            fprintf(strcat(string(sum(filtered))," candidates\n"))
        case 2
            xplot = 65:90;
            tideGroup = tideGroup2;
            xlabelName = 'Inclination [deg]';
            ylims = [0,8];
            name = "aliasingInclination";
            filtered = orbits(:,3)<90;
            fprintf(strcat(string(sum(filtered))," candidates after removing retrograde orbits\n"))
        case 3
            xplot = 9:35;
            tideGroup = tideGroup4;
            xlabelName = 'Cycle [days]';
            ylims = [0,16];
            name = "aliasingCycles";
            filtered = filtered & all(aliasCPYs(:,tideGroup3)>2,2);
            fprintf(strcat(string(sum(filtered))," candidates after applying 2 cpy threshold based on inclination (except K1)\n"))
            fprintf(strcat(string(max(orbits(filtered,3))),"º max inclination\n"))
            filtered = filtered & aliasCPYs(:,3)>1;
            fprintf(strcat(string(sum(filtered))," candidates after applying 1 cpy threshold for K1 exclusively\n"))
            fprintf(strcat(string(max(orbits(filtered,3))),"º max inclination\n"))
        case 4
            xplot = 10:33;
            tideGroup = tideGroup4;
            xlabelName = 'Cycle [days]';
            ylims = [2,16];
            name = "aliasingCyclesOver2";
            filtered = filtered & all(aliasCPYs(:,tideGroup4)>2,2);
            fprintf(strcat(string(sum(filtered))," candidates after applying 2 cpy threshold based on cycles\n"))
            uniqueCycles = unique(orbits(filtered,2));
            for i=1:length(uniqueCycles)
                fprintf(strcat(string(uniqueCycles(i)),'\t'))
            end
            fprintf(" available repeat cycles\n")
            fprintf(strcat(string(max(orbits(filtered,3))),"º max inclination\n"))
        case 5
            xplot = 10:22;
            tideGroup = tideGroup4;
            xlabelName = 'Cycle [days]';
            ylims = [2,16];
            name = "aliasingWithSeparability";
            filtered = filtered & tidalSeparability<5;
            fprintf(strcat(string(sum(filtered))," candidates after requiring tide separability of 5 years or less\n"))
            uniqueCycles = unique(orbits(filtered,2));
            for i=1:length(uniqueCycles)
                fprintf(strcat(string(uniqueCycles(i)),'\t'))
            end
            fprintf(" available repeat cycles\n")
            fprintf(strcat(string(max(orbits(filtered,3))),"º max inclination\n"))
    end
    
    figure;
    hold on
    for m=1:length(tideGroup)
        j = tideGroup(m);
        yplot = zeros(1,length(xplot));
        width = ones(1,length(xplot));
        for i=1:length(xplot)
            switch option
                case {1,2}
                    selected = filtered & orbits(:,3)==xplot(i);
                otherwise
                    selected = filtered & orbits(:,2)==xplot(i);
            end
            results = aliasCPYs(selected,j);
            if isempty(results)
                yplot(i) = NaN;
                continue
            end
            yplot(i) = (max(results)+min(results))/2;
            width(i) = max(max(results)-yplot(i),0.05);
        end
        color = eval(['C' num2str(m+1)]);
        
        errorbar(xplot, yplot, width, '.', 'Color',color,'DisplayName',tideNames(j),'LineWidth',6,'CapSize',0)
    end
    
    yline(2,'Color','red','DisplayName','Threshold')
    
    legend('Location','best');
    grid on
    xlabel(xlabelName)
    ylabel('aliasing freq. [cpy]')
    xlim([min(xplot)*0.99,max(xplot)*1.01])
    ylim(ylims)
    fontsize(16,"points")
    
    print(gcf, strcat('output/',string(name),'.png'), '-dpng', '-r200');
end

%% Remaining orbits
orbits = orbits(filtered,:);
tidalSeparability = tidalSeparability(filtered);

figure;
scatter(orbits(:,4), orbits(:,2), 16, tidalSeparability,'filled')
grid on
xlabel('Altitude [km]')
ylabel('Cycle [days]')
fontsize(16,"points")
colorbar; colormap jet
title("Tidal separability of remaining orbits [years]")
print(gcf, 'output/remainingCycles.png', '-dpng', '-r200');

figure;
scatter(orbits(:,4), orbits(:,3), 12, tidalSeparability,'filled')
grid on
xlabel('Altitude [km]')
ylabel('Inclination [deg]')
fontsize(16,"points")
colorbar; colormap jet
title("Tidal separability of remaining orbits [years]")
print(gcf, 'output/remainingInclination.png', '-dpng', '-r200');

%%
RearthAverage = 6371010E-3;
Rearth = 6378137.0E-3;
J2 = 0.0010826261738522227;         % EGM2008
J4 = -1.6198975999169731E-6;        % EGM2008
Tearth = 86164.0905;
muEarth = astroConstants(13);

orbitsParams = zeros(length(orbits),7);
% Norbits[-], cycle[days], inclination[deg], altitude[km], exact repeat cycle[days]
% wEffective[km], w[km], half scan angle[deg], iEffective[deg], coverage[-], precession[], latchange
for i=1:length(orbits)
    a = Rearth+orbits(i,4); e = orbits(i,6);

    inc = orbits(i,3); p = a*(1-e^2); RoverP = Rearth/p; eFactor = sqrt(1-e^2);
    n0 = sqrt(muEarth/a^3);
    n = meanMotion(a, inc, J2, J4, muEarth, RoverP^2, RoverP^4, eFactor, e);
    RAANsecular = -3*pi*J2*RoverP^2*cosd(inc) ...
        +3*pi/16*J2^2*RoverP^4*cosd(inc)*(-36-4*e^2+48*eFactor+(40-5*e^2-72*eFactor)*sind(inc)^2) ...
        +n0/n*15*pi/16*J4*RoverP^4*cosd(inc)*(8+12*e^2-(14+21*e^2)*sind(inc)^2);
    precession = 1/abs((-360/365.2422+rad2deg(RAANsecular*n/2/pi)*86400)/360); % days for 24h RAAN change
    latchange = abs(-4*pi^2/n/Tearth + RAANsecular)/2/pi*360;

    wEffective = (360-latchange)/360*2*pi*Rearth/orbits(i,1);
    w = wEffective*sind(orbits(i,3));
    scanAngle = 0.5*atan2d(w,a-sqrt(Rearth^2-w^2));
    iEffective = orbits(i,3)+asind(w/Rearth);
    coverage = sind(iEffective);

    orbitsParams(i,:) = [wEffective,w,scanAngle,iEffective,coverage,precession,latchange];
end

figure;
scatter(orbits(:,3), orbitsParams(:,2), 16,tidalSeparability,'filled')
grid on
xlabel('Inclination [deg]')
ylabel('Minimum swath [km]')
fontsize(16,"points")
colorbar; colormap jet
print(gcf, 'output/remainingInclinationSwath.png', '-dpng', '-r200');

figure;
scatter(orbits(:,4), orbitsParams(:,5)*100, 16,tidalSeparability,'filled')
grid on
xlabel('Altitude [km]')
ylabel('Coverage [%]')
fontsize(16,"points")
colorbar; colormap jet
print(gcf, 'output/remainingCoverage.png', '-dpng', '-r200');

filtered = orbitsParams(:,2) < 150;
fprintf(strcat(string(sum(filtered))," candidates after limiting the minimum swath width to 150 km\n"))
filtered = orbitsParams(:,2) < 150 & tidalSeparability < 3;
fprintf(strcat(string(sum(filtered))," candidates after limiting tide separability to 3 years\n"))
orbits = orbits(filtered,:); orbitsParams = orbitsParams(filtered,:); tidalSeparability = tidalSeparability(filtered);

%% Final selection
filtered = orbits(:,4) < 920 & orbits(:,4) > 870;
orbits = orbits(filtered,:); orbitsParams = orbitsParams(filtered,:); tidalSeparability = tidalSeparability(filtered);

figure;
scatter(orbits(:,4), orbits(:,3), 32, tidalSeparability,'filled')
grid on
xlabel('Altitude [km]')
ylabel('Inclination [deg]')
fontsize(16,"points")
colorbar; colormap jet
print(gcf, 'output/finalClusters.png', '-dpng', '-r200');

%%
options = [length(orbits),1];
figure;
for k=1:2
    i = options(k);
    a = Rearth+orbits(i,4); e = orbits(i,6); inc = orbits(i,3);
    
    norbit = 0:1:orbits(i,1);
    distances = zeros(1,length(norbit));
    p = a*(1-e^2); RoverP = Rearth/p; eFactor = sqrt(1-e^2);
    n0 = sqrt(muEarth/a^3);
    n = meanMotion(a, inc, J2, J4, muEarth, RoverP^2, RoverP^4, eFactor, e);
    drift = abs(-4*pi^2/n/Tearth ...
        -3*pi*J2*RoverP^2*cosd(inc) ...
        +3*pi/16*J2^2*RoverP^4*cosd(inc)*(-36-4*e^2+48*eFactor+(40-5*e^2-72*eFactor)*sind(inc)^2) ...
        +n0/n*15*pi/16*J4*RoverP^4*cosd(inc)*(8+12*e^2-(14+21*e^2)*sind(inc)^2) );
    
    wEarth = 2*pi/Tearth;

    for j=1:length(norbit)
        distance = drift*norbit(j)/2/pi;
        distance = abs(distance-round(distance))*360;
        distances(j) = distance;
    end
    
    xplot = norbit/orbits(i,1)*orbits(i,2);
    approaches = isolateRevisits(distances, 0.2);
    approaches2 = isolateRevisits(distances(approaches), 0.02)-1;
    
    plot(xplot, distances, 'LineWidth',4/k,'DisplayName',strcat('h=',string(orbits(i,4)),' km'))
    hold on
end
grid on
xlabel('Time [days]')
xlim([0,max(xplot)])
legend
ylabel('Angular distance [deg]')
fontsize(14,"points")
print(gcf, 'output/subCycles.png', '-dpng', '-r200');

subcycles = zeros(length(orbits),1);
for i=1:length(orbits)
    a = Rearth+orbits(i,4); e = orbits(i,6); inc = orbits(i,3);
    
    norbit = 0:1:orbits(i,1);
    distances = zeros(1,length(norbit));
    p = a*(1-e^2); RoverP = Rearth/p; eFactor = sqrt(1-e^2);
    n0 = sqrt(muEarth/a^3);
    n = meanMotion(a, inc, J2, J4, muEarth, RoverP^2, RoverP^4, eFactor, e);
    drift = abs(-4*pi^2/n/Tearth ...
        -3*pi*J2*RoverP^2*cosd(inc) ...
        +3*pi/16*J2^2*RoverP^4*cosd(inc)*(-36-4*e^2+48*eFactor+(40-5*e^2-72*eFactor)*sind(inc)^2) ...
        +n0/n*15*pi/16*J4*RoverP^4*cosd(inc)*(8+12*e^2-(14+21*e^2)*sind(inc)^2) );
    
    wEarth = 2*pi/Tearth;
    
    for j=1:length(norbit)
        distance = drift*norbit(j)/2/pi;
        distance = abs(distance-round(distance));
        distances(j) = distance;
    end

    distances2 = 0:0.5:orbits(i,1);
    for j=1:length(distances2)
        distance = drift*distances2(j);
        distance = abs(distance-floor(distance))*360;
        distances2(j) = distance;
    end
    longs = 0:0.001:360;
    maxSwathAngle = 0;
    for j=1:length(longs)
        minDistance = 360;
        for m=1:length(distances2)
            distance = abs(longs(j)-distances2(m));
            if distance > 180
                distance = 360-distance;
            end
            if distance < minDistance
                minDistance = distance;
            end
        end
        if minDistance > maxSwathAngle
            maxSwathAngle = minDistance;
        end
    end
    maxLon = maxSwathAngle*2*pi/180;
    w = 1.1*0.5*Rearth*asin(sind(orbits(i,3))*sin(maxLon));
    scanAngle = atan2d(w,a-sqrt(Rearth^2-w^2));
    iEffective = orbits(i,3)+asind(w/Rearth);
    coverage = sind(iEffective);
    orbitsParams(i,1:5) = [maxLon*Rearth/2,w,scanAngle,iEffective,coverage];


    approaches = isolateRevisits(distances, 0.2);
    approaches2 = isolateRevisits(distances(approaches), 0.02)-1;
    subcycles(i) = length(approaches2)-2;
end

filtered = subcycles > 0;
orbits = orbits(filtered,:); orbitsParams = orbitsParams(filtered,:); tidalSeparability = tidalSeparability(filtered);

[~,chosen] = min(tidalSeparability);

%% Description
orbits(chosen,:)
orbitsParams(chosen,:)
tidalSeparability(chosen)

%% Functions
function aliasedPeriod = aliasingPeriods(revisitPeriod, tidePeriod)
    freqRevisit = 1/revisitPeriod;
    freqTide = 1/tidePeriod;
    
    freqNyquist = freqRevisit/2;
    ratio = floor(freqTide/freqNyquist);
    
    freqAliasing = freqTide - ratio*freqNyquist;
    if (rem(ratio,2)~=0)
       freqAliasing = freqNyquist-freqAliasing;
    end
    aliasedPeriod = 1/freqAliasing;
end

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

function n = meanMotion(a, i, J2, J4, muEarth, RoverPpower2, RoverPpower4, eFactor, e)
    n0 = sqrt(muEarth/a^3);
    n = n0*(1+...
        +3/4*J2*RoverPpower2*eFactor*(2-3*sind(i)^2) ...
        +3/128*J2^2*RoverPpower4*eFactor*(120+64*eFactor-40*eFactor^2+(-240-192*eFactor+40*eFactor^2)*sind(i)^2+(105+144*eFactor+25*eFactor^2)*sind(i)^4) ...
        -45/128*J4*RoverPpower4*eFactor*e^2*(-8+40*sind(i)^2-35*sind(i)^4));
end