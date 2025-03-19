%% test
clear
close all
set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [50 100 1000 500]);
addpath(genpath('../shared'))

%%
muEarth = astroConstants(13);
Rearth = astroConstants(23);
Tearth = 86164.0905;
wEarth = 2*pi/Tearth;

% Static parameters
a = Rearth + 892.09;
e = 0.00099976;
inc = 77;
sOmega = 90;

Norbits = 292;
k = 21;
exactRepeat = 20.861;         % days

Npoints = 10000;


%% Calculation
p = a*(1-e^2);
RoverPpower2 = (Rearth/p)^2; RoverPpower4 = RoverPpower2^2;
n0 = sqrt(muEarth/a^3);
Torb = 2*pi/n0;
J2 = 0.0010826261738522227;         % EGM2008
J4 = -1.6198975999169731E-6;        % EGM2008
eFactor = sqrt(1-e^2);

bOmegadot = (-3*pi*J2*RoverPpower2*cosd(i) ...                                                                   % J2
            +3*pi/16*J2^2*RoverPpower4*cosd(i)*(-36-4*e^2+48*eFactor+(40-5*e^2-72*eFactor)*sind(i)^2) ...       % J2^2
            +15*pi/16*J4*RoverPpower4*cosd(i)*(8+12*e^2-(14+21*e^2)*sind(i)^2))/Torb;

output = zeros(Npoints,5);
for i=1:Npoints
    output(i,1) = i-1;

    theta = wrapTo360(Norbits*360/Npoints*i);
    
    time = exactRepeat*86400/Npoints*i;
    bOmega = time*rad2deg(bOmegadot);
    %bOmega = 45;

    [r, v] = kep2car( a, e, inc, bOmega, sOmega, theta, muEarth, 'deg');
    [~,~,~,bOmega,~,~] = car2kep( r, v, muEarth, 'deg' )
    output(i,2) = r(1); output(i,3) = r(2); output(i,4) = r(3);

    output(i,5) = wrapTo2Pi(wEarth*time);
end

figure;
scatter(output(:,2),output(:,3),10,output(:,1))
hold on
t = 0:0.1:360;
plot(Rearth*sind(t),Rearth*cosd(t))
axis equal;

rad2deg(output(end,5))

%% Output
headers = {'frame', 'x', 'y', 'z'};
filename = 'trajectoryData.csv';
writecell(headers, filename);
writematrix(output, filename, 'WriteMode', 'append');