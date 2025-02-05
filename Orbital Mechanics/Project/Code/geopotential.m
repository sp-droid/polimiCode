clear
close all

% EGM96 from https://cddis.nasa.gov/926/egm96/getit.html
% https://www.orbiter-forum.com/threads/j4-j5-perturbations.36042/page-2
% http://mitgcm.org/~mlosch/geoidcookbook/node11.html
% https://people.sc.fsu.edu/~lb13f/projects/space_environment/egm96.php
CS = load('egm96/egm96_to360.ascii', '-ascii');

R = astroConstants(23);
mu = astroConstants(13);
nmax = 360;
[A,B] = legendreAB(nmax);

%% Undulation geopotential
ngrid = 400;
rnorm = R;%7625.3;
lat = linspace(-90,90,ngrid);
long = linspace(-180,180,ngrid);
[lat,long] = meshgrid(lat,long);
N = zeros(ngrid);
for i=1:ngrid
    for j=1:ngrid
        z = rnorm*sind(lat(i,j));
        x = cosd(long(i,j))*sqrt(rnorm^2-z^2);
        y = sind(long(i,j))*sqrt(rnorm^2-z^2);
        r = [x;y;z];
        N(i,j) = egm96undulation(r,R,mu,nmax,CS,A,B)*1000;
    end
end

figure;
img = imread('earth2Doutline','jpg');
image('CData',img,'XData',[-180 180],'YData',[90,-90]);
hold on
LevelList = linspace(min(min(N)),max(max(N)),40);
contourf(long,lat,N,'LevelList',LevelList,'EdgeColor','none','FaceAlpha',0.7)
xlabel('Longitude'); ylabel('Latitude');
xlim([-180,180]);
xticks([-180,-120,-60,0,60,120,180])
ylim([-90,90]);
yticks([-90,-60,-30,0,30,60,90])
colorbar;
colormap(jet)
grid on;
hold off

%% Geopotential's force
r = [ -4431.86444818923;338.703956088788;-6195.86273753096 ];
egm96(r, 0, R, mu, nmax, CS, A, B)

