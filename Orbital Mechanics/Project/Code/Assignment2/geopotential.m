clear
close all

% EGM96 from https://cddis.nasa.gov/926/egm96/getit.html
% https://www.orbiter-forum.com/threads/j4-j5-perturbations.36042/page-2
% http://mitgcm.org/~mlosch/geoidcookbook/node11.html
% https://people.sc.fsu.edu/~lb13f/projects/space_environment/egm96.php
coefs = load('egm96/egm96_to360.ascii', '-ascii');

R = astroConstants(23);
mu = astroConstants(13);
nmax = 360;

anmM = zeros(nmax+1); bnmM = zeros(nmax+1);
for n=2:nmax
    for m=0:n-1
        anmM(n+1,m+1) = sqrt((2*n-1)*(2*n+1)/(n-m)/(n+m));
        bnmM(n+1,m+1) = sqrt((2*n+1)*(n+m-1)*(n-m-1)/(n-m)/(n+m)/(2*n-3));
    end
end

%% test
ngrid = 100;
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
        N(i,j) = egm96undulation(r,R,mu,nmax,coefs,anmM,bnmM)*1000;
    end
end

figure;
img = imread('earth2Doutline','png');
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

%% test
r = [ -4431.86444818923;338.703956088788;-6195.86273753096 ];
rnorm = norm(r); x = r(1); y = r(2); z = r(3);
long = atan2(y,x); ro = sqrt(x^2 + y^2);
xi = z/rnorm; u = sqrt(1-xi^2); Roverr = R/rnorm;

tic
a = [0;0;0];
PnmM = zeros(nmax+1); PnmM(1,1) = 1; PnmM(2,1) = xi*sqrt(3); PnmM(2,2) = u*sqrt(3);
PnmdM = zeros(nmax+1); PnmdM(1,1) = 0; PnmdM(2,1) = sqrt(3); PnmdM(2,2) = -xi/u*sqrt(3);
index = 0;
for n=2:nmax
    aN = [0;0;0];
    for m=0:n
        if (m==n)
            Pnm = sqrt((2*n+1)/2/n)*u*PnmM(n,n);
            Pnmd = sqrt((2*n+1)/2/n)*(u*PnmdM(n,n)-xi/u*PnmM(n,n));
            PnmM(n+1,n+1) = Pnm; PnmdM(n+1,n+1) = Pnmd;
        else
            anm = anmM(n+1,m+1);
            bnm = bnmM(n+1,m+1);
            Pnm = anm*xi*PnmM(n,m+1)-bnm*PnmM(n-1,m+1);
            Pnmd = anm*PnmM(n,m+1)+anm*xi*PnmdM(n,m+1)-bnm*PnmdM(n-1,m+1);
            PnmM(n+1,m+1) = Pnm; PnmdM(n+1,m+1) = Pnmd;
        end
        aC = [Pnm*rnorm*((1+n)*x*cos(m*long)-m*(rnorm/ro)^2*y*sin(m*long))+Pnmd*x*z*cos(m*long)
              Pnm*rnorm*((1+n)*y*cos(m*long)+m*(rnorm/ro)^2*x*sin(m*long))+Pnmd*y*z*cos(m*long)
              Pnm*rnorm*(1+n)*z*cos(m*long)-Pnmd*(x^2+y^2)*cos(m*long)];
        aS = [Pnm*rnorm*((1+n)*x*sin(m*long)+m*(rnorm/ro)^2*y*cos(m*long))+Pnmd*x*z*sin(m*long)
              Pnm*rnorm*((1+n)*y*sin(m*long)-m*(rnorm/ro)^2*x*cos(m*long))+Pnmd*y*z*sin(m*long)
              Pnm*rnorm*(1+n)*z*sin(m*long)-Pnmd*(x^2+y^2)*sin(m*long)];

        index = index + 1;
        aN = aN - aC*coefs(index,3) - aS*coefs(index,4);
    end
    a = a + aN*Roverr^n;
end
a = a*mu/rnorm^4;
toc

%% funcs
function N = egm96undulation(r, R, mu, nmax, coefs, anmM, bnmM)
rnorm = norm(r); x = r(1); y = r(2); z = r(3);
long = atan2(y,x); xi = z/rnorm; u = sqrt(1-xi^2); Roverr = R/rnorm;

% Perturbing gravitational potential and legendre terms
U = 0;
PnmM = zeros(nmax+1); PnmM(1,1) = 1; PnmM(2,1) = xi*sqrt(3); PnmM(2,2) = u*sqrt(3);
index = 0;
for n=2:nmax
    aN = 0;
    for m=0:n
        if (m==n)
            Pnm = sqrt((2*n+1)/2/n)*u*PnmM(n,n);
            PnmM(n+1,n+1) = Pnm;
        else
            anm = anmM(n+1,m+1);
            bnm = bnmM(n+1,m+1);
            Pnm = anm*xi*PnmM(n,m+1)-bnm*PnmM(n-1,m+1);
            PnmM(n+1,m+1) = Pnm;
        end
        index = index + 1;
        aN = aN + Pnm*(coefs(index,3)*cos(m*long)+coefs(index,4)*sin(m*long));
    end
    U = U + aN*Roverr^n;
end
U = mu/rnorm*(1+U);

% Ellipsoid normal potential
V = 0;
index = 0;
for n=2:nmax
    for m=0:n
        index = index + 1;
        if (~mod(n,2)) && (m==0)
            Pnm = PnmM(n+1,1);
            Jn = coefs(index,3);
            
            V = V + Roverr^n*Jn*Pnm;
        end
    end
end
V = mu/rnorm*(1+V);

% Disturbing potential
T = U-V;

% Gamma or normal gravity on the surface of the ellipsoid
gamma = mu/R^2;

% Undulation given by Brun's formula
N = T/gamma;

end