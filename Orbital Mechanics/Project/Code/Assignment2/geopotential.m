clear
close all

% EGM96 from https://cddis.nasa.gov/926/egm96/getit.html
% https://www.orbiter-forum.com/threads/j4-j5-perturbations.36042/page-2
% http://mitgcm.org/~mlosch/geoidcookbook/node11.html
coefs = load('egm96/egm96_to360.ascii', '-ascii');

r = [ -4431.86444818923;338.703956088788;-6195.86273753096 ];
R = astroConstants(23);
mu = astroConstants(13);
nmax = 22;

rnorm = norm(r); x = r(1); y = r(2); z = r(3);
long = atan2(y,x); ro = sqrt(x^2 + y^2);
xi = z/rnorm; u = sqrt(1-xi^2); Roverr = R/rnorm;

anmM = zeros(nmax+1); bnmM = zeros(nmax+1);
for n=2:nmax
    for m=0:n-1
        anmM(n+1,m+1) = sqrt((2*n-1)*(2*n+1)/(n-m)/(n+m));
        bnmM(n+1,m+1) = sqrt((2*n+1)*(n+m-1)*(n-m-1)/(n-m)/(n+m)/(2*n-3));
    end
end

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
            PnmM(n+1,n+1) = Pnm; PnmdM(n+1,m+1) = Pnmd;
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