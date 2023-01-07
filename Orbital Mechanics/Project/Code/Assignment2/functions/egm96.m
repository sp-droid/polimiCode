function a = egm96(r, thetaEarth, Rearth, muEarth, nmax, CS, A, B)
% Earth perturbing acceleration at distance r
%
% PROTOTYPE
% a = egm96acc( r, Rearth, muEarth, nmax, CS, A, B )
%
% INPUTS:
% r[3x1]		- Position vector [L]
% thetaEarth[1] - Angle relative of x axis relative to the Greenwich meridian [rad]
% Rearth[1]     - Earth's radius [L]
% muEarth[1]    - Earth's gravitational parameter [L^3/T^2]
% nmax[1]       - Number of terms used. Egm96 is complete at nmax = 360
% CS[65338x6]   - Matrix of coefficients in the egm96 model
% A[nmax+1]     - Matrix of anm coefficients, calculated through legendreAB.m
% B[nmax+1]     - Matrix of bnm coefficients, calculated through legendreAB.m
%
% OUTPUTS:
% a[3x1]        - Undulation potential [ L/T^2 ]
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
% -------------------------------------------------------------------------

% This assumes r comes from a fixed frame while the Earth is rotating, so we need thetaEarth
rnorm = norm(r); x = r(1); y = r(2); z = r(3);
long = atan2(y,x)-thetaEarth; ro = sqrt(x^2 + y^2);
xi = z/rnorm; u = sqrt(1-xi^2); Roverr = Rearth/rnorm;

a = [0;0;0];
P = zeros(nmax+1); P(1,1) = 1; P(2,1) = xi*sqrt(3); P(2,2) = u*sqrt(3);
PD = zeros(nmax+1); PD(1,1) = 0; PD(2,1) = sqrt(3); PD(2,2) = -xi/u*sqrt(3);
index = 0;
for n=2:nmax
    aN = [0;0;0];
    for m=0:n
        if (m==n)
            Pnm = sqrt((2*n+1)/2/n)*u*P(n,n);
            Pnmd = sqrt((2*n+1)/2/n)*(u*PD(n,n)-xi/u*P(n,n));
            P(n+1,n+1) = Pnm; PD(n+1,n+1) = Pnmd;
        else
            anm = A(n+1,m+1);
            bnm = B(n+1,m+1);
            Pnm = anm*xi*P(n,m+1)-bnm*P(n-1,m+1);
            Pnmd = anm*P(n,m+1)+anm*xi*PD(n,m+1)-bnm*PD(n-1,m+1);
            P(n+1,m+1) = Pnm; PD(n+1,m+1) = Pnmd;
        end
        aC = [Pnm*rnorm*((1+n)*x*cos(m*long)-m*(rnorm/ro)^2*y*sin(m*long))+Pnmd*x*z*cos(m*long)
              Pnm*rnorm*((1+n)*y*cos(m*long)+m*(rnorm/ro)^2*x*sin(m*long))+Pnmd*y*z*cos(m*long)
              Pnm*rnorm*(1+n)*z*cos(m*long)-Pnmd*(x^2+y^2)*cos(m*long)];
        aS = [Pnm*rnorm*((1+n)*x*sin(m*long)+m*(rnorm/ro)^2*y*cos(m*long))+Pnmd*x*z*sin(m*long)
              Pnm*rnorm*((1+n)*y*sin(m*long)-m*(rnorm/ro)^2*x*cos(m*long))+Pnmd*y*z*sin(m*long)
              Pnm*rnorm*(1+n)*z*sin(m*long)-Pnmd*(x^2+y^2)*sin(m*long)];

        index = index + 1;
        aN = aN - aC*CS(index,3) - aS*CS(index,4);
    end
    a = a + aN*Roverr^n;
end
a = a*muEarth/rnorm^4;
end