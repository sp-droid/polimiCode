function N = egm96undulation(r, Rearth, muEarth, nmax, CS, A, B)
% Earth undulation at distance r
%
% PROTOTYPE
% N = egm96undulation( r, Rearth, muEarth, nmax, CS, A, B )
%
% INPUTS:
% r[3x1]		- Position vector [L]
% Rearth[1]     - Earth's radius [L]
% muEarth[1]    - Earth's gravitational parameter [L^3/T^2]
% nmax[1]       - Number of terms used. Egm96 is complete at nmax = 360
% CS[65338x6]   - Matrix of coefficients in the egm96 model
% A[nmax+1]     - Matrix of anm coefficients, calculated through legendreAB.m
% B[nmax+1]     - Matrix of bnm coefficients, calculated through legendreAB.m
%
% OUTPUTS:
% N[1]	        - Undulation potential [ L*? ]
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
% -------------------------------------------------------------------------
rnorm = norm(r); x = r(1); y = r(2); z = r(3);
long = atan2(y,x); xi = z/rnorm; u = sqrt(1-xi^2); Roverr = Rearth/rnorm;

% Perturbing gravitational potential and legendre terms
U = 0;
P = zeros(nmax+1); P(1,1) = 1; P(2,1) = xi*sqrt(3); P(2,2) = u*sqrt(3);
index = 0;
for n=2:nmax
    aN = 0;
    for m=0:n
        if (m==n)
            Pnm = sqrt((2*n+1)/2/n)*u*P(n,n);
            P(n+1,n+1) = Pnm;
        else
            Anm = A(n+1,m+1);
            Bnm = B(n+1,m+1);
            Pnm = Anm*xi*P(n,m+1)-Bnm*P(n-1,m+1);
            P(n+1,m+1) = Pnm;
        end
        index = index + 1;
        aN = aN + Pnm*(CS(index,3)*cos(m*long)+CS(index,4)*sin(m*long));
    end
    U = U + aN*Roverr^n;
end
U = muEarth/rnorm*(1+U);

% Ellipsoid normal potential
V = 0;
index = 0;
for n=2:nmax
    for m=0:n
        index = index + 1;
        if (~mod(n,2)) && (m==0)
            Pnm = P(n+1,1);
            Jn = CS(index,3);
            
            V = V + Roverr^n*Jn*Pnm;
        end
    end
end
V = muEarth/rnorm*(1+V);

% Disturbing potential
T = U-V;

% Gamma or normal gravity on the surface of the ellipsoid
gamma = muEarth/Rearth^2;

% Undulation given by Brun's formula
N = T/gamma;
end