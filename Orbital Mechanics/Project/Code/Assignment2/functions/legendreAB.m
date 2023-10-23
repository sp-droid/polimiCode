function [A,B] = legendreAB(nmax)
% Coefficient matrix for legendre polynomials, to avoid unnecessary computations
%
% PROTOTYPE
% [A,B] = legendreAB( r, Rearth, muEarth, nmax, CS, A, B )
%
% INPUTS:
% nmax[1]       - Number of terms used. Egm96 is complete at nmax = 360
%
% OUTPUTS:
% A[nmax+1]     - Matrix of anm coefficients, calculated through legendreAB.m
% B[nmax+1]     - Matrix of bnm coefficients, calculated through legendreAB.m
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
% -------------------------------------------------------------------------
A = zeros(nmax+1); B = zeros(nmax+1);
for n=2:nmax
    for m=0:n-1
        A(n+1,m+1) = sqrt((2*n-1)*(2*n+1)/(n-m)/(n+m));
        B(n+1,m+1) = sqrt((2*n+1)*(n+m-1)*(n-m-1)/(n-m)/(n+m)/(2*n-3));
    end
end
end