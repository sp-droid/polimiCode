function rho = densitySimplified(h)
% INPUTS:
% I had no time to write the descriptions, i'm leaving them like this for now. The scripts on Code/aasignment2 should be self explanatory, in particular the realcomparison.m. It's 100% updated
%
% OUTPUT:
% rho[1]	Density 
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% -------------------------------------------------------------------------

% References used, rho @122km ~= 2e-8, rho @400 ~= 3e-12
h0 = 122;
rho0 = 2e-8;
H0 = -(400-122)/log(3e-12/2e-8);

rho = rho0*exp(-(h-h0)/H0);
end