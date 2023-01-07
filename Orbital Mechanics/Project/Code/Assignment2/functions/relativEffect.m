function a = relativEffect(r, v, mu)
% INPUTS:
% I had no time to write the descriptions, i'm leaving them like this for now. The scripts on Code/aasignment2 should be self explanatory, in particular the realcomparison.m. It's 100% updated
%
% OUTPUT:
% a[3x1] Perturbing acceleration
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% -------------------------------------------------------------------------

c = 299792.458; % Speed of light in km/s

rnorm = norm(r);

a = - mu/rnorm^3*(((2/c)^2*mu/rnorm-dot(v,v)/c^2)*r+(2/c)^2*dot(r,v)*v);
end