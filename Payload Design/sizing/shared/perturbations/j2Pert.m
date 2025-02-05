function a = j2Pert(r, J2, R, mu)
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

rnorm = norm(r);

a = 5*(r(3)/rnorm)^2;
a = [   r(1)/rnorm*(a-1)
		r(2)/rnorm*(a-1)
		r(3)/rnorm*(a-3)];
a = a*1.5*J2*mu*(R/rnorm^2)^2;
end