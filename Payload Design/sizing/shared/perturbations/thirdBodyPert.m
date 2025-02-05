function a = thirdBodyPert(r, rThirdBody, muThirdBody)
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

a = muThirdBody*((rThirdBody-r)/norm(rThirdBody-r)^3-rThirdBody/norm(rThirdBody)^3);
end