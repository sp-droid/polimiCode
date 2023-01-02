function R1 = rotRx( angle )
%Returns a rotation matrix around the first axis
%
% PROTOTYPE
% A = rotx( angle )
%
% INPUT:
% angle[1] Rotation value
%
% OUTPUT:
% A[3x3] Rotation matrix of around X axis
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% VERSIONS
% 2022-10-13: v1
%
% -------------------------------------------------------------------------
R1 = [1 0 0;
    0 cos(angle) sin(angle);
    0 -sin(angle) cos(angle)];
end