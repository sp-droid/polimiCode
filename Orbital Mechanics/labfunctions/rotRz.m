function R3 = rotRz( angle )
%Returns a rotation matrix around the third axis
%
% PROTOTYPE
% A = rotx( angle )
%
% INPUT:
% angle[1] Rotation value
%
% OUTPUT:
% A[3x3] Rotation matrix of around Z axis
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% VERSIONS
% 2022-10-13: v1
%
% -------------------------------------------------------------------------
R3 = [cos(angle) sin(angle) 0;
      -sin(angle) cos(angle) 0;
      0 0 1];
end