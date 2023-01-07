function R2 = rotRy( angle )
%Returns a rotation matrix around the second axis
%
% PROTOTYPE
% A = rotx( angle )
%
% INPUT:
% angle[1] Rotation value
%
% OUTPUT:
% A[3x3] Rotation matrix of around Y axis
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% -------------------------------------------------------------------------
R2 = [cos(angle) 0 sin(angle);
      0 1 0;
      -sin(angle) 0 cos(angle)];
end