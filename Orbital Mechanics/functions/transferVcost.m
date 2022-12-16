function [vCost, vc1, vc2, eflag] = transferVcost( t1, t2, body1, body2, mu )
% Simple orbit determination, for fast plotting (low accuracy)
%
% PROTOTYPE
% vCost = transferT( t1, t2, body1, body2, mu )
%
% INPUT:
% r[3x1] Position ( rx, ry, rz ) [ L ]
% v[3x1] Velocity ( vx, vy, vz ) [ L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
% ngrid[1] Number of points (low number is recommended)
% time[1] Length of time plotted. If not specified, it reverts to 1 orbital period
%
% OUTPUT:
% vCost[1] Total deltaV
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% VERSIONS
% 2018-09-26: First version
%
% -------------------------------------------------------------------------
[kep, ~] = uplanet(t1, body1);
[r1, v1] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6), mu, 'rad');
[kep, ~] = uplanet(t2, body2);
[r2, v2] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6), mu, 'rad');

deltaT = (t2-t1)*24*3600;
if (deltaT < 1)
	vc1 = []; vc2 = []; vCost = [];
	eflag = 5;
	return
end
[~,~,~,eflag,vt1,vt2,~,~] = lambertMR( r1, r2, deltaT, mu, 0, 0, 0 );

vc1 = vecnorm(vt1'-v1);
vc2 = vecnorm(vt2'-v2);
vCost = vc1+vc2;

end