function dv = lambert_dv(date, ID1, ID2, mu, min_type)
% 
% This function computes the total dv for a mission from body 1 to body 2
% using a Lambert solver.
% 
% INPUTS:
%   date    [2x1]   Column vector containing departure and arrival dates in
%                   mjd2000.
%   ID1             Natural number representing departure body 1
%   ID2             Natural number represeting arrival body 2
%   mu              Planetary constant
%   min_type        Variable defining whether the efficiency of the mission in terms of:
%                       - 0: Both injection manoeuvre (from orbit 1 to transfer arc)
%                            and arrival manoeuvre (from trasnfer arc to orbit 2)
%                       - 1: Only injection manoeuvre (from orbit 1 to transfer arc 1) 
%                       - 2: Only arrival manoeuvre (from trasnfer arc to orbit 2) 
% 
% OUTPUT:
%   dv      [1x1]   Mission cost
% 
% VERSIONS:
%   2022-11-27
% 
% CONTRIBUTORS:
%   Alessandro Michelazzi

date1 = date(1,:);
date2 = date(2,:);

t1 = date1 * (24*60*60);
t2 = date2 * (24*60*60);
tof = t2 - t1;

if ID1 < 12
    [kep1,~] = uplanet(date1, ID1);
else
    [kep1,~,~] = ephNEO(date1,ID1);
end
if ID2 < 12
    [kep2,~] = uplanet(date2, ID2);
else
    [kep2,~,~] = ephNEO(date2, ID2);
end

[r1, v1] = kep2car(kep1, mu);
[r2, v2] = kep2car(kep2, mu);

[~,~,~,ERR,v1_t,v2_t,~,~] = lambertMR(r1, r2, tof, mu, 0, 0, 0, 2);
dv_1 = v1_t' - v1;
dv_2 = v2 - v2_t';

dv = NaN;
if ERR == 0
    switch min_type
        case 0
            dv = abs(norm(dv_1)) + abs(norm(dv_2));
        case 1
            dv = abs(norm(dv_1));
        case 2
            dv = abs(norm(dv_2));
    end
end
