function [newlong,newlat] = removeLonLines( long, lat )
%Adds NaN values every time longitude wraps around the planet for better line plots
%
% PROTOTYPE
% long = rotx( long )
%
% INPUT:
% long[nx1] Longitude in degrees
% lat[nx1] Latitude in any unit
%
% OUTPUT:
% newlong[(n+k)x1] Longitude with intermediate NaN values
% newlat[nx1] Latitude with intermediate NaN values
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
% -------------------------------------------------------------------------
newlong = [];
newlat = [];
for j=1:length(long)-1
    if long(j)>150 && long(j+1)<30
        newlong = [newlong; long(j); NaN];
        newlat = [newlat; lat(j); NaN];
    elseif long(j)<30 && long(j+1)>150
        newlong = [newlong; long(j); NaN];
        newlat = [newlat; lat(j); NaN];
    else
        newlong = [newlong; long(j)];
        newlat = [newlat; lat(j)];
    end
end
newlong = [newlong; long(end)];
newlat = [newlat; lat(j)];
end