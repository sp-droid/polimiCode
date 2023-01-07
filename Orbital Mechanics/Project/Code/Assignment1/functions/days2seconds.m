function [seconds] = days2seconds (days)
% 
% Seconds from days.
% 
% INPUT:
%   days        [day]
% 
% OUTPUT:
%   seconds     [s]
% 
% VERSIONS:
%   2022-11-29
% 
% CONTRIBUTORS:
%   Alessandro Michelazzi

seconds = days * (24*60*60);