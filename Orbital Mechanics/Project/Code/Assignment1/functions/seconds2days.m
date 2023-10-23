function [days] = seconds2days (seconds)
% 
% Days from seconds.
% 
% INPUT:
%   seconds     [s]
% 
% OUTPUT:
%   days        [day]
% 
% VERSIONS:
%   2022-11-29
% 
% CONTRIBUTORS:
%   Alessandro Michelazzi

days = seconds / (24*60*60);