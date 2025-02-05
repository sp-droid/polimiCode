function [Y, tname] = timescaling( T )
% Timescaling for relevant time information in plots
%
% PROTOTYPE
% T, tname = timescaling( T )
%
% INPUT:
% T[nx1] Time vector
%
% OUTPUT:
% Y[nx1] Scaled time vector
% tname String of the relevant time scale
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% -------------------------------------------------------------------------
maxT = max(T);
if maxT<181
    Y = T;
    tname = 's';
elseif maxT<10860
    Y = T/60;
    tname = 'min';
elseif maxT<259200
    Y = T/3600;
    tname = 'h';
elseif maxT<7776000
    Y = T/86400;
    tname = 'days';
elseif maxT<62208000
    Y = T/2592000;
    tname = 'months';
else
    Y = T/31557600;
    tname = 'years';
end
end