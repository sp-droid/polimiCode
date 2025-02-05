function f107 = f107estimation(t, MJD0)
% INPUTS:
% t[1]      Time from the start, in seconds.
% MJD0[1]   Start time in MJD2000
% OUTPUT:
% f107[1]	f107 in solar flux units
%
% CONTRIBUTORS:
% Pablo Arbelo Cabrera
%
% -------------------------------------------------------------------------

MJDcurrent = MJD0+t/86400;
MJDdic2019 = 7287.5; % Solar cycle 25 start date
ellapsed = (MJDcurrent-MJDdic2019)/365.25;

f107 = 140 - 70*cos(2*pi/11*ellapsed);

end