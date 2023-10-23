function tpar_min = min_parabolic_time(dep_win,fb_win,arr_win,ID1,ID2,ID3,mu)
% 
% This function computes the minimum parbolic flight time for two transfer arcs: 
% from body 1 to body 2 and from body 2 to body 3 for the given departure, flyby 
% and arrival time windows. The output is the smaller time between the two computed.
% 
% INPUTS:
%   dep_win     [1xl]   Vector of possible days (mjd2000) for departure
%   fb_win      [1xm]   Vector of possible days (mjd2000) for fly-by
%   arr_win     [1xn]   Vector of possible days  (mjd2000) for arrival
%   ID1         [1x1]   ID of departure body (Natural number)
%   ID2         [1x1]   ID of flyby body (Natural number)
%   ID3         [1x1]   ID of arrival body (Natural number)
%   mu          [1x1]   Planetary constant [km^3/s^2]
% 
% OUTPUT:
%   tpar_min    [1x1]   Smallest parabolic flight time found [s]
% 
% VERSIONS:
%   2023-01-03: First version
% 
% CONTRIBUTORS:
%   Alessandro Michelazzi

l = length(dep_win); m = length(fb_win); n = length(arr_win);
tspan1 = dep_win * (24*60*60); % [s]
tspan2 = fb_win * (24*60*60); % [s]
tspan3 = arr_win * (24*60*60); % [s]

% Velocity initialization matrix
tpar_matrix = NaN(l,m,n);

% Cycle on departure dates
for i = 1:l
    if ID1 < 12
        [kep1_i, ~] = uplanet(dep_win(i),ID1);
    else
        [kep1_i,~,~] = ephNEO(dep_win(i),ID1);
    end

    % Cycle on fly-by dates
    for j = 1:m
        if tspan2(j) < tspan1(i) % Fly-by date is before departure date            
            continue
        end
        if ID2 < 12
            [kep2_j, ~] = uplanet(fb_win(j),ID2);
        else
            [kep2_j,~,~] = ephNEO(fb_win(j),ID2);
        end 
        [r1_i, ~] = kep2car(kep1_i, mu);
        [r2_j, ~] = kep2car(kep2_j, mu); 

        % First transfer arc
        tof1 = tspan2(j) - tspan1(i);         
        [~,~,~,ERR1,~,~,tpar1,~] = lambertMR(r1_i, r2_j, tof1, mu, 0, 0, 0, 2);

        % Cycle on arrival dates
        for k = 1:n
            if tspan3(k) < tspan2(j) % Arrival date is before fly-by date
                continue
            end
            if ID3 < 12
                [kep3_k, ~] = uplanet(arr_win(k),ID3);
            else
                [kep3_k,~,~] = ephNEO(arr_win(k),ID3);
            end 
            [r3_k, ~] = kep2car(kep3_k, mu);
            
            % Second transfer arc
            tof2 = tspan3(k) - tspan2(j);        
            [~,~,~,ERR2,~,~,tpar2,~] = lambertMR(r2_j, r3_k, tof2, mu, 0, 0, 0, 2); 
            
            if ~ all([ERR1, ERR2])
                tpar_matrix(i,j,k) = min(tpar1,tpar2); 
            end
        end
    end
end

tpar_min = min(min(min(tpar_matrix)));