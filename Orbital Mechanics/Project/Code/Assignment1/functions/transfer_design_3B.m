function [dv_min,date] = transfer_design_3B(dep_win,fb_win,arr_win,IDs,v_inf,h,opts)
% 
% This function computes the optimal transfer from body 1 to body 3, with a powered gravity assist 
% fly-by at body 2.
% A first coarse search is executed. For the best regions (the one with the
% lowest dv) a second rifined grid search is executed.
% 
% INPUT:
%   dep_win     [1xl]   Vector of possible days (mjd2000) for departure
%                       l MUST be an odd number
%   fb_win      [1xm]   Vector of possible days (mjd2000) for fly-by
%                       l MUST be an odd number
%   arr_win     [1xn]   Vector of possible days  (mjd2000) for arrival
%                       l MUST be an odd number
%   IDs         [1x3]   Vector of natural numbers representing in order departure body,
%                       fly-by body and arrival body:
%                           1: Mercury          6: Saturn
%                           2: Venus            7: Uranus
%                           3: Earth            8: Neptune
%                           4: Mars             9: Pluto
%                           5: Jupiter          10: Moon
%                           or ID > 12 for NEOs
%   v_inf               Maximum excess velocity
%   h                   Atmosphere scale height of fly-by body
%   opts                Struct containing function options
%       - orbitType1    Logical variable defining whether first transfer is
%                           0: direct transfer from R1 to R2 (counterclockwise)
%                           1: retrograde transfer from R1 to R2 (clockwise)
%       - orbitType2    Logical variable defining whether second transfer is
%                           0: direct transfer from R1 to R2 (counterclockwise)
%                           1: retrograde transfer from R1 to R2 (clockwise)
%       - RegionsNumber Number of best regions for search. For more than 2 best regions
%                       time discretisation is lowered to reduce computational cost.
%                       Default number is 1: function continues the refined grid search only for
%                       the minimum dv fund in first grid search.
%       - ToFmax1       Maximum Time of Flight desired for first transfer arc
%       - ToFmax2       Maximum Time of Flight desired for second transfer arc
%       - tpar_min      Minimum parabolic flight time between ID1-ID2 or ID2-ID3
% 
% OUTPUT:
%   dv_min      [1x1]   Minimum cost for the given mission
%   date        [3x1]   Column vecotor containing departure, fly-by and
%                       arrival dates expressed in mjd2000
% 
% VERSIONS
%   2022-12-23: First version
% 
% CONTRIBUTORS:
%   Alessandro Michelazzi


l = length(dep_win); m = length(fb_win); n = length(arr_win);
tspan1 = dep_win * (24*60*60); % [s]
tspan2 = fb_win * (24*60*60); % [s]
tspan3 = arr_win * (24*60*60); % [s]

if all([IDs(1) > 0, IDs(1) ~= 11])
    ID1 = IDs(1);
else
    error('ID1 not valid.')
end
if all([IDs(2) > 0, IDs(2) ~= 11])
    ID2 = IDs(2);
else
    error('ID2 not valid.')
end
if all([IDs(3) > 0, IDs(3) ~= 11])
    ID3 = IDs(3);
else
    error('ID3 not valid.')
end

if any([nargin < 7, ~isfield(opts,'tpar_min')])
    tpar_min = 0;
else
    tpar_min = opts.tpar_min;
end

mu = astroConstants(4); % Sun planetary constant

if any([nargin < 7, ~isfield(opts,'orbitType1')])
    orbitType1 = 0;
else
    orbitType1 = opts.orbitType1;
end
if any([nargin < 7, ~isfield(opts,'orbitType2')])
    orbitType2 = 0;
else
    orbitType2 = opts.orbitType2;
end

Refine = 1;
if any([nargin < 7, ~isfield(opts,'ToFmax1')])
    ToFmax1 = tspan3(end) - tspan1(1);
else
    ToFmax1 = opts.ToFmax1;
    Refine = 0;
end
if any([nargin < 7, ~isfield(opts,'ToFmax2')])
    ToFmax2 = tspan3(end) - tspan1(1);
else
    ToFmax2 = opts.ToFmax2;
    Refine = 0;
end

% Velocity initialization matrix
deltaV_matrix = NaN(l,m,n);

% Search for an approximate solution
% Cycle on departure dates
for i = 1:l
    if ID1 < 12
        [kep1_i, ~] = uplanet(dep_win(i),ID1);
    else
        [kep1_i,~,~] = ephNEO(dep_win(i),ID1);
    end
    % Cycle on fly-by dates
    for j = 1:m
        tof1 = tspan2(j) - tspan1(i);
        if any([tof1 < tpar_min, tof1 > ToFmax1])
            % Fly-by date either before departure date or before the minimum ToF possible 
            continue
        end
        if ID2 < 12
            [kep2_j, ~] = uplanet(fb_win(j),ID2);
        else
            [kep2_j,~,~] = ephNEO(fb_win(j),ID2);
        end 
        [r1_i, v1_i] = kep2car(kep1_i, mu);
        [r2_j, v2_j] = kep2car(kep2_j, mu); 
        % First transfer arc                 
        [~,~,~,ERR1,v1_t1,v2_t1,~,~] = lambertMR(r1_i, r2_j, tof1, mu, orbitType1, 0, 0, 2);        
        v1_t1 = v1_t1'; v2_t1 = v2_t1';
        deltaV_1 = norm(v1_t1 - v1_i); % Cost first manoeuvre
        % Cycle on arrival dates
        for k = 1:n
            tof2 = tspan3(k) - tspan2(j); 
            if any([tof2 < tpar_min, tof2 > ToFmax2]) 
                % Arrival date is before fly-by date or before the minimum ToF possible
                continue
            end
            if ID3 < 12
                [kep3_k, ~] = uplanet(arr_win(k),ID3);
            else
                [kep3_k,~,~] = ephNEO(arr_win(k),ID3);
            end 
            [r3_k, v3_k] = kep2car(kep3_k, mu);          
            % Second transfer arc       
            [~,~,~,ERR2,v1_t2,v2_t2,~,~] = lambertMR(r2_j, r3_k, tof2, mu, orbitType2, 0, 0, 2);            
            v1_t2 = v1_t2'; v2_t2 = v2_t2';
            deltaV_2 = norm(v3_k - v2_t2); % Cost second manoeuvre
            % Fly-by
            v_inf_m = v2_t1 - v2_j;
            v_inf_p = v1_t2 - v2_j;
            [r_p, ~, delta_v, ~, ~, ~] = power_gravity_assist(v_inf_m, v_inf_p, ID2, h);        
            
            if ~ all([ERR1, ERR2])
                if all([deltaV_1 < v_inf, deltaV_2 < v_inf, isfinite(r_p')])
                     deltaV_matrix(i,j,k) = abs(deltaV_1) + abs(deltaV_2) + abs(delta_v);                   
                end
            end
        end
    end
end

if any([nargin < 7, ~isfield(opts,'RegionsNumber')])
    p_max = 1;
else
    p_max = opts.RegionsNumber;
end
% Initializing parameters for second grid search
deltaV_matrix_copy = deltaV_matrix;
deltaV_min_vec = NaN(1,p_max);
cord = NaN(3,p_max);
for p = 1:p_max 
    deltaV_min = min(min(min(deltaV_matrix_copy)));
    deltaV_max = max(max(max(deltaV_matrix_copy)));
    % Finding coordinates of the deltaV_min found
    for i = 1:l
        for j = 1:m
            for k = 1:n
                if deltaV_matrix_copy(i,j,k) == deltaV_min
                    deltaV_min_vec(p) = deltaV_min; % Vector containing the first-p minimum dvs
                    cord(:,p) = [i; j; k];
                    deltaV_matrix_copy(i,j,k) = deltaV_max;
                end
            end
        end
    end
end 

% Search for a more accurate solution
mjd2000_matrix = NaN(3,p_max); true_date_matrix = NaN(3,p_max);
if p_max > 2
    if rem(ceil(l/1.5),2)
        l = ceil(l/1.5);
    else
        l = floor(l/1.5);
    end 
    if rem(ceil(m/1.5),2)
        m = ceil(m/1.5);
    else
        m = floor(m/1.5);
    end 
    if rem(ceil(n/1.5),2)
        n = ceil(n/1.5);
    else
        n = floor(n/1.5);
    end
    % To decrease computational cost
end
    
for p = 1:p_max
    dv_vec = [0, deltaV_min];
    tspan1_p = tspan1; tspan2_p = tspan2; tspan3_p = tspan3;
    for q = 1
        % Initializing deltaV matrix
        deltaV_matrix_q = NaN(l,m,n);
        dt1 = round(tspan1_p(end)-tspan1_p(end-8)); 
        dt2 = round(tspan2_p(end)-tspan2_p(end-4)); dt3 = round(tspan3_p(end)-tspan3_p(end-4));
        % Decreasing time range around minimum solution at each iteration
        tspan1_p = linspace(tspan1_p(cord(1,p)) - dt1, tspan1_p(cord(1,p)) + dt1, l);
        tspan2_p = linspace(tspan2_p(cord(2,p)) - dt2, tspan2_p(cord(2,p)) + dt2, m);
        tspan3_p = linspace(tspan3_p(cord(3,p)) - dt3, tspan3_p(cord(3,p)) + dt3, n);
        if all([p_max > 2, q == 1])
            for ip = floor(l/2)-1:ceil(l/2)+1
                if tspan1_p(ip) == tspan1(cord(1,p))
                    cord(1,p) = ip;
                end
            end
            for jp = floor(m/2)-1:ceil(m/2)+1
                if tspan2_p(jp) == tspan2(cord(2,p))
                    cord(2,p) = jp;
                end
            end            
            for kp = floor(n/2)-1:ceil(n/2)+1
                if tspan3_p(kp) == tspan3(cord(3,p))
                    cord(3,p) = kp;
                end
            end
        end
        mjd2000_dep = tspan1_p / (24*60*60);
        mjd2000_GA = tspan2_p / (24*60*60);
        mjd2000_arr = tspan3_p / (24*60*60);
        % Cycle on departure dates
        for i = 1:length(tspan1_p)
            if ID1 < 12
                [kep1_i, ~] = uplanet(mjd2000_dep(i),ID1);
            else
                [kep1_i,~,~] = ephNEO(mjd2000_dep(i),ID1);
            end    
            % Cycle on fly-by dates
            for j = 1:length(tspan2_p)
                tof1 = tspan2_p(j) - tspan1_p(i);
                if any([tof1 < tpar_min, tof1 > ToFmax1])
                    % Fly-by date either before departure date or before the minimum ToF possible            
                    continue
                end
                if ID2 < 12
                    [kep2_j, ~] = uplanet(mjd2000_GA(j),ID2);
                else
                    [kep2_j,~,~] = ephNEO(mjd2000_GA(j),ID2);
                end      
                [r1_i, v1_i] = kep2car(kep1_i, mu);
                [r2_j, v2_j] = kep2car(kep2_j, mu);
                % First transfer arc       
                [~,~,~,ERR1,v1_t1,v2_t1,~,~] = lambertMR(r1_i, r2_j, tof1, mu, orbitType1, 0, 0, 2);            
                v1_t1 = v1_t1'; v2_t1 = v2_t1';
                deltaV_1 = norm(v1_t1 - v1_i);             
                % Cycle on arrival dates
                for k = 1:length(tspan3_p)
                    tof2 = tspan3_p(k) - tspan2_p(j); 
                    if any([tof2 < tpar_min, tof2 > ToFmax2])
                        % Arrival date is before fly-by date or before the minimum ToF possible
                        continue
                    end
                    if ID3 < 12
                        [kep3_k, ~] = uplanet(mjd2000_arr(k),ID3);
                    else
                        [kep3_k,~,~] = ephNEO(mjd2000_arr(k),ID3);
                    end 
                    [r3_k, v3_k] = kep2car(kep3_k, mu);
                    % Second transfer arc      
                    [~,~,~,ERR2,v1_t2,v2_t2,~,~] = lambertMR(r2_j, r3_k, tof2, mu, orbitType2, 0, 0, 2);                
                    v1_t2 = v1_t2'; v2_t2 = v2_t2';
                    deltaV_2 = norm(v3_k - v2_t2);
                    % Fly-by
                    v_inf_m = v2_t1 - v2_j;
                    v_inf_p = v1_t2 - v2_j;
                    [r_p, ~, delta_v, ~, ~, ~] = power_gravity_assist(v_inf_m, v_inf_p, ID2, h);
    
                    if ~ all([ERR1, ERR2])
                        if all([deltaV_1 < v_inf, deltaV_2 < v_inf, isfinite(r_p')])
                             deltaV_matrix_q(i,j,k) = abs(deltaV_1) + abs(deltaV_2) + abs(delta_v);                         
                        end
                    end
                end
            end
        end
        dv_q = min(min(min(deltaV_matrix_q))); % Minimum manoeuvre cost at iteration q
        if dv_q <= dv_vec(2)
            dv_vec(1) = dv_vec(2);
            dv_vec(2) = dv_q;
        else
            continue
        end
        % Retrive position of minimum dv in matrix
        for i = 1:length(tspan1_p)
            for j = 1:length(tspan2_p)
                for k = 1:length(tspan3_p)
                    if deltaV_matrix_q(i,j,k) == dv_q
                        cord(:,p) = [i; j; k];
                    end
                end
            end
        end
        % Saving date for minimum dv
        mjd2000_matrix(:,p) = [mjd2000_dep(cord(1,p));mjd2000_GA(cord(2,p)); mjd2000_arr(cord(3,p))];        
        % Comparing dv_min found in the last two iterations
        diff = dv_vec(1) - dv_vec(2);
        toll = 0.1;
        if all([diff < toll, diff >= 0])
            break
        end
    end

    deltaV_min_vec(p) = dv_q; % Vector containing minimum manouvre cost for each region

    if Refine
        % Refining the solution
        dep_approx = mjd2000_matrix(1,p);
        GA_approx = mjd2000_matrix(2,p);
        arr_approx = mjd2000_matrix(3,p);
        dvmin = @(date_mjd2000) lambert_dv_3B (date_mjd2000, ID1, ID2, ID3, mu, opts, h);
        date0 = [dep_approx; GA_approx; arr_approx];
        
        it = 0;
        exitflag = 5;
        while all([exitflag==5, it < 5])
            [true_date, delta_v, exitflag] = fminunc(dvmin, date0);
            date0 = true_date;
            it = it + 1;
        end 
        true_date_matrix(:,p) = true_date; deltaV_min_vec(p) = delta_v;
    else
        true_date_matrix(:,p) = mjd2000_matrix(:,p);
    end
end

dv_min = min(deltaV_min_vec);
for p = 1:p_max
    if deltaV_min_vec(p) == dv_min
        date = true_date_matrix(:,p);
    end
end













