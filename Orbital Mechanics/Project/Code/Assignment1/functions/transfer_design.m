function [dep_date, arr_date, delta_v, delta_t, r1, v1, r2, v2, dv1, dv2] = transfer_design(dep_window, arr_window, ID1, ID2, mu, min_type, v_infty)
% 
% Given a departure window on planet 1 and an arrival window on planet 2,
% the function computes the most efficient transfer in terms of cost.
% WARNING: this function was implemented for Earth-Saturn and
% Saturn-NEO85 transfers only.
% 
% INPUT:
%   dep_window      Departure window on planet 1. It's an [2x1] array containing
%                   on row 1 the earliest departure date and on row 2 the
%                   latest departure date expressed in mjd2000.
%   arr_window      Arrival window on planet 2. It's an [2x1] array containing
%                   on row 1 the earliest arrival date and on row 2 the
%                   latest arrival date expressed in mjd2000.
%   ID1             ID of planet 1 for function uplanet(mjd2000, ibody).
%   ID2             ID of planet 2 for function uplanet(mjd2000, ibody).
%   mu              Planetary constant (mu = mass * G) [km^3/s^2].
%   min_type        Variable defining wheter the efficiency of the mission in terms of:
%                       - 0: Both injection manoeuvre (from orbit 1 to transfer arc)
%                            and arrival manoeuvre (from trasnfer arc to orbit 2)
%                       - 1: Only injection manoeuvre (from orbit 1 to transfer arc 1) 
%                       - 2: Only arrival manoeuvre (from trasnfer arc to orbit 2) 
%   v_infty         Maximum velocity that can be given by launcher.  [km/s]
% 
% OUTPUT:
%   dep_date        Optimal departure date in mjd2000.
%   arr_date        Optimal arrival date in mjd2000.
%   delta_v         Minimum delta_velocity for the given problem.
%   delta_t         Time of Flight required for the given problem.
%   r1              Position vecotr of the 1st manoeuvre.
%   v1              Velocity vector on transfer orbit after 1st manouvre.
%   r2              Position vector of the 2nd manoeuvre.
%   v2              Velocity vector on transfer orbit before 2nd manouvre.
%   dv1             Cost injection manoeuvre
%   dv2             Cost injection manoeuvre
% 
% VERSIONS:
%   2022-11-27
% 
% CONTRIBUTORS:
%   Alessandro Michelazzi

if nargin < 7
    v_infty = 100;
end

n = 1000;
mjd2000_vec1 = linspace(dep_window(1), dep_window(2), n);
mjd2000_vec2 = linspace(arr_window(1), arr_window(2), n);
tspan1 = mjd2000_vec1 * (24*60*60); 
tspan2 = mjd2000_vec2 * (24*60*60);

% Departure planet
if ID1 < 12
    [kep1, ~] = uplanet(mjd2000_vec1(1),ID1);
else
    [kep1,~,~] = ephNEO(mjd2000_vec1(1),ID1);
end
T1 = 2*pi * sqrt(kep1(1)^3 / mu); T1_mjd = seconds2days(T1); % Period departure orbit
% Arrival planet
if ID2 < 12
    [kep2, ~] = uplanet(mjd2000_vec2(1),ID2);
else
    [kep2,~,~] = ephNEO(mjd2000_vec2(1),ID2);
end
T2 = 2*pi * sqrt(kep2(1)^3 / mu); T2_mjd = seconds2days(T2); % Period arrival orbit

% Initialising matrix
deltaV_matrix = NaN(n,n);

% Search for an approximate solution
for i = 1:n
    if ID1 < 12
        [kep1_i, ~] = uplanet(mjd2000_vec1(i),ID1);
    else
        [kep1_i,~,~] = ephNEO(mjd2000_vec1(i),ID1);
    end
    for j = 1:n
        if ID2 < 12
            [kep2_i, ~] = uplanet(mjd2000_vec2(j),ID2);
        else
            [kep2_i,~,~] = ephNEO(mjd2000_vec2(j),ID2);
        end 
        [r1_i, v1_i] = kep2car(kep1_i, mu);
        [r2_j, v2_j] = kep2car(kep2_i, mu);

        tof = tspan2(j) - tspan1(i);
        if tof > 0
            [~,~,~,ERR,v1_t,v2_t,~,~] = lambertMR(r1_i, r2_j, tof, mu, 0, 0, 0, 2);
        else
            continue
        end
        v1_t = v1_t'; v2_t = v2_t';
        delta_v1 = abs(norm((v1_t - v1_i)));
        delta_v2 = abs(norm((v2_j - v2_t)));
        deltaV_tot = delta_v1 + delta_v2;
        switch min_type
            case 0
                if all([ERR == 0, deltaV_tot <= v_infty])
                    deltaV_matrix(j,i) = deltaV_tot;
                end
            case 1
                if all([ERR == 0, delta_v1 <= v_infty])
                    deltaV_matrix(j,i) = delta_v1;
                end
            case 2
                if all([ERR == 0, delta_v2 <= v_infty])
                    deltaV_matrix(j,i) = delta_v2;
                end
        end
    end
end

v_min = min(min(deltaV_matrix));

% For porkchop plots Earth-Saturn or Saturn-NEO85
if ID1 == 6
    x_contour = (mjd2000_vec1 - mjd2000_vec1(1)) .* 360/T1_mjd + rad2deg(kep1(6));
else
    x_contour = (mjd2000_vec1 - mjd2000_vec1(1)) ./ T1_mjd;
end
if ID2 == 6
    y_contour = (mjd2000_vec2 - mjd2000_vec2(1)) .* 360/T2_mjd + rad2deg(kep2(6));
else
    y_contour = (mjd2000_vec2 - mjd2000_vec2(1)) ./ T2_mjd;
end

% Porkchop plot Earth-Saturn or Saturn-NEO85
figure()
hold on
contourf(x_contour, y_contour, deltaV_matrix)
colorbar
if ID1 == 6
    xlabel('Saturn true anomaly  [deg]');
else
    xlabel('Number of periods of Earth from earliest departure  [-]');
end
if ID2 == 6
    ylabel('Saturn true anomaly  [deg]');
else
    ylabel('Number of periods of NEO 85 from earliest departure  [-]');
end
grid on

% Search for precise solution
for i = 1:n
    for j = 1:n
        if deltaV_matrix(j,i) == v_min
            dep_approx = mjd2000_vec1(i);
            arr_approx = mjd2000_vec2(j);
            delta_t = tspan2(j) - tspan1(i);
        end
    end
end

dvmin = @(date_mjd2000) lambert_dv(date_mjd2000, ID1, ID2, mu, min_type);

date0 = [dep_approx; arr_approx];

lb = [dep_window(1); arr_window(1)];
ub = [dep_window(2); arr_window(2)];
[true_date, delta_v] = fmincon(dvmin, date0, [], [], [], [], lb, ub);
tof = days2seconds(true_date(2) - true_date(1));

% Plotting minimum found
if ID1 == 6
    xplot = (true_date(1) - mjd2000_vec1(1)) * 360/T1_mjd + rad2deg(kep1(6));
else
    xplot = (true_date(1) - mjd2000_vec1(1)) / T1_mjd;
end
if ID2 == 6
    yplot = (true_date(2) - mjd2000_vec2(1)) * 360/T2_mjd + rad2deg(kep2(6));
else
    yplot = (true_date(2) - mjd2000_vec2(1)) / T2_mjd;
end
plot(xplot,yplot,'r','LineWidth', 5, 'Marker','o')
legend('','Most efficient transfer','Location','best')

% Computing r1, v1, r2, v2, dv1, dv2 for solution found
if ID1 < 12
    [kep1,~] = uplanet(true_date(1), ID1);
else
    [kep1,~,~] = ephNEO(true_date(1),ID1);
end
if ID2 < 12
    [kep2,~] = uplanet(true_date(2), ID2);
else
    [kep2,~,~] = ephNEO(true_date(2), ID2);
end

[r1, v_1] = kep2car(kep1, mu);
[r2, v_2] = kep2car(kep2, mu);

[~,~,~,~,v1_t,v2_t,~,~] = lambertMR(r1, r2, tof, mu, 0, 0, 0, 2);
v1 = v1_t';
v2 = v2_t';
dv1 = norm(v1 - v_1);
dv2 = norm(v_2 - v2);

dep_date = true_date(1);
arr_date = true_date(2);
