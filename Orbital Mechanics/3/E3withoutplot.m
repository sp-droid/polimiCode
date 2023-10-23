clear
close all

%% Input data
t1Window = [date2mjd2000([2023;11;1;0;0;0]); date2mjd2000([2025;1;1;0;0;0])];
t2Window = [date2mjd2000([2024;4;1;0;0;0]); date2mjd2000([2025;3;1;0;0;0])];
window1 = (t1Window(2)-t1Window(1))*24*3600;
window2 = (t2Window(2)-t2Window(1))*24*3600;

muSun = astroConstants(4);

p1 = 3;
p2 = 1;
porkchopUL = 50;

%% Porkchop plot and optimization
ngrid = 200;
[T1,T2] = meshgrid(linspace(t1Window(1),t1Window(2),ngrid),...
                   linspace(t2Window(1),t2Window(2),ngrid));
vCostgrid = ones(ngrid,ngrid)*porkchopUL;
deltaTgrid = zeros(ngrid,ngrid);
for i=1:ngrid
    for j=1:ngrid
        deltaTgrid(i,j) = T2(i,j)-T1(i,j);

        [vCost,vc1,vc2,eflag] = transferVcost(T1(i,j),T2(i,j),p1,p2,muSun);
        if (eflag == 0 && vc1<=7)
            vCostgrid(i,j) = min(vCost, porkchopUL);
        end
    end
end
% Best grid point
[Mp,Ip] = min(vCostgrid);
[M,I] = min(Mp);
point = [T1(Ip(I),I); T2(Ip(I),I)];

% Optimization
% Define function to solve for
func = @(x) transferVcost(x(1),x(2),p1,p2,muSun);

% Define constraint (times inside the departure and arrival spans)
lb = [t1Window(1); t2Window(1)];
ub = [t1Window(2); t2Window(2)];

% Solve
options = optimoptions('fmincon', 'Display', 'off');
[point, minVcost] = fmincon(func, point, [], [], [], [], lb, ub, [], options);
t1 = point(1); t2 = point(2);

deltaT = (t2-t1)*24*3600; %deltaT or time of flight (ToF) in [s]
[scaledT, Tname] = timescaling(deltaT);

% Porkchop
for i=1:ngrid
    for j=1:ngrid
        date = round(mjd20002date(T1(i,j)));
        T1(i,j) = datenum(date(1),date(2),date(3),date(4),date(5),date(6));
        date = round(mjd20002date(T2(i,j)));
        T2(i,j) = datenum(date(1),date(2),date(3),date(4),date(5),date(6));
    end
end
date = round(mjd20002date(t1));
point(1) = datenum(date(1),date(2),date(3),date(4),date(5),date(6));
date = round(mjd20002date(t2));
point(2) = datenum(date(1),date(2),date(3),date(4),date(5),date(6));

figure;
% Contour of lines with equal times of flight
contour(T1,T2,deltaTgrid,'ShowText','on', 'LevelStep',60,'EdgeColor','black')
hold on
% Contour of curves with equal deltaV
contour(T1,T2,vCostgrid,'LineWidth',2)
% Lowest deltaV in the grid
scatter(point(1),point(2),30,'filled','red')
title(strcat('Lowest deltaV:', {' '}, num2str(minVcost), '[km/s]', {' '},...
    '(', num2str(scaledT), {' '}, Tname, ')'))

xlabel('Departure date'); ylabel('Arrival date');
cbar = colorbar;
clim([min(vCostgrid,[],'all'), porkchopUL]);
cbar.Title.String = strcat('deltaV [km/s]');
datetick('x',1,'keeplimits')
datetick('y',1,'keeplimits')
grid on;
hold off