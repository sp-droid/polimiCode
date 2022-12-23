clear
close all

ngrid = 86400;
muEarth = astroConstants(13);
r0 = [ 7495.3;0.0000;0.0000 ];
v0 = [ 0.0000;0.2686;-7.3239 ];
y0 = [ r0; v0 ];

tic
[yode113, t] = timed2BP(y0, muEarth, [], ngrid, 86400);
yode113 = yode113';
toc

deltaT = t(2)-t(1);

% Symbolic
tic
% Gravitational potential and derivatives
syms x(ti) y(ti) z(ti) vx(ti) vy(ti) vz(ti) xt yt zt vxt vyt vzt rnorm
fx = -muEarth*x/(x^2+y^2+z^2)^(3/2); Ux = fx;
fy = -muEarth*y/(x^2+y^2+z^2)^(3/2); Uy = fx;
fz = -muEarth*z/(x^2+y^2+z^2)^(3/2); Uz = fx;
%feval(derivs{1,j},1,1,1,1,1,1,0)
norder = 8;
derivs = cell(1,norder);
deriv = subs(Ux,...
    [x(ti),y(ti),z(ti),vx(ti),vy(ti),vz(ti)], ...
    [xt,yt,zt,vxt,vyt,vzt]);
derivs{1,1} = matlabFunction(deriv,'Vars',{xt,yt,zt,vxt,vyt,vzt,ti});
for j=2:norder
    Ux = subs(...
        diff(Ux,ti,1), ...
        [diff(x,ti),diff(y,ti),diff(z,ti),diff(vx,ti),diff(vy,ti),diff(vz,ti)], ...
        [vx, vy, vz, fx, fy, fz]);
    deriv = subs(Ux,...
        [x(ti),y(ti),z(ti),vx(ti),vy(ti),vz(ti)], ...
        [xt,yt,zt,vxt,vyt,vzt]);
    derivs{1,j} = matlabFunction(deriv,'Vars',{xt,yt,zt,vxt,vyt,vzt,ti});
end
toc
% Taylor expansion substitution
%https://www.esa.int/gsp/ACT/projects/taylorpropagation/
%https://www.researchgate.net/publication/264235964_Numerical_solution_of_Kepler's_planetary_equation_for_orbit_propagation_Taylor_series_approach
yN = zeros(6,ngrid); yN(:,1) = y0;
tic
for i=1:ngrid-1
    yN(1:3,i+1) = yN(1:3,i)+deltaT*yN(4:6,i);
    yN(4:6,i+1) = yN(4:6,i);

    x = yN(1,i); y = yN(2,i); z = yN(3,i);
    vx = yN(4,i); vy = yN(5,i); vz = yN(6,i);
    for j=1:norder
        U = derivs{1,j};
        term = [U(x,y,z,vx,vy,vz,0); U(y,x,z,vy,vx,vz,0); U(z,y,x,vz,vy,vx,0)];

        yN(1:3,i+1) = yN(1:3,i+1)+deltaT^(j+1)/factorial(j+1)*term;
        yN(4:6,i+1) = yN(4:6,i+1)+deltaT^(j)/factorial(j)*term;
    end
end
toc

% Analytical
yT = zeros(6,ngrid); yT(:,1) = y0;

tic
[a,e,i,bOmega,sOmega,theta0] = car2kep(r0,v0,muEarth,'rad');
E0 = 2*atan2(sqrt((1-e)/(1+e))*sin(theta0/2),cos(theta0/2));
t0 = t(1);
for j=1:ngrid-1
    E = keplerSolver(t(j+1), e, a, muEarth, t0, E0);
    theta = 2*atan2(sqrt((1+e)/(1-e))*sin(E/2),cos(E/2));
    [yT(1:3,j+1), yT(4:6,j+1)] = kep2car(a,e,i,bOmega,sOmega,theta,muEarth,'rad');
end
toc

figure;
plot(t,yT(1,:),'LineWidth',2)
hold on
plot(t,yode113(1,:),'LineWidth',1)
plot(t,yN(1,:),'LineWidth',2)
legend('Analytical','Ode113','Taylor')
grid on
hold off

figure;
plot(t,yT(1,:)-yode113(1,:),'LineWidth',1)
hold on
plot(t,yT(1,:)-yN(1,:),'LineWidth',1)
legend('Ode113','Taylor')
grid on
hold off