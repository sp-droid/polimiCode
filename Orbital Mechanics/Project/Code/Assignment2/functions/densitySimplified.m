function rho = densitySimplified(h)

% References used, rho @122km ~= 2e-8, rho @400 ~= 3e-12
h0 = 122;
rho0 = 2e-8;
H0 = -(400-122)/log(3e-12/2e-8);

rho = rho0*exp(-(h-h0)/H0);
end