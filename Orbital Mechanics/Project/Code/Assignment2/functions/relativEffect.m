function a = relativEffect(r, v, mu)

c = 299792.458; % Speed of light in km/s

rnorm = norm(r);

a = - mu/rnorm^3*(((2/c)^2*mu/rnorm-dot(v,v)/c^2)*r+(2/c)^2*dot(r,v)*v);
end