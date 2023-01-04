function a = j2Pert(r, J2, R, mu)

rnorm = norm(r);

a = 5*(r(3)/rnorm)^2;
a = [   r(1)/rnorm*(a-1)
		r(2)/rnorm*(a-1)
		r(3)/rnorm*(a-3)];
a = a*1.5*J2*mu*(R/rnorm^2)^2;
end