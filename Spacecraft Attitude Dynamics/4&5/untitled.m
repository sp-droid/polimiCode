I = diag([1,1,1]);
X0 = I;

options=optimset('Display','iter');   % Option to display output
[x,fval] = fsolve(@funcc,X0,options)

function f0 = funcc(X)
I = diag([1,1,1]);
Inertia = [0.07;0.005;0.0011];
f0 = [norm(X'*X-I);
      norm(X*Inertia-[1;0;0])];
end