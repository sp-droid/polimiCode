function Cd0_capsule=zerolift_drag_capsule(Mach)
% Parameters
q=1.304737*10^5;
% ref at 20000 ft

% Shape of the launcher
x=[0.308754 2.174555]; %[m]
l=x(2);
d=x(1);

% Shape of the nose %[m]
a=[x(1) 0.922519];
b=a;
d_nose=a(1);
ln=a(2);

% Reference surface (cross section of the launcher)
AREF=pi*(x(2)/2)^2; % has to be the same for all
Sur_ref=AREF;


    if Mach>1
        Cd0_capsule_wave=(1.586+1.834/(Mach^2))*(atan(0.5/(ln/d)))^(1.69);
        Cd0_capsule_base=0.25/Mach;        
    else 
        Cd0_capsule_wave=0;
        Cd0_capsule_base=0.12+0.13*Mach^2;
    end

Cd0_capsule_friction=0.053*(ln^0.8/d)*(Mach.^-0.2)/(q)^(0.2); %ln = l nose= l capsule, à modifier si besoin

% Total contribution of the body 
Cd0_capsule=Cd0_capsule_wave+Cd0_capsule_base+Cd0_capsule_friction;