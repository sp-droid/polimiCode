
function Cd0_tot=zerolift_drag_coast(Mach)

% Parameters
q=1.304737*10^5;


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

% Nozzle exit area A VERIFIER
Ae=0.27854;

% Dimensions of the strakes and the fins  %(1)=strakes (2)=fins
g=0.55*[1.9 0.96]; % base dimension [m]
p=0.55*[1.85 0.75];  % tip dimension 
h=0.55*[0.08 0.8];   % Height dimension A VERIFIER

% Number of strakes and fins
nb_strakes=4; % ou 2?
nb_fins=4;

% Surfaces of the strakes
Sur=h(1)*(p(1)+g(1))/2;

% Surface of the fins
Sur_2=h(2)*(p(2)+g(2))/2;


nose_type='C';

% Parameters for the wave "wings" contribution ------------------

 % for the strakes
t=p(1)/g(1);    %taper ratio
cmax=g(1)*2/3*((1+t+t^2)/(1+t));
ALE=deg2rad(45);
Mle=Mach*cos(ALE);
deltaLE=deg2rad(10); % A VOIR
tmac=0.55*0.02;
bwing=h(1); % A VERIFIER

  % for the fins
t2=p(2)/g(2);    %taper ratio
cmax2=g(2)*2/3*((1+t2+t2^2)/(1+t2));
ALE2=deg2rad(38.4);
Mle2=Mach*cos(ALE2);
tmac2=0.55*0.003;
bwing2=h(2); % A VERIFIER

    if Mach>1
        Cd0_body_wave=(1.586+1.834/(Mach^2))*(atan(0.5/(ln/d)))^(1.69);
        Cd0_body_base=0.25/Mach;        
    else 
        Cd0_body_wave=0;
        Cd0_body_base=0.12+0.13*Mach^2;
    end

Cd0_body_friction=0.053*(l^0.8/d)*(Mach.^-0.2)/(q)^(0.2); % ATTENTION -0.2 POUR AVOIR LA BONNE COURBE....   

% Total contribution of the body 
Cd0_body=Cd0_body_wave+Cd0_body_base+Cd0_body_friction;


% Contribution of the strakes
 if Mle>1
     Cd0_strakes_wave=nb_strakes*(1.429/(Mle^2))*((1.2*Mle^2)^(3.5)*(2.4/(2.8*Mle^2-0.4))^(2.5)-1)*((sin(deltaLE))^2*cos(ALE)*3.28*tmac*bwing/AREF);
 else
     Cd0_strakes_wave=0;
 end

    Cd0_strakes_friction=nb_strakes*0.01333*(Mach/(q*cmax))^(-0.2)*2*(Sur/AREF);

% Contribution of the fins

if Mle>1
     Cd0_fins_wave=nb_fins*(1.429/(Mle2^2))*((1.2*Mle2^2)^(3.5)*(2.4/(2.8*Mle2^2-0.4))^(2.5)-1)*((sin(deltaLE))^2*cos(ALE2)*3.28*tmac2*bwing2/AREF);
 else
     Cd0_fins_wave=0;
 end

    Cd0_fins_friction=nb_fins*0.01333*(Mach/(q*cmax2))^(-0.2)*2*(Sur_2/AREF);

% Total contribution of the wings (strakes + fins)

    Cd0_wings=Cd0_strakes_wave+Cd0_strakes_friction+Cd0_fins_wave+ Cd0_fins_friction;

% TOTAL
    
    Cd0_tot=Cd0_body+Cd0_wings;
end

