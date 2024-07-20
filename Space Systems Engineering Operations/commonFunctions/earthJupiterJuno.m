function [Earth, Jupiter, Juno, time, distance, relVelocity] = earthJupiterJuno()
    AU = 149597870700;%m
    ephemeridsPath = strcat(fileparts(mfilename('fullpath')),"\ephemerids");
    addpath(ephemeridsPath)

    load('ephemerids/Earth.mat');
    load('ephemerids/Jupiter.mat');
    load('ephemerids/Juno.mat');
    
    time.year = Juno(:,1);%years
    Juno = Juno(:,2:7);
    distance.EarthToSun = vecnorm(Earth(:,1:3), 2, 2);%AU
    distance.JupiterToSun = vecnorm(Jupiter(:,1:3), 2, 2);
    distance.JunoToSun = vecnorm(Juno(:,1:3), 2, 2);
    distance.JunoToEarth = distanceObjects(Earth, Juno);
    distance.JunoToJupiter = distanceObjects(Jupiter, Juno);
    relVelocity.JunoEarth = relVelocityObjects(Earth, Juno)*AU/1000/24/3600;%km/s
    relVelocity.JunoJupiter = relVelocityObjects(Jupiter, Juno)*AU/1000/24/3600;
end