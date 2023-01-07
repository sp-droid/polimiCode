function celestial_body (ID, x, y, z, r)
%
% The function plots in the 3D space the space object requested. 
% If x,y,z not specified, the default option is [0,0,0]. 
% If r not specified, the default option is r = 1.
% 
% INPUTS:
%       ID       Body ID
%       x       x-coordinate of body centre
%       y       y-coordinate of body centre
%       z       z-coordinate of body centre
%       r       radius scale (R*r)
% 
% Celestial objects ID:
%   [1]     Mercury
%   [2]     Venus
%   [3]     Earth
%   [4]     Mars
%   [5]     Jupiter
%   [6]     Saturn
%   [7]     Uranus
%   [8]     Neptune
%   [10]    Moon
%   [11]    Sun
% 
% VERSIONS:
%   2022-09-25
% 
% CONTRIBUTORS:
%   Alessandro Michelazzi

if nargin < 5 
    r = 1;
    if nargin < 4
        if nargin < 2
            x = 0; y = 0; z = 0;
        else
            error ('Not enough inputs!')
        end
    end
end

switch ID
    case 1
        object = imread('mercury.jpg');
        R = astroConstants(21) * r; % [km]
    case 2
        object = imread('venus.jpg');
        R = astroConstants(22) * r; % [km]
    case 3
        object = imread('earth.jpg');
        R = astroConstants(23) * r; % [km]
    case 4
        object = imread('mars.jpg');
        R = astroConstants(24) * r; % [km]
    case 5
        object = imread('jupiter.jpg');
        R = astroConstants(25) * r; % [km]
    case 6
        object = imread('saturn.jpg');
        R = astroConstants(26) * r; % [km]
    case 7
        object = imread('uranus.jpg');
        R = astroConstants(27) * r; % [km]
    case 8
        object = imread('neptune.jpg');
        R = astroConstants(28) * r; % [km]  
    case 10
        object = imread('moon.jpg');
        R = astroConstants(30) * r; % [km]
    case 11
        object = imread('sun.jpg');
        R = astroConstants(3) * r; % [km]
    otherwise
           error('Space object identifier %d is not defined!',ID);
end

[X,Y,Z] = sphere(100);
X = X*R + x; Y = Y*R + y; Z = Z*R - z;

surface(X,Y,-Z,'CData', object,'FaceColor','texturemap','EdgeColor','none')
axis equal
