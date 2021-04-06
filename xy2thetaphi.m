function [theta,phi] = xy2thetaphi(xy,v)
%XY2THETAPHI Converts the xy coordinates in the equirectangular frame to
%   polar angle theta and azimuthal angle phi
%
%   xy is a list of positions (Nx2, subsequent colums are not considered)
%   v informs about size of equirectangular frames; it is either: 
%       - a movie structure 
%       - an array containing [width height]
%   
%
% Raphael Sarfati, 03/2020
% Peleg Lab, University of Colorado Boulder

if nargin == 1
    error('You need to specify the size of equirectangular frames.')
else
    if isobject(v)
        w = v.Width;
        h = v.Height;
    elseif isnumeric(v)
        w = v(1);
        h = v(2);
        if w~=2*h
            error('Frames dimensions are not consistent; width must be twice height.')
        end
    end
end

theta = xy(:,1)*(2*pi)/w;
phi = xy(:,2)*pi/h;

end

