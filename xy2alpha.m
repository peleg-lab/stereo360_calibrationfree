function [alpha] = xy2alpha(xy,v)
%XY2ALPHA Converts the xy (Nx2) coordinates in the equirectangular frame to
%   spherical projections alpha (Nx3)
%
%   xy is list of positions (Nx(2+p), subsequent p columns added at the end)
%   v informs about size of equirectangular frames; it is either: 
%       - a movie structure 
%       - an array containing [width height]
%
%   alpha is Nx(3+p) array (alpha projections are first 3 columns)
%   
%
% Raphael Sarfati, 03/2020
% Peleg Lab, University of Colorado Boulder

[theta,phi] = xy2thetaphi(xy,v);

alpha = horzcat(cos(theta).*sin(phi),sin(theta).*sin(phi),cos(phi),xy(:,3:end));

end

