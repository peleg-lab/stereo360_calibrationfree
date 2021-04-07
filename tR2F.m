function [F] = tR2F(t,R)
%TR2F Calculates the fundamental matrix F from t and R
%   t is a vector 
%   R should be rotation matrix but can also be a vector of
%   the three rotation angles around the x, y, and z axes respectively.
%   NB: rotx, roty, rotz require Phased Array System Toolbox; if unable to
%   get, uncomment function at the end of this file.
%
% RS, 10/2019

Tx = vcross(t);

if ismatrix(R)
    F = Tx*R';
    
elseif isvector(R)
    phixyz = R;
    R = rotx(phixyz(1))*roty(phixyz(2))*rotz(phixyz(3));
    F = Tx*R';
    
else
    error('R must be either a rotation matrix or a vector of three angles.');
    
end

end

function [Vx] = vcross(v)
%VCROSS Calculate matrix Vx such that Vx*u = cross(v,u).

Vx = [0 -v(3) v(2) ; v(3) 0 -v(1); -v(2) v(1) 0]; 

end

%% Manual definition of rotx, roty, rotz if unable to get Phased Array System Toolbox 
% Uncomment if necessary.
% See https://www.mathworks.com/help/phased/ref/rotx.html for details.
% 
function Rx = rotx(phi)

Rx = [1 0 0 ; 0 cos(phi) -sin(phi) ; 0 sin(phi) cos(phi)];

end

function Ry = roty(phi)

Ry = [cos(phi) 0 sin(phi) ; 0 1 0 ; -sin(phi) 0 cos(phi)];

end

function Rz = rotz(phi)

Rz = [cos(phi) -sin(phi) 0 ; sin(phi) cos(phi) 0 ; 0 0 1];

end
