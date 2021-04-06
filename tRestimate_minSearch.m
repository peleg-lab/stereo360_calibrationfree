function [out] = tRestimate_minSearch(alpha1,alpha2,x0)
%TRESTIMATE_MINSEARCH Calculates t and R by optimization. 
%   Minimizes sum of distances given t and R, represented as a 6D vector
%   x = {tx,ty,tz,phix,phiy,phiz}, with phi the 3 rotation angles.
%   alpha1, alpha2 are matching sets of projected points
%   x0 (optional) is the initial guess vector
%
%   Uses fmincon, since we have the constraint |t| = 1.
%   The cost function is defined below.
%
% Raphael Sarfati, 03/2020
% Peleg Lab, University of Colorado Boulder

% if guess vector not provided, assumes t = tx and R = Id
if nargin == 2
    x0 = [1 0 0 0 0 0]';
end

% arguments for fmincon
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = @fcon;
options = optimset('MaxFunEvals',1e6,'MaxIter',1e6,'Display','off');

% optimization, applying fmincon to the cost function defined separately
% and with the constraint |t| = 1
xmin = fmincon(@(x) tRcostFunction(alpha1,alpha2,x),x0,...
    A,b,Aeq,beq,lb,ub,nonlcon,options);


out.t = xmin(1:3);
out.phi = xmin(4:6);
out.R = rotx(xmin(4))*roty(xmin(5))*rotz(xmin(6));
out.F = tR2F(out.t,out.R);

end

function [c,ceq] = fcon(x)
%FCON Constraint for fmincon.
%   t(1)^2 + t(2)^2 + t(3)^2 = 1

c = [];
ceq = x(1)^2 + x(2)^2 +x(3)^2 - 1;

end


%% Cost Function

function [c] = tRcostFunction(alpha1,alpha2,x)
%TRCOSTFUNCTION Cost (sum of squared distances) given t and R.
%   Calculates sum of squared distances between points in camera 1 and
%   camera 2 given the 6 parameters in x (translation and rotation
%   coordinates).
%   Given the projected points alpha1, alpha2, and a test t, R, there is an
%   optimum set (r1,r2) to mimimize || r1*a1 - (t + R'*r2*a2) ||, 
%   so (r1,r2) don't need to be put into the test vector x.
%
% RS, 10/2019

% number of points
N = size(alpha1,1);

% t and R from test vector x
t = [x(1) x(2) x(3)];
R = rotx(x(4))*roty(x(5))*rotz(x(6));

% calculate optimal r1 and r2
r1 = zeros(N,1);
r2 = zeros(N,1);

for i=1:N
    
    rr = rOptimum(alpha1(i,:),alpha2(i,:),t,R);
    r1(i) = rr(1);
    r2(i) = rr(2);
    
end

% separations between X1 and t+R'*X2 
s = r1.*alpha1 - r2.*(R'*alpha2')' - t;

% sum of distances
c = sum(vecnorm(s,2,2),1);

end


function [rr] = rOptimum(a1,a2,t,R)
%ROPTIMUM Calculates r1 and r2 that minimize N = || r1*a1 - t - R'*r2*a2 ||
% Matrix method based on Ma 2015, but equivalent to setting dN/dr1 = 0 and
% dN/dr2 = 0.

a1 = a1(:);
a2 = a2(:);
t = t(:);

C = [a1 , -R'*a2];

%rr = (C'*C)\C'*t;

A = -eye(2);
b = zeros(2,1);
Aeq = [];
beq = [];
lb = [];
ub = [];
x0 = [];
options = optimoptions(@lsqlin,'Display','off');
rr = lsqlin(C,t,A,b,Aeq,beq,lb,ub,x0,options);

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

