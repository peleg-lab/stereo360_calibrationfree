function [out] = F2tR(F,tguess,Rguess)
%F2TR Returns best estimate for t and R, given F (= Tx*R')
%   F is the fundamental matrix, given by estimateFundamentalMatrix
%   tguess, Rguess are estimates useful to resolve ambiguities on t, R
%   Method is based on single-value decomposition (SVD) and W, Z matrices.
%   See link: http://www.maths.lth.se/matematiklth/personal/calle/datorseende13/notes/forelas6.pdf
%
% RS, 10/2019

% if no estimate for t, R, assume t = tx and R = Id.
if nargin == 1
    tguess = [1 0 0];
    Rguess = eye(3);
elseif nargin == 2
    Rguess = eye(3);
end

% make tguess vertical
tguess = tguess(:);

% Decomposition using SVD and W, Z as explained in link above.
[U,~,V] = svd(F);

W = [0 -1  0 ; 
     1  0  0 ;
     0  0  1]; 
 
Z = [0 1 0; 
    -1 0 0; 
     0 0 0];
 
S1 = -U*Z*U';
S2 = U*Z*U';

R1 = U*W'*V';
R2 = U*W*V';

e1 = null(S1);
out.e1 = e1;

% corrects to get real rotations
if det(R1)<0
    out.R1 = -R1';
else
    out.R1 = R1';
end

if det(R2)<0
    out.R2 = -R2';
else
    out.R2 = R2';
end

if norm(out.R1-Rguess) < norm(out.R2-Rguess)
    out.R = out.R1;
else
    out.R = out.R2;
end

if norm(e1-tguess) < norm(-e1-tguess)
    out.t = e1;
else
    out.t = -e1;
end


end

