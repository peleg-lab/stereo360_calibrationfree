function stereo360Params = tRestimate_fundamentalMatrix(alpha1,alpha2,estMethod)
%TRESTIMATE_FUNDAMENTALMATRIX Estimates t and R from calculation of F.
%   Uses built-in function estimateFundamentalMatrix to estimate F.
%   Calculates t and R from F.
%   Different algorithms can be used to estimate F.
%
% Raphael Sarfati, 03/2020
% Peleg Lab, University of Colorado Boulder


% set the default estimation method to RANSAC
% other options are 'Norm8Point', 'LMedS', 'MSAC', 'LTS'
% see doc for more details
if nargin == 2
    estMethod = 'RANSAC';
end

% renormalize to get third coordinate equal to 1 (see Matlab doc)
points1 = alpha1(:,1:2)./alpha1(:,3);
points2 = alpha2(:,1:2)./alpha2(:,3);

% calculate fundamental matrix
% for consistency, points of camera 2 enter as first argument
F = estimateFundamentalMatrix(points2,points1,'Method',estMethod);

% estimate t and R from F
tR = F2tR(F);

stereo360Params.t = tR.t;
stereo360Params.R = tR.R;
stereo360Params.F = F;
stereo360Params.Method = estMethod;

end


