function [xyz,err] = triangulate360(matchedAlpha1,matchedAlpha2,stereo360Params)
%TRIANGULATE360 Triangulates a set of matched pairs of spherical
%projections given the camera pose (t,R) contained in stereo360Params.
%
%   matchedAlpha1, matchedAlpha2 (Nx3) can be obtained from the matchPoints
%   function
%   stereo360Params can be obtained from estimate360CameraParameters
%
%   xyz (Nx3) is the corresponding list of triangulated positions
%   err is an error estimation
%
%   Note: we use R' instead of inv(R) for speed since they are equivalent.
%
%   NB: xyz coordinates are not in real world units, they are based on
%   camera separation = 1;
%   NB2: xyz coordinates may need to be permutated (and/or reversed) to
%   match intuitive coordinate system
%
%
% Raphael Sarfati, 03/2020
% raphael.sarfati@aya.yale.edu
% Peleg Lab, University of Colorado Boulder

nPoints = size(matchedAlpha1,1);
xyz = NaN(nPoints,3);
err = NaN(nPoints,1);

for i=1:nPoints
    
    alpha1 = matchedAlpha1(i,1:3);
    alpha2 = matchedAlpha2(i,1:3);
    beta2 = (stereo360Params.R'*alpha2')';
    
    % estimates optimal r1, r2 (see Ma 2015 for details)
    C = [alpha1' -beta2'];
    rr = lsqnonneg(C,stereo360Params.t);
    
    r1 = rr(1);
    r2 = rr(2);
    
    P1 = r1*alpha1;
    P2 = r2.*(stereo360Params.R'*alpha2')' + stereo360Params.t(:)';
    
    xyz(i,:) = (P1+P2)/2;
    err(i) = vecnorm(P1-P2);
    
    %w = waitbar(i/nPoints);
    
end

%close(w)

end

