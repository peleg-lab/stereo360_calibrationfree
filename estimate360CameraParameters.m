function stereo360Params = estimate360CameraParameters(matchedAlpha1,matchedAlpha2,estMethod)
%ESTIMATE360CAMERAPARAMETERS Estimates camera pose (t,R) from set of matched
%spherical projections alpha.
% Estimation methods estMethod are: 
%   - RANSAC, MSAC, LTS, Norem8Point, LMedS for Fundamental Matrix
%   - or minSearch for optimization search
% (see corresponding functions for details)
%
%   
% Raphael Sarfati, 03/2020
% raphael.sarfati@aya.yale.edu
% Peleg Lab, University of Colorado Boulder

% default method is fundamental matrix estimation using RANSAC
if nargin == 2
    estMethod = 'RANSAC';
end

if strcmp(estMethod,'RANSAC')
    stereo360Params = tRestimate_fundamentalMatrix(matchedAlpha1,matchedAlpha2,'RANSAC');
    
elseif strcmp(estMethod,'MSAC')
    stereo360Params = tRestimate_fundamentalMatrix(matchedAlpha1,matchedAlpha2,'MSAC');
    
elseif strcmp(estMethod,'LTS')
    stereo360Params = tRestimate_fundamentalMatrix(matchedAlpha1,matchedAlpha2,'LTS');
    
elseif strcmp(estMethod,'Norm8Point')
    stereo360Params = tRestimate_fundamentalMatrix(matchedAlpha1,matchedAlpha2,'Norm8Point');
    
elseif strcmp(estMethod,'LMedS')
    stereo360Params = tRestimate_fundamentalMatrix(matchedAlpha1,matchedAlpha2,'LMedS');
    
elseif strcmp(estMethod,'minSearch')
    stereo360Params = tRestimate_minSearch(matchedAlpha1,matchedAlpha2);
    
end

end

