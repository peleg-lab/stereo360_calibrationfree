function [matchedAlpha1t,matchedAlpha2t] = matchPoints(xyt1,xyt2,v,stereo360Params,dThreshold)
%MATCHPOINTS Matches pairs of x,y coordinates in equirectangular frames given t and R.
%
%   xyt1 (N1x3) and xyt2 (N2x3) contain list of xy coordinates and
%   corresponding frame number t; they need not have same size
%   v informs about frame size, see xy2alpha
%   stereo360Params contains t and R
%   dThreshold is cost threshold, see matchAlphaPoints 
%   
%
% Raphael Sarfati, 03/2020
% Peleg Lab, University of Colorado Boulder

matchedAlpha1t = [];
matchedAlpha2t = [];

for t = unique(xyt1(:,3))'
    alpha1 = xy2alpha(xyt1(xyt1(:,3)==t,1:2),v);
    alpha2 = xy2alpha(xyt2(xyt2(:,3)==t,1:2),v);
    
    if ~isempty(alpha1) && ~isempty(alpha2)
        
        [ma1,ma2] = matchAlphaPoints(alpha1,alpha2,stereo360Params,dThreshold);
        matchedAlpha1t = vertcat(matchedAlpha1t,[ma1 repmat(t,size(ma1,1),1)]);
        matchedAlpha2t = vertcat(matchedAlpha2t,[ma2 repmat(t,size(ma1,1),1)]);
        
    end
end



end

