function [xyzt,swo] = reconstruct360v02(xyt1,xyt2,frameWidth,scaleFactor,stereo360Params)
%RECONSTRUCT360 Wrapper for calibration-free 3D reconstruction from 360-degree cameras.
%   No calibration (spatial or temporal) needed:
%       - uses temporal cross-correlation to estimate frame delay;
%       - uses unique-flash frames for stereo calibration.
%
%   INPUTS
%   xyt1 -- (N1x3 array) flash positions in equirectangular frames of
%       camera1; t is time in frame number (positive integers)
%   xyt2 -- (N2x3 array) same for camera 2
%   frameWidth -- number of pixels for the width of equirectangular frames
%   scaleFactor (optional) -- actual distance between the two cameras
%   stereo360Params (optional) -- if stereo360Params already calculated; if
%                               empty, only does xytCleaning and dk
%
%   OUTPUTS
%   xyzt -- 3D reconstructed positions; xy is the horizontal plane, z is
%       vertical up direction; t is time in frame number
%   swo -- swarm 360 structure collecting all intermediate variables, if more
%       details are needed
%
%   reconstruct360 v2
%
%   Raphael Sarfati, 03/2021
%   raphael.sarfati@aya.yale.edu
%   Peleg Lab, University of Colorado Boulder


%% initialize
swo.xyt1in = xyt1;
swo.xyt2in = xyt2;
swo.frameWidth = frameWidth;

if nargin > 3
    swo.scaleFactor = scaleFactor;
else
    scaleFactor = 1;
    warning('no scale factor specified; space dimensions will be unitless')
end

% saves processing date and code
swo.dateProcessed = datestr(now,31);
swo.code = fileread([mfilename('fullpath') '.m']);


%% clean xyt
xyt1 = xytClean(xyt1,frameWidth);
xyt2 = xytClean(xyt2,frameWidth);

swo.xyt1 = xyt1;
swo.xyt2 = xyt2;


%% calculate dk
nObjCam1 = Ntimeseries(xyt1);
nObjCam2 = Ntimeseries(xyt2);

[dk,dkRes] = dkRobust(nObjCam1,nObjCam2);

disp(['frame delay between cameras: dk = ' num2str(dk)])
disp(datetime('now'))

swo.dk = dk;
swo.dkRes = dkRes;


%% calculate stereo360Params
if nargin == 5
    if isempty(stereo360Params)
        xyzt = [];
        disp('no spatial calibration')
        return
    else
        disp('using stereo360Params input')
        disp(datetime('now'))
    end
else

[calTraj1,calTraj2] = extractCalibrationTrajectories(xyt1,xyt2,dk);

swo.calTraj1 = calTraj1;
swo.calTraj2 = calTraj2;

frameHeight = frameWidth/2;
calAlpha1 = xy2alpha(calTraj1,[frameWidth frameHeight]);
calAlpha2 = xy2alpha(calTraj2,[frameWidth frameHeight]);

disp('starting calibration; this may take a while (up to ~1hr)...')
disp(datetime('now'))

stereo360ParamsRANSAC = estimate360CameraParameters(calAlpha1(:,1:3),calAlpha2(:,1:3),'RANSAC');
stereo360Params = estimate360CameraParameters(calAlpha1(:,1:3),calAlpha2(:,1:3),'minSearch');

% save to workspace -- just in case
assignin('base','stereo360Params',stereo360Params)

swo.stereo360ParamsRANSAC = stereo360ParamsRANSAC;
swo.stereo360Params = stereo360Params;

disp('calibration complete, yay!')
disp(datetime('now'))

end


%% match and triangulate
disp('starting match and triangulate...')
disp(datetime('now'))
xyzt = mt360(xyt1,xyt2,frameWidth,dk,stereo360Params);

xyzt(:,1:3) = xyzt(:,1:3)*scaleFactor;

swo.xyzt = xyzt;
successRate = size(xyzt,1)/min(size(xyt1,1),size(xyt2,1));
swo.successRate = successRate;

disp(['3D reconstruction completed; triangulation success: ' num2str(successRate,2)])
disp(datetime('now'))

end

function n = Ntimeseries(xyt)

t = xyt(:,3);
tmin = 1;
tmax = max(t);
binEdges = (tmin-0.5):(tmax+0.5);
n = histcounts(t,binEdges);

end

function [calTraj1,calTraj2] = extractCalibrationTrajectories(xyt1,xyt2,dk)
%EXTRACTCALIBRATIONTRAJECTORIES Extracts calibration trajectories by
% finding synchronized frames with only one flash detected.

rp = regionprops(xyt1(:,3));
a1 = vertcat(rp.Area);

rp = regionprops(xyt2(:,3)-dk);
a2 = vertcat(rp.Area);

f1 = find(a1==1);
f2 = find(a2==1);

f = intersect(f1,f2);

idx1 = ismember(xyt1(:,3),f);
idx2 = ismember(xyt2(:,3)-dk,f);

calTraj1 = xyt1(idx1,:);
calTraj2 = xyt2(idx2,:);

% limit to first 5000 points to avoid too long calibration
if size(calTraj1,1)>5000
    calTraj1 = calTraj1(1:5000,:);
    calTraj2 = calTraj2(1:5000,:);
end

end


function xyzt = mt360(xyt1,xyt2,frameWidth,dk,stereo360Params)
%MT360 Match and Triangulate

xyt2(:,3) = xyt2(:,3)-dk;

dThreshold = 0.2;

% match
[matchedAlpha1t,matchedAlpha2t] = matchPoints(xyt1,xyt2,[frameWidth frameWidth/2],stereo360Params,dThreshold);

% triangulate
xyz = triangulate360(matchedAlpha1t,matchedAlpha2t,stereo360Params);
xyzt = [xyz matchedAlpha1t(:,4)];


end


