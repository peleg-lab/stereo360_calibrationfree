function [trajCal] = extractCalibrationTrajectories(xyt1,xyt2,dk)
%EXTRACTCALIBRATIONTRAJECTORIES Extracts calibration trajectories by
% finding synchronized frames with only one flash detected.
% dk returned by dkRobust
%   
% RS, 7/2020

rp = regionprops(xyt1(:,3));
a1 = vertcat(rp.Area);

rp = regionprops(xyt2(:,3)-dk);
a2 = vertcat(rp.Area);

f1 = find(a1==1);
f2 = find(a2==1);

f = intersect(f1,f2);

idx1 = ismember(xyt1(:,3),f);
idx2 = ismember(xyt2(:,3)-dk,f);

trajCal.j1 = xyt1(idx1,:);
trajCal.j2 = xyt2(idx2,:);

end

