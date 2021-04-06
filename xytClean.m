function [xytOut] = xytClean(xytIn,frameWidth)
%XYTCLEAN Remove start and end point clusters and persistent objects.
%   
% RS, 03/2021

%% time
t = xytIn(:,3);

%% remove beginning and end 
tFirstZero = find(~ismember(1:max(t),t),1);
frameFirstZero = find(t>tFirstZero,1);

tLastZero = find(~ismember(1:max(t),t),1,'last');
frameLastZero = find(t<tLastZero,1,'last');

xytOut = xytIn(frameFirstZero:frameLastZero,:);

%% remove persistent objects
xytOut = removePersistentObjects(xytOut,frameWidth);



end
