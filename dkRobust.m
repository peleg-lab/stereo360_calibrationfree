function [dk,res] = dkRobust(N1,N2)
%DKROBUST Calculates dk from raw N (number of objects) time series
%   Removes outliers due to artificial light pollution,
%   then cross-correlates cleaned traces over random intervals.
%
% Raphael Sarfati, 02/2021
% raphael.sarfati@aya.yale.edu
% Peleg Lab, University of Colorado Boulder

%% removes trailing zeros at the beginning
startFrame1 = find(N1>0,1);
startFrame2 = find(N2>0,1);
startFrame = min(startFrame1,startFrame2);

N1c = N1(startFrame:end);
N2c = N2(startFrame:end);

res.N1c = N1c;
res.N2c = N2c;
res.startFrame = startFrame;

%% removes potential punctual noise
N1cc = removeInnerNoise(N1c);
N2cc = removeInnerNoise(N2c);
res.N1cc = N1cc;
res.N2cc = N2cc;

%% cross-correlation of random N intervals
[dk,out] = ransacNxcorr(N1cc,N2cc);
res.numTrials = out.numTrials;
res.dks = out.dks;
res.gof = out.gof;

end

function [Ncc] = removeInnerNoise(Nc)
%REMOVEINNERNOISE Removes sporadic outliers
%   Fits strictly positive N values to exponential distribution;
%   Sets outliers to 0.

Ncc = Nc;
p = Nc(Nc>0 & Nc<10);
exppdf = fitdist(p(:),'exponential');

outliers = cdf(exppdf,Ncc) > 1-1/length(Ncc);
Ncc(outliers) = 0;

end

function [dk,out] = ransacNxcorr(N1cc,N2cc)
%RANSACNXCORR Cross-correlation of cleaned traces over random intervals.

numTrials = 100;
maxIntervalSize = max(length(N1cc),length(N2cc));

for i=1:numTrials
    windowSize = randi(round([maxIntervalSize/2 maxIntervalSize]));
    k = randi(maxIntervalSize-windowSize);
    s1 = circshift(N1cc,k);
    s2 = circshift(N2cc,k);
    dks(i) = finddelay(s1,s2);    
end

[dk,f] = mode(dks);

out.numTrials = numTrials;
out.dks = dks;
out.gof = f/numTrials;

if out.gof < 0.8
    warning(['dk estimation seems unreliable; gof =' num2str(out.gof)])
end

end



