function [valuesCorr] = monotonicCorrTEST(values)
% ensure Lambda1 and Lambda2 for CAMP correction are monotonically
% decreasing within each decile
%*NOTE: 'values' must be inputted in the order they are meant to be plotted
%(first to last)
%
% SYNTAX:
% function monotonicCorrTEST(values)
%
% INPUT:
% obs (scalar) 1xn An array of values of points on a line
% 
% OUTPUT:
% valuesCorr (scalar) 1xn An array of inputted values corrected to always
% be monotonically increasing from the middle point (*floor if n is even) to the last and
% monotonically decreasing from the middle point (*floor if n is even)
% testing CAMP monotonic forcing function
%
% EXAMPLES:
%with random integers
% r = randi(10,10,1);
% [rCorr] = monotonicCorrTEST(r);
% rCenterline=floor(length(r)/2);
% figure;
% hold on;
% plot(r);
% plot(rCorr);
% xline(rCenterline,'k:');
% legend('Random Values','Monotonic Correction','Point of Monotonic Forcing');
% 
% %with a sine curve
% bInc = rand(1);
% b=sin(-pi:bInc:pi);
% [bCorr] = monotonicCorrTEST(b);
% bCenterline=floor(length(b)/2);
% figure;
% hold on;
% plot(b);
% plot(bCorr);
% xline(bCenterline,'k:');
% legend('Random Values','Monotonic Correction','Point of Monotonic Forcing');
% 
% %with a polynomial
% p = randi([-5,5],10,1)';
% x = randi(10,10,1)';
% xSort = sort(x);
% y=polyval(p,xSort);
% [yCorr] = monotonicCorrTEST(y);
% yCenterline=floor(length(y)/2);
% figure;
% hold on;
% plot(y);
% plot(yCorr);
% xline(yCenterline,'k:');
% legend('Random Values','Monotonic Correction','Point of Monotonic Forcing');

valuesCorr = [];

valuesCorr(floor(length(values)/2)) = values(floor(length(values)/2));

for i = floor(length(values)/2)+1:length(values)
    %make monotonically increasing - fix for loop
    if values(i) < valuesCorr(i-1) %if last decile is less than penultimate, replace last with penultimate
        valuesCorr(i) = valuesCorr(i-1);
    else
        valuesCorr(i) = values(i);
    end
end

for i = floor(length(values)/2)-1:-1:1
    %make monotonically decreasing - fix for loop
    if values(i) > valuesCorr(i+1) %if decile is less than previous, replace decile with next decile value
        valuesCorr(i) = valuesCorr(i+1);
    else
        valuesCorr(i) = values(i);
    end
end
