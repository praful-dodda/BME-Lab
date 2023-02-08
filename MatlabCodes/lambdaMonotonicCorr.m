function [lambdaCorr] = lambdaMonotonicCorr(lambda)
% ensure Lambda1 and Lambda2 for CAMP correction are monotonically
% decreasing within each decile

lambdaCorr = [];

% monotonically decreasing from last decile
% for i = length(lambda):-1:1
%     if i == length(lambda)
%         lambdaCorr(i) = lambda(i);
%     elseif i == 1
%             if lambda(i) > lambdaCorr(i+1) %if last decile is less than penultimate, replace last with penultimate
%                 lambdaCorr(i) = lambdaCorr(i+1);
%             else 
%                 lambdaCorr(i) = lambda(i);
%             end
%     else
%         if lambda(i) >= lambda(i-1) & lambda(i) <= lambdaCorr(i+1) %if decile is less than previous, replace decile with next decile value
%             lambdaCorr(i) = lambda(i);
%         else
%             lambdaCorr(i) = lambdaCorr(i+1);
%         end
%     end
% end

lambdaCorr(floor(length(lambda)/2)) = lambda(floor(length(lambda)/2));

for i = floor(length(lambda)/2)+1:length(lambda)
    %make monotonically increasing - fix for loop
    if lambda(i) < lambdaCorr(i-1) %if last decile is less than penultimate, replace last with penultimate
        lambdaCorr(i) = lambdaCorr(i-1);
    else
        lambdaCorr(i) = lambda(i);
    end
end

for i = floor(length(lambda)/2)-1:-1:1
    %make monotonically decreasing - fix for loop
    if lambda(i) > lambdaCorr(i+1) %if decile is less than previous, replace decile with next decile value
        lambdaCorr(i) = lambdaCorr(i+1);
    else
        lambdaCorr(i) = lambda(i);
    end
end
