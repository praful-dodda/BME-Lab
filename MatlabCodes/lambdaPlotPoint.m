function [plotVals] = lambdaPlotPoint(lambda1,lambda2,decileAvgModelVal,pairedObs)
% extrapolate average slope over all deciles to find the first and last
% lambda1/lambda2 values; determine 12 plot points for curve

%find average slope for extrapolation of lambda mins/max
lambda1slope=(lambda1(:,10)-lambda1(:,1))/(decileAvgModelVal(:,10)-decileAvgModelVal(:,1));
lambda2slope=(lambda2(:,10)-lambda2(:,1))/(decileAvgModelVal(:,10)-decileAvgModelVal(:,1));

lambda1intercept=lambda1(:,1)-decileAvgModelVal(:,1)*lambda1slope;
obsValMin=lambda1slope*min(pairedObs.modelVal)+lambda1intercept;
obsValMax=lambda1slope*max(pairedObs.modelVal)+lambda1intercept;
lambda2intercept=lambda2(:,1)-decileAvgModelVal(:,1)*lambda2slope;
obsVarMin=lambda2slope*min(pairedObs.modelVal)+lambda2intercept;
obsVarMax=lambda2slope*max(pairedObs.modelVal)+lambda2intercept;

% plot points
plotVals(1,:)=[obsValMin min(pairedObs.modelVal) obsVarMin];
plotVals(12,:)=[obsValMax max(pairedObs.modelVal) obsVarMax];
for i=1:10
    plotVals(i+1,:)=[lambda1(i) decileAvgModelVal(i) lambda2(i)];
end
end