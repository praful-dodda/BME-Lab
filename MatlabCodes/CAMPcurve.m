function [lambda1,lambda2,decileModelVal,idxDecile,modelVal,obsVal,modelLoc,tME]=CAMPcurve(pairedObs,figName,lt,plotCAMPcurve,modelVar,backLog);

if nargin<4, plotCAMPcurve=1; end;
if nargin<5, modelVar=NaN; end;
if nargin<6, backLog=1; end;   % 0:mean, 1:median

%sort into deciles
idxNaN=find(~isnan(pairedObs.obsVal)&~isnan(pairedObs.modelVal));
modelVal=pairedObs.modelVal(idxNaN);
obsVal=pairedObs.obsVal(idxNaN);
obsVar=pairedObs.obsVar(idxNaN);
modelLoc=pairedObs.modelLoc((idxNaN),:);
tME(1,:)=pairedObs.tME(idxNaN);
[~,idxSort]=sort(modelVal,'ascend');
decileSize=length(idxSort)*0.1;
for i=1:10
    idxDecile(:,i)=idxSort((i-1)*decileSize+1:decileSize*i);
    decileModelVal(:,i)=modelVal(idxDecile(:,i));
    decileObsVal(:,i)=obsVal(idxDecile(:,i));
    decileObsVar(:,i)=obsVar(idxDecile(:,i));
end

decileAvgModelVal=mean(decileModelVal,"omitnan"); %mean model value
decileAvgObsVal=mean(decileObsVal,"omitnan"); %mean obs value
decileVarModel=var(decileModelVal,"omitnan"); %var model value
decileStdModel=sqrt(decileVarModel); %std model value
decileVarObs=mean(decileObsVar,"omitnan");  %mean obs value variance    %change this to something else besides mean if using
decileVar=var(decileObsVal,"omitnan"); %var obs value
decileStd=sqrt(decileVar); %std obs val

% Lambda 1 and 2 unadjusted plot points
lambda1Raw=decileAvgObsVal;
lambda2=decileStd;

[plotValsUnadj] = lambdaPlotPoint(lambda1Raw,lambda2,decileAvgModelVal,pairedObs);

% monotonically decreasing lambda1 correction and plot points
[lambda1] = lambdaMonotonicCorr(lambda1Raw);
[plotVals] = lambdaPlotPoint(lambda1,lambda2,decileAvgModelVal,pairedObs);

% monotonically decreasing lambda1 and lambda2 correction and plot points
[lambda2Adj] = lambdaMonotonicCorr(lambda2);
[plotValsAdj] = lambdaPlotPoint(lambda1,lambda2Adj,decileAvgModelVal,pairedObs);

%add coverage function

%plot deciles/ramp curve
if plotCAMPcurve>0

    %if starting with log-t values
    if lt == 1

        %log-t figure
        figure;
        hold on;

        yyaxis left
        % 1:1 line
        plot(-10:3,-10:3)
        % paired values
        plot(pairedObs.modelVal,pairedObs.obsVal,'.','Color', [0.4660 0.6740 0.1880]);
        % monotonically adjusted lambda1
        plot(plotVals(:,2),plotVals(:,1),'-*','LineWidth',2);
%        % Unadjusted lambda1 plot values
%        plot(plotValsUnadj(:,2),plotValsUnadj(:,1),'k--.','LineWidth',1);
        ylabel('Observed Values (ln-ppb)')

        yyaxis right
        % lambda2
        plot(plotVals(:,2),plotVals(:,3).^2,'-*','LineWidth',2);
        % monotonically adjusted lambda2 plot values
        plot(plotValsAdj(:,2),plotValsAdj(:,3).^2,'k--.','LineWidth',1);
        ylabel('\lambda_2 Variance (ln-ppb)^2')

        % plot decile lines
        xline(decileModelVal(1,:));

        xlim([min(modelVal) max(modelVal)]);
        xlabel('Model Values (ln-ppb)');
        legend('1:1 Line','Paired Observed-Model Values','Monotonic \lambda_1','Unadjusted \lambda_2','Monotonic \lambda_2');
        title('Model Value vs. Observed Value and Std. at each Decile (ln-ppb)');
        figName1=sprintf('%s_CAMP_logt',figName);
        figName1=(['..\analysis\8CAMx\figs\' figName1]);
        print(figName1, '-dpng');

        %non log-t figure
        modPlotVals=pairedObs.modelVal;
        obsPlotVals=pairedObs.obsVal;
        CAMPmodelPlotVals=plotVals(:,2);
        lambda1PlotVals=plotVals(:,1);
        lambda2PlotVals=plotVals(:,3);
        decilePlotVals=decileModelVal(1,:);
        %backlog to mean >> ASK MARC!!
        if backLog==0
            disp("Error.")
%             %model values
%             modPlotVals=exp(modelValsPlot+(modelVar/2));
%             %observed values
%             obsPlotVals=exp(modelValsPlot+(modelVar/2));
%             %CAMP curve model values
%             CAMPmod_mZv=exp(modelValsPlot+(modelVar/2));
%             %lambda1
%             lambda1_mZv=exp(modelVals+(modelVar/2));
%             %lambda2
%             lambda2_mZv=exp(modelVals+(modelVar/2));
%             %decile model values
%             decile_mZv=exp(modelVals+(modelVar/2));
        %backlog to median
        elseif backLog==1
            %model values
            mod_mZv=exp(modPlotVals);
            %obs values
            obs_mZv=exp(obsPlotVals);
            %CAMP curve model values
            CAMPmod_mZv=exp(CAMPmodelPlotVals);
            %lambda1
            lambda1_mZv=exp(lambda1PlotVals);
            %lambda2
            lambda2_mZv=exp(lambda2PlotVals);
            %decile plot values
            decile_mZv=exp(decilePlotVals);
        end
        figure;
        hold on;
        plot(mod_mZv,obs_mZv,'.');
        % monotonically adjusted lambda1 and unadjusted lambda2 plot values
        plot(CAMPmod_mZv,lambda1_mZv,'-*','LineWidth',2);
        plot(CAMPmod_mZv,lambda2_mZv,'-*','LineWidth',2);
        % plot decile lines
        xline(decile_mZv);
        % 1:1 line
        plot(-10:3,-10:3)
        xlim([min(mod_mZv) max(mod_mZv)]);
        xlabel('Model Values (ppb)');
        ylabel('Obs Values (ppb)');
        legend('Model Value','Monotonic \lambda_1','Unadjusted \lambda_2^1^/^2','1:1 Line');
        title('Model Value vs. Observed Value and Std. at each Decile (ppb)');
        figName2=sprintf('%s_CAMP',figName);
        figName2=(['..\analysis\8CAMx\figs\' figName2]);
        print(figName2, '-dpng');

    elseif lt == 0
        %non-logt figure
        figure;
        hold on;
        plot(pairedObs.modelVal,pairedObs.obsVal,'.');
        % monotonically adjusted lambda1 and unadjusted lambda2 plot values
        plot(plotVals(:,2),plotVals(:,1),'-*','LineWidth',2);
        plot(plotVals(:,2),plotVals(:,3),'-*','LineWidth',2); 
        % Unadjusted lambda1 plot values
        plot(plotValsUnadj(:,2),plotValsUnadj(:,1),'k--.','LineWidth',1);
        % monotonically adjusted lambda2 plot values
        plot(plotValsAdj(:,2),plotValsAdj(:,3),'k--.','LineWidth',1);
        % plot decile lines
        xline(decileModelVal(1,:));
        % 1:1 line
        plot(-10:3,-10:3)
        xlim([min(modelVal) max(modelVal)]);
        xlabel('Model Values (ppb)');
        ylabel('Obs Values (ppb)');
        legend('Model Value','Monotonic \lambda_1','Unadjusted lambda_2^1^/^2','Unadjusted lambda_1','Monotonic lambda_2^1^/^2','1:1 Line');
        title('Model Value vs. Observed Value and Std. at each Decile (ppb)');
        figName1=sprintf('%s_CAMP',figName);
        figName1=(['..\analysis\8CAMx\figs\' figName1]);
        print(figName1, '-dpng');

        %logt figure
        %model values
        modPlotVals=pairedObs.modelVal;
        logModPlotVals=log(modPlotVals);
        %obs values
        obsPlotVals=pairedObs.obsVal;
        logObsPlotVals=log(obsPlotVals);
        %CAMP curve model values
        CAMPmodelPlotVals=plotVals(:,2);
        logCAMPmodelPlotVals=log(CAMPmodelPlotVals);
        %lambda1 values
        lambda1PlotVals=plotVals(:,1);
        logLambda1PlotVals=log(lambda1PlotVals);
        %lambda2 values
        lambda2PlotVals=plotVals(:,3);
        logLambda2PlotVals=log(lambda2PlotVals);
        %decile plot values
        decilePlotVals=decileModelVal(1,:);
        logDecilePlotVals=log(decilePlotVals);

        figure;
        hold on;
        plot(logModPlotVals,logObsPlotVals,'.');
        % monotonically adjusted lambda1 and unadjusted lambda2 plot values
        plot(logCAMPmodelPlotVals,logLambda1PlotVals,'-*','LineWidth',2);
        plot(logCAMPmodelPlotVals,logLambda2PlotVals,'-*','LineWidth',2);
        % plot decile lines
        xline(logDecilePlotVals);
        % 1:1 line
        plot(-10:3,-10:3)
        xlim([min(logModPlotVals) max(logModPlotVals)]);
        xlabel('Model Values (ln-ppb)');
        ylabel('Obs Values (ln-ppb)');
        legend('Model Value','Monotonic \lambda_1','Unadjusted lambda_2^1^/^2','1:1 Line');
        title('Model Value vs. Observed Value and Std. at each Decile (log-ppb)');
        figName2=sprintf('%s_CAMP_logt',figName);
        figName2=(['..\analysis\8CAMx\figs\' figName2]);
        print(figName2, '-dpng');
    else
        disp("Error. Input valid 'lt' value");
    end  
end


lambda1=plotVals(:,1);
lambda2=((plotVals(:,3)).^2);
decileModelVal=plotVals(:,2);






