function [CAMxCC]=CAMPcorrection(pairedObs,lambda1,lambda2,decileModelVal,modelVal,modelLoc,tME,filename)

CAMxDir='..\analysis\8CAMx';
if exist(CAMxDir)~=7
    mkdir(CAMxDir)
end

CCval = interp1(decileModelVal,lambda1,modelVal,'linear','extrap');
CCstd = interp1(decileModelVal,lambda2,modelVal,'linear','extrap');

%[CCval,CCloc,CCtME,nanratio]=valstv2stg(ch,CCval);
%[CCstd,CCloc,CCtME,nanratio]=valstv2stg(ch,CCstd);

CAMxCC.lambda1=lambda1;
CAMxCC.lambda2=lambda2;
CAMxCC.decileModelVal=decileModelVal;
CAMxCC.loc=modelLoc;
CAMxCC.tME=tME;
CAMxCC.val=CCval;
CAMxCC.std=CCstd;
CAMxCC.var=CCstd.*CCstd;
CAMxCC.nameShort=pairedObs.nameShort;
CAMxCC.unit=pairedObs.unit;
CAMxCC.nameAndUnit=pairedObs.nameAndUnit;
CAMxCC.grid=pairedObs.grid;
CAMxCC.year=pairedObs.year;
CAMxCC.zeroType=pairedObs.zeroType;
CAMxCC.lt=pairedObs.lt;
CAMxCC.timeUnit=pairedObs.timeUnit;
CAMxCC.spaceUnit=pairedObs.spaceUnit;
CAMxCC.dateZero=pairedObs.dateZero;
CAMxCC.dateVec=pairedObs.dateVec;

calculateLambda2Coverage(lambda1,lambda2,decileModelVal,pairedObs)

save(['..\analysis\8CAMx\' filename], 'CAMxCC')
save ([CAMxDir '\' filename]);




