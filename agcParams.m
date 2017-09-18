function [ yDot, ACE, PMeasured, OmegaMeasured] = agcParams(omega,y, v,theta, pg)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% system constants [these do not change]
global OMEGAS Sbase N G L NodeSet GenSet LoadSet YMat GMat BMat Cg...
    YffVec YftVec  YtfVec YttVec

%  indices [these  do not change]
global deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  fIdx prefIdx  SlackIdx GenSlackIdx NonSlackSet GenNonSlackSet


% machine [these do not change]
global  TauVec XdVec XqVec XprimeVec DVec MVec TchVec FreqRVec...
    

% dynamical simulations 
global TFinal TPert FSample NSamples NPertSamples Mass...
PertSet PPertValues QPertValues NoiseVarianceSet 

% system constants [these do not change]
global y0 y0Plus yS yDot0Plus yDotS yIdx...
    ParticipationFactors NumberOfAreas AreaSet TieLineFromSet TieLineToSet...
    ACE0Plus PScheduledS GensPerArea BusesPerArea...
    KI KACE KPG KPflow KSumPG KThetaSlack



% next time-slot conditions
global xS omegaS deltaS eS mS...
    aS vS thetaS pgS qgS...
    uS prefS fS...
    pdS qdS...
    vgS thetagS...
    zS NetworkS




     PMeasured=zeros(NumberOfAreas,1);
     OmegaMeasured=zeros(NumberOfAreas,1);

ACE=zeros(NumberOfAreas,1);
yDot=zeros(G,1);
    for ii=1:NumberOfAreas
        [pFrom,~]=determineLineFlows(TieLineFromSet{ii,1},v,theta);
        [~,pTo]=determineLineFlows(TieLineToSet{ii,1},v,theta);
        PMeasured(ii)=sum(pFrom)-sum(pTo);
        OmegaMeasured(ii)=mean(omega(GensPerArea{ii,1})-OMEGAS);
        ACE(ii)=KPflow*( PMeasured(ii)-PScheduledS(ii,1))+ sum((1./FreqRVec(GensPerArea{ii,1})+...
            DVec(GensPerArea{ii,1}))).*mean(omega(GensPerArea{ii,1})-OMEGAS)./(2*pi);
yDot(GensPerArea{ii,1})=KI*(-y(GensPerArea{ii,1})-KACE.*ACE(ii)+...
    KSumPG.*sum(pgS(GensPerArea{ii,1}))+KPG.*(pgS(GensPerArea{ii,1})-pg(GensPerArea{ii,1})));

       


    end
end

