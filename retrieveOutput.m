function [deltaVec, omegaVec, eVec, mVec,...
    thetaVec, vVec, pgVec, qgVec, ...
    prefVec,fVec, ...
    pdVec, qdVec,yVec,...
    deltaDotVec, omegaDotVec, eDotVec, mDotVec,ACEVec] = retrieveOutput( t, ZNEW , NoiseVector)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global ControlMode
% system constants
global  N G  GenSet 

% indices
global deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  


    

% dynamical simulations
global TFinal TPert FSample NSamples NPertSamples Mass...
PertSet PPertValues QPertValues NoiseVarianceSet 




global pd0 qd0


if strcmp(ControlMode,'AGC')
global y0 y0Plus yS yDot0Plus yDotS yIdx...
    ParticipationFactors NumberOfAreas AreaSet TieLineFromSet TieLineToSet...
    ACE0Plus PScheduledS GensPerArea BusesPerArea...
    KI KACE KPG KPflow KSumPG KThetaSlack

end


disp('Retrieving the output');
% retrieving time-varying quantities: 
deltaVec=ZNEW(deltaIdx,:); 
omegaVec=ZNEW(omegaIdx,:); 
eVec=ZNEW(eIdx,:); 
mVec=ZNEW(mIdx,:);
thetaVec=ZNEW(mIdx(end)+thetaIdx,:); 
vVec=ZNEW(mIdx(end)+vIdx,:); 
pgVec=ZNEW(mIdx(end)+pgIdx,:); 
qgVec=ZNEW(mIdx(end)+qgIdx,:); 

yVec=[];
if strcmp(ControlMode,'AGC')
    yVec=ZNEW(mIdx(end)+qgIdx(end)+yIdx,:);
end



% retrieving load 
pdVec=zeros(N,NSamples); 
qdVec=zeros(N,NSamples);





prefVec=zeros(G,NSamples); 
fVec=zeros(G,NSamples);


if strcmp(ControlMode,'LQR')
   for tt=1:NSamples
 [prefVec(:,tt), fVec(:,tt)]=   control_law('LQR',deltaVec(:,tt),omegaVec(:,tt),eVec(:,tt),mVec(:,tt),...
                    vVec(:,tt), thetaVec(:,tt), pgVec(:,tt),qgVec(:,tt),[]);
          
    
   end
else
    
      for tt=1:NSamples
 [prefVec(:,tt), fVec(:,tt)]=   control_law('LQR',deltaVec(:,tt),omegaVec(:,tt),eVec(:,tt),mVec(:,tt),...
                    vVec(:,tt), thetaVec(:,tt), pgVec(:,tt),qgVec(:,tt),yVec(:,tt));
          
    
   end
end



% computing xdot:
deltaDotVec=zeros(G,NSamples); 
omegaDotVec=zeros(G,NSamples); 
eDotVec=zeros(G,NSamples); 
mDotVec=zeros(G,NSamples); 

if strcmp(ControlMode,'LQR')
for tt=1:NSamples
[ deltaDotVec(:,tt), omegaDotVec(:,tt), eDotVec(:,tt) , mDotVec(:,tt)] = gTildeFunctionVectorized( ...
    deltaVec(:,tt), omegaVec(:,tt), eVec(:,tt),mVec(:,tt),...
     vVec(:,tt), thetaVec(:,tt),pgVec(:,tt), qgVec(:,tt),[]);
end
else
    for tt=1:NSamples
[ deltaDotVec(:,tt), omegaDotVec(:,tt), eDotVec(:,tt) , mDotVec(:,tt)] = gTildeFunctionVectorized( ...
    deltaVec(:,tt), omegaVec(:,tt), eVec(:,tt),mVec(:,tt),...
     vVec(:,tt), thetaVec(:,tt),pgVec(:,tt), qgVec(:,tt),yVec(:,tt));
    end
end
    

for tt=1:NSamples
    [pdVec(:,tt),qdVec(:,tt)]=loadPert('Transient', t(tt),pd0,qd0,PertSet,PPertValues, QPertValues,TPert,TFinal,NoiseVector) ;   
end



%% Reconstruct the ACE 
ACEVec=[];
if strcmp(ControlMode,'AGC')
ACEVec=zeros(NumberOfAreas,NSamples); 
yDotVec=zeros(G,NSamples);
PMeasuredVec=zeros(NumberOfAreas,NSamples); 
OmegaMeasuredVec=zeros(NumberOfAreas,NSamples);
for tt=1:NSamples

     [yDotVec(:,tt),ACEVec(:,tt), PMeasuredVec(:,tt), OmegaMeasuredVec(:,tt)]=...
         agcParams(omegaVec(:,tt),yVec(:,tt), vVec(:,tt),thetaVec(:,tt), pgVec(:,tt));
     

end
end



end

