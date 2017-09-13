function [deltaVec, omegaVec, eVec, mVec,...
    thetaVec, vVec, pgVec, qgVec, ...
    prefVec,fVec, ...
    pdVec, qdVec,...
    deltaDotVec, omegaDotVec, eDotVec, mDotVec] = retrieveOutput( t, ZNEW , NoiseVector)
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



% retrieving load 
pdVec=zeros(N,NSamples); 
qdVec=zeros(N,NSamples);





prefVec=zeros(G,NSamples); 
fVec=zeros(G,NSamples);



   for tt=1:NSamples
 [prefVec(:,tt), fVec(:,tt)]=   control_law('LQR',deltaVec(:,tt),omegaVec(:,tt),eVec(:,tt),mVec(:,tt),...
                    vVec(:,tt), thetaVec(:,tt), pgVec(:,tt),qgVec(:,tt));
    
end



% computing xdot:
deltaDotVec=zeros(G,NSamples); 
omegaDotVec=zeros(G,NSamples); 
eDotVec=zeros(G,NSamples); 
mDotVec=zeros(G,NSamples); 

for tt=1:NSamples
[ deltaDotVec(:,tt), omegaDotVec(:,tt), eDotVec(:,tt) , mDotVec(:,tt)] = gFunctionVectorized( ...
    deltaVec(:,tt), omegaVec(:,tt), eVec(:,tt),mVec(:,tt),...
     vVec(GenSet,tt), thetaVec(GenSet,tt),pgVec(:,tt), prefVec(:,tt), fVec(:,tt));
end

for tt=1:NSamples
    [pdVec(:,tt),qdVec(:,tt)]=loadPert('Transient', t(tt),pd0,qd0,PertSet,PPertValues, QPertValues,TPert,TFinal,NoiseVector) ;   
end







end

