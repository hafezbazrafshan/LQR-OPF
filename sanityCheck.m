function [ Sanitycheck1,SanityCheck2,Sanitycheck3 , Success] = sanityCheck(...
    deltaVec, omegaVec, eVec, mVec, ...
    thetaVec, vVec, pgVec, qgVec,...
    prefVec, fVec, ...
    ploadVec,qloadVec,...
    deltaDotVec, omegaDotVec, eDotVec, mDotVec)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

global  N G L  GenSet LoadSet... 
 FSample NSamples...


disp('Performing sanity checks on dynamical simulation'); 
% sanity check 1: checking power flow equatios at each time-step:
% check power flow equations
CheckpfVec=zeros(1,NSamples); 
CheckEqsVec=zeros(2*N,NSamples); 
RealGenCheckVec=zeros(G,NSamples); 
ReactiveGenCheckVec=zeros(G,NSamples); 
RealLoadCheckVec=zeros(L,NSamples); 
ReactiveLoadCheckVec=zeros(L,NSamples); 
for tt=1:NSamples
    
pload=ploadVec(:,tt); 
qload=qloadVec(:,tt);
[CheckpfVec(:,tt), CheckEqsVec(:,tt),RealGenCheckVec(:,tt), ReactiveGenCheckVec(:,tt), ...
    RealLoadCheckVec(:,tt),ReactiveLoadCheckVec(:,tt)]=...
   checkPowerFlows(vVec(:,tt),thetaVec(:,tt),pgVec(:,tt),qgVec(:,tt), pload,qload);


end

if all(CheckpfVec==1)
    disp(' Sanity check 1: in transient simulations, power flow equations were all satisfied at every time step.'); 
    Sanitycheck1=1;
else 
    disp(' Sanity check 1: FAILED, power flow equaitons were not satisfied at every time step.'); 
        Sanitycheck1=0;
end

% Sanity check 2:  checking the entire algebraic equations (kind of
% redundant for power flow equations that are already checked but OK)
h1Idx=1:G;
h2Idx=G+1:2*G;
h3Idx=2*G+1:3*G;
h4Idx=3*G+1:4*G;
h5Idx=4*G+1:4*G+L;
h6Idx=4*G+L+1:4*G+2*L;
d=zeros(h6Idx(end),1);

h1Vec=zeros(G,NSamples); 
h2Vec=zeros(G,NSamples); 
h3Vec=zeros(G,NSamples); 
h4Vec=zeros(G,NSamples); 
h5Vec=zeros(L,NSamples); 
h6Vec=zeros(L,NSamples); 
Sanity2Vec=zeros(2*N+2*G,NSamples); 

for tt=1:NSamples

pload=ploadVec(:,tt); 
qload=qloadVec(:,tt);
ploadg=pload(GenSet);
qloadg=qload(GenSet);
ploadl=pload(LoadSet);
qloadl=qload(LoadSet);
d(h3Idx)=-ploadg;
d(h4Idx)=-qloadg;
d(h5Idx)=-ploadl;
d(h6Idx)=-qloadl;   
[Sanity2Vec(:,tt),~]=hTildeAlgebraicFunctionVectorized(deltaVec(:,tt),omegaVec(:,tt), eVec(:,tt),mVec(:,tt),...
    vVec(:,tt),thetaVec(:,tt),pgVec(:,tt),qgVec(:,tt),...
    d);


end

if max(max(abs(Sanity2Vec)))<1e-3
    disp('Sanity check 2: in transient simulations, all algebraic equations were satisfied.'); 
    SanityCheck2=1;
else 
    disp('Sanity check 2: FAILED, some algebraic equations were not satisfied'); 
    SanityCheck2=0;
end

% sanity check 3: comparison of calculated state derivative and numerical
% derivative:

Err1=abs(deltaDotVec(:,1:end-1) - numericalDerivative(deltaVec,1/FSample)); 
Err2=abs(omegaDotVec(:,1:end-1)- numericalDerivative(omegaVec,1/FSample));
Err3=abs(eDotVec(:,1:end-1) - numericalDerivative(eVec,1/FSample)); 
 Err4=abs(mDotVec(:,1:end-1)- numericalDerivative(mVec, 1/FSample)); 

Err1=Err1(:,2:end); 
Err2=Err2(:,2:end); 
Err3=Err3(:,2:end);
Err4=Err4(:,2:end);

if sum(max([Err1;Err2;Err3;Err4])>1e-1)<0.2*NSamples
    disp('Sanity check 3: numerical and computed time derivatives for states match'); 
    Sanitycheck3=1;
else 
      disp('Sanity check 3: FAILED, numerical and computed time derivatives for states DO NOT match'); 
      Sanitycheck3=0;
end


if sum([Sanitycheck1,SanityCheck2,Sanitycheck3])==3
    Success=1; 
else 
    Success=0;
end

end

