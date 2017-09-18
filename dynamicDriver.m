function [ out ] = dynamicDriver(SteadyStateData, LfControl, varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

load(SteadyStateData);
YesPlots=false;
if length(varargin)>0
    YesPlots=true;
end

global ControlMode
ControlMode=LfControl;

global SteadyStateMode
SteadyStateMode=SsControl;


%% Defining some global variables
% Some of these variables refer to the network. For example, Sbase, Ymat
% (the bus admittance matrix), etc.  Some others are specific indices
% within a vector.  Some are LQR parameters.  
% These variables are declared global to avoid extra function
% arguments. 



% dynamical simulations 
global TFinal TPert FSample NSamples NPertSamples Mass...
 NoiseVarianceSet 





% 0plus conditions
global x0Plus omega0Plus delta0Plus e0Plus m0Plus...
    a0Plus v0Plus theta0Plus pg0Plus qg0Plus...
    u0Plus pref0Plus f0Plus...
    pd0Plus qd0Plus...
    vg0Plus thetag0Plus...
    z0Plus...
    deltaDot0Plus omegaDot0Plus eDot0Plus mDot0Plus






if strcmp(ControlMode,'AGC')
global y0 y0Plus yS yDot0Plus yDotS yIdx...
    ParticipationFactors NumberOfAreas AreaSet TieLineFromSet TieLineToSet...
    ACE0Plus PScheduledS GensPerArea BusesPerArea...
    KI KACE KPG KPflow KSumPG KThetaSlack

end
%% Areas:
if strcmp(ControlMode,'AGC')
AreaSet=unique(Network.bus(:,7)); 
NumberOfAreas=length(AreaSet); 
BusesPerArea=cell(NumberOfAreas,1); 
GensPerArea=cell(NumberOfAreas,1);
TieLineFromSet=cell(NumberOfAreas,1);
TieLineToSet=cell(NumberOfAreas,1);
for ii=1:NumberOfAreas
    BusesPerArea{ii,1}=find(Network.bus(:,7)==AreaSet(ii)); 
    GensPerArea{ii,1}=find(Network.bus(GenSet,7)==AreaSet(ii));
%     TieLineFromSet{ii,1}= find(and(Network.bus(Network.branch(:,1),7)==AreaSet(ii),Network.bus(Network.branch(:,2),7)~=AreaSet(ii)));
%     TieLineToSet{ii,1}=find(and(Network.bus(Network.branch(:,1),7)~=AreaSet(ii),Network.bus(Network.branch(:,2),7)==AreaSet(ii)));
    TieLineFromSet{ii,1}= find(and(Network.bus(getNodeNumbersFromLabels(Network.branch(:,1) ),7)==AreaSet(ii),...
        Network.bus(getNodeNumbersFromLabels(Network.branch(:,2) ),7)~=AreaSet(ii)));
    TieLineToSet{ii,1}=find(and(Network.bus(getNodeNumbersFromLabels(Network.branch(:,1) ),7)~=AreaSet(ii),...
        Network.bus(getNodeNumbersFromLabels(Network.branch(:,2) ),7)==AreaSet(ii)));
end
end



%% Dynamical simulation section
if strcmp(ControlMode,'AGC')
yIdx=(1:G).';
end

if strcmp(ControlMode,'AGC')
y0=zeros(G,1); 
end

%  12. Define the MASS matrix (The E matrix in $E\dot{x}$ descriptor systems)
if strcmp(ControlMode,'LQR')
    Mass=zeros(length(x0)+length(a0), length(x0)+length(a0));
Mass(sub2ind(size(Mass), [deltaIdx,omegaIdx,eIdx,mIdx],[deltaIdx,omegaIdx,eIdx,mIdx]))=1;
elseif strcmp(ControlMode,'AGC')
Mass=zeros(length(x0)+length(a0)+length(y0), length(x0)+ length(a0)+length(y0));
Mass(sub2ind(size(Mass), [deltaIdx,omegaIdx,eIdx,mIdx],[deltaIdx,omegaIdx,eIdx,mIdx]))=1;
Mass(sub2ind(size(Mass), [length(x0)+length(a0)+yIdx],[length(x0)+length(a0)+yIdx]))=1;
end
 %% 13. Set dynamical simulation parameters:
DynamicSolverOptions = odeset('Mass',Mass,'MassSingular', 'yes','MStateDependence','none', ...
    'RelTol',1e-7,'AbsTol',1e-6,'Stats','off');

TFinal=25;
TPert=0; % NEEDS TO BE CLOSE TO t=0 for increased accuracy 
FSample = 100;
NSamples = TFinal * FSample+1;
NPertSamples = max(TPert,0) * FSample+1;
t = 0:1/FSample:TFinal;
NoiseVarianceSet=0*Network.bus(:,3)/Sbase;
NoiseVector=repmat(NoiseVarianceSet,1,NSamples).*randn(length(NoiseVarianceSet),NSamples);

 
 %% Simulation intial conditions at zero plus
display('Configuring 0plus intial conditions due to disturbance');
[pd0Plus,qd0Plus]=loadPert('Transient',0,pd0,qd0,PertSet, PPertValues,QPertValues, TPert,TFinal,NoiseVector);
delta0Plus=delta0;
omega0Plus=omega0;
e0Plus=e0;
m0Plus=m0;
h1Idx=(1:G).';
h2Idx=(G+1:2*G).';
h3Idx=(2*G+1:3*G).';
h4Idx=(3*G+1:4*G).';
h5Idx=(4*G+1:4*G+L).';
h6Idx=(4*G+L+1:4*G+2*L).';
d0Plus=zeros(h6Idx(end),1);
pdg0Plus=pd0Plus(GenSet);
qdg0Plus=qd0Plus(GenSet);
pdl0Plus=pd0Plus(LoadSet);
qdl0Plus=qd0Plus(LoadSet);
d0Plus(h3Idx)=-pdg0Plus;
d0Plus(h4Idx)=-qdg0Plus;
d0Plus(h5Idx)=-pdl0Plus;
d0Plus(h6Idx)=-qdl0Plus;








% % from this point on we would like to have a closed loop system using
% % KLQRstep and the optimal controller set-points 
% % 
% % First
    disp('Solving for the 0plus initial conditions using levenberg-marquardt');
InitialOptions= optimoptions('fsolve','Display','Iter','Algorithm','levenberg-marquardt','InitDamping',0.5, 'ScaleProblem','jacobian',...
    'SpecifyObjectiveGradient',true,'MaxIterations',100,'MaxFunctionEvaluations',200,'OptimalityTolerance',1e-6);
% 
[a0Plus,Res0Plus, exitflag]=fsolve(@( a) hTildeAlgebraicFunctionVectorized(...
    delta0Plus,omega0Plus,e0Plus,m0Plus, ...
    a(vIdx), a(thetaIdx), a(pgIdx), a(qgIdx),...
  d0Plus), [v0;theta0;pg0;qg0],InitialOptions);

% [a10Plus,Res0Plus,exitflag]=fsolve(@(a1) hAlgebraic(delta0Plus,omega0Plus,e0Plus,m0Plus,...
%     a1(1:N), a1(N+1:2*N), pd0Plus,qd0Plus), [v0;theta0], InitialOptions);


% v0Plus=a10Plus(1:N); 
% theta0Plus=a10Plus(N+1:2*N); 
% 
% [pg0Plus,qg0Plus]=computepgqg(delta0Plus,e0Plus,v0Plus,theta0Plus);


v0Plus=a0Plus(vIdx);
theta0Plus=a0Plus(thetaIdx); 
pg0Plus=a0Plus(pgIdx); 
qg0Plus=a0Plus(qgIdx);

if exitflag
   disp([PreSuccessStr,'Successful']);
   pause(PauseTime);
else
        disp([PreSuccessStr,'Failed!!']);
pause;
end
% InitialOptions= optimoptions('fsolve','Display','Iter','Algorithm','trust-region','InitDamping',0.1, 'ScaleProblem','jacobian',...
%     'SpecifyObjectiveGradient',true,'MaxIterations',100000,'MaxFunctionEvaluations',50000,'OptimalityTolerance',1e-6,...
%     'PlotFcn',@optimplotfval);
% [a0Plus]=fsolve(@( a) hTildeAlgebraicFunctionVectorized(...
%     delta0Plus,omega0Plus,e0Plus,m0Plus, ...
%     a(vIdx), a(thetaIdx), a(pgIdx), a(qgIdx),...
%   d0Plus), a0Plus,InitialOptions);
% v0Plus=a0Plus(vIdx);
% theta0Plus=a0Plus(thetaIdx); 
% pg0Plus=a0Plus(pgIdx); 
% qg0Plus=a0Plus(qgIdx);


% [ deltaDot0Plus, omegaDot0Plus, eDot0Plus, mDot0Plus ] = gTildeFunctionVectorized(...
%     delta0Plus, omega0Plus, e0Plus,m0Plus,...
%      v0Plus,theta0Plus,pg0Plus,qg0Plus);




if strcmp(ControlMode,'AGC')
     ParticipationFactors=cell(NumberOfAreas,1);
    y0Plus=zeros(G,1);
     PScheduledS=zeros(NumberOfAreas,1);
     KI=1000;
     KACE=	1;
     KPG=0;
     KSumPG=1;
     KPflow=1;
     KThetaSlack=0;
     for ii=1:NumberOfAreas
          [PFromS,~]=determineLineFlows(TieLineFromSet{ii,1},vS,thetaS);
               [~,PToS]=determineLineFlows(TieLineToSet{ii,1},vS,thetaS); 
               PScheduledS(ii,1)=sum(PFromS)-sum(PToS);
ParticipationFactors{ii,1}=pgS(GensPerArea{ii,1})./sum(pgS(GensPerArea{ii,1}));

     end

[yDot0Plus, ACE0Plus, PMeasured0plus, OmegaMeasured0plus] = agcParams(omega0Plus,y0Plus, v0Plus,theta0Plus, pg0Plus);
end



%%

    [ deltaDot0Plus, omegaDot0Plus, eDot0Plus, mDot0Plus ] = gFunctionVectorized(...
    delta0Plus, omega0Plus, e0Plus,m0Plus,...
     v0Plus(GenSet),theta0Plus(GenSet),pg0Plus,pref0,f0);
%  
% 
disp('Running dynamical simulations');
if strcmp(ControlMode,'LQR')
DynamicSolverOptions.Jacobian=@dynamicsJacobian; % this is for the ode solver
end
Znew0=[delta0Plus; omega0Plus; e0Plus; m0Plus; v0Plus; theta0Plus; pg0Plus; qg0Plus];
if strcmp(ControlMode,'AGC')
    Znew0=[Znew0;y0Plus];
end
ZDot0Plus=zeros(size(Znew0));
ZDot0Plus([deltaIdx;omegaIdx;eIdx;mIdx])=[deltaDot0Plus;omegaDot0Plus;eDot0Plus;mDot0Plus];
if strcmp(ControlMode,'AGC')
    ZDot0Plus(mIdx(end)+qgIdx(end)+yIdx)=yDot0Plus;
end
DynamicSolverOptions.InitialSlope=ZDot0Plus;
% 
% 
% 
% %%


[~,ZNEW]=ode15s(@(t,znew)...
    runDynamics(t,znew,NoiseVector), t, Znew0, DynamicSolverOptions);
ZNEW = transpose(ZNEW);
  disp([PreSuccessStr,'Successful']);
   pause(PauseTime);
 
% 
% 
% 
% 
% 
% 
% 
%  
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% 
% 
disp('Retrieving output states, algebraic variables, and controls as a function of time');
[deltaVec, omegaVec, eVec, mVec,...
    thetaVec, vVec, pgVec, qgVec, ...
    prefVec,fVec, ...
    ploadVec, qloadVec,yVec,...
    deltaDotVec, omegaDotVec, eDotVec, mDotVec,ACEVec] =...
    retrieveOutput( t, ZNEW , NoiseVector); 

  disp([PreSuccessStr,'Successful']);
   pause(PauseTime);
% 
disp('Validating results of the dynamical simulation') 
[ SanityCheck1,SanityCheck2,SanityCheck3 , Success] = sanityCheck(...
    deltaVec, omegaVec, eVec, mVec, ...
    thetaVec, vVec, pgVec, qgVec,...
    prefVec, fVec, ...
    ploadVec,qloadVec,...
    deltaDotVec, omegaDotVec, eDotVec, mDotVec);
if Success==1
  disp([PreSuccessStr,'Successful']);
   pause(PauseTime);
else 
    disp([PreSuccessStr,'Failed!']); 
    disp('Unfortunately, dynamical results are not reliable'); 
    pause;
end
    
%
disp('Evaluating dynamical costs');
[ TrCost] = calculateTrCostUsingIntegration(pgS, qgS, Alpha,...
    deltaVec, omegaVec, eVec, mVec, prefVec, fVec, ...
    deltaS, omegaS, eS, mS, prefS, fS,Tlqr);
  disp([PreSuccessStr,'Successful']);
   pause(PauseTime);



TotalCost=SsCost+TrCost;

MaxFreqDev=max(max(abs(omegaVec-OMEGAS)))./(2*pi);
MaxVoltDev=max(max(abs(vVec-repmat(vS,1,length(t)))));


% 
% 

% 
cd(SteadyStatePath);

 
if exist(LfControl)~=7
    mkdir(LfControl);
end
% 
cd(LfControl)
DynamicFileName=[SteadyStateFileName,'_',LfControl];
DynamicPath=pwd;
save(DynamicFileName); 
out=load(DynamicFileName);
cd(CurrentDirectory);

if (YesPlots) 
    plotsForSolvedCase(out); 
end
end

