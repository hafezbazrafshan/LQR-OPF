function [ yDot, ACE, Pmeasured, OmegaMeasured] = agcParams(omega,y, v,theta, pg)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



% system constants [these do not change]
global OMEGA_S Sbase N G L node_set gen_set load_set Ymat Gmat Bmat Cg...
    yff_vec yft_vec  ytf_vec ytt_vec

%  indices [these  do not change]
global deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  fIdx prefIdx   zIdx 

% machine [these do not change]
global  tau_vec xd_vec xq_vec xprime_vec d_vec m_vec Tch_vec freqR_vec ...
    

% dynamical simulations 
global Tfinal Tpert fsample n_samples n_pertSamples Mass...
pertSet pPertValues qPertValues NoiseVarianceSet 





% next time-slot conditions
global xS omegaS deltaS eS mS...
    aS vS thetaS pgS qgS...
    uS prefS fS...
    pdS qdS...
    vgS thetagS...
    zS networkS
 
 % global 
global  KLQRstep


global y0 y0plus yS yDot0plus yDotS yIdx...
    ParticipationFactors NumberOfAreas AreaSet TieLineFromSet TieLineToSet...
    ACE0plus PscheduledS GensPerArea BusesPerArea...
    K_I K_ACE K_PG K_Pflow K_sumPG K_ThetaSlack
slackIdx=find(networkS.bus(:,2)==3);

     Pmeasured=zeros(NumberOfAreas,1);
     OmegaMeasured=zeros(NumberOfAreas,1);

ACE=zeros(NumberOfAreas,1);
yDot=zeros(G,1);
    for ii=1:NumberOfAreas
        [pFrom,~]=determineLineFlows(TieLineFromSet{ii,1},v,theta);
        [~,pTo]=determineLineFlows(TieLineToSet{ii,1},v,theta);
        Pmeasured(ii)=sum(pFrom)-sum(pTo);
        OmegaMeasured(ii)=mean(omega(GensPerArea{ii,1})-OMEGA_S);
        ACE(ii)=K_Pflow*( Pmeasured(ii)-PscheduledS(ii,1))+ sum((1./freqR_vec(GensPerArea{ii,1})+...
            d_vec(GensPerArea{ii,1}))).*mean(omega(GensPerArea{ii,1})-OMEGA_S)./(2*pi);
yDot(GensPerArea{ii,1})=K_I*(-y(GensPerArea{ii,1})-K_ACE.*ACE(ii)+...
    K_sumPG.*sum(pgS(GensPerArea{ii,1}))+K_PG.*(pgS(GensPerArea{ii,1})-pg(GensPerArea{ii,1})))+...
    K_ThetaSlack*(thetaS(slackIdx)-theta(slackIdx));
       


    end
end

