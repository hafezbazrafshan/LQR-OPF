function [Pload,Qload, NoiseVector]=loadPert(TimeMod, t,Pload0,Qload0,PertSet, PPertValues, QPertValues,Tpert,Tfinal,NoiseVector)
%LOADPERT produces a new load starting at Tpert
%   [pload,qload]=loadPert(timeMod, t,pload0,qload0,nodepert_set, pertvalue_set, Tpert)
%   produces a new  real and reactive power load at Tpert 
% 
% Description of Outputs: 
% 1. pload: the new real power load of size(N,1);
% 2. qload: the new reactive power load of size(N,1);
% 
% Description of Inputs:
% 1. timeMod: a switch option for 'Steady-State' or 'Transient'. 
% 2. t: the time instant (a scalar). If timeMod is set to 'Steady-State' is not used.
% 3. pload0: the initial real power load of size(N,1); 
% 4. qload0: the intial reactive power load of size(N,1); 
% 5. nodepert_set: set of nodes in which we would like the load to be
% perturbed. 
% 6. pertvalue_set: set of pertubration values corresponding to
% nodepert_set
% 7. Tpert: time of load perturbation. 
%
% Required:
% 1. Extend to Tpert vector

global NoiseVarianceSet FSample

switch TimeMod
    
    case 'Steady-State'
 N=length(Pload0);
pstep=zeros(N,1);
pstep(PertSet)=PPertValues;
qstep=zeros(N,1);
qstep(PertSet)=QPertValues;

Pload= Pload0 +pstep;
Qload=Qload0+qstep;
        
    case 'Transient'
% N=length(pload0);
% step=zeros(N,1);
% step(pertSet)=pPertValues;
% Tramp=Tfinal/2;
% 
noise=zeros(size(NoiseVector(:,1)));
if t>0
n1=floor(t*FSample);
n2=ceil(t*FSample);
alpha=n2-t*FSample;
noise=alpha*NoiseVector(:,n1+1)+(1-alpha)*NoiseVector(:,n2+1);
end
% 
% if t<=Tramp
% pload=pload0+step.*((t-Tpert)./(Tramp-Tpert)).*(t>Tpert)'+noise.*(t>Tpert)';
% else 
%     pload=pload0+step.*(t>Tpert).'+noise.*(t>Tpert)';
% %     qload=qload0+step.*(t>Tpert).'+noise.*(t>Tpert)';
% end
% 
% qload=qload0;


 N=length(Pload0);
pstep=zeros(N,1);
pstep(PertSet)=PPertValues;
qstep=zeros(N,1);
qstep(PertSet)=QPertValues;

Pload= Pload0 +pstep+noise;
Qload=Qload0+qstep;
        
end



end

