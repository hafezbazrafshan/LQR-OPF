function [ fz] = runDynamics(t,z,NoiseVector)
%odeFunction15i puts together the equations of the descriptor system based
%on the requirements of the MATLAB ode15i solver. 
%  [ fz] = odeFunction15i(t,znew,zDot, K,z0, zS, pload0, qload0, Tpert,
%  nodepert_set,pertvalue_set) puts the DAE equations of $\dot{x}=g(x,a,u)$ and
%  $d=h(x,a,u)$ of CDC 2016 equations (4) and (5) in the following form:
% $0 = fz= f(t,z, \dot{z}) = [E\dot{z} -gz; hz-d]$
% 
% Description of Outputs: 
% 1. fz: a vector of size(2N+2*G,1) at time t fz= f(t,z, zDot) = [MASS*zDot -gz; hz-d]$
% 
% Description of Inputs: 
% 1. t: the time instant for the functions to be evaluated
% 2. z: the variable z, size(2*N+5*G,1).
% 3. zDot: the derivative of variable z, size(2*N+5*G,1).
% 4. K: the linear feedback gain that relates the controls to states,
% size(2*G, 2*G).
% 5. z0: the initial conditions to start the simulations, where variable z
% starts from
% 6. zS: the desired condition to end the simulations, where we hope that
% the controller u=Kx will take variable z to. 
% 7. pload0: the initial real load condition, size(G,1). 
% 8. qload0; the initial reactive load condition, size(G,1).
% 9. Tpert: the perturbation time instant
% 10. nodepert_set: the set of nodes where load pertubration occurs
% 12. pertvalue_set: the set of load perturbation values corresponding to
% nodepert_set
% 
% Required:
%

% system constants [these do not change]
global ControlMode

global  G L  GenSet LoadSet 


%  indices [these  do not change]
global deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  yIdx

    

% dynamical simulations 
global TFinal TPert ...
    PertSet PPertValues QPertValues 


% initial conditions
global pd0 qd0...
    









delta=z(deltaIdx);
omega=z(omegaIdx);
e=z(eIdx);
m=z(mIdx); 
theta=z(mIdx(end)+thetaIdx);
v=z(mIdx(end)+vIdx);
pg=z(mIdx(end)+pgIdx);
qg=z(mIdx(end)+qgIdx);

y=[];
if strcmp(ControlMode,'AGC')
y=z(mIdx(end)+qgIdx(end)+yIdx);
end



        















h1Idx=(1:G).';
h2Idx=(G+1:2*G).';
h3Idx=(2*G+1:3*G).';
h4Idx=(3*G+1:4*G).';
h5Idx=(4*G+1:4*G+L).';
h6Idx=(4*G+L+1:4*G+2*L).';
d=zeros(h6Idx(end),1);

[pload,qload]=loadPert('Transient',t,pd0,qd0,PertSet, PPertValues,QPertValues, TPert,TFinal,NoiseVector);
ploadg=pload(GenSet);
qloadg=qload(GenSet);
ploadl=pload(LoadSet);
qloadl=qload(LoadSet);
d(h3Idx)=-ploadg;
d(h4Idx)=-qloadg;
d(h5Idx)=-ploadl;
d(h6Idx)=-qloadl;

[ deltaDot, omegaDot, eDot, mDot] =gTildeFunctionVectorized( ...
    delta, omega, e,m,...
    v,theta, pg,qg,y);



gz=[deltaDot; omegaDot; eDot;mDot];


[h1,h2,h3,h4,h5,h6] = hFunctionVectorized(delta,omega,e,m,...
    v,theta,...
    pg,qg);

hz=[h1;h2;h3;h4;h5;h6];

if strcmp(ControlMode,'AGC')

 [ yDot, ACE, PMeasured, OmegaMeasured] = agcParams(omega,y, v,theta, pg);
fz=[gz; hz-d; yDot];
else 
    fz=[gz;hz-d];
end

end

