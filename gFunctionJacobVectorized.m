function [ gx,ga,gu ] = gFunctionJacobVectorized(delta, omega, e, m, v, theta, pg, qg)
%GFUNCTIONJACOBVECTORIZED calculates the jacobian of g(x,a,u); 
% [gx,ga,gu]= gFunctionJacobVectorized(z, m_vec, d_vec, tau_vec, xd_vec,xprime_vec)
% calculates the jacobian of g(x,a,u). 
% 
% Description of Outputs: 
% 1. gx: the jacobian of g with respect to x, size(4*G,4*G); 
% 2. ga: the jacobian of g with respect to a, size(4*G, 2*N+2*G)
% 3. gu: the jacobian of g with respect to u, size(4*G, 2*G)
% 
% Description of Inputs: 
% 1. z: vector of z=(x,a,u) combining states, algebraic and control
% variables
% 2. m_vec: vector of generator inertia constants in pu, size(G,1). 
% 3. d_vec: vector of damping coefficients, size(G,1)
% 4. tau_vec: vector of direct axis transient open-circuit time constant,
% size(G,1). 
% 5. xd_vec:  vector of direct axis synchronous reactance pu, size(G,1).
% 6. xprime_vec: direct axis transient reactance pu, size(G,1).
% 7. Tch_vec
% 8. freqR_vec
% See also gFunctionJacob

global OMEGAS Sbase N G L NodeSet GenSet LoadSet YMat GMat BMat Cg...
    YffVec YftVec  YtfVec YttVec

%  indices [these  do not change]
global  deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  fIdx prefIdx  

% machine [these do not change]
global  TauVec XdVec XqVec XprimeVec DVec MVec TchVec FreqRVec...

% converting to matrices
M=diag(sparse(MVec));
D=diag(sparse(DVec));
T=diag(sparse(TauVec));
Xprime=diag(sparse(XprimeVec)); 
Xd=diag(sparse(XdVec));
Tch=diag(sparse(TchVec)); 
FreqR=diag(sparse(FreqRVec)); 




Vg=v(sparse(GenSet));
V_g=diag(sparse(Vg));
thetag=theta(sparse(GenSet));
VgIdx=vIdx(sparse(GenSet));
thetagIdx=thetaIdx(sparse(GenSet));





gx=sparse(length([delta; omega;e; m]), length([delta; omega;e; m])); 
ga=sparse(length([delta; omega;e; m]), length([v;theta;pg;qg]));
gu=sparse(length([delta; omega;e; m]),length([prefIdx;fIdx]));



g1Idx=(1:G).'; 
g2Idx=(g1Idx(end)+1:g1Idx(end)+G).';
g3Idx=(g2Idx(end)+1:g2Idx(end)+G).';
g4Idx=(g3Idx(end)+1:g3Idx(end)+G).';

gx(g1Idx, omegaIdx) =speye(G);
gx( g2Idx, omegaIdx)=- mldivide(M,speye(size(M)))*D;
gx(g2Idx,mIdx)=mldivide(M,speye(size(M))); 
gx(g3Idx,deltaIdx)= -mldivide(T,speye(size(T)))*mldivide(M,speye(size(M)))*(Xd-Xprime)*V_g*diag(sin(delta-thetag)); 
gx(g3Idx,eIdx)=-mldivide(T,speye(size(T)))*mldivide(Xprime,speye(size(Xprime)))*Xd;
gx(g4Idx,mIdx)=-mldivide(Tch,speye(size(Tch)));
gx(g4Idx,omegaIdx)=-mldivide(Tch,speye(size(Tch)))*mldivide(FreqR,speye(size(FreqR)));



ga(g2Idx,pgIdx)=-mldivide(M,speye(size(M)));
ga(g3Idx,VgIdx)=mldivide(T,speye(size(T)))*mldivide(Xprime,speye(size(Xprime)))*(Xd-Xprime)*diag(cos(delta-thetag));
ga(g3Idx,thetagIdx) =mldivide(T,speye(size(T)))*mldivide(Xprime,speye(size(Xprime)))*(Xd-Xprime)*V_g*diag(sin(delta-thetag)); 


gu(g3Idx, fIdx)=mldivide(T,speye(size(T)));
gu(g4Idx,prefIdx)=mldivide(Tch,speye(size(Tch)));






end

