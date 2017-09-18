function [ pg,qg] =computepgqg(delta,e,v,theta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% system constants [these do not change]
global OMEGAS Sbase N G L NodeSet GenSet LoadSet YMat GMat BMat Cg...
    YffVec YftVec  YtfVec YttVec

%  indices [these  do not change]
global deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  fIdx prefIdx  SlackIdx

% machine [these do not change]
global  TauVec XdVec XqVec XprimeVec DVec MVec TchVec FreqRVec...
    

% dynamical simulations 
global TFinal TPert FSample NSamples NPertSamples Mass...
PertSet PPertValues QPertValues NoiseVarianceSet 



Xprime=diag(sparse(XprimeVec)); 
Xq=diag(sparse(XqVec));
XqInv=mldivide(Xq,speye(size(Xq)));
XprimeInv=mldivide(Xprime,speye(size(Xprime)));

vg=v(GenSet);
thetag=theta(GenSet);

EMat=diag(sparse(e));
V_g=diag(sparse(vg));
VMat=diag(sparse(v));
CosMat=diag( cos(sparse(theta)));
SinMat=diag( sin(sparse(theta)));



pg=XprimeInv*EMat*diag(sparse(sin(delta-thetag)))*vg+0.5*XqInv*XprimeInv*(Xprime-Xq)*diag(sparse(sin(2*(delta-thetag))))*V_g*vg;


                
 qg=XprimeInv*EMat* diag(sparse(cos(delta-thetag)))*vg-...
     0.5*XqInv*XprimeInv*(Xprime+Xq)*V_g*vg+0.5*XqInv*XprimeInv*(Xprime-Xq)*diag(sparse(cos(2*(delta-thetag))))*V_g*vg;
 
 

end

