function [ dfdz ] = dynamicsJacobian( t,z)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
global deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx 

global Mass

global KLQRstep


delta=z(deltaIdx); 
omega=z(omegaIdx); 
e=z(eIdx);
m=z(mIdx); 
theta=z(mIdx(end)+thetaIdx);
v=z(mIdx(end)+vIdx);
pg=z(mIdx(end)+pgIdx);
qg=z(mIdx(end)+qgIdx);

[ gx,ga,gu] = gFunctionJacobVectorized(...
     delta, omega, e,m,...
     v,theta,pg,qg);
 

 
 [ hx,ha ] = hFunctionJacobVectorized( ...
    delta, omega, e , m , v, theta, pg, qg );

dfdz=[[gx+gu*KLQRstep, ga];[hx,ha]];


end

