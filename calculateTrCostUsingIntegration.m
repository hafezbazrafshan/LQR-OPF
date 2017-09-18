function [ TrCost] = calculateTrCostUsingIntegration(pgS, qgS, alpha,...
    deltaVec, omegaVec, eVec, mVec, prefVec, fVec, ...
    deltaS, omegaS, eS, mS, prefS, fS, Tlqr)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
global G NSamples NPertSamples deltaIdx omegaIdx eIdx mIdx prefIdx fIdx NetworkS


[Qinv,Rinv]=QinvRinv(pgS,qgS,alpha, zeros(4*G), zeros(2*G),NetworkS);
Q=inv(Qinv); 
R=inv(Rinv); 

xMat=[deltaVec;omegaVec;eVec;mVec];
uMat=[prefVec;fVec];
xs=[deltaS;omegaS;eS;mS];
us=[prefS;fS];
xsMat=repmat(xs,1, NSamples-NPertSamples+1); 
usMat=repmat(us, 1, NSamples-NPertSamples+1);
TrCost=Tlqr*calculateTrCostNumeric(xMat,uMat,Q,R, xsMat, usMat); 

end

