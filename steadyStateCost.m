function [ SsCost ] = steadyStateCost(pg, Network)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global G Sbase


c2k=Network.gencost(:,5).*Sbase.^2;
c1k=Network.gencost(:,6).*Sbase; 
c0k=Network.gencost(:,7);


SsCost=quad_form(pg,diag(c2k)) + c1k.'*pg+c0k.'*ones(G,1);


end

