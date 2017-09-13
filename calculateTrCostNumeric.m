function [ TrCost] = calculateTrCostNumeric(x, u, Q, R , xs, us)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global FSample
x=x-xs; 
u=u-us;

TrCost=(1./2/FSample)*sum(diag(x.'*Q*x+u.'*R*u)); 

 
end