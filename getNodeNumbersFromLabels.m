function [OutNumbers ] = getNodeNumbersFromLabels(InLabels )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

global NodeLabels

L=length(InLabels); 
OutNumbers=zeros(L,1);
for ll=1:L
    OutNumbers(ll)=find(NodeLabels==InLabels(ll));
end

end

