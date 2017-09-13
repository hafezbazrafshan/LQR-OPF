function [hTilde,JHTilde]=hTildeAlgebraicFunctionVectorized(delta,omega, e,m,...
    v,theta,pg,qg,...
    d)









[h1,h2,h3,h4,h5,h6] = hFunctionVectorized(...
    delta,omega, e,m,...
v,theta, pg,qg);

hTilde=[[h1;h2;h3;h4;h5;h6]-d];

[~,JHTilde]=hFunctionJacobVectorized(delta,omega,e,m,v,theta,pg,qg);

