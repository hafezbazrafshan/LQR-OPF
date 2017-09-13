function [ JHTilde] = hTildeAlgebraicFunctionJacobVectorized( ...
    delta, omega, e , m , v, theta, pg, qg )
%HFUNCTIONJACOBVECTORIZED calculates the jacobian of g(x,a,u); 
% [ hx, ha,hu ] = hFunctionJacob( z,...
 %    xq_vec, xprime_vec ) calculates the jacobian of h. 
% 
% Description of Outputs: 
% 1. hx: the jacobian of h with respect to x, size(2*N+2*G,4*G); 
% 2. ha: the jacobian of h with respect to a, size(2*N+2*G, 2*N+2*G)
% 3. hu: the jacobian of h with respect to u, size(2*N+2*G, 2*G)
% 
% Description of Inputs: 
% 1. z: vector of z=(x,a,u) combining states, algebraic and control
% variables, size(2*N+7*G,1).
% 5. xq_vec:  vector of quadrature axis synchronous reactance (pu) size(G,1)
% 6. xprime_vec: direct axis transient reactance pu, size(G,1).
%
% See also hFunctionJacob






end

