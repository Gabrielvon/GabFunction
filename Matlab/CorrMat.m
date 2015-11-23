function [Xa,CorrMat] = CorrMat(mu,sigma,rho,distribution,simpaths)
%**************************************************************************
% Generation Correlation Matrix and simulate the matrix according
% specific distribution from simulations
%
%   CorrMat(mu,sigma,rho,distribution,simpaths)
%
%==========================================================================
% INPUTS:     
%   mu    - average vector
%   sigma - variance vector
%   rho   - correlation coefficient matrix
%   distribution - 'Normal', 'Unique';
%   simpaths = simulation paths.
%
%==========================================================================
% OUTPUTS:
%
%   Xa   - Simulation numbers
%   CorrMat   - Correlation Matrix
%==========================================================================
% EXAMPLE:
%
%       Xa = CorrMat(0.01,0.3,[1.0000,0.5000,0.3500;0,-0.5000,1.0000;...
%               0,0,0.8000],'Normal',10000)
%
%**************************************************************************

    rho=rho+rho';
    CorrMat=diag(sigma)*rho*diag(sigma);
    % something about Jordan Decomposition
    [P,D]=eig(CorrMat);
    K=D.^(.5);
    N=P*K/P ;
    
    %% simulation
    n=simpaths; 
    Y=random(distribution,0,1,n,size(N,1));
    X=Y*N;
    Mu=ones(size(X))*diag(mu);
    Xa=X+Mu;

end