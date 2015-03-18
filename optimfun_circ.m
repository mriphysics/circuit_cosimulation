% Optimisation cost function for calculating optimal capcitor values as 
% part of circuit co-simulation framework.
%
% Currently only supports optimisation of capacitor values but can easily
% be modified to optimise all lumped element parameters.
%
% Created by Arian Beqiri, King's College London, November 2014.
% Email: arian.beqiri@kcl.ac.uk
%
% This code is free under the terms of the MIT license. Please cite Beqiri et al. 2014, 
% DOI: 10.1002/mrm.25504 in any published work

function [capValsOpt] = optimfun_circ(capValsInit,portCapInds,s_ind,y0mat,R,L,lambda,nPorts,nLumped,selFreq)
    
    % Create impedance vector for optimisation capacitor values
    capvec = R + 1./(1i*2*pi*selFreq*capValsInit) + (1i*2*pi*selFreq*L);
    
    % Allocate to correct ports
    capvec_full = zeros(nLumped,1);
    for ii=1:length(capvec)
        capvec_full(portCapInds{ii}) = capvec(ii);
    end

    % Impedance matrix (z-parameters)
    z_cap = diag(capvec_full);
    
    % S-parameters from z-parameters and characteristic admittances
    sigma = (inv(y0mat*z_cap*y0mat + eye(nLumped)))*(y0mat*z_cap*y0mat - eye(nLumped));
    
    % Calculate S-matrix for cost-function
    cftmp = ((s_ind.ports + s_ind.right*sigma*(inv(eye(nLumped) -...
            s_ind.large*sigma))*s_ind.left));

    % Evaluate cost-function
    capValsOpt = lambda*max(abs(cftmp(~logical(eye(nPorts))))) +...
                    norm(diag(cftmp));

end