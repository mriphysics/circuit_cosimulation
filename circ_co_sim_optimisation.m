% Circuit Co-Simulation Optimisation function - Top-level script to attempt
% to find optimal capacitor values for a given EM simulation using circuit
% simulation at a given frequency
%
%     circ_co_sim_optimisation(S,freqs,nPorts,nLumped,capValsInit,
%                                portCapInds,lambda,varargin)
%
%     S           = matrix variable of all S-parameter matrices for all
%                   ports with size (nPorts +nLumped) x (nPorts +nLumped) x 
%                   (number of frequency points)
%
%     freqs       = vector of frequencies in Hz
%     nPorts      = number of driving ports
%     nLumped     = number of lumped element ports
%     capValsInit = vector of initial set of capacitor values. This is not
%                   likely to be of length nLumped due to symmetries within
%                   simulation models - the number of different capacitor
%                   values will be smaller
%
%     portCapInds = cell of size 1 x (length(capValsInit)). Each cell entry
%                   contains indices of lumped element ports corresponding 
%                   to the different capacitor values.
%                   N.B. these indices are numbered from 1 to nLumped and 
%                   should be conisdered as such when the indices are 
%                   defined but the actual EM simulation S-paramter ports 
%                   they correspond to will have index = index + nPorts
% 
%     lambda      = tuning paramter for optimisation. As lambda is
%                   increased, optimisation favours minimisation of 
%                   off-diagonal elements of S-matrix (i.e. minimise 
%                   inter-element coupling)
%
%     Default characteristic admittance is set as z0 = 50 Ohm. Add argument
%     'z0' to modify this
%     Default selected operational frequency is set as 128MHz. Add argument
%     'selfreq' to modify this
%     Default resistance of lumped element ports is set to R = 0.3 Ohms. 
%     Add argument 'R' to modify this
%     Default inductance of lumped element ports is set to L = 0 H. 
%     Add argument 'L' to modify this
%
% Created by Arian Beqiri, King's College London, October 2014.
% Email: arian.beqiri@kcl.ac.uk

% This code is free under the terms of the MIT license. Please cite Beqiri et al. 2014, 
% DOI: 10.1002/mrm.25504 in any published work

function [lumpedVals,S_out] = circ_co_sim_optimisation(S,freqs,nPorts,nLumped,capValsInit,portCapInds,lambda,varargin)

% Default characteristic impedance
z0 = 50;

% Default calculation frequency
selFreq = 128e6;

% Default serial resistance and inductance of lumped element ports
R = 0.3;
L = 0;

% Check args
for ii=1:length(varargin)
    
    % Define characteristic impedance of all ports
    if strcmpi(varargin{ii},'z0')
        z0 = varargin{ii+1};
    end
    
    % Define selected frequency to run calculation at
    if strcmpi(varargin{ii},'selFreq')
        selFreq = varargin{ii+1};
    end
    
    % Define resistance of all lumped elements
    if strcmpi(varargin{ii},'R')
        R = varargin{ii+1};
    end
    
    % Define inductance of all lumped elements
    if strcmpi(varargin{ii},'L')
        L = varargin{ii+1};
    end
end

% Separate S-Matrix into separate sub-matrices
sm.ports = S(1:nPorts,1:nPorts,:);
sm.large = S((nPorts+1):(nPorts+nLumped),(nPorts+1):(nPorts+nLumped),:);
sm.left = S((nPorts+1):(nPorts+nLumped),1:nPorts,:);
sm.right = S(1:nPorts,(nPorts+1):(nPorts+nLumped),:);

% Characteristic Admittances
y0mat = inv(diag(sqrt(z0)*ones(size(sm.large,1),1)));

%% ============= Find frequency index and squeeze matrices ================

% Find index of correct frequency
[~,ind] = min(abs(freqs-selFreq));

% Creat sub-matrices at resonant frequency
s_ind.ports = sm.ports(:,:,ind);
s_ind.large = sm.large(:,:,ind);
s_ind.left = sm.left(:,:,ind);
s_ind.right = sm.right(:,:,ind);

%% ======= Run optimisation for capacitor values ===============

% Cost function
CostFnHandle = @(Cvals) optimfun_circ(Cvals,...
    portCapInds,s_ind,y0mat,R,L,lambda,nPorts,nLumped,selFreq);

%%% fminsearch (not bounded)
options = optimset('Display','iter','MaxFunEvals',1e10,'TolFun',1e-12,...
        'MaxIter',2000,'TolX',1e-12);
[capValsOpt,~] = fminsearch(CostFnHandle,capValsInit,options);



%% ========== Generate post-optimisation lumped element array =============

% Preallocate and assign first column to resistance
lumpedVals = ones(nLumped,3)*R;

% Assign second column to capacitor values from optimisation
for ii=1:length(capValsOpt)
    lumpedVals(portCapInds{ii},2) = capValsOpt(ii);
end

% Assign third column to inductances of elements
lumpedVals(:,3) = L;

%% ============== Assess Final S-Parameter response =======================

% Output S-parameter curves over frequency range
S_out = gen_sparam_response(S,freqs,nPorts,nLumped,lumpedVals,varargin);

end