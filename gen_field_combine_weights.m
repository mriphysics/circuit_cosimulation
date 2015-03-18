function a_K = gen_field_combine_weights(S,freqs,nPorts,nLumped,lumpedVals,varargin)
% Function to generate weights a_K to apply to s-parameter ports and
% generate EM fields for a given drive at a given frequency
% 
%     gen_sparam_response(S,freqs,y0mat,nPorts,nLumped,lumpedVals)
% 
%     S          = matrix variable of all s-parameter matrices for all
%                  ports with size (nPorts +nLumped) x (nPorts +nLumped) x 
%                  (number of frequency points)
%
%     freqs      = vector of frequencies in Hz
%     y0mat      = characteristic admittance matrix
%     nPorts     = number of driving ports
%     nLumped    = number of lumped element ports
%     lumpedVals = [R C L] = nLumped x 3 array - units are in Ohms, Farads 
%                  and Henries respectively
%
%     Default characteristic admittance is set as z0 = 50 Ohm. Add argument
%     'z0' to modify this
%     Default selected operational frequency is set as 128MHz. Add argument
%     'selfreq' to modify this
%
% Created by Arian Beqiri, King's College London, (C) October 2014.
% Email: arian.beqiri@kcl.ac.uk
%
% This code is free for use but please cite Beqiri et al. 2014, 
% DOI: 10.1002/mrm.25504 if used for any published work

%% =========== Set up variables and S-parameter sub-matrices ==============

% Assign default caracteristic impedance of ports as 50 Ohms
z0 = 50;

% Assign default selected frequency as 128 MHz
selFreq = 128e6;

% Check args
for ii=1:length(varargin)
    
    % Allocate characteristic impedance of all ports
    if strcmpi(varargin{ii},'z0')
        z0 = varargin{ii+1};
    end
    % Define selected frequency to run calculation at
    if strcmpi(varargin{ii},'selFreq')
        selFreq = varargin{ii+1};
    end
end

% Create necessary sub-matrices from S-parameter matrices set
sm.large = S((nPorts+1):(nPorts+nLumped),(nPorts+1):(nPorts+nLumped),:);
sm.left = S((nPorts+1):(nPorts+nLumped),1:nPorts,:);

% Characteristic Admittances
y0mat = inv(diag(sqrt(z0)*ones(size(sm.large,1),1)));

%% ============= Find frequency index and squeeze matrices ================

% Find index of correct frequency
[~,ind] = min(abs(freqs-selFreq));

% Creat sub-matrices at resonant frequency
s_large_ind = sm.large(:,:,ind);
s_left_ind = sm.left(:,:,ind);

%% =================== Define Impedance Matrices ==========================

% Impedance matrix (z-parameters) at selected frequency - calculates serial
% impedance of lumped elements
z_cap = diag(lumpedVals(:,1) + 1./(1i*2*pi*freqs(ind)*lumpedVals(:,2))...
        + 1i*2*pi*freqs(ind)*lumpedVals(:,3));

% S-parameters from z-parameters and characteristic admittances
sigma = (y0mat*z_cap*y0mat + eye(nLumped))\...
        (y0mat*z_cap*y0mat - eye(nLumped));

%% ============ Calculate Field Combination weights, a_k ==================

% Preallocate for final matrix
a_K = zeros(nLumped + nPorts, nPorts);

% Preallocate for port drives
a = [1; zeros((nPorts-1),1)];

% Calculate all port weights
for ii=1:nPorts 
    a_port = circshift(a,(ii-1));
    a_K(:,ii) = [a_port; sigma*((eye(nLumped) - s_large_ind*sigma)\...
                (s_left_ind*a_port))];    
end

end
