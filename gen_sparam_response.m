function s_final = gen_sparam_response(S,freqs,nPorts,nLumped,lumpedVals,varargin)
% Function to output S-parameter response for given lumped element set
% across whole frequency range
% 
%     gen_sparam_response(S,freqs,y0mat,nPorts,nLumped,lumpedVals)
% 
%     S          = matrix variable of all S-parameter matrices for all
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
%
% Created by Arian Beqiri, King's College London, October 2014.
% Email: arian.beqiri@kcl.ac.uk
%
% This code is free under the terms of the MIT license. Please cite Beqiri et al. 2014, 
% DOI: 10.1002/mrm.25504 in any published work

%% =========== Set up variables and S-parameter sub-matrices ==============

% Assign default caracteristic impedance of ports as 50 Ohms
z0 = 50;

% Check args
for ii=1:length(varargin)
    
    % Allocate characteristic impedance of all ports
    if strcmpi(varargin{ii},'z0')
        z0 = varargin{ii+1};
    end
end

% Create sub-matrices from S-parameter matrices set
sm.ports = S(1:nPorts,1:nPorts,:);
sm.large = S((nPorts+1):(nPorts+nLumped),(nPorts+1):(nPorts+nLumped),:);
sm.left = S((nPorts+1):(nPorts+nLumped),1:nPorts,:);
sm.right = S(1:nPorts,(nPorts+1):(nPorts+nLumped),:);

% Characteristic Admittances
y0mat = inv(diag(sqrt(z0)*ones(size(sm.large,1),1)));

%% ================== Calculate Final S-matrix Set ========================

% Preallocate final s_matrix
s_final = zeros(nPorts,nPorts,length(freqs));

% Loop over all frequencies
parfor ii=1:length(freqs)

    % Impedance matrix (z-parameters) at each frequency - calculates
    % serial impedance of lumped element with resistance, capacitance and
    % induactance accounted for
    z_cap = diag(lumpedVals(:,1) + 1./(1i*2*pi*freqs(ii)*...
            lumpedVals(:,2)) + 1i*2*pi*freqs(ii)*lumpedVals(:,3));
    
    % S-parameters from z-parameters and characteristic admittances
    sigma = (y0mat*z_cap*y0mat + eye(nLumped))\...
            (y0mat*z_cap*y0mat - eye(nLumped));
    
    % Squeeze sub-matrices to value at frequency
    sp = squeeze(sm.ports(:,:,ii)); sr = squeeze(sm.right(:,:,ii));
    slar = squeeze(sm.large(:,:,ii)); sl = squeeze(sm.left(:,:,ii));

    % Perform S-matrix calculation
    s_final(:,:,ii) = sp + sr*sigma*(inv(eye(nLumped) - slar*sigma))*sl;

end

%% =================== Plot S-parameter Response ==========================
figure
for ii=1:nPorts
    for jj=1:nPorts
        if ii==jj;c=[1 0 0];else c=[0 0 1];end;
        plot(freqs*1e-6,20*log10(squeeze(abs(s_final(ii,jj,:)))),'color',c)
        hold on
    end
end
legend('Sii','Sij','Location','sw');grid on
set(gca,'fontsize',20)
xlabel('Frequency (MHz)','FontSize',20)
ylabel('Magnitude (dB)','FontSize',20)

end
