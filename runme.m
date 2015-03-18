% This test script contains a simple example based on a simulation of the
% Duke phantom placed within an 8-element body coil at 128MHz. The device
% has 8 input ports with an additional 120 lumped elements in the structure
% simulated as ports. This code shows how to take the resulting 128x128 
% S-matrix that exists for each frequency, and return a combined 8x8 
% S-matrix for the principal input ports replacing the lumped elements by
% any complex impedance. 

% Created by Arian Beqiri, King's College London, March 2015.
% Email: arian.beqiri@kcl.ac.uk
%
% This code is free under the terms of the MIT license. Please cite Beqiri et al. 2014, 
% DOI: 10.1002/mrm.25504 in any published work
%
% Note also that the binary data provided for use with this test script are
% to be used as an example test case only. This data set is NOT one of
% those described in the publication listed above and MAY NOT be used for
% any further research purposes without explicit permission from the
% authors.

%%% Load in the data and example lumped element values.
load S-Params-allports.mat % <- Download this example set from the releases page
load bin/LumpedValsTest


%% Generate combined S-parameters and plot them
S_reduced = gen_sparam_response(S,Freqs,8,120,lumpedVals);


%% Optimization
% If good impedance values are unknown, an optimization may be run to
% optimize matching and decoupling at a sepcific frequency. Details are
% given in DOI:10.1002/mrm.25504


% Load some initial capacitance values and indices of relevant ports (this
% latter variable allows multiple ports to be treated as having the same
% capacitance such that the user doesn't need to optimise every port
% independently)
load bin/OptInitVals

% Perform optimization
lambda = 2; %<- sets trade off between decoupling and matching
[lumpOpt,Sopt] = circ_co_sim_optimisation(S,Freqs,8,120,capValsInit,portCapInds,lambda);

