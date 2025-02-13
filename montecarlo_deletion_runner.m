%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   script: montecarlo_deletion_runner                   %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Runs a Monte Carlo simulation with the desired protocol                 %
%                                                                         %
% Parameters:                                                             %
% -L:           the number of steps to simulate per episode [scalar]      %
% -N:           the number of nodes [scalar]                              %
% -E:           the number of episodes to simulate [scalar]               %
% -lambda:      the generation rate for all nodes [scalar]                %
% -epsilon:     the wireless channel error probability [scalar]           %
% -algo:        the number of slots resolved in each BT round [scalar]    %
% -K:           number of cleared slots in BT [scalar]                    %
% -p1:          alpha for ZW/GZW/LZW [scalar]                             %
% -p2:          beta for GZW/LZW [scalar]                                 %
% -M:           the maximum AoII [scalar]                                 %
% -omegas:      the feedback error probabilities                          %
% -setting:     0 for simulating over sigma, 1 for simulating over nu     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clearvars

%%% PARAMETERS
N = 20;
L = 1e5 + 1000;
E = 20;
K = 50;
rho = 0.5;
epsilon = 0.05;
algo = 'delta+';
M = 100;
p1 = 0.29;
p2 = 0.13;

omegas = 0 : 0.02 : 0.2;
   
% Auxiliary vectors for CDFs
aois = zeros(length(omegas), M + 1);
max_aois = zeros(length(omegas), M + 1);
aoiis = zeros(length(omegas), M + 1);

for is = 1 : length(omegas)
    omega = omegas(is)
    for e = 1 : E
        % Run Monte Carlo simulation and compute CDF (skipping first 1000 steps)
        lambdas = ones(1, N) * rho / N;
        [aoi, aoii] = montecarlo_deletion(L, N, lambdas, epsilon, omega, algo, K, p1, p2);
        aoi = aoi(:, 1001 : L);
        aoii = aoii(:, 1001 : L);
        [aoi_dist, ~] = hist(aoi(:), 0 : M);
        [max_aoi_dist, ~] = hist(max(aoi, [], 1), 0 : M);
        [aoii_dist, ~] = hist(aoii(:), 0 : M);
        aoi_cdf = cumsum(aoi_dist) / sum(aoi_dist);                
        max_aoi_cdf = cumsum(max_aoi_dist) / sum(max_aoi_dist);                
        aoii_cdf = cumsum(aoii_dist) / sum(aoi_dist);
        aois(is, :) = aois(is, :) + aoi_cdf / E;
        max_aois(is, :) = max_aois(is, :) + max_aoi_cdf / E;
        aoiis(is, :) = aoiis(is, :) + aoii_cdf / E;
    end
end