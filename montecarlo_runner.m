%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        script: montecarlo_runner                        %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Runs a Monte Carlo simulation with the desired protocol                 %
%                                                                         %
% Parameters:                                                             %
% -L:           the number of steps to simulate [scalar]                  %
% -N:           the number of nodes [scalar]                              %
% -lambda:      the generation rate for all nodes [scalar]                %
% -epsilon:     the wireless channel error probability [scalar]           %
% -algo:        the number of slots resolved in each BT round [scalar]    %
% -K:           number of cleared slots in BT [scalar]                    %
% -p1:          alpha for ZW/GZW/LZW [scalar]                             %
% -p2:          beta for GZW/LZW [scalar]                                 %
% -M:           the maximum AoII [scalar]                                 %
% -nodes:       0 for simulating over lambda, 1 for simulating over N     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clearvars

%%% PARAMETERS
epsilons = 0.05;
algo = 'delta';
M = 100;
L = 1e5 + 1000;
K = 50;
p1 = 1;
p2 = 0.2;

nodes = 0;
if (nodes == 1)
    Ns = 4 : 2 : 50;
    lambdas = 0.025;
else
    lambdas = 0.001 : 0.001 : 0.04;
    Ns = 20;
end

% Auxiliary vectors for CDFs
aois = zeros(length(Ns), length(lambdas), length(epsilons), M + 1);
max_aois = zeros(length(Ns), length(lambdas), length(epsilons), M + 1);
aoiis = zeros(length(Ns), length(lambdas), length(epsilons), M + 1);

% Loop over parameters
for in = 1 : length(Ns)
    N = Ns(in);
    for il = 1 : length(lambdas)
        lambda = lambdas(il) / N * 20
        for ie = 1 : length(epsilons)
            epsilon = epsilons(ie);
                % Run Monte Carlo and compute CDFs (skipping first 1000 steps)
                [aoi, aoii] = montecarlo(L, N, ones(1, N) * lambda, epsilon, algo, K, p1, p2);
                aoi = aoi(:, 1001 : L);
                aoii = aoii(:, 1001 : L);
                [aoi_dist, ~] = hist(aoi(:), 0 : M);
                [max_aoi_dist, ~] = hist(max(aoi, [], 1), 0 : M);
                [aoii_dist, ~] = hist(aoii(:), 0 : M);
                aoi_cdf = cumsum(aoi_dist) / sum(aoi_dist);                
                max_aoi_cdf = cumsum(max_aoi_dist) / sum(max_aoi_dist);                
                aoii_cdf = cumsum(aoii_dist) / sum(aoi_dist);
                aois(in, il, ie, :) = aoi_cdf;
                max_aois(in, il, ie, :) = max_aoi_cdf;
                aoiis(in, il, ie, :) = aoii_cdf;
        end
    end
end

aois = squeeze(aois);
max_aois = squeeze(max_aois);
aoiis = squeeze(aoiis);