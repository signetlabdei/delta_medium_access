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
% -setting:     0 for simulating over rho, 1 over N, 2 over nu            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clearvars

%%% PARAMETERS
epsilons = 0.05;
algo = 'delta';
M = 100;
L = 5e6 + 1000;
Kf = 5 / 2;
p1 = 0.17;
p2 = 0.13;

setting = 1;

% Iterate over rho
if (setting == 0)
    rhos = 0.02 : 0.02 : 0.6;
    nus = 0;
    Ns = 20;
end

% Iterate over N
if (setting == 1)
    Ns = 4 : 2 : 50;
    nus = zeros(1, length(Ns));
    rhos = 0.5;
end

% Iterate over nu
if (setting == 2)
    nus = 0 : 0.05 : 0.5;
    Ns = ones(1, length(nus)) * 20;
    rhos = 0.5;
end

% Auxiliary vectors for CDFs
aois = zeros(length(Ns), length(rhos), length(epsilons), M + 1);
max_aois = zeros(length(Ns), length(rhos), length(epsilons), M + 1);
aoiis = zeros(length(Ns), length(rhos), length(epsilons), M + 1);

% Loop over parameters
for in = 1 : length(Ns)
    N = Ns(in)
    nu = nus(in)
    for il = 1 : length(rhos)
        lambda = rhos(il) / N
        K = N * Kf;
        K = best_pessimistic(in);
        % p1 = p1s(in);
        % p2 = p2s(in);
        for ie = 1 : length(epsilons)
            epsilon = epsilons(ie);
                % Run Monte Carlo and compute CDFs (skipping first 1000 steps)
                lambdas = (ones(1, N)  + nu * 2 * (rand(1, N) - 0.5)) * lambda;
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