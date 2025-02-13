%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           script: optimize_K                            %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Runs a Monte Carlo simulation with the desired protocol                 %
%                                                                         %
% Parameters:                                                             %
% -epsilon:     the wireless channel error probability [scalar]           %
% -M:           the maximum AoII [scalar]                                 %
% -nodes:       0 for simulating over lambda, 1 for simulating over N     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clearvars

%%% PARAMETERS
Ks = 4 : 2 : 120;
epsilon = 0.05;
M = 100;
nodes = 1;
if (nodes == 1)
    Ns = 4 : 2 : 50;
    lambdas = 0.015;
else
    lambdas = 0.001 : 0.001 : 0.04;
    Ns = 20;
end



% Auxiliary variables
best_optimistic = zeros(length(Ns), length(lambdas));
best_optimistic_prob = zeros(length(Ns), length(lambdas));
best_pessimistic = zeros(length(Ns), length(lambdas));
best_pessimistic_prob = zeros(length(Ns), length(lambdas));
optimistic_probs = zeros(length(Ns), length(lambdas), length(Ks));
pessimistic_probs = zeros(length(Ns), length(lambdas), length(Ks));

% Iterate over parameters
for in = 1 : length(Ns)
    N = Ns(in)
    for il = 1 : length(lambdas)
        lambda = lambdas(il) / N * 20;
        for id = 1 : length(Ks)
            K = Ks(id);
            if (K < 10 * N && K > N)
                % Optimistic evaluation
                pi_opt = analytical_pi(N, lambda, epsilon, K, M, 1);
                % Pessimistic evaluation
                pi_pess = analytical_pi(N, lambda, epsilon, K, M, 0);
                optimistic_probs(in, il, id) = pi_opt(1);
                pessimistic_probs(in, il, id) = pi_pess(1);
                % Save improved values
                if (pi_opt(1) > best_optimistic_prob(in, il))
                    best_optimistic_prob(in, il) = pi_opt(1);
                    best_optimistic(in, il) = K;
                end
                if (pi_pess(1) > best_pessimistic_prob(in, il))
                    best_pessimistic_prob(in, il) = pi_pess(1);
                    best_pessimistic(in, il) = K;
                end
            end
        end
    end
end

best_optimistic = squeeze(best_optimistic);
best_optimistic_prob = squeeze(best_optimistic_prob);
best_pessimistic = squeeze(best_pessimistic);
best_pessimistic_prob = squeeze(best_pessimistic_prob);

figure;
hold on
plot(best_optimistic)
plot(best_pessimistic)

figure;
hold on
plot(best_optimistic_prob)
plot(best_pessimistic_prob)