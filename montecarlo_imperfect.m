function [psi, theta] = montecarlo_imperfect(L, N, lambda, epsilon, feedback_err, tx_method, K, p1, p2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     function: montecarlo_imperfect                      %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Runs a Monte Carlo simulation with the desired protocol                 %
%                                                                         %
% Inputs:                                                                 %
% -L:               the number of steps to simulate [scalar]              %
% -N:               the number of nodes [scalar]                          %
% -lambda:          the generation rate for each node [1 x N]             %
% -epsilon:         the wireless channel error probability [scalar]       %
% -feedback_err:    the feedback channel error probability [scalar]       %
% -tx_method:       the selected protocol [string]                        %
% -K:               number of cleared slots in BT [scalar]                %
% -p1:              alpha for ZW/GZW/LZW [scalar]                         %
% -p2:              beta for GZW/LZW [scalar]                             %
%                                                                         %
% Outputs:                                                                %
% -psi:             the maximum AoII for all nodes, step by step [N x L]  %
% -theta:           the real AoII for all nodes, step by step [N x L]     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize auxiliary variables
psi = zeros(N, L);
theta = zeros(N, L);
state = ones(1, N);
coll = 0;
coll_sequence = 0;
colliders = zeros(1, N);
psi_belief = zeros(N, N);
theta_belief = zeros(1, N);
phase_belief = ones(1, N);

% Compute threshold
threshold = exp(K * log(1 - mean(lambda)));

% Pre-compute CR transmission probabilities
tx_probs = ones(1, N);
for j = 1 : N - 1
    active = 1 - (1 - mean(lambda)) ^ (floor(K / N));
    tx_probs(j) = optimize_cr(active, epsilon, N - j + 1, 0.0001);
end

for l = 1 : L
    % Update beliefs on max AoII and private AoII
    if (l > 1)
        psi(:, l) = psi(:, l - 1) + 1;
        psi_belief = psi_belief + 1;
        theta(:, l) = theta(:, l - 1) + (state' == 2);
        for n = 1 : N
            if (theta_belief(n)) > 0 || state(n) == 2
                theta_belief(n) = theta_belief(n) + 1;
            end
        end
    end
    % Update system state (anomaly generation)
    state = min(2, state + (rand(1, N) < lambda));
    % Simulate feedback channel
    ack = rand(1, N) > feedback_err;
    % Transmission
    tx_ind = [];
    if (strcmp(tx_method, 'zero_wait'))
        % Run ZW algorithm
        tx = theta_belief > 0;
        tx = tx' .* (rand(N, 1) < p1);
        tx_ind = find(tx);
    end
    if (strcmp(tx_method, 'zero_wait_local'))
        % Run LZW algorithm
        tx = theta_belief > 0;
        for n = 1 : N
            p = p1;
            if(coll_belief(n) > 0)
                p = p2;
            end
            tx(n) = tx(n) * (rand < p);
        end
        tx_ind = find(tx);
    end
    if (strcmp(tx_method, 'zero_wait_global'))
        % Run GZW algorithm
        tx = theta_belief > 0;
        for n = 1 : N
            p = p1;
            if (coll_belief(n) == 1)
                p = p2;
            end
            tx(n) = tx(n) .* (rand < p);
        end
        tx_ind = find(tx);
    end
    if (strcmp(tx_method, 'delta'))
        % Run DELTA
        tx = zeros(1, N);
        for n = 1 : N
            if (phase_belief(n) == 1)
                if (max(psi_belief(n, :)) == 1)
                    % ZW phase
                    tx(n) = theta_belief(n) > 0;
                else
                    % BT phase
                    p_tx = 0;
                    % Compute belief that node has the highest AoII
                    if (theta_belief(n) > 0)
                        p_tx = 1;
                        for j = 1 : N
                            if (j ~= n && psi_belief(n, j) >= theta_belief(n))
                                p_tx = p_tx * (1 - lambda(j)) ^ (psi_belief(n, j) - theta_belief(n) + 1);
                            end
                        end
                    end
                    tx(n) = p_tx > threshold;
                end
            else
                if (phase_belief(n) == 3)
                    %CE phase
                    tx(n) = colliders(n);
                else
                    % CR phase
                    tx(n) = colliders(n) * (rand < tx_probs(coll_sequence + 1));
                end
            end
        end
        tx_ind = find(tx);
        % Correct psi (real)
        if (coll == 0)
            if (sum(psi(:, l)) <= K)
                psi(:, l) = zeros(N, 1);
            else
                maxage = max(psi(:, l));
                new_age = psi(:, l);
                while (sum(psi(:, l)) - sum(new_age) < K)
                    maxage = maxage - 1;
                    new_age = min(psi(:, l), maxage);
                end
                psi(:, l) = min(psi(:, l), maxage + 1);
            end
        end
        % Correct psi (believed)
        for n = 1 : N
            if (phase_belief(n) < 2)
                if (sum(psi_belief(n, :)) <= K)
                    psi_belief(n, :) = zeros(N, 1);
                else
                    maxage = max(psi_belief(n, :));
                    new_age = psi_belief(n, :);
                    while (sum(psi_belief(n, :)) - sum(new_age) < K)
                        maxage = maxage - 1;
                        new_age = min(psi_belief(n, :), maxage);
                    end
                    psi_belief(n, :) = min(psi_belief(n, :), maxage + 1);
                end
            end
            if (theta_belief(n) > psi_belief(n, n))
                theta_belief(n) = 0;
            end
            if (ack(n))
                psi_belief(n, :) = min(psi_belief(n, :), max(psi(:, l)));
            end
        end
    end
    if (strcmp(tx_method, 'max_age'))
        % Run MAF algorithm
        if (rand > feedback_err)
            [~, tx] = max(psi(:, l));
            tx_ind = tx(1);
        else
            tx_ind = [];
        end
    end
    if (strcmp(tx_method, 'round_robin'))
        % Run RR algorithm
        tx_ind = 1 + mod(l, N);
    end
    % No transmission: successful CE
    if (coll == 1 && isempty(tx_ind))
        coll = 0;
        coll_sequence = 0;
    end
    % One node transmits
    if (isscalar(tx_ind))
        if (rand > epsilon)
            if (ack (tx_ind))
                % ACKed tx: the node resets its AoII
                theta_belief(tx_ind) = 0;
            end
            % Reset real state and AoII
            state(tx_ind) = 1;
            psi(tx_ind, l) = 0;
            % Check ACK
            for n = 1 : N
                if (ack(n))                
                    coll_belief(n) = 0;
                    psi_belief(n, tx_ind) = 0;
                end
            end
            % Update real values
            theta(tx_ind, l) = 0;
            coll = max(0, coll - 1);
            colliders(tx_ind) = 0;
            coll_sequence = 0;
        else
            % Wireless channel error
            colliders(tx_ind) = 1;
            if (coll == 1)
                coll_sequence = coll_sequence + 1;
            end
            coll = 2;
        end
    end
    % Collision
    if (length(tx_ind) > 1)
        colliders(tx_ind) = 1;
        coll_belief = colliders;
        if (coll == 1)
            coll_sequence = coll_sequence + 1;
        end
        coll = 2;
    end
    if (strcmp(tx_method, 'delta'))
        % Phase belief: ACK from previous step
        for n = 1 : N
            if (ack(n))
                if (coll == 0)
                    phase_belief(n) = 1;
                else
                    if (coll == 2)
                        phase_belief(n) = 2;
                    else
                        phase_belief(n) = 3;
                    end
                end
            else
                if (phase_belief(n) == 3)
                    if (isempty(tx_ind))
                        phase_belief(n) = 1;
                    else
                        phase_belief(n) = 2;
                    end
                end
            end
        end
    end
end