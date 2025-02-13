function [psi, theta] = montecarlo_noisy(L, N, lambda, epsilon, sigma, tx_method, K, p1, p2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       function: montecarlo_noisy                        %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Runs a Monte Carlo simulation with the desired protocol, using a noisy  %
% feedback channel model.                                                 %
%                                                                         %
% Inputs:                                                                 %
% -L:               the number of steps to simulate [scalar]              %
% -N:               the number of nodes [scalar]                          %
% -lambda:          the generation rate for each node [1 x N]             %
% -epsilon:         the wireless channel error probability [scalar]       %
% -sigma:           the ACK error variance [scalar]                       %
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
collider_belief = zeros(N);

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
            if (theta_belief(n) > 0 || state(n) == 2)
                theta_belief(n) = theta_belief(n) + 1;
            end
        end
    end
    % Simulate feedback channel
    ack_error = randn(1, N) * sigma;
    % Transmission
    tx_ind = [];
    if (strcmp(tx_method, 'zero_wait'))
        % Run ZW algorithm
        tx = theta(:, l) > 0;
        tx = tx .* (rand(N, 1) < p1);
        tx_ind = find(tx);
    end
    if (strcmp(tx_method, 'zero_wait_local'))
        % Run LZW algorithm
        tx = theta(:, l) > 0;
        for n = 1 : N
            p = p1;
            if(colliders(n) > 0)
                p = p2;
            end
            tx(n) = tx(n) * (rand < p);
        end
        tx_ind = find(tx);
    end
    if (strcmp(tx_method, 'zero_wait_global'))
        % Run GZW algorithm
        tx = theta(:, l) > 0;
        for n = 1 : N
            p = p1;
            if (coll == 2)
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
            if (coll == 0)
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
                if (coll == 1)
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
            if (coll == 0)
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
        end
    end
    if (strcmp(tx_method, 'delta+'))
        % Run DELTA
        tx = zeros(1, N);
        for n = 1 : N
            if (coll == 0)
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
                if (coll == 1)
                    %CE phase
                    tx(n) = colliders(n);
                else
                    % CR phase
                    tx_probs(n) = optimize_cr_belief(collider_belief(n, :), 0.001);
                    tx(n) = colliders(n) * (rand < tx_probs(n));
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
            if (coll == 0)
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
        end
    end
    if (strcmp(tx_method, 'max_age'))
        % Run MAF algorithm
        [~, legit] = max(psi(:, l));
        tx_ind = [];
        for n = 1 : N
            ack_id = min(max(legit + round(ack_error(n)), 1), N);
            if (ack_id == n)
                tx_ind = [tx_ind, n];
            end
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
    % No transmission in CR: update beliefs
    if (coll == 2 && isempty(tx_ind) && strcmp(tx_method, 'delta+'))
        for n = 1 : N
            collider_belief(n, :) = update_belief(collider_belief(n, :), 2, tx_probs(n), epsilon);
        end
    end

    % One node transmits
    if (isscalar(tx_ind))
        if (rand > epsilon)
            % ACKed tx: the node resets its state and real and believed AoII
            state(tx_ind) = 1;
            psi(tx_ind, l) = 0;
            theta_belief(tx_ind) = 0;
            psi_belief(tx_ind, tx_ind) = 0;
            theta(tx_ind, l) = 0;
            coll = max(0, coll - 1);
            colliders(tx_ind) = 0;
            coll_sequence = 0;
            % Check ACK for other nodes
            for n = 1 : N
                if (n ~= tx_ind)
                    ack_id = min(max(tx_ind + round(ack_error(n)), 1), N);
                    if (ack_id ~= n)
                        psi_belief(n, ack_id) = 0;
                    end
                end
            end
            if (coll == 2)
                % Update belief over number of colliders
                if (strcmp(tx_method, 'delta+'))
                    for n = 1 : N
                        collider_belief(n, :) = update_belief(collider_belief(n, :), 0, tx_probs(n), epsilon);
                    end
                end
            end
        else
            % Wireless channel error: NACK
            colliders(tx_ind) = 1;
            if (coll == 0)
                % Initialize belief over number of colliders
                if (strcmp(tx_method, 'delta+'))
                    for n = 1 : N
                        % Initialize belief
                        if (max(psi_belief(n, :)) == 1)
                            collider_belief(n, 1) = epsilon * N * mean(lambda) * (1 - mean(lambda)) ^ (N - 1);
                            for c = 2 : N
                                collider_belief(n, c) = nchoosek(N, c) * mean(lambda) ^ c * (1 - mean(lambda)) ^ (N - c);
                            end
                            collider_belief(n, :) = collider_belief(n, :) / sum(collider_belief(n, :));
                        else
                            % Count possible collision steps
                            maxage = max(psi_belief(n, :));
                            new_age = psi_belief(n, :);
                            while (sum(psi_belief(n, :)) - sum(new_age) < K)
                                maxage = maxage - 1;
                                new_age = min(psi_belief(n, :), maxage);
                            end
                            collision_steps = psi_belief(n, :) - new_age;
                            % Evaluate transmission probability in BT
                            active = 1 - (1 - mean(lambda)) ^ (max(collision_steps));
                            Na = length(find(psi_belief(n, :) > max(psi_belief(n, :)) - collision_steps));
                            collider_belief(n, 1) = epsilon * Na * mean(active) * (1 - mean(active)) ^ (N - 1);
                            for c = 2 : Na
                                collider_belief(n, c) = nchoosek(Na, c) * mean(active) ^ c * (1 - mean(active)) ^ (Na - c);
                            end
                            collider_belief(n, :) = collider_belief(n, :) / sum(collider_belief(n, :));
                        end
                    end
                end
            end
            if (coll == 1)
                coll_sequence = coll_sequence + 1;
                % Update belief over number of colliders
                if (strcmp(tx_method, 'delta+'))
                    for n = 1 : N
                        collider_belief(n, :) = update_belief(collider_belief(n, :), 1, 1, epsilon);
                    end
                end
            end
            if (coll == 2)
                % Update belief over number of colliders
                if (strcmp(tx_method, 'delta+'))
                    for n = 1 : N
                        collider_belief(n, :) = update_belief(collider_belief(n, :), 1, tx_probs(n), epsilon);
                    end
                end
            end
            coll = 2;
        end
    end
    % Collision
    if (length(tx_ind) > 1)        
        colliders(tx_ind) = 1;
        if (coll == 0)
            % Initialize belief over number of colliders
            if (strcmp(tx_method, 'delta+'))
                for n = 1 : N
                    % Initialize belief
                    if (max(psi_belief(n, :)) == 1)
                        collider_belief(n, 1) = epsilon * N * mean(lambda) * (1 - mean(lambda)) ^ (N - 1);
                        for c = 2 : N
                            collider_belief(n, c) = nchoosek(N, c) * mean(lambda) ^ c * (1 - mean(lambda)) ^ (N - c);
                        end
                        collider_belief(n, :) = collider_belief(n, :) / sum(collider_belief(n, :));
                    else
                        % Count possible collision steps
                        maxage = max(psi_belief(n, :));
                        new_age = psi_belief(n, :);
                        while (sum(psi_belief(n, :)) - sum(new_age) < K)
                            maxage = maxage - 1;
                            new_age = min(psi_belief(n, :), maxage);
                        end
                        collision_steps = psi_belief(n, :) - new_age;
                        % Evaluate transmission probability in BT
                        active = 1 - (1 - mean(lambda)) ^ (max(collision_steps));
                        Na = length(find(psi_belief(n, :) > max(psi_belief(n, :)) - collision_steps));
                        collider_belief(n, 1) = epsilon * Na * mean(active) * (1 - mean(active)) ^ (N - 1);
                        for c = 2 : Na
                            collider_belief(n, c) = nchoosek(Na, c) * mean(active) ^ c * (1 - mean(active)) ^ (Na - c);
                        end
                        collider_belief(n, :) = collider_belief(n, :) / sum(collider_belief(n, :));
                    end
                end
            end
        end
        if (coll == 1)
            coll_sequence = coll_sequence + 1;
            % Update belief over number of colliders
            if (strcmp(tx_method, 'delta+'))
                for n = 1 : N
                    collider_belief(n, :) = update_belief(collider_belief(n, :), 1, 1, epsilon);
                end
            end
        end
        if (coll == 2)
            % Update belief over number of colliders
            if (strcmp(tx_method, 'delta+'))
                for n = 1 : N
                    collider_belief(n, :) = update_belief(collider_belief(n, :), 1, tx_probs(n), epsilon);
                end
            end
        end
        coll = 2;
    end
    % Update system state (anomaly generation)
    state = min(2, state + (rand(1, N) < lambda));
end

end