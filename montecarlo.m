function [psi, theta] = montecarlo(L, N, lambda, epsilon, tx_method, K, p1, p2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          function: montecarlo                           %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Runs a Monte Carlo simulation with the desired protocol                 %
%                                                                         %
% Inputs:                                                                 %
% -L:           the number of steps to simulate [scalar]                  %
% -N:           the number of nodes [scalar]                              %
% -lambda:      the generation rate for each node [1 x N]                 %
% -epsilon:     the wireless channel error probability [scalar]           %
% -tx_method:   the selected protocol [string]                            %
% -K:           number of cleared slots in BT [scalar]                    %
% -p1:          alpha for ZW/GZW/LZW [scalar]                             %
% -p2:          beta for GZW/LZW [scalar]                                 %
%                                                                         %
% Outputs:                                                                %
% -psi:         the maximum AoII for all nodes, step by step [N x L]      %
% -theta:       the real AoII for all nodes, step by step [N x L]         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize auxiliary variables
psi = zeros(N, L);
theta = zeros(N, L);
state = ones(1, N);
coll = 0;
coll_sequence = 0;
colliders = zeros(1, N);
collider_belief = zeros(1, N);
collision_steps = zeros(1, N);
tx_prob = 0;

% Compute threshold
threshold = exp(K * log(1 - mean(lambda)));

for l = 1 : L
    % Update public max AoII and private AoII
    if (l > 1)
        psi(:, l) = psi(:, l - 1) + 1;
        theta(:, l) = theta(:, l - 1) + (state' == 2);
    end
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
        p = p1;
        if (coll == 2)
            p = p2;
        end
        tx = tx .* (rand(N, 1) < p);
        tx_ind = find(tx);
    end
    if (strcmp(tx_method, 'delta'))
        % Run DELTA
        if (l < 100)
            % First 100 steps: random transmission (hot start)
            tx = rand(1, N) < 1 / N;
            if (sum(tx) == 1)
                tx_ind = find(tx, 1);
            end
        else
            if (coll == 0)
                if (max(psi(:, l)) == 1)
                    % ZW phase
                    collision_steps = ones(1, N);
                    tx = state - 1;
                else
                    % BT phase
                    p_tx = zeros(N, 1);
                    for n = 1 : N
                        if (theta(n, l) > 0)
                            p_tx(n) = 1;
                            % Compute belief
                            for j = 1 : N
                                if (j ~= n && psi(j, l) >= theta(n, l))
                                    p_tx(n) = p_tx(n) * (1 - lambda(j)) ^ (psi(j, l) - theta(n, l) + 1);
                                end
                            end
                        end
                    end
                    if (sum(psi(:, l)) <= K)
                        % Reset to ZW phase
                        collision_steps = ones(1, N) * floor(K / N);
                    else
                        % Count possible collision steps
                        maxage = max(psi(:, l));
                        new_age = psi(:, l);
                        while (sum(psi(:, l)) - sum(new_age) < K)
                            maxage = maxage - 1;
                            new_age = min(psi(:, l), maxage);
                        end
                        collision_steps = psi(:, l) - new_age;
                    end
                    tx = p_tx > threshold;
                end
            else
                if (coll == 1)
                    %CE phase
                    tx = colliders';
                else
                    % CR phase
                    active = 1 - (1 - mean(lambda)) ^ (max(collision_steps));
                    tx_prob = optimize_cr(active, epsilon, N - coll_sequence, 0.0001);
                    tx = colliders' .* (rand(N, 1) < tx_prob);
                end
            end
            tx_ind = find(tx);
            % Correct maximum AoII
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
        end
    end
    if (strcmp(tx_method, 'delta+'))
        % Run DELTA+
        if (l < 100)
            % First 100 steps: random transmission (hot start)
            tx = rand(1, N) < 1 / N;
            if (sum(tx) == 1)
                tx_ind = find(tx, 1);
            end
        else
            if (coll == 0)
                if (max(psi(:, l)) == 1)
                    % ZW phase
                    collision_steps = ones(1, N);
                    tx = state - 1;
                else
                    % BT phase
                    p_tx = zeros(N, 1);
                    for n = 1 : N
                        if (theta(n, l) > 0)
                            p_tx(n) = 1;
                            % Compute belief
                            for j = 1 : N
                                if (j ~= n && psi(j, l) >= theta(n, l))
                                    p_tx(n) = p_tx(n) * (1 - lambda(j)) ^ (psi(j, l) - theta(n, l) + 1);
                                end
                            end
                        end
                    end
                    if (sum(psi(:, l)) <= K)
                        % Reset to ZW phase
                        collision_steps = ones(1, N) * floor(K / N);
                    else
                        % Count possible collision steps
                        maxage = max(psi(:, l));
                        new_age = psi(:, l);
                        while (sum(psi(:, l)) - sum(new_age) < K)
                            maxage = maxage - 1;
                            new_age = min(psi(:, l), maxage);
                        end
                        collision_steps = psi(:, l) - new_age;
                    end
                    tx = p_tx > threshold;
                end
            else
                if (coll == 1)
                    %CE phase
                    tx = colliders';
                else
                    % CR phase
                    tx_prob = optimize_cr_belief(collider_belief, 0.0001);
                    tx = colliders' .* (rand(N, 1) < tx_prob);
                end
            end
            tx_ind = find(tx);
            % Correct maximum AoII
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
        end
    end
    if (strcmp(tx_method, 'max_age'))
        % Run MAF algorithm
        [~, tx] = max(psi(:, l));
        tx_ind = tx(1);
    end
    if (strcmp(tx_method, 'round_robin'))
        % Run RR algorithm
        tx_ind = 1 + mod(l, N);
    end
    % CE successful: no transmission
    if (coll == 1 && isempty(tx_ind))
        coll = 0;
        coll_sequence = 0;
    end
    % CR: no transmission
    if (coll == 2 && isempty(tx_ind))
        % Update belief over number of colliders
        if (strcmp(tx_method, 'delta+'))
            collider_belief = update_belief(collider_belief, 2, tx_prob, epsilon);
        end
    end

    % Only one transmitter
    if (isscalar(tx_ind))
        if (rand > epsilon)
            % Successful transmission
            state(tx_ind) = 1;
            psi(tx_ind, l) = 0;
            theta(tx_ind, l) = 0;
            coll = max(0, coll - 1);
            colliders(tx_ind) = 0;
            coll_sequence = 0;
            if (coll == 2)
                % Update belief over number of colliders
                if (strcmp(tx_method, 'delta+'))
                    collider_belief = update_belief(collider_belief, 0, tx_prob, epsilon);
                end
            end
        else
            % Wireless channel error: NACK
            colliders(tx_ind) = 1;
            if (coll == 0)
                % Initialize belief over number of colliders
                if (strcmp(tx_method, 'delta+'))
                    if (max(psi(:, l)) == 1)
                        collider_belief(1) = epsilon * N * mean(lambda) * (1 - mean(lambda)) ^ (N - 1);
                        for c = 2 : N
                            collider_belief(c) = nchoosek(N, c) * mean(lambda) ^ c * (1 - mean(lambda)) ^ (N - c);
                        end
                        collider_belief = collider_belief / sum(collider_belief);
                    else
                        % Count possible collision steps
                        maxage = max(psi(:, l));
                        new_age = psi(:, l);
                        while (sum(psi(:, l)) - sum(new_age) < K)
                            maxage = maxage - 1;
                            new_age = min(psi(:, l), maxage);
                        end
                        collision_steps = psi(:, l) - new_age;
                        % Evaluate transmission probability in BT
                        active = 1 - (1 - mean(lambda)) ^ (max(collision_steps));
                        Na = length(find(psi(:, l) > max(psi(:, l)) - collision_steps));
                        collider_belief(1) = epsilon * Na * mean(active) * (1 - mean(active)) ^ (N - 1);
                        for c = 2 : Na
                            collider_belief(c) = nchoosek(Na, c) * mean(active) ^ c * (1 - mean(active)) ^ (Na - c);
                        end
                        collider_belief = collider_belief / sum(collider_belief);
                    end
                end
            end
            if (coll == 1)
                coll_sequence = coll_sequence + 1;
                % Update belief over number of colliders
                if (strcmp(tx_method, 'delta+'))
                    collider_belief = update_belief(collider_belief, 1, 1, epsilon);
                end
            end
            if (coll == 2)
                % Update belief over number of colliders
                if (strcmp(tx_method, 'delta+'))
                    collider_belief = update_belief(collider_belief, 1, tx_prob, epsilon);
                end
            end
            coll = 2;
        end
    end
    % More than one transmitter: NACK
    if (length(tx_ind) > 1)
        colliders(tx_ind) = 1;
        if (coll == 0)
            % Initialize belief over number of colliders
            if (strcmp(tx_method, 'delta+'))
                if (max(psi(:, l)) == 1)
                    collider_belief(1) = epsilon * N * mean(lambda) * (1 - mean(lambda)) ^ (N - 1);
                    for c = 2 : N
                        collider_belief(c) = nchoosek(N, c) * mean(lambda) ^ c * (1 - mean(lambda)) ^ (N - c);
                    end
                    collider_belief = collider_belief / sum(collider_belief);
                else
                    % Count possible collision steps
                    maxage = max(psi(:, l));
                    new_age = psi(:, l);
                    while (sum(psi(:, l)) - sum(new_age) < K)
                        maxage = maxage - 1;
                        new_age = min(psi(:, l), maxage);
                    end
                    collision_steps = psi(:, l) - new_age;
                    % Evaluate transmission probability in BT
                    active = 1 - (1 - mean(lambda)) ^ (max(collision_steps));
                    Na = length(find(psi(:, l) > max(psi(:, l)) - collision_steps));
                    collider_belief(1) = epsilon * Na * mean(active) * (1 - mean(active)) ^ (N - 1);
                    for c = 2 : Na
                        collider_belief(c) = nchoosek(Na, c) * mean(active) ^ c * (1 - mean(active)) ^ (Na - c);
                    end
                    collider_belief = collider_belief / sum(collider_belief);
                end
            end
        end
        if (coll == 1)
            coll_sequence = coll_sequence + 1;
            % Update belief over number of colliders
            if (strcmp(tx_method, 'delta+'))
                collider_belief = update_belief(collider_belief, 1, 1, epsilon);
            end
        end
        if (coll == 2)
            % Update belief over number of colliders
            if (strcmp(tx_method, 'delta+'))
                collider_belief = update_belief(collider_belief, 1, tx_prob, epsilon);
            end
        end
        coll = 2;
    end
    % Update system state (anomaly generation)
    state = min(2, state + (rand(1, N) < lambda));
end