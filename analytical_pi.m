function [pi] = analytical_pi(N, lambda, epsilon, delta, M, optimistic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         function: analytical_pi                         %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Computes the steady-state distribution of the semi-Markov model         %
%                                                                         %
% Inputs:                                                                 %
% -N:       the number of nodes [scalar]                                  %
% -lambda:  the generation rate for all nodes [scalar]                    %
% -epsilon: the wireless channel error probability [scalar]               %
% -delta:   the number of slots resolved in each BT round [scalar]        %
% -M:       the maximum AoII [scalar]                                     %
%                                                                         %
% Outputs:                                                                %
% -pi: the steady-state distribution [1 x N]                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auxiliary vectors and parameters
K = delta + 1;
ps = zeros(K + 1, N);
xi = zeros(M + 1, N);
nu = zeros(1, N);

for k = 1 : K + 1
    % Optimize transmission probabilities
    for n = 1 : N
        ps(k, n) = optimize_cr(1 - (1 - lambda) ^ k, epsilon, N - n + 1, 0.00001);
    end
    % Compute collision probabilities
    active = 1 - (1 - lambda) ^ k;
    for n = 2 : N
        xi(k, n) = 1 - (1 - active) ^ n - n * active * (1 - active) ^ (n - 1);
    end
    nu(1) = delta;
    for n = 2 : N
        nu(n) = floor(delta / n);
    end
end

% Compute collision resolution transition matrix
Pcs = {};
for k = 1 : K + 1
    for c = 2 : N
        Pc = zeros(c);
        Pc(c, c) = 1;
        for i = 1 : c - 1
            Pc(i, i + 1) = (1 - epsilon) * (c - i + 1) * ps(k, i) * (1 - ps(k, i)) ^ (c - i);
            Pc(i, i) = 1 - Pc(i, i + 1);
        end
        Pcs{k, c} = Pc;
    end
end

% Compute distribution of collision resolution time
zetas = zeros(K + 1, M);
for k = 1 : K + 1
    active = 1 - (1 - lambda) ^ k;
    coeff = 1 - (1 - active) ^ N - (1 - epsilon) * N * (1 - active) ^ (N - 1) * active;
    % Singleton collision
    prob = epsilon * N * (1 - active) ^ (N - 1) * active / coeff;
    for m = 1 : M - 1
        zetas(k, m + 1) = prob * (1 - (1 - (1 - epsilon) * ps(k, 1)) ^ m);
    end
    for c = 2 : N
        prob = nchoosek(N, c) * (1 - active) ^ (N - c) * active ^ c / coeff;
        for m = 1 : M - c + 1
            Pc = (Pcs{k, c}) ^ m;
            p_cr = prob * Pc(1, c);
            zetas(k, m + c - 1) = zetas(k, m + c - 1) + p_cr * (1 - epsilon);
            for ell = 1 : M - m - c
                p_geo = (1 - (1 - epsilon) * ps(k, c)) ^ (ell - 1) * (1 - epsilon) * ps(k, c);
                zetas(k, m + c + ell) = zetas(k, m + c + ell) + p_cr * epsilon * p_geo; 
            end
        end
    end
end
for k = 1 : K + 1
    zetas(k, 2 : M) = diff(zetas(k, :));
    zetas(k, M) = 1 - sum(zetas(k, 1 : M - 1));
end

% Compute number of contending nodes (optimistic estimate)
eta = zeros(1, M);
if (optimistic == 1)
    etas = ones(M, N);
    for c = 1 : N
        prob = nchoosek(N, c) * lambda ^ c * (1 - lambda) ^ (N - c);
        if (c == 1)
            prob = prob * epsilon;
        end
        sigma = (1 - epsilon) * c * ps(1, 1) * (1 - ps(1, 1)) ^ (c - 1);
        for m = 2 : M
            pm = sigma * (1 - sigma) ^ (m - 1);
            etas(m, c) = prob * pm / zetas(1, m);
        end
    end
    for m = 2 : M
        etas(m, :) = etas(m, :) / sum(etas(m, :));
        eta(m) = sum(etas(m, :) .* (1 : N));
    end
end

% Transition and soujourn matrices
P = zeros(2 * M + 2);
T = ones(2 * M + 2);
% Transition from state ZW
P(1, M + 2) = 1;
T(1, M + 2) = min(1 / xi(1, N - floor(eta(M))), 1e9);

% Transitions from state BT(i)
for i = 1 : M
    Kp = nu(N - floor(eta(i)));
    j = max(0, i - Kp + 1);
    if (Kp == 0)
        Kp = 1;
    end
    if (j > M)
        j = M;
    end
    % Failure
    P(1 + i, j + M + 2) = xi(Kp, N - floor(eta(i)));
    P(i + 1, j + 1) = 1 - xi(Kp, N - floor(eta(i)));
end

% Transitions from state CR(j)
for i = 1 : M - 1
    for j = 0 : i - 1
        Kp = 1;
        if (j > 0)
            Kp = nu(N - floor(eta(j)));
        end
        if (Kp == 0)
            Kp = 1;
        end
        P(j + M + 2, i + 1) = zetas(min(Kp, j + 1), i - j);
        T(j + M + 2, i + 1) = i - j;
    end
end
for j = 0 : M
    P(j + M + 2, M + 1) = 1 - sum(P(j + M + 2, 1 : 2 * M + 1));
end

% Compute steady-state distribution
[~, D, W] = eig(P);
[~, m_idx] = min(abs(1 - diag(D)));
alphas = W(:, m_idx) / sum(W(:, m_idx));
pi = zeros(2 * M + 2, 1);
for s = 1 : 2 * M + 2
    for sn = 1 : 2 * M + 2
        g = alphas(s) * T(s, sn) * P(s, sn);
        pi(s) = pi(s) + alphas(s) * T(s, sn) * P(s, sn);
    end
end
pi = pi / sum(pi);

end