function [belief] = update_belief(belief, outcome, p_tx, epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         function: update_belief                         %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Update nodes' beliefs over the number of colliders in a CR/CE cycle     %
%                                                                         %
% Inputs:                                                                 %
% -belief:          the previous belief PMF [1 x N]                       %
% -outcome:         0 for ACK, 1 for NACK, 2 for silence [scalar]         %
% -p_tx:            the transmission probability for colliders [scalar]   %
% -epsilon:         the channel error probability [scalar]                %
%                                                                         %
% Outputs:                                                                %
% -belief:          the updated belief distribution PMF [1 x N]           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(belief);

% ACK
if (outcome == 0)
    % Compute the probability of an ACK if there are c colliders
    p_ack = zeros(1, N);
    for c = 1 : N
        p_ack(c) = c * p_tx * (1 - epsilon) * (1 - p_tx) ^ (c - 1);
    end
    overall_ack = belief * p_ack';
    % Compute the new belief
    new_belief = zeros(1, N);
    % We account for the case in which the last node transmits
    new_belief(1) = belief(1) * p_ack(1) / overall_ack;
    for c = 1 : N - 1
        % The ACK removes a collider
        new_belief(c) = new_belief(c) + belief(c + 1) * p_ack(c + 1) / overall_ack;
    end
end

% NACK
if (outcome == 1)
    % Compute the probability of a NACK if there are c colliders
    p_nack = zeros(1, N);
    for c = 1 : N
        p_nack(c) = 1 - (1 - p_tx) ^ c - c * p_tx * (1 - epsilon) * (1 - p_tx) ^ (c - 1);
    end
    overall_nack = belief * p_nack';
    % Compute the new belief
    new_belief = zeros(1, N);
    for c = 1 : N
        new_belief(c) = belief(c) * p_nack(c) / overall_nack;
    end
end

% Silence
if (outcome == 2)
    % Compute the probability of a silent slot if there are c colliders
    p_sil = zeros(1, N);
    for c = 1 : N
        p_sil(c) = (1 - p_tx) ^ c;
    end
    overall_sil = belief * p_sil';
    % Compute the new belief
    new_belief = zeros(1, N);
    for c = 1 : N
        new_belief(c) = belief(c) * p_sil(c) / overall_sil;
    end
end

belief = new_belief / sum(new_belief);
