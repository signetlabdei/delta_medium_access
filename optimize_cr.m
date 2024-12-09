function [p] = optimize_cr(lambda, epsilon, N, precision)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          function: optimize_cr                          %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Optimize CR transmission probability through bisection search           %
%                                                                         %
% Inputs:                                                                 %
% -lambda:          the generation rate [scalar]                          %
% -epsilon:         the wireless channel error probability [scalar]       %
% -N:               the number of nodes [scalar]                          %
% -precision:       the required precision on the output [scalar]         %
%                                                                         %
% Outputs:                                                                %
% -p:               the transmission probability [scalar]                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


p_low = 0;
p_high = 1;

if (N == 1)
    p = 1;
    return;
end

% Bisection search over the expected resolution time
while (p_low < p_high - precision)
    p = (p_high + p_low) / 2;
    val = N * lambda * (1- lambda) ^ (N - 1) * epsilon / p;
    for c = 2 : N
        if (N > 50)
            prob = binopdf(c, N, lambda);
        else
            prob = nchoosek(N, c) * lambda ^ c * (1- lambda) ^ (N - c);
        end
        val = val + prob * (1 - c * p) / (c * p ^ 2 * (1 - p) ^ c);
    end
    if (val > 0)
        p_low = p;
    else
        p_high = p;
    end
end


end