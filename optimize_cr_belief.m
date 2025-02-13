function [p] = optimize_cr_belief(belief, precision)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       function: optimize_cr_belief                      %
%           author: Federico Chiariotti (chiariot@dei.unipd.it)           %
%                             license: GPLv3                              %
%                                                                         %
%                                                                         %
%                                                                         %
% Optimize CR transmission probability through bisection search using the %
% belief PMF over the number of colliders                                 %
%                                                                         %
% Inputs:                                                                 %
% -belief:          the collider number belief PMF [1 x N]                %
% -precision:       the required precision on the output [scalar]         %
%                                                                         %
% Outputs:                                                                %
% -p:               the transmission probability [scalar]                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


p_low = 0;
p_high = 1;

N = length(belief);

if (N == 1)
    p = 1;
    return;
end

% Bisection search over the expected resolution time
while (p_low < p_high - precision)
    p = (p_high + p_low) / 2;
    val = belief(1) / p ^ 2;
    for c = 2 : N
        val = val + belief(c) * (1 - c * p) / (c * p ^ 2 * (1 - p) ^ c);
    end
    if (val > 0)
        p_low = p;
    else
        p_high = p;
    end
end


end