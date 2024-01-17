
function [pi_evolution] = matrixEvolution(pi_0,P)
% Compute the Limiting distribution starting from an initial state
% distribution.

%% Initialisation.

k = 1000000;

%% Treatment.

pi_evolution = pi_0;

for n = 1:k
    pi_evolution = pi_evolution*P;
end
end