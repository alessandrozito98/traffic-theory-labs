function [P_steady_state] = steadyStateDistribution(P)
    % Provide steady-state distribution
    % P: Transition matrix
    % Perform the necessary computations and return the steady-state distribution
    [eigenvector, ~] = eig(P', "vector");
    P_steady_state = eigenvector(:, 1) / sum(eigenvector(:, 1));
end