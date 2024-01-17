%% Initialization.

% Different initial state distributions for P1.

P1_01 = [0.1 0.3 0.4 0.2];
P1_02 = [0.8 0.2 0 0];
P1_03 = [0.75 0.1 0.1 0.05];
P1_04 = [0.4 0.4 0.2 0];

% Transition probability matrix P1.

P1 = [0.4 0.6 0 0; 0.7 0.3 0 0; 0 0 0.5 0.5; 0 0 0.5 0.5];

% Different initial state distributions for P2.

P2_01 = [0.2 0.4 0.3 0.1 0];
P2_02 = [0.2 0.5 0.1 0.2 0];
P2_03 = [0 1 0 0 0];
P2_04 = [0.4 0.4 0 0.1 0.1];

% Transition probability matrix P2.

P2 = [0.4 0.6 0 0 0; 0.1 0.7 0.2 0 0; 0 0 1 0 0; 0 0 0.3 0.1 0.6; 
    0 0 0 0.5 0.5];

%% Treatment P1

% Analyze the limit. Check if there is a convergence.

num_iterations = 1000;
limit_P1 = P1;

for i = 1:num_iterations
    limit_P1 = limit_P1 * limit_P1^(i-1);
end

% Analyze the Limiting Distribution starting from different initial
% state distributions

% From PO1

p1_evolution1 = (matrixEvolution(P1_01,P1));

% From PO2

p1_evolution2 = (matrixEvolution(P1_02,P1));

% From PO3

p1_evolution3 = (matrixEvolution(P1_03,P1));

% From PO4

p1_evolution4 = (matrixEvolution(P1_04,P1));

% Since all the evolutions don't converge, there is no Steady State
% Distribution.

% Since P1 is not irreducible because we are blocked into state {1,2} or 
% {3,4}, we don't have a unique Stationary Distribution.



%% P2.


% Analyze the limit. Check if there is a convergence.

num_iterations = 1000;
limit_P2 = P2;

for i = 1:num_iterations
    limit_P2 = limit_P2 * limit_P2^(i-1);
end

% Therefore, the limit exists which means that we have a limiting
% distribution. Let's check for the limiting distribution starting from
% initial distributions.

% Analyze the Limiting Distribution starting from different initial
% state distributions

% From PO1

p2_evolution1 = (matrixEvolution(P2_01,P2));

% From PO2

p2_evolution2 = (matrixEvolution(P2_02,P2));

% From PO3

p2_evolution3 = (matrixEvolution(P2_03,P2));

% From PO4

p2_evolution4 = (matrixEvolution(P2_04,P2));


% For every initial state probability distribution evolves to the same
% state probability distribution. Since the limiting distribution
% converges, we have a Steady State Dsitribution.

p2_steady_state_distribution = steadyStateDistribution(P2);

% Since P2 is irreducible because we always go into the third state, we
% have a unique Stationary Distribution. 

% Compute the Stationary Distribution by solving GBE in Matlab.

% It's a reducible MC, therefore we can't compute a single unique
% stationary probability distribution pi. Let's find one with the Global
% Balance Equation (GBE).

A = (P2 - eye(size(P2)))';
b = zeros(1,length(P2))';

An = [A(1:end-1,:); ones(1,length(P2))];
bn = [b(1:end-1); 1];

pi_transp = An \ bn;
pi = pi_transp';
