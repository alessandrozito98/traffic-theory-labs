clear;
close all;
rng(9); % random number generator with seed = 9
lambda = 2; % arrival rate expressed in [arrivals/s]
nb_dropped_points = 1000;
nb_instances = 100000; 

% generating nb_instances of a Poisson Process with lambda = 2
inter_arr=exprnd((1/lambda).*ones(nb_dropped_points,nb_instances));
arr_epochs = cumsum(inter_arr,1);

% two intervals expressed in [s]
interval_1 = [0.5, 20.5];
interval_2 = [15.5, 45.5];

% number of arrivals in each interval
nb_arrivals_int_1 = sum(arr_epochs >= interval_1(1) & arr_epochs <= interval_1(2),1);
nb_arrivals_int_2 = sum(arr_epochs >= interval_2(1) & arr_epochs <= interval_2(2),1);

% average number of arrivals in each interval
m_1 = mean(nb_arrivals_int_1)
m_2 = mean(nb_arrivals_int_2)

hist_edges = [0:100];

% figure 1 referring to interval 1
figure
hold on
histo_1 = histogram(nb_arrivals_int_1,hist_edges,'Normalization','probability');
theoretical_1 = poisspdf(hist_edges(1:end-1),lambda.*(interval_1(2)-interval_1(1)));
plot(hist_edges(1:end-1)+0.5,theoretical_1,'LineWidth',2)
legend('Empirical','Theoretical')
title('Distribution of the number of arrivals in interval 1')
xlabel('Number of arrivals')
ylabel('Probability')
hold off

% figure 2 referring to interval 2
figure
hold on
histo_2 = histogram(nb_arrivals_int_2,hist_edges,'Normalization','probability');
theoretical_2 = poisspdf(hist_edges(1:end-1),lambda.*(interval_2(2)-interval_2(1)));
plot(hist_edges(1:end-1)+0.5,theoretical_2,'LineWidth',2)
legend('Empirical','Theoretical')
title('Distribution of the number of arrivals in interval 2')
xlabel('Number of arrivals')
ylabel('Probability')
hold off

% joint pmf distribution
joint = zeros(length(hist_edges)-1);    % empty matrix 100x100
% populating the matrix
for i=1:nb_instances
    joint(nb_arrivals_int_1(i)+1, nb_arrivals_int_2(i)+1) = joint(nb_arrivals_int_1(i)+1, nb_arrivals_int_2(i)+1) + 1;
end
joint_pmf = joint ./ sum(joint(:)); % normalization

% plotting the joint pmf
figure
imagesc(joint_pmf)
colorbar
title('Joint pmf of the two intervals')
xlabel('# arrivals in interval 2')
ylabel('# arrivals in interval 1')


% product of marginals pmfs
% retrieving values from the histograms previously computed
pmf_int_1 = histo_1.Values; 
pmf_int_2 = histo_2.Values;

product_marg_pmf = pmf_int_1' .* pmf_int_2;

% plotting the product of the two marginal pmfs
figure
imagesc(product_marg_pmf)
colorbar
title('Product of the two marginal pmfs')
xlabel('# arrivals in interval 2')
ylabel('# arrivals in interval 1')

% Mean-squared error between joint pmf and the product of the two marginal pmfs
mse = immse(joint_pmf,product_marg_pmf);
