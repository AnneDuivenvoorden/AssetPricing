clear
close all

data_market  = readtable('dataMP.xlsx');
data_returns = readtable('dataConverted.xlsx');
variance_import = readtable('varianceSizeMomentum.xlsx');
sigma = table2array(variance_import);
rf = 0.0015;

gross_returns = table2array(data_returns(:,2:26));
gross_returns = gross_returns + rf;


one_vector = ones(25,1);
mu = transpose(mean(gross_returns));
mu_excess = mu - rf;
inv_sigma = inv(sigma);

A = transpose(mu)*inv_sigma*mu;
B = transpose(mu)*inv_sigma*one_vector;
C = transpose(one_vector)*inv_sigma*one_vector;

portfolio_gmv = (1/C)*inv_sigma*one_vector;
portfolio_mu = (1/B)*inv_sigma*mu;

mean_gmv = B/C;
vol_gmv = sqrt(1/C);
mean_mu = transpose(mu)*portfolio_mu;
vol_mu = sqrt(transpose(portfolio_mu)*sigma*portfolio_mu);

start = -1;
dx = 0.01;
N = 400;
mu_target = transpose(start + (0:N-1)*dx);
nr_inef = floor((mean_gmv-start)/0.01) + 1;
nr_ef = N - nr_inef;

lambda = zeros(N,1);
vol_lambda = zeros(N,1);
for i = 1:N
    lambda(i,1) = (B*C*mu_target(i,1)-B^2)/(A*C - B^2);
    vol_lambda(i,1) = sqrt((A-2*B*mu_target(i,1)+C*mu_target(i,1)^2)/(A*C-B^2));
end

portfolio_tang = 1/(transpose(one_vector)*inv_sigma*mu_excess)*inv_sigma*mu_excess;
mean_tang = transpose(mu)*portfolio_tang;
vol_tang = sqrt(transpose(portfolio_tang)*sigma*portfolio_tang);

riskfree_target = [-1; rf; 2.75];
vol_riskfree = zeros(3,1);
for i = 1:3
    vol_riskfree(i,1) = abs(riskfree_target(i,1)-rf)/(sqrt(transpose(mu_excess)*inv_sigma*mu_excess));
end

mu_test = mu_target(nr_inef:N, 1);
%plot
figure
%plot(vol_lambda_inef, mu_target(1:nr_inef, 1), 'b');
plot(vol_lambda(1:nr_inef,1), mu_target(1:nr_inef,1), 'b', 'DisplayName', 'Inefficient part, risky assets');
axis([0 3 -1 2]);
xlabel('Volatility');
ylabel('Gross return');
hold on
plot(vol_lambda(nr_inef:N-1,1), mu_target(nr_inef:N-1, 1), 'color', '#EDB120', 'DisplayName', 'Efficient part, risky assets');
plot(vol_riskfree(1:2,1), riskfree_target(1:2,1),'b--', 'DisplayName', 'Inefficient part');
plot(vol_riskfree(2:3,1), riskfree_target(2:3,1),'--', 'color', '#EDB120', 'DisplayName', 'Efficient part');
plot(vol_gmv, mean_gmv, '*', 'color', '#00841a', 'MarkerFaceColor', '#00841a', 'MarkerSize', 8, 'DisplayName', 'GMV portfolio');
plot(vol_mu, mean_mu, '*r', 'MarkerFaceColor', 'r','MarkerSize', 8, 'DisplayName', 'Portfolio pi-mu');
plot(vol_tang, mean_tang, 'db', 'MarkerSize', 8, 'DisplayName', 'Tangency portfolio');
legend('Location','northwest', 'Orientation', 'horizontal');
lgd = legend;
lgd.NumColumns = 1;


