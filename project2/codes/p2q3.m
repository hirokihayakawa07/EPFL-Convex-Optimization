%%%%%%%%%%%%%%%%%%%%%% MGT-418 Convex Optimization %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Project 2 / Question 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Solve the SVM problem using the Gaussian kernel function %% 
     
clc;
clearvars;
load('p2data2.mat');

rho = 0.0001; %regularization parameter
sigma = 3; %bandwidth of the Gaussian kernel
m = length(y);

%%
%Solve the dual problem (4) with the Gaussian kernel 
%Denote the dual decision variables by lambda

%Type your code here...
% define problem parameters
d = size(x,2);

% define decision variables
lambda = sdpvar(1,m);

% define objective function
% yalmip assumes min problem. So add minus on the objective function of the dual problem
% dual summation part of the objective function is calculated in function
% calc
objective = - ( sum(lambda - m/2 * lambda.^2) - 1/(2*rho) * calc(sigma, m, lambda, x, y) );

% initialize constraints
constraints = [];

% add equality constraints
constraints = [constraints, lambda * y == 0];

% add inequality constraints
for i=1:m
    constraints = [constraints, 0 <= lambda(i)];
    constraints = [constraints, lambda(i) <= 1/m];
end

% specify solver settings
opt_settings = sdpsettings('solver', 'gurobi', 'verbose', 0);

% run solver
diagnosis = optimize(constraints, objective, opt_settings);

% display solver report
disp('solver report:');
disp(diagnosis);

% retrieve and display optimal objective value
disp('optimal objective value:');
opt_objective = value(objective);
disp(opt_objective);

% retrieve and display optimal solution values
disp('optimal solution values:');
opt_lambda = value(lambda);
disp(opt_lambda);

%%
%Compute optimal b (denote by b_opt) using the optimal dual solution

%Type your code here...
%Compute optimal b from the equation shown in quenstion 3.4
b_opt = 0;
for i = 1:m
    b_opt = b_opt + 1/rho * opt_lambda(i) * y(i) * exp(- (x(i,:) - x(1,:))*(x(i,:) - x(1,:)).' / (2*sigma^2));
end
b_opt = b_opt + y(1)*(m*opt_lambda(1) - 1);

%%
%Discretize each feature range to 100 discretization points to get 100^d
%total number of discretization points in the feature space
%Construct a feature matrix (denote by feature) of discretization points
%Specifically, feature will be a matrix in R^((100^d) x d), where each row
%represents a distinct feature vector
%Compute the label of each discrete point by using optimal w and b
%Construct a label vector (denoted by label) containing the respective labels
%Specifically, label will be a vector in R^((100^d) x 1)

%Type your code here...
%Take max and min value of x
mx = max(x);
mn = min(x);
%Create 100 discretization points between max x and min x
f1 = [mn(1):(mx(1) - mn(1))/99:mx(1)];
f2 = [mn(2):(mx(2) - mn(2))/99:mx(2)];
%Construct feature matrix
feature = combvec(f1,f2).';
%Compute the label of each point in the feature matrix
label = zeros(size(feature,1),1);
for i = 1:size(feature,1)
    for j = 1:m
        label(i) = label(i) + opt_lambda(j) * y(j) * exp(-(x(j,:) - feature(i,:)) * (x(j,:) - feature(i,:)).' /(2*sigma^2));
    end
end
label = 1/rho*label - b_opt;
% label = label > 0;
% label = (label*2 -1).';

%%
%Visualization - feel free to comment out and construct your own plots.

figure;
mf = length(feature);
colorMap = [zeros(mf , 1), zeros(mf , 1), ones(mf ,1)];
for k = 1 : mf
    if label(k) >= 0
        if label(k) >= 1
            colorMap(k, :) = [1,0.3,0]; % Light red
        else
            colorMap(k, :) = [1,0,0]; % Dark red
        end
    else
        if label(k) <= -1
            colorMap(k, :) = [0,0.3,1]; % Light blue
        else
            colorMap(k, :) = [0,0,1]; % Dark blue
        end
    end
end
scatter(feature(:,1), feature(:,2),23* ones(mf, 1), colorMap, '.', 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.4); grid on; hold on;

colorMap = [zeros(m, 1), zeros(m, 1), ones(m,1)];
for k = 1 : m
  if y(k) >= 1
    colorMap(k, :) = [1,0,0]; % Red
  else
    colorMap(k, :) = [0,0,1]; % Blue
  end
end
scatter(x(:,1), x(:,2),35* ones(m, 1), colorMap, 'filled'); grid on; hold on;