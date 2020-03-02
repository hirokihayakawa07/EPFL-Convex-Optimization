%%%%%%%%%%%%%%%%%%%%% MGT-418 Convex Optimization %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Project 2 / Answer 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       %% Solve the SVM problem %%

clear;
%load('p2data1.mat');
load('p2data2.mat'); %uncomment to solve for the second data set

rho = 0.0001; %regularization parameter
m = length(y);

%%
%Solve SVM problem with smooth Hinge loss to compute the SVM coefficients 
%w and b (denote them by w and b and describe w as a column vector)

%Type your code here...
% define problem parameters
d = size(x,2);

% define decision variables
w = sdpvar(d,1);
b = sdpvar(1,1);
t = sdpvar(1,m);
s = sdpvar(1,m);

% define objective function
objective = 1/m * sum(s) + rho * w.'* w;

% initialize constraints
constraints = [];

% add constraints
for i=1:m
    constraints = [constraints, 1/2*t(i)^2+1-y(i)*(w.'* x(i,:).' - b) + t(i) <= s(i)];
    constraints = [constraints, 1/2*t(i)^2 <= s(i)];
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
opt_w = value(w);
disp(opt_w);
opt_b = value(b);
disp(opt_b);
opt_t = value(t);
disp(opt_t);
opt_s = value(s);
disp(opt_s);

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
label = opt_w.' * feature.' - opt_b;
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
