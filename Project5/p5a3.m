%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPFL | MGT-418: Convex Optimization | Project 5, Solution 3 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;
clc;

%% Load observed review matrix (either '20c50m' or '50c200m')
load 20c50m;
%load 50c200m; % !! long runtime - use only once code is stable !!
[m,n] = size(R_masked);
tic;

%% Question 3.1

% Decision variables
X = sdpvar(m,n);        
Lambda1 = sdpvar(n,n);
Lambda2 = sdpvar(m,m);

% Objective function
obj = trace(Lambda1) + trace(Lambda2);

% Constraints enforcing observed ratings
constr = [];
for i = 1:m
    for j = 1:n
        if R_masked(i,j) > 0
            constr = [constr; X(i,j) == R_masked(i,j)];
        end
    end
end

% Constraints enforcing range beween 1 and 10
constr = [constr; 1 <= X(:); X(:) <= 10];

% Semidefinite constraints
constr = [constr; [Lambda1, -0.5*X'; -0.5*X, Lambda2] >= 0];

%% Run the solver and retrieve rounded optimal solution

% Specify solver settings and run solver
ops = sdpsettings('solver', 'mosek', 'verbose',0);
diagnosis = optimize(constr, obj, ops);
runtime = toc;

% Retrieve and round optimal solution
X_star = round(value(X));

%% Load true review matrix (either '20c50m_truth' or '50c200m_truth')
load 20c50m_truth;
%load 50c200m_truth; % !! long runtime - use only once code is stable !!

%% Question 3.2

% Compute various performance metrics
f_obs = sum(sum(R_masked > 0)) / (m*n);
f_pm0 = sum(sum(X_star == R_truth)) / (m*n);
f_pm1 = f_pm0 + (sum(sum(X_star == R_truth-1)) + sum(sum(X_star == R_truth+1)) ) / (m*n);
f_pm2 = f_pm1 + (sum(sum(X_star == R_truth-2)) + sum(sum(X_star == R_truth+2)) ) / (m*n);
