%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPFL | MGT-418: Convex Optimization | Project 5, Question 3 %
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
% X = ...
% Lambda1 = ...
% Lambda2 = ...

% Objective function
% obj = ...

% Constraints enforcing observed ratings
% constr = ...

% Constraints enforcing range beween 1 and 10
% ...

% Semidefinite constraints
% ...

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
%f_obs = ...
%f_pm0 = ...
%f_pm1 = ...
%f_pm2 = ...
