clear; close all;

% Load data
load(fullfile('.', 'data', 'diag-bench-500-1.000e-01-opt.mat'));
X = data.X;

% Setup problem
param.ptype = 'S';
% param.method = "SDP"; % Use SDP to solve the problem
param.method = "SINF"; % Use Semi-inf solver to solve the problem
% Add identity and Jacobian preconditioner
n = size(X, 1);
param.subspace = [ones(n, 1), full(diag(X))];
DSprob = getoptprob(X, param);

% Solve problem
DSopt = optprecond(DSprob);

% Get performance
fprintf("Condition number of X                            : %5.2e \n", cond(full(X)));
fprintf("Condition number of DX for subspace constrained D: %5.2e \n", cond(full(DSopt.pX)));