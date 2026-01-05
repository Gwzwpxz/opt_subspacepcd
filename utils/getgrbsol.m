function [z, tau, t, pi] = getgrbsol(A1s, A2s, Ds, rhs1, rhs2, mode)

if nargin < 6
    mode = 1;
end % End if

k = size(A1s, 2);
n = size(Ds, 1);
m1 = size(A1s, 1);
m2 = size(A2s, 1);

model.A = sparse([A1s; A2s; Ds, zeros(n, 1)]);
model.rhs = [rhs1; rhs2; zeros(n, 1)];
model.obj = zeros(k, 1);
model.obj(end) = 1;
model.modelsense = 'max';
model.sense = [repelem('>', m1, 1); repelem('<', m2, 1); repelem('>', n, 1)];
model.lb = -inf(k, 1);
model.lb(end) = 0.0;

if mode == 2
    param.Presolve = 0;
    param.Method = 0;
    param.FeasibilityTol = 1e-09;
    param.OptimalityTol = 1e-09;
end % End if

model.obj = model.obj;
param.ObjScale = 1e-03;
param.Presolve = 0;
param.Method = 2;
param.OutputFlag = 0;
param.Crossover = 0;
% param.BarHomogeneous = 1;
param.BarConvTol = 1e-10;
param.NumericFocus = 3;

sol = gurobi(model, param);

try
    z = sol.x(1:end - 1);
catch
    param.Method = 0;
    sol = gurobi(model, param);
    z = sol.x(1:end - 1);
end % End try

t = sol.runtime;
tau = sol.x(end);
pi = sol.pi(1:m1 + m2);

end % End function